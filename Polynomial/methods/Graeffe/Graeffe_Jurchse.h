// Модификация Юркша метода Греффе
// Нормализация коэффициентов на каждой итерации для предотвращения переполнения.
// Полностью лог-доменная реализация шага Грефе.
// Модули корней: |x_j| = exp((L[j] - L[j-1]) / 2^r)
// Реализация по: Graeffe.pdf, раздел 3
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <algorithm>
#include <cmath>
#include <vector>

#include "NumericConcepts.h"
#include "NumericConstants.h"
#include "mathUtils.h"

namespace jurchse_detail {

// log(|a| + |b|) при известных log|a|, log|b| (оба одного знака)
template <RealNumber T>
T log_add(T log_a, T log_b, const T &NEG_INF) {
	if (log_a == NEG_INF)
		return log_b;
	if (log_b == NEG_INF)
		return log_a;
	T mx = max_val(log_a, log_b);
	T mn = (log_a + log_b) - mx; // min without branch
	T diff = mn - mx;
	if (diff < T(-50))
		return mx;
	return mx + log_val(T(1) + exp_val(diff));
}

// log(|a| - |b|) при log|a| > log|b| (оба положительные, a > b)
template <RealNumber T>
T log_sub(T log_a, T log_b, const T &NEG_INF) {
	if (log_b == NEG_INF)
		return log_a;
	if (log_a == NEG_INF)
		return NEG_INF;
	T diff = log_b - log_a;
	if (diff < T(-50))
		return log_a;
	T arg = T(1) - exp_val(diff);
	if (arg <= T(0))
		return NEG_INF;
	return log_a + log_val(arg);
}

// Шаг Грефе полностью в лог-домене (descending коэффициенты)
// L[i] = log|a[i]|, S[i] = sign(a[i])
// Классическая формула:
// b[k] = (-1)^k * [a[k]^2 + 2*sum_{j=1}^{min(k,n-k)} (-1)^j * a[k-j]*a[k+j]]
template <RealNumber T>
void graeffe_step_log(std::vector<T> &L, std::vector<int> &S, int n,
                      const T &NEG_INF) {
	std::vector<T> newL(n + 1, NEG_INF);
	std::vector<int> newS(n + 1, 1);
	T LOG2 = log_val(T(2));

	for (int k = 0; k <= n; ++k) {
		// Собираем все слагаемые: term_0 = a[k]^2 (sign: +1)
		// term_j = 2*(-1)^j * a[k-j]*a[k+j]  (j = 1..min(k,n-k))

		// Группируем по знакам
		std::vector<T> pos_logs, neg_logs;

		// Слагаемое a[k]^2 — всегда положительное
		if (L[k] != NEG_INF)
			pos_logs.push_back(T(2) * L[k]);

		int lim = std::min(k, n - k);
		for (int j = 1; j <= lim; ++j) {
			if (L[k - j] == NEG_INF || L[k + j] == NEG_INF)
				continue;
			T log_term = LOG2 + L[k - j] + L[k + j];
			int cross_sign = S[k - j] * S[k + j];
			int factor_sign = (j % 2 == 0) ? 1 : -1;
			int total_sign = cross_sign * factor_sign;

			if (total_sign > 0)
				pos_logs.push_back(log_term);
			else
				neg_logs.push_back(log_term);
		}

		// log-sum-exp для каждой группы
		T log_pos = NEG_INF;
		for (auto &lp : pos_logs)
			log_pos = log_add(log_pos, lp, NEG_INF);

		T log_neg = NEG_INF;
		for (auto &ln : neg_logs)
			log_neg = log_add(log_neg, ln, NEG_INF);

		// Вычисляем log|pos - neg| и знак
		T result_log;
		int result_sign;

		if (log_pos == NEG_INF && log_neg == NEG_INF) {
			result_log = NEG_INF;
			result_sign = 1;
		} else if (log_neg == NEG_INF) {
			result_log = log_pos;
			result_sign = 1;
		} else if (log_pos == NEG_INF) {
			result_log = log_neg;
			result_sign = -1;
		} else if (log_pos > log_neg) {
			result_log = log_sub(log_pos, log_neg, NEG_INF);
			result_sign = 1;
		} else if (log_neg > log_pos) {
			result_log = log_sub(log_neg, log_pos, NEG_INF);
			result_sign = -1;
		} else {
			result_log = NEG_INF;
			result_sign = 1;
		}

		// Применяем (-1)^k
		if (k % 2 == 1)
			result_sign = -result_sign;

		newL[k] = result_log;
		newS[k] = result_sign;
	}

	L = newL;
	S = newS;
}

} // namespace jurchse_detail

/**
 * @brief Вычисление модулей корней методом Юркша.
 *
 * @tparam T Тип коэффициентов (RealNumber: double, long double, float_precision).
 * @param coeffs_desc Коэффициенты (descending: старшая степень первой).
 * @param epsilon Требуемая точность.
 * @param maxIter Максимальное число итераций Грефе.
 * @return Вектор модулей корней, отсортированный по убыванию.
 */
template <RealNumber T>
std::vector<T> find_moduli_by_Jurchse(
    const std::vector<T> &coeffs_desc,
    T epsilon = numeric_constants::adaptive_epsilon<T>(
        numeric_constants::EPSILON_SCALE_PRECISE),
    int maxIter = 40) {

	const T NEG_INF = neg_inf_val<T>();

	auto safe_log = [&](const T &x) -> T {
		T v = abs_val(x);
		if (is_zero_val(v))
			return NEG_INF;
		T result = log_val(v);
		return isfinite_val(result) ? result : NEG_INF;
	};

	// Пропуск ведущих нулей
	size_t first_nz = 0;
	while (first_nz < coeffs_desc.size() && is_zero_val(coeffs_desc[first_nz]))
		++first_nz;
	if (first_nz >= coeffs_desc.size())
		return {};

	std::vector<T> raw(coeffs_desc.begin() + first_nz, coeffs_desc.end());
	int n = int(raw.size()) - 1;
	if (n <= 0)
		return {};

	// Нормализация до мониального полинома
	T lead = raw[0];
	for (auto &c : raw)
		c /= lead;

	// Инициализация лог-коэффициентов (descending: raw[0]=1)
	std::vector<T> L(n + 1);
	std::vector<int> S(n + 1);
	for (int j = 0; j <= n; ++j) {
		L[j] = safe_log(raw[j]);
		S[j] = (raw[j] >= T(0)) ? 1 : -1;
	}

	// Аппроксимации модулей
	std::vector<T> logN(n, T(0));
	std::vector<T> prevLogN(n, T(0));

	int iter_done = 0;

	for (int r = 0; r < maxIter; ++r) {
		jurchse_detail::graeffe_step_log(L, S, n, NEG_INF);

		// Нормализация Юркша: L[0] -> 0 (мониальный полином)
		if (L[0] != NEG_INF) {
			T lead_log = L[0];
			for (int j = 0; j <= n; ++j) {
				if (L[j] != NEG_INF)
					L[j] -= lead_log;
			}
		}

		iter_done = r + 1;
		T pow2k = pow_val(T(2), T(iter_done));

		prevLogN = logN;
		for (int j = 0; j < n; ++j) {
			if (L[j + 1] != NEG_INF && L[j] != NEG_INF)
				logN[j] = (L[j + 1] - L[j]) / pow2k;
			else if (L[j + 1] != NEG_INF)
				logN[j] = L[j + 1] / pow2k;
		}

		// Проверка сходимости
		if (r >= 3) {
			int converged = 0;
			for (int j = 0; j < n; ++j) {
				T diff = abs_val(logN[j] - prevLogN[j]);
				T scale = max_val(abs_val(logN[j]), T(1));
				if (diff < epsilon * scale)
					++converged;
			}
			if (converged == n)
				break;
		}
	}

	std::vector<T> moduli;
	moduli.reserve(n);
	for (int j = 0; j < n; ++j) {
		T mod = exp_val(logN[j]);
		if (isfinite_val(mod) && mod > T(0))
			moduli.push_back(mod);
	}

	std::sort(moduli.begin(), moduli.end(),
	          [](const T &a, const T &b) { return a > b; });

	return moduli;
}
