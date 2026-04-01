// Модификация Сан-Хуана метода Греффе
// После k итераций Грефе корни малого модуля экспоненциально подавляются,
// позволяя аппроксимировать m доминирующих корней через усечённый полином.
//
// Алгоритм (Graeffe.pdf, раздел 8):
// 1. Выполнить k итераций Грефе
// 2. Извлечь модули из нормализованных коэффициентов
// 3. Для доминирующих m корней: из усечённого полинома степени m

#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <algorithm>
#include <complex>
#include <vector>

#include "NumericConcepts.h"
#include "NumericConstants.h"
#include "mathUtils.h"

namespace san_juan_detail {

// log(|a| + |b|) (log-sum-exp для двух значений)
template <RealNumber T>
T log_add(T log_a, T log_b, const T &NEG_INF) {
	if (log_a == NEG_INF)
		return log_b;
	if (log_b == NEG_INF)
		return log_a;
	T mx = max_val(log_a, log_b);
	T diff = (log_a + log_b - mx) - mx;
	if (diff < T(-50))
		return mx;
	return mx + log_val(T(1) + exp_val(diff));
}

// log(|a| - |b|), assumes log_a >= log_b
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

// Шаг Грефе в лог-домене (descending коэффициенты)
template <RealNumber T>
void graeffe_step_log(std::vector<T> &L, std::vector<int> &S, int n,
                      const T &NEG_INF) {
	std::vector<T> newL(n + 1, NEG_INF);
	std::vector<int> newS(n + 1, 1);
	T LOG2 = log_val(T(2));

	for (int k = 0; k <= n; ++k) {
		std::vector<T> pos_logs, neg_logs;

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

		T log_pos = NEG_INF;
		for (auto &lp : pos_logs)
			log_pos = log_add(log_pos, lp, NEG_INF);
		T log_neg = NEG_INF;
		for (auto &ln : neg_logs)
			log_neg = log_add(log_neg, ln, NEG_INF);

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

		if (k % 2 == 1)
			result_sign = -result_sign;

		newL[k] = result_log;
		newS[k] = result_sign;
	}

	L = newL;
	S = newS;
}

} // namespace san_juan_detail

/**
 * @brief Модификация Сан-Хуана: нахождение m доминирующих корней.
 *
 * @tparam T Вещественный тип (double, long double, float_precision).
 * @param coeffs_desc Коэффициенты (descending: старшая степень первой).
 * @param m Количество доминирующих корней (1 <= m < n). 0 = автоопределение.
 * @param graeffe_iters Количество итераций Грефе.
 * @param epsilon Точность.
 * @return Вектор модулей m доминирующих корней (по убыванию).
 */
template <RealNumber T>
std::vector<T> find_dominant_moduli_by_SanJuan(
    const std::vector<T> &coeffs_desc, int m = 0, int graeffe_iters = 8,
    T epsilon = numeric_constants::adaptive_epsilon<T>(
        numeric_constants::EPSILON_SCALE_PRECISE)) {

	if (coeffs_desc.size() < 2)
		return {};

	const T NEG_INF = neg_inf_val<T>();

	auto safe_log = [&](const T &x) -> T {
		T v = abs_val(x);
		if (is_zero_val(v))
			return NEG_INF;
		T result = log_val(v);
		return isfinite_val(result) ? result : NEG_INF;
	};

	std::vector<T> poly = coeffs_desc;
	size_t first_nz = 0;
	while (first_nz < poly.size() && is_zero_val(poly[first_nz]))
		++first_nz;
	if (first_nz >= poly.size())
		return {};
	poly.assign(poly.begin() + first_nz, poly.end());

	int n = int(poly.size()) - 1;
	if (n <= 0)
		return {};
	if (m <= 0 || m >= n)
		m = (n > 2) ? n / 2 : 1;

	// Нормализация
	T lead = poly[0];
	for (auto &c : poly)
		c /= lead;

	// Лог-коэффициенты
	std::vector<T> L(n + 1);
	std::vector<int> S(n + 1);
	for (int j = 0; j <= n; ++j) {
		L[j] = safe_log(poly[j]);
		S[j] = (poly[j] >= T(0)) ? 1 : -1;
	}

	// Грефе-итерации
	for (int r = 0; r < graeffe_iters; ++r) {
		san_juan_detail::graeffe_step_log(L, S, n, NEG_INF);
		if (L[0] != NEG_INF) {
			T lead_log = L[0];
			for (int j = 0; j <= n; ++j)
				if (L[j] != NEG_INF)
					L[j] -= lead_log;
		}
	}

	T pow2k = pow_val(T(2), T(graeffe_iters));

	// Извлечение m доминирующих модулей из первых m+1 коэффициентов
	std::vector<T> moduli;
	for (int j = 1; j <= m; ++j) {
		if (L[j] != NEG_INF && L[j - 1] != NEG_INF) {
			T log_mod = (L[j] - L[j - 1]) / pow2k;
			T mod = exp_val(log_mod);
			if (isfinite_val(mod) && mod > T(0))
				moduli.push_back(mod);
		}
	}

	std::sort(moduli.begin(), moduli.end(),
	          [](const T &a, const T &b) { return a > b; });

	return moduli;
}

/**
 * @brief Полная версия Сан-Хуана: находит ВСЕ модули корней.
 *
 * @return Вектор всех модулей корней (по убыванию).
 */
template <RealNumber T>
std::vector<T> find_all_moduli_by_SanJuan(
    const std::vector<T> &coeffs_desc, int graeffe_iters = 10,
    T epsilon = numeric_constants::adaptive_epsilon<T>(
        numeric_constants::EPSILON_SCALE_PRECISE)) {

	if (coeffs_desc.size() < 2)
		return {};

	const T NEG_INF = neg_inf_val<T>();

	auto safe_log = [&](const T &x) -> T {
		T v = abs_val(x);
		if (is_zero_val(v))
			return NEG_INF;
		T result = log_val(v);
		return isfinite_val(result) ? result : NEG_INF;
	};

	std::vector<T> poly = coeffs_desc;
	size_t first_nz = 0;
	while (first_nz < poly.size() && is_zero_val(poly[first_nz]))
		++first_nz;
	if (first_nz >= poly.size())
		return {};
	poly.assign(poly.begin() + first_nz, poly.end());

	int n = int(poly.size()) - 1;
	if (n <= 0)
		return {};

	T lead = poly[0];
	for (auto &c : poly)
		c /= lead;

	std::vector<T> L(n + 1);
	std::vector<int> S(n + 1);
	for (int j = 0; j <= n; ++j) {
		L[j] = safe_log(poly[j]);
		S[j] = (poly[j] >= T(0)) ? 1 : -1;
	}

	for (int r = 0; r < graeffe_iters; ++r) {
		san_juan_detail::graeffe_step_log(L, S, n, NEG_INF);
		if (L[0] != NEG_INF) {
			T lead_log = L[0];
			for (int j = 0; j <= n; ++j)
				if (L[j] != NEG_INF)
					L[j] -= lead_log;
		}
	}

	T pow2k = pow_val(T(2), T(graeffe_iters));

	std::vector<T> moduli;
	for (int j = 1; j <= n; ++j) {
		if (L[j] != NEG_INF && L[j - 1] != NEG_INF) {
			T log_mod = (L[j] - L[j - 1]) / pow2k;
			T mod = exp_val(log_mod);
			if (isfinite_val(mod) && mod > T(0))
				moduli.push_back(mod);
		}
	}

	std::sort(moduli.begin(), moduli.end(),
	          [](const T &a, const T &b) { return a > b; });

	return moduli;
}
