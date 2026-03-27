// Вспомогательные функции для метода цепных дробей
// Теория + основная реализация: Шаймарданов Санджар, КМБО-07-23
// Интеграция в фреймворк + исправление ошибок + настройка шаблонности для
// пользовательских типов: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#include "NumericConstants.h"
#include "mathUtils.h"
#include "polynomialUtils.h"
#include <algorithm>
#include <cmath>
#include <vector>

// Сдвиг полинома: P(x) -> P(x + a) по схеме Горнера (ascending порядок)
template <typename T>
std::vector<T> cf_taylor_shift(const std::vector<T> &coeffs, T a) {
	std::vector<T> res = coeffs;
	int n = int(res.size());
	for (int i = n - 2; i >= 0; --i) {
		for (int j = i; j < n - 1; ++j)
			res[j] += a * res[j + 1];
	}
	return res;
}

// Обратный полином: x^n * P(1/x) = реверс коэффициентов
template <typename T>
std::vector<T> cf_reciprocal(const std::vector<T> &coeffs) {
	std::vector<T> res(coeffs.rbegin(), coeffs.rend());
	return res;
}

// Нижняя граница наименьшего положительного корня (Local-Max-Quadratic bound)
// Возвращает floor этой границы (>= 1 если есть положительные корни)
template <typename T> T compute_partial_quotient(const std::vector<T> &coeffs) {
	if (coeffs.size() < 2)
		return T(1);

	int n = int(coeffs.size()) - 1;
	T threshold = numeric_constants::adaptive_epsilon<T>(
		numeric_constants::EPSILON_SCALE_PRECISE);

	// Старший коэффициент
	T lead = coeffs[n];
	if (abs_val(lead) < threshold)
		return T(1);

	// Граница Коши для обратного полинома (нижняя граница положительных корней)
	// LB = 1 / upper_bound(x^n * P(1/x))
	// Простая оценка: a_k = наименьшая положительная нижняя граница
	// Используем формулу Коши: max(|a_i/a_n|)^{1/(n-i)} для отрицательных коэфф.
	T best = T(0);
	for (int i = 0; i < n; ++i) {
		if ((lead > T(0) && coeffs[i] < T(0)) ||
		    (lead < T(0) && coeffs[i] > T(0))) {
			T ratio = abs_val(coeffs[i]) / abs_val(lead);
			T exponent = T(1) / T(n - i);
			T bound = pow_val(ratio, exponent);
			if (bound > best)
				best = bound;
		}
	}

	if (best < T(1))
		return T(1);

	T result = floor_val(best);
	if (result < T(1))
		result = T(1);
	return result;
}

// Вычисление подходящей дроби по частичным частным цепной дроби
template <typename T>
T build_convergent(const std::vector<T> &a, T &approx, T eps) {
	if (a.empty()) {
		approx = T(0);
		return approx;
	}

	T Pm2 = T(1), Pm1 = a[0];
	T Qm2 = T(0), Qm1 = T(1);
	T prev_approx = Pm1 / Qm1;

	for (size_t k = 1; k < a.size(); ++k) {
		T P = a[k] * Pm1 + Pm2;
		T Q = a[k] * Qm1 + Qm2;

		if (abs_val(Q) < eps) {
			approx = prev_approx;
			return approx;
		}

		approx = P / Q;

		if (abs_val(approx - prev_approx) < eps)
			break;
		prev_approx = approx;

		Pm2 = Pm1;
		Pm1 = P;
		Qm2 = Qm1;
		Qm1 = Q;
	}

	return approx;
}

// Удаление рациональных корней (целые кандидаты в [-100, 100])
template <typename T>
std::vector<T>
remove_rational_roots(std::vector<T> &coeffs, std::vector<T> &rational_roots,
                      T eps = numeric_constants::adaptive_epsilon<T>(
                          numeric_constants::EPSILON_SCALE_STANDARD)) {
	std::vector<T> result = coeffs;
	for (int i = -100; i <= 100; ++i) {
		T r = T(i);
		int deflate_count = 0;
		while (abs_val(eval_poly(result, r)) < eps) {
			rational_roots.push_back(r);
			result = deflate_poly(result, r);
			if (++deflate_count > 20)
				break;
		}
	}
	return result;
}
