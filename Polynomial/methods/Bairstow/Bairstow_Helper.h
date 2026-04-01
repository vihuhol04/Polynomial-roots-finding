// Вспомогательные функции для метода Берстоу
// Синтетическое деление на квадратичный множитель (x^2 - r*x - s)
// и решение линейной системы 2x2

#pragma once

#include <complex>
#include <cstddef>
#include <vector>

#include "NumericConcepts.h"
#include "mathUtils.h"

/**
 * @brief Синтетическое деление P(x) на (x^2 - r*x - s).
 *
 * Для P(x) = a_0 + a_1*x + ... + a_n*x^n вычисляет массив b[0..n]:
 *   b[n] = a[n]
 *   b[n-1] = a[n-1] + r * b[n]
 *   b[k] = a[k] + r * b[k+1] + s * b[k+2],  k = n-2, ..., 0
 *
 * Остаток: b[1] (при x) и b[0] (свободный).
 * Частное: b[2], b[3], ..., b[n] (ascending, степень n-2).
 *
 * @param[out] b Полный массив b-коэффициентов (размер = coeffs.size()).
 */
template <RealNumber T>
void bairstow_synth_div(const std::vector<T> &coeffs, T r, T s,
                        std::vector<T> &b) {
	int n = int(coeffs.size()) - 1;
	b.resize(coeffs.size());

	b[n] = coeffs[n];
	if (n >= 1)
		b[n - 1] = coeffs[n - 1] + r * b[n];
	for (int k = n - 2; k >= 0; --k)
		b[k] = coeffs[k] + r * b[k + 1] + s * b[k + 2];
}

/**
 * @brief Решение линейной системы 2x2 методом Крамера.
 *
 * | a11 a12 | |x1|   |rhs1|
 * | a21 a22 | |x2| = |rhs2|
 *
 * @return true если система решена успешно (определитель != 0).
 */
template <RealNumber T>
bool solve_2x2(T a11, T a12, T a21, T a22, T rhs1, T rhs2, T &x1, T &x2) {
	T det = a11 * a22 - a12 * a21;
	T threshold = numeric_constants::adaptive_epsilon<T>(
		numeric_constants::EPSILON_SCALE_COARSE);
	if (abs_val(det) < threshold)
		return false;

	x1 = (rhs1 * a22 - rhs2 * a12) / det;
	x2 = (a11 * rhs2 - a21 * rhs1) / det;
	return true;
}

/**
 * @brief Извлечение корней из x^2 - r*x - s = 0.
 *
 * Дискриминант: D = r^2 + 4*s.
 * @return Пара комплексных корней.
 */
template <RealNumber T>
std::pair<std::complex<T>, std::complex<T>> solve_quadratic_rs(T r, T s) {
	T disc = r * r + T(4) * s;
	if (disc >= T(0)) {
		T sd = sqrt_val(disc);
		return {{(r + sd) / T(2), T(0)}, {(r - sd) / T(2), T(0)}};
	} else {
		T sd = sqrt_val(-disc);
		return {{r / T(2), sd / T(2)}, {r / T(2), -sd / T(2)}};
	}
}
