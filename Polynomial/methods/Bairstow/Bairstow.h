// Метод Берстоу для нахождения всех корней полинома с вещественными коэффициентами
// Итеративно извлекает квадратичные множители (x^2 - r*x - s) через
// двойное синтетическое деление и ньютоновскую коррекцию системы 2x2
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <complex>
#include <limits>
#include <vector>

#include "Bairstow_Helper.h"
#include "Helper_for_all_methods.h"
#include "NumericConstants.h"
#include "polynomialUtils.h"

/**
 * @brief Поиск всех корней полинома методом Берстоу.
 *
 * Стандартная формулировка: делим на (x^2 - r*x - s), вычисляем
 * b-коэффициенты, затем d-коэффициенты (та же рекурсия по b),
 * и решаем линейную систему:
 *
 *   | d[2]  d[3] | |dr|   | -b[1] |
 *   | d[1]  d[2] | |ds| = | -b[0] |
 *
 * @tparam T Тип коэффициентов, должен удовлетворять RealNumber.
 * @param input_coeffs Коэффициенты по возрастанию степеней.
 * @param eps Требуемая точность.
 * @return Вектор комплексных корней.
 */
template <RealNumber T>
std::vector<std::complex<T>>
find_roots_by_Bairstow(const std::vector<T> &input_coeffs,
                        T eps = numeric_constants::adaptive_epsilon<T>(
                            numeric_constants::EPSILON_SCALE_PRECISE)) {
	std::vector<std::complex<T>> all_roots;

	if (input_coeffs.size() < 2)
		return all_roots;

	std::vector<T> coeffs = input_coeffs;
	while (coeffs.size() > 1 && abs_val(coeffs.back()) < eps)
		coeffs.pop_back();

	if (coeffs.size() < 2)
		return all_roots;

	const int max_iter = 1000;
	const int max_restarts = 15;

	while (coeffs.size() > 3) {
		int n = int(coeffs.size()) - 1;

		// Начальное приближение r, s
		T r = T(0), s = T(0);
		if (abs_val(coeffs[n]) > eps) {
			r = coeffs[n - 1] / coeffs[n];
			s = (n >= 2) ? coeffs[n - 2] / coeffs[n] : T(0);
		}

		bool converged = false;

		for (int restart = 0; restart < max_restarts && !converged; ++restart) {
			for (int iter = 0; iter < max_iter; ++iter) {
				// Первое синтетическое деление: a -> b
				std::vector<T> b;
				bairstow_synth_div(coeffs, r, s, b);

				// Второе синтетическое деление (та же рекурсия): b -> d
				std::vector<T> d;
				bairstow_synth_div(b, r, s, d);

				// Якобиан (стандартная формула Burden & Faires):
				// | d[2]  d[3] | |dr|   | -b[1] |
				// | d[1]  d[2] | |ds| = | -b[0] |
				T d1 = (n >= 1) ? d[1] : T(0);
				T d2 = (n >= 2) ? d[2] : T(0);
				T d3 = (n >= 3) ? d[3] : T(0);

				T dr, ds;
				if (!solve_2x2(d2, d3, d1, d2, -b[1], -b[0], dr, ds)) {
					r += T(0.5) / T(restart + 1);
					s -= T(0.3) / T(restart + 1);
					break;
				}

				r += dr;
				s += ds;

				if (abs_val(dr) + abs_val(ds) < eps) {
					converged = true;
					break;
				}
			}

			if (!converged) {
				// Рестарт с другим начальным приближением
				r = T(restart + 1) * T(0.3) - T(1);
				s = -T(restart + 1) * T(0.2) + T(0.5);
			}
		}

		auto [r1, r2] = solve_quadratic_rs(r, s);
		all_roots.push_back(r1);
		all_roots.push_back(r2);

		// Дефляция: частное = b[2], b[3], ..., b[n]
		std::vector<T> b;
		bairstow_synth_div(coeffs, r, s, b);
		std::vector<T> quotient(b.begin() + 2, b.end());
		coeffs = quotient;
	}

	// Остаток: степень 2 или 1
	if (coeffs.size() == 3) {
		T a2 = coeffs[2], a1 = coeffs[1], a0 = coeffs[0];
		// x^2 + (a1/a2)*x + (a0/a2) = 0  ->  x^2 - r*x - s = 0 с r = -a1/a2, s = -a0/a2
		T r_last = -a1 / a2;
		T s_last = -a0 / a2;
		auto [r1, r2] = solve_quadratic_rs(r_last, s_last);
		all_roots.push_back(r1);
		all_roots.push_back(r2);
	} else if (coeffs.size() == 2) {
		T root = -coeffs[0] / coeffs[1];
		all_roots.push_back({root, T(0)});
	}

	// Финальная полировка Ньютоном по исходному полиному
	for (auto &root : all_roots) {
		for (int iter = 0; iter < 50; ++iter) {
			std::complex<T> pz(T(0), T(0));
			std::complex<T> dpz(T(0), T(0));
			int n = int(input_coeffs.size()) - 1;

			for (int k = n; k >= 0; --k) {
				dpz = dpz * root + pz;
				pz = pz * root + std::complex<T>(input_coeffs[k], T(0));
			}

			if (std::abs(dpz) < eps)
				break;

			std::complex<T> correction = pz / dpz;
			root -= correction;

			if (std::abs(correction) < eps * std::max(T(1), std::abs(root)))
				break;
		}
	}

	return all_roots;
}
