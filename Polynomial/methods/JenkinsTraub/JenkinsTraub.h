// Метод Дженкинса-Трауба (CPOLY) для нахождения всех корней
// комплексного полинома. Трёхэтапный алгоритм с H-полиномами:
//   Stage 1: Без сдвига (No-shift) — разделение нулей по модулю
//   Stage 2: Фиксированный сдвиг (Fixed-shift) — сходимость к конкретному нулю
//   Stage 3: Переменный сдвиг (Variable-shift) — квадратичная сходимость
//
// Поддержка произвольной точности через mathUtils (float_precision и др.)
//
// Реализация по: C-Poly.pdf (Метод Дженкинса-Трауба: поиск корней
// комплексных полиномов, 2 февраля 2026 г.)

#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <algorithm>
#include <complex>
#include <limits>
#include <vector>

#include "JenkinsTraub_Helper.h"
#include "NumericConcepts.h"
#include "NumericConstants.h"
#include "mathUtils.h"

namespace cpoly_impl {

template <typename T>
T jt_machine_eps() {
	if constexpr (std::numeric_limits<T>::is_specialized)
		return std::numeric_limits<T>::epsilon();
	else
		return T(1) / T(numeric_constants::FACTOR_1B);
}

template <typename T>
T jt_pi() {
	return atan2_val(T(0), T(-1));
}

template <typename T>
std::vector<std::complex<T>>
update_H(const std::vector<std::complex<T>> &H,
         const std::vector<std::complex<T>> &P,
         const std::complex<T> &shift) {
	using Complex = std::complex<T>;

	Complex H_at_s = cpoly_detail::horner_eval(H, shift);
	Complex P_at_s = cpoly_detail::horner_eval(P, shift);

	T eps100 = jt_machine_eps<T>() * T(100);
	Complex c;
	if (abs_val(H_at_s) > eps100)
		c = -P_at_s / H_at_s;
	else
		c = Complex(T(0), T(0));

	int deg_H = int(H.size()) - 1;
	int deg_P = int(P.size()) - 1;
	int max_deg = (deg_H > deg_P) ? deg_H : deg_P;
	std::vector<Complex> T_poly(max_deg + 1, Complex(T(0), T(0)));

	for (int k = 0; k <= deg_H; ++k)
		T_poly[k] += H[k];
	for (int k = 0; k <= deg_P; ++k)
		T_poly[k] += c * P[k];

	return cpoly_detail::complex_deflate(T_poly, shift);
}

template <typename T>
std::vector<std::complex<T>>
stage1_no_shift(const std::vector<std::complex<T>> &P, int M) {
	using Complex = std::complex<T>;
	int n = int(P.size()) - 1;
	if (n < 1)
		return {};

	std::vector<Complex> H(n);
	for (int k = 0; k < n; ++k)
		H[k] = P[k + 1] * Complex(T(k + 1) / T(n), T(0));

	Complex P0 = P[0];
	Complex zero(T(0), T(0));
	T eps_mach = jt_machine_eps<T>();

	for (int l = 0; l < M; ++l) {
		if (abs_val(P0) < eps_mach)
			break;
		H = update_H(H, P, zero);
	}

	return H;
}

template <typename T>
bool stage23_find_root(const std::vector<std::complex<T>> &P,
                       std::vector<std::complex<T>> H,
                       std::complex<T> shift,
                       std::complex<T> &alpha,
                       T eps) {
	using Complex = std::complex<T>;
	const T eps_mach = jt_machine_eps<T>();
	const T eps100 = eps_mach * T(100);
	const T eps10 = eps_mach * T(10);

	const int stage2_max = 10;
	Complex prev_t(T(0), T(0));
	bool stage2_converged = false;

	for (int l = 0; l < stage2_max; ++l) {
		Complex H_at_s = cpoly_detail::horner_eval(H, shift);
		Complex P_at_s = cpoly_detail::horner_eval(P, shift);

		if (abs_val(P_at_s) < eps100) {
			alpha = shift;
			return true;
		}

		Complex t(T(0), T(0));
		if (abs_val(H_at_s) > eps100)
			t = -P_at_s / H_at_s;

		if (l > 0 && abs_val(t) > eps_mach) {
			T rel = abs_val(t - prev_t) / abs_val(t);
			if (rel < T(0.5))
				stage2_converged = true;
		}
		prev_t = t;

		H = update_H(H, P, shift);

		if (stage2_converged)
			break;
	}

	Complex s = shift;
	if (stage2_converged) {
		Complex H_at_s = cpoly_detail::horner_eval(H, shift);
		Complex P_at_s = cpoly_detail::horner_eval(P, shift);
		if (abs_val(H_at_s) > eps_mach)
			s = shift + P_at_s / H_at_s;
	}

	const int stage3_max = 50;
	for (int l = 0; l < stage3_max; ++l) {
		auto [P_at_s, dP_at_s] = cpoly_detail::horner_eval_and_deriv(P, s);

		if (abs_val(P_at_s) < eps100) {
			alpha = s;
			return true;
		}

		Complex H_at_s = cpoly_detail::horner_eval(H, s);
		Complex t(T(0), T(0));

		if (abs_val(H_at_s) > eps10)
			t = -P_at_s / H_at_s;
		else if (abs_val(dP_at_s) > eps_mach)
			t = P_at_s / dP_at_s;

		Complex new_s = s - t;

		if (abs_val(t) < eps * max_val(T(1), abs_val(s))) {
			alpha = new_s;
			return true;
		}

		H = update_H(H, P, s);
		s = new_s;
	}

	for (int l = 0; l < 100; ++l) {
		auto [pz, dpz] = cpoly_detail::horner_eval_and_deriv(P, s);
		if (abs_val(pz) < eps100) {
			alpha = s;
			return true;
		}
		if (abs_val(dpz) < eps_mach)
			break;
		Complex t = pz / dpz;
		s -= t;
		if (abs_val(t) < eps * max_val(T(1), abs_val(s))) {
			alpha = s;
			return true;
		}
	}

	alpha = s;
	return abs_val(cpoly_detail::horner_eval(P, s)) <
	       sqrt_val(eps) * max_val(T(1), abs_val(s));
}

} // namespace cpoly_impl

/**
 * @brief Поиск всех корней комплексного полинома методом Дженкинса-Трауба.
 *
 * Поддерживает стандартные типы (float, double, long double) и типы
 * повышенной точности (float_precision) через обёртки mathUtils.
 *
 * @tparam T Базовый вещественный тип, должен удовлетворять RealNumber.
 * @param input_coeffs Комплексные коэффициенты (ascending: input_coeffs[0] = a_0).
 * @param eps Требуемая точность.
 * @return Вектор комплексных корней.
 */
template <RealNumber T>
std::vector<std::complex<T>>
find_roots_by_JenkinsTraub(const std::vector<std::complex<T>> &input_coeffs,
                            T eps = numeric_constants::adaptive_epsilon<T>(
                                numeric_constants::EPSILON_SCALE_PRECISE)) {
	using Complex = std::complex<T>;
	using CPolynomial = std::vector<Complex>;

	std::vector<Complex> all_roots;

	if (input_coeffs.size() < 2)
		return all_roots;

	CPolynomial poly = input_coeffs;
	while (poly.size() > 1 && abs_val(poly.back()) < eps)
		poly.pop_back();

	if (poly.size() < 2)
		return all_roots;

	{
		Complex lead = poly.back();
		if (abs_val(lead) < eps)
			return all_roots;
		for (auto &c : poly)
			c /= lead;
	}

	const int degree = int(poly.size()) - 1;
	const T pi = cpoly_impl::jt_pi<T>();
	const T eps_mach = cpoly_impl::jt_machine_eps<T>();

	CPolynomial current_poly = poly;
	int remaining = degree;

	while (remaining > 1) {
		int n = int(current_poly.size()) - 1;

		T cauchy_r = T(0);
		for (int k = 0; k < n; ++k) {
			T ratio = pow_val(abs_val(current_poly[k] / current_poly[n]),
			                  T(1) / T(n - k));
			cauchy_r = max_val(cauchy_r, ratio);
		}
		if (cauchy_r < eps_mach)
			cauchy_r = T(1);

		CPolynomial H = cpoly_impl::stage1_no_shift<T>(current_poly, 5);

		bool root_found = false;
		Complex alpha;
		int num_shifts = 5;
		if (n > 5) num_shifts = n;

		for (int attempt = 0; attempt < num_shifts && !root_found; ++attempt) {
			T theta = T(2) * pi * (T(attempt) + T(0.7)) / T(num_shifts);
			Complex shift = Complex(cauchy_r * cos_val(theta),
			                        cauchy_r * sin_val(theta));

			root_found = cpoly_impl::stage23_find_root(
			    current_poly, H, shift, alpha, eps);
		}

		if (!root_found) {
			for (int attempt = 0; attempt < 20 && !root_found; ++attempt) {
				T theta = T(2) * pi * T(attempt) / T(20);
				T r = cauchy_r * (T(0.5) + T(attempt) * T(0.1));
				Complex start =
				    Complex(r * cos_val(theta), r * sin_val(theta));

				for (int iter = 0; iter < 200; ++iter) {
					auto [pz, dpz] =
					    cpoly_detail::horner_eval_and_deriv(current_poly, start);
					if (abs_val(pz) < eps_mach * T(100)) {
						alpha = start;
						root_found = true;
						break;
					}
					if (abs_val(dpz) < eps_mach)
						break;
					Complex t = pz / dpz;
					start -= t;
					if (abs_val(t) < eps * max_val(T(1), abs_val(start))) {
						alpha = start;
						root_found = true;
						break;
					}
				}
			}
		}

		if (!root_found) {
			alpha = Complex(cauchy_r * T(0.5), cauchy_r * T(0.3));
		}

		all_roots.push_back(alpha);
		current_poly = cpoly_detail::complex_deflate(current_poly, alpha);
		--remaining;
	}

	if (current_poly.size() == 2 && abs_val(current_poly[1]) > eps_mach) {
		all_roots.push_back(-current_poly[0] / current_poly[1]);
	}

	// Финальная полировка Ньютоном по исходному нормализованному полиному
	for (auto &root : all_roots) {
		for (int iter = 0; iter < 50; ++iter) {
			auto [pz, dpz] = cpoly_detail::horner_eval_and_deriv(poly, root);

			if (abs_val(pz) < eps_mach * T(100))
				break;
			if (abs_val(dpz) < eps_mach)
				break;

			Complex t = pz / dpz;
			root -= t;

			if (abs_val(t) < eps_mach * max_val(T(1), abs_val(root)))
				break;
		}
	}

	return all_roots;
}

/**
 * @brief Перегрузка для полиномов с вещественными коэффициентами.
 */
template <RealNumber T>
std::vector<std::complex<T>>
find_roots_by_JenkinsTraub(const std::vector<T> &real_coeffs,
                            T eps = numeric_constants::adaptive_epsilon<T>(
                                numeric_constants::EPSILON_SCALE_PRECISE)) {
	return find_roots_by_JenkinsTraub(cpoly_detail::to_complex(real_coeffs), eps);
}
