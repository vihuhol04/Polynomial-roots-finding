// Метод Дженкинса-Трауба (CPOLY) для нахождения всех корней
// комплексного полинома. Трёхэтапный алгоритм с H-полиномами:
//   Stage 1: Без сдвига (No-shift) — разделение нулей по модулю
//   Stage 2: Фиксированный сдвиг (Fixed-shift) — сходимость к конкретному нулю
//   Stage 3: Переменный сдвиг (Variable-shift) — квадратичная сходимость
//
// Реализация по: C-Poly.pdf (Метод Дженкинса-Трауба: поиск корней
// комплексных полиномов, 2 февраля 2026 г.)
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "JenkinsTraub_Helper.h"
#include "NumericConcepts.h"
#include "NumericConstants.h"

namespace cpoly_impl {

/**
 * @brief Обновление H-полинома: H_new = [H(z) + c*P(z)] / (z - shift)
 *
 * Где c = -P(shift)/H(shift). Деление на (z - shift) выполняется
 * синтетическим делением (Горнер). Результат — полином степени n-1.
 */
template <typename T>
std::vector<std::complex<T>>
update_H(const std::vector<std::complex<T>> &H,
         const std::vector<std::complex<T>> &P,
         const std::complex<T> &shift) {
	using Complex = std::complex<T>;

	Complex H_at_s = cpoly_detail::horner_eval(H, shift);
	Complex P_at_s = cpoly_detail::horner_eval(P, shift);

	Complex c;
	if (std::abs(H_at_s) > std::numeric_limits<T>::epsilon() * T(100))
		c = -P_at_s / H_at_s;
	else
		c = Complex(T(0), T(0));

	// T(z) = H(z) + c * P(z)
	int deg_H = int(H.size()) - 1;
	int deg_P = int(P.size()) - 1;
	int max_deg = std::max(deg_H, deg_P);
	std::vector<Complex> T_poly(max_deg + 1, Complex(T(0), T(0)));

	for (int k = 0; k <= deg_H; ++k)
		T_poly[k] += H[k];
	for (int k = 0; k <= deg_P; ++k)
		T_poly[k] += c * P[k];

	// Деление T(z) / (z - shift) синтетическим делением
	return cpoly_detail::complex_deflate(T_poly, shift);
}

/**
 * @brief Stage 1: No-shift (M итераций).
 *
 * H_0(z) = P'(z) / n
 * H_{l+1}(z) = [H_l(z) + c_l * P(z)] / z, где c_l = -H_l(0) / P(0)
 */
template <typename T>
std::vector<std::complex<T>>
stage1_no_shift(const std::vector<std::complex<T>> &P, int M) {
	using Complex = std::complex<T>;
	int n = int(P.size()) - 1;
	if (n < 1)
		return {};

	// H_0 = P'(z) / n
	std::vector<Complex> H(n);
	for (int k = 0; k < n; ++k)
		H[k] = P[k + 1] * Complex(T(k + 1) / T(n), T(0));

	Complex P0 = P[0];
	Complex zero(T(0), T(0));

	for (int l = 0; l < M; ++l) {
		if (std::abs(P0) < std::numeric_limits<T>::epsilon())
			break;
		H = update_H(H, P, zero);
	}

	return H;
}

/**
 * @brief Stage 2+3: Fixed-shift, затем variable-shift для одного корня.
 *
 * @return true если корень найден (записан в alpha).
 */
template <typename T>
bool stage23_find_root(const std::vector<std::complex<T>> &P,
                       std::vector<std::complex<T>> H,
                       std::complex<T> shift,
                       std::complex<T> &alpha,
                       T eps) {
	using Complex = std::complex<T>;
	const T eps_mach = std::numeric_limits<T>::epsilon();

	// Stage 2: Fixed-shift (до 10 итераций, проверяем сходимость t_l)
	const int stage2_max = 10;
	Complex prev_t(T(0), T(0));
	bool stage2_converged = false;

	for (int l = 0; l < stage2_max; ++l) {
		Complex H_at_s = cpoly_detail::horner_eval(H, shift);
		Complex P_at_s = cpoly_detail::horner_eval(P, shift);

		if (std::abs(P_at_s) < eps_mach * T(100)) {
			alpha = shift;
			return true;
		}

		Complex t(T(0), T(0));
		if (std::abs(H_at_s) > eps_mach * T(100))
			t = -P_at_s / H_at_s;

		// Проверка сходимости t
		if (l > 0 && std::abs(t) > eps_mach) {
			T rel = std::abs(t - prev_t) / std::abs(t);
			if (rel < T(0.5))
				stage2_converged = true;
		}
		prev_t = t;

		H = update_H(H, P, shift);

		if (stage2_converged)
			break;
	}

	// Stage 3: Variable-shift (Newton с H-полиномом)
	Complex s = shift;
	if (stage2_converged) {
		Complex H_at_s = cpoly_detail::horner_eval(H, shift);
		Complex P_at_s = cpoly_detail::horner_eval(P, shift);
		if (std::abs(H_at_s) > eps_mach)
			s = shift + P_at_s / H_at_s;
	}

	const int stage3_max = 50;
	for (int l = 0; l < stage3_max; ++l) {
		auto [P_at_s, dP_at_s] = cpoly_detail::horner_eval_and_deriv(P, s);

		if (std::abs(P_at_s) < eps_mach * T(100)) {
			alpha = s;
			return true;
		}

		Complex H_at_s = cpoly_detail::horner_eval(H, s);
		Complex t(T(0), T(0));

		if (std::abs(H_at_s) > eps_mach * T(10))
			t = -P_at_s / H_at_s;
		else if (std::abs(dP_at_s) > eps_mach)
			t = P_at_s / dP_at_s;

		Complex new_s = s - t;

		if (std::abs(t) < eps * std::max(T(1), std::abs(s))) {
			alpha = new_s;
			return true;
		}

		H = update_H(H, P, s);
		s = new_s;
	}

	// Резервный путь: классический Ньютон
	for (int l = 0; l < 100; ++l) {
		auto [pz, dpz] = cpoly_detail::horner_eval_and_deriv(P, s);
		if (std::abs(pz) < eps_mach * T(100)) {
			alpha = s;
			return true;
		}
		if (std::abs(dpz) < eps_mach)
			break;
		Complex t = pz / dpz;
		s -= t;
		if (std::abs(t) < eps * std::max(T(1), std::abs(s))) {
			alpha = s;
			return true;
		}
	}

	alpha = s;
	return std::abs(cpoly_detail::horner_eval(P, s)) <
	       std::sqrt(eps) * std::max(T(1), std::abs(s));
}

} // namespace cpoly_impl

/**
 * @brief Поиск всех корней комплексного полинома методом Дженкинса-Трауба.
 *
 * @tparam T Базовый вещественный тип (double, float, long double).
 * @param input_coeffs Комплексные коэффициенты (ascending: input_coeffs[0] = a_0).
 * @param eps Требуемая точность.
 * @return Вектор комплексных корней.
 */
template <typename T>
std::vector<std::complex<T>>
find_roots_by_JenkinsTraub(const std::vector<std::complex<T>> &input_coeffs,
                            T eps = std::numeric_limits<T>::epsilon() * T(1e6)) {
	static_assert(std::is_floating_point_v<T>,
		"Jenkins-Traub требует float, double или long double");

	using Complex = std::complex<T>;
	using CPolynomial = std::vector<Complex>;

	std::vector<Complex> all_roots;

	if (input_coeffs.size() < 2)
		return all_roots;

	CPolynomial poly = input_coeffs;
	while (poly.size() > 1 && std::abs(poly.back()) < eps)
		poly.pop_back();

	if (poly.size() < 2)
		return all_roots;

	// Нормализация (ведущий коэффициент = 1)
	{
		Complex lead = poly.back();
		if (std::abs(lead) < eps)
			return all_roots;
		for (auto &c : poly)
			c /= lead;
	}

	const int degree = int(poly.size()) - 1;
	const T pi = T(3.14159265358979323846L);
	const T eps_mach = std::numeric_limits<T>::epsilon();

	CPolynomial current_poly = poly;
	int remaining = degree;

	while (remaining > 1) {
		int n = int(current_poly.size()) - 1;

		// Радиус Коши для начальных сдвигов
		T cauchy_r = T(0);
		for (int k = 0; k < n; ++k) {
			T ratio = std::pow(std::abs(current_poly[k] / current_poly[n]),
			                   T(1) / T(n - k));
			cauchy_r = std::max(cauchy_r, ratio);
		}
		if (cauchy_r < eps_mach)
			cauchy_r = T(1);

		// Stage 1: No-shift (5 итераций)
		CPolynomial H = cpoly_impl::stage1_no_shift<T>(current_poly, 5);

		// Пробуем несколько сдвигов на окружности радиуса cauchy_r
		bool root_found = false;
		Complex alpha;
		const int num_shifts = std::max(5, n);

		for (int attempt = 0; attempt < num_shifts && !root_found; ++attempt) {
			T theta = T(2) * pi * (T(attempt) + T(0.7)) / T(num_shifts);
			Complex shift = cauchy_r * Complex(std::cos(theta), std::sin(theta));

			root_found =
			    cpoly_impl::stage23_find_root(current_poly, H, shift, alpha, eps);
		}

		if (!root_found) {
			// Резервный путь: Ньютон от нескольких стартовых точек
			for (int attempt = 0; attempt < 20 && !root_found; ++attempt) {
				T theta = T(2) * pi * T(attempt) / T(20);
				T r = cauchy_r * (T(0.5) + T(attempt) * T(0.1));
				Complex start = r * Complex(std::cos(theta), std::sin(theta));

				for (int iter = 0; iter < 200; ++iter) {
					auto [pz, dpz] =
					    cpoly_detail::horner_eval_and_deriv(current_poly, start);
					if (std::abs(pz) < eps_mach * T(100)) {
						alpha = start;
						root_found = true;
						break;
					}
					if (std::abs(dpz) < eps_mach)
						break;
					Complex t = pz / dpz;
					start -= t;
					if (std::abs(t) < eps * std::max(T(1), std::abs(start))) {
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

	// Последний линейный множитель
	if (current_poly.size() == 2 && std::abs(current_poly[1]) > eps_mach) {
		all_roots.push_back(-current_poly[0] / current_poly[1]);
	}

	// Stage финальная: Ньютон-полировка по исходному (ненормализованному) полиному
	for (auto &root : all_roots) {
		for (int iter = 0; iter < 50; ++iter) {
			auto [pz, dpz] = cpoly_detail::horner_eval_and_deriv(poly, root);

			if (std::abs(pz) < eps_mach * T(100))
				break;
			if (std::abs(dpz) < eps_mach)
				break;

			Complex t = pz / dpz;
			root -= t;

			if (std::abs(t) < eps_mach * std::max(T(1), std::abs(root)))
				break;
		}
	}

	return all_roots;
}

/**
 * @brief Перегрузка для полиномов с вещественными коэффициентами.
 */
template <typename T>
std::vector<std::complex<T>>
find_roots_by_JenkinsTraub(const std::vector<T> &real_coeffs,
                            T eps = std::numeric_limits<T>::epsilon() * T(1e6)) {
	static_assert(std::is_floating_point_v<T>,
		"Jenkins-Traub требует float, double или long double");

	return find_roots_by_JenkinsTraub(cpoly_detail::to_complex(real_coeffs), eps);
}
