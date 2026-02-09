// Вспомогательные функции для метода CIsolate
// Разделение для уменьшения размера основного файла
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#include "NumericConcepts.h"
#include "mathUtils.h"
#include <complex>
#include <vector>

template <typename T> using cplx = std::complex<T>;

namespace detail_sagralov {
	/**
	 * @brief Шаблонный факториал для произвольных типов.
	 *
	 * Работает с float, double, long double и arbitrary precision типами.
	 * Использует мемоизацию для эффективности.
	 *
	 * @tparam T Вещественный тип.
	 */
	template <RealNumber T> T factorial_t(int k) {
		// Для arbitrary precision — вычисляем напрямую через T
		T result = T(1);
		for (int i = 2; i <= k; ++i)
			result = result * T(i);
		return result;
	}

	// Специализация для стандартных типов с мемоизацией
	inline long double factorial(int k) {
		static std::vector<long double> f = { 1.0L };
		while (static_cast<int>(f.size()) <= k)
			f.push_back(f.back() * static_cast<long double>(f.size()));
		return (k >= 0 && k < static_cast<int>(f.size())) ? f[k] : 1.0L;
	}

	/**
	 * @brief Вычисление многочлена методом Горнера.
	 */
	template <RealNumber T>
	cplx<T> eval_poly(const std::vector<cplx<T>>& a, const cplx<T>& z) {
		cplx<T> result(0, 0);
		for (const auto& coef : a)
			result = result * z + coef;
		return result;
	}

	/**
	 * @brief Производная многочлена.
	 */
	template <RealNumber T>
	std::vector<cplx<T>> derivative(const std::vector<cplx<T>>& a) {
		int n = static_cast<int>(a.size()) - 1;
		if (n < 1)
			return { cplx<T>(0, 0) };
		std::vector<cplx<T>> d(n);
		for (int i = 0; i < n; ++i)
			d[i] = a[i] * cplx<T>(static_cast<T>(n - i), T(0));
		return d;
	}

	/**
	 * @brief k-я производная многочлена.
	 */
	template <RealNumber T>
	std::vector<cplx<T>> derivative_k(const std::vector<cplx<T>>& a, int k) {
		std::vector<cplx<T>> cur = a;
		for (int iter = 0; iter < k; ++iter)
			cur = derivative(cur);
		return cur;
	}

	/**
	 * @brief Построение F_delta для теста Пелле.
	 */
	template <RealNumber T>
	std::vector<cplx<T>> build_Fdelta(const std::vector<cplx<T>>& a,
		const cplx<T>& m, T r) {
		int n = static_cast<int>(a.size()) - 1;
		if (n < 0)
			return {};

		std::vector<cplx<T>> c(n + 1);
		std::vector<cplx<T>> curr_poly = a;

		T r_pow = T(1);
		for (int i = 0; i <= n; ++i) {
			cplx<T> val = eval_poly(curr_poly, m);

			// Используем шаблонный факториал для совместимости с arbitrary precision
			T fac = factorial_t<T>(i);

			c[i] = val * cplx<T>(r_pow / fac);

			if (i < n) {
				curr_poly = derivative(curr_poly);
				r_pow *= r;
			}
		}
		return c;
	}

	/**
	 * @brief Итерация Греффе: F^[1](x) = (-1)^n (F_e^2 - x*F_o^2).
	 *
	 * Корни нового полинома — квадраты корней исходного.
	 */
	template <RealNumber T>
	std::vector<cplx<T>> graeffe_iteration(const std::vector<cplx<T>>& a) {
		int n = static_cast<int>(a.size()) - 1;
		if (n == 0)
			return a;

		// Разделяем на чётную и нечётную части
		std::vector<cplx<T>> Fe(n + 1, cplx<T>(0, 0));
		std::vector<cplx<T>> Fo(n + 1, cplx<T>(0, 0));

		for (int i = 0; i <= n; ++i) {
			int power = n - i;
			Fe[i] = a[i] * cplx<T>(T((power % 2) == 0));
			Fo[i] = a[i] * cplx<T>(T((power % 2) == 1));
		}

		// Возведение в квадрат
		auto sq = [](const std::vector<cplx<T>>& p) -> std::vector<cplx<T>> {
			int d = static_cast<int>(p.size()) - 1;
			std::vector<cplx<T>> r(2 * d + 1, cplx<T>(0, 0));
			for (int i = 0; i <= d; ++i)
				for (int j = 0; j <= d; ++j)
					r[i + j] += p[i] * p[j];
			return r;
			};

		auto Fe2 = sq(Fe);
		auto Fo2 = sq(Fo);

		// Fe2 - x*Fo2
		std::vector<cplx<T>> full(2 * n + 2, cplx<T>(0, 0));
		for (size_t i = 0; i < Fe2.size(); ++i)
			full[i] += Fe2[i];
		for (size_t i = 0; i < Fo2.size(); ++i)
			full[i + 1] -= Fo2[i];

		// Множитель (-1)^n
		if (n % 2 == 1)
			for (auto& c : full)
				c = -c;

		// Сжатие: берём только чётные степени
		std::vector<cplx<T>> result(n + 1);
		for (int j = 0; j <= n; ++j)
			result[j] = full[2 * j];

		return result;
	}

	/**
	 * @brief Нормализация старшего коэффициента: 1/4 < |a_n| <= 1.
	 */
	template <RealNumber T>
	std::vector<cplx<T>> normalize_polynomial(const std::vector<cplx<T>>& a) {
		std::vector<cplx<T>> result = a;
		if (result.empty())
			return result;

		T leading = abs_val(result[0]);
		if (is_zero_val(leading))
			return result;

		// Нормализация степенями двойки (устойчиво)
		T scale = T(1);
		T normalized = leading;

		while (normalized > T(1)) {
			normalized /= T(2);
			scale /= T(2);
		}
		while (normalized <= T(0.25)) {
			normalized *= T(2);
			scale *= T(2);
		}

		if (scale != T(1)) {
			for (auto& c : result)
				c *= scale;
		}
		return result;
	}

} // namespace detail_sagralov