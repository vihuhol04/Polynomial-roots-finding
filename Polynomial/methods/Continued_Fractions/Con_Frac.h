// Метод цепных дробей (теория и реализация)
// Автор основной логики: Шаймарданов Санджар, КМБО-07-23
// Интеграция в фреймворк + исправление ошибок + настройка шаблонности для
// пользовательских типов: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <algorithm>
#include <limits>
#include <vector>

#include "Con_Frac_Helper.h"
#include "Helper_for_all_methods.h"
#include "NumericConstants.h"

/**
 * @brief Поиск вещественных корней полинома методом цепных дробей.
 *
 * Алгоритм на основе теоремы Винсента: последовательное преобразование
 * полинома (сдвиг + обращение) порождает цепную дробь, сходящуюся
 * к положительному вещественному корню.
 *
 * @tparam T Тип коэффициентов, должен удовлетворять RealNumber.
 * @param input_coeffs Коэффициенты по возрастанию степеней.
 * @param eps Требуемая точность.
 * @return Вектор найденных вещественных корней.
 */
template <typename T>
std::vector<T> find_roots_by_Continued_Fractions(
    const std::vector<T> &input_coeffs,
    T eps = numeric_constants::adaptive_epsilon<T>(
        numeric_constants::EPSILON_SCALE_PRECISE)) {
	std::vector<T> coeffs = input_coeffs;
	std::vector<T> all_roots;

	std::vector<T> rational_roots;
	coeffs = remove_rational_roots(coeffs, rational_roots, eps);
	all_roots.insert(all_roots.end(), rational_roots.begin(),
	                 rational_roots.end());

	// Ищем положительные корни текущего полинома
	auto find_positive_root = [&](std::vector<T> &poly) -> bool {
		std::vector<T> current = poly;
		std::vector<T> deriv = derivative(poly);
		std::vector<T> partials;
		int iter = 0;
		const int max_iter = numeric_constants::DEFAULT_MAX_ITERATIONS;
		T best_approx = T(0);
		T min_error = std::numeric_limits<T>::max();

		while (iter < max_iter) {
			int sc = sign_changes(current);
			if (sc == 0)
				break;

			if (sc == 1) {
				// Ровно один положительный корень — уточняем через подходящую дробь
				T a = compute_partial_quotient(current);
				partials.push_back(a);

				T approx = T(0);
				build_convergent(partials, approx, eps);
				approx = newton_method(poly, deriv, approx, eps);

				T val = abs_val(eval_poly(poly, approx));
				if (val < numeric_constants::scale_epsilon(
				               eps, numeric_constants::RELAXED_TOLERANCE)) {
					all_roots.push_back(approx);
					poly = deflate_poly(poly, approx);
					return true;
				}
				if (val < min_error) {
					min_error = val;
					best_approx = approx;
				}
			}

			T a = compute_partial_quotient(current);
			partials.push_back(a);

			// P(x) -> P(x + a)
			current = cf_taylor_shift(current, a);
			// P(x) -> x^n * P(1/x) (реверс коэффициентов)
			current = cf_reciprocal(current);
			++iter;

			T approx = T(0);
			build_convergent(partials, approx, eps);
			approx = newton_method(poly, deriv, approx, eps);

			T val = abs_val(eval_poly(poly, approx));
			if (val < numeric_constants::scale_epsilon(
			               eps, numeric_constants::RELAXED_TOLERANCE)) {
				all_roots.push_back(approx);
				poly = deflate_poly(poly, approx);
				return true;
			}

			if (val < min_error) {
				min_error = val;
				best_approx = approx;
			}
		}

		if (min_error < sqrt_val(eps)) {
			best_approx = newton_method(poly, deriv, best_approx, eps);
			T val = abs_val(eval_poly(poly, best_approx));
			if (val < numeric_constants::scale_epsilon(
			               eps, numeric_constants::RELAXED_TOLERANCE)) {
				all_roots.push_back(best_approx);
				poly = deflate_poly(poly, best_approx);
				return true;
			}
		}
		return false;
	};

	// Ищем положительные корни
	while (coeffs.size() > 2) {
		if (!find_positive_root(coeffs))
			break;
	}

	// Ищем отрицательные корни: подставляем x -> -x и ищем положительные корни
	std::vector<T> neg_coeffs = coeffs;
	for (size_t i = 0; i < neg_coeffs.size(); ++i) {
		if (i % 2 == 1)
			neg_coeffs[i] = -neg_coeffs[i];
	}

	while (neg_coeffs.size() > 2) {
		size_t old_size = all_roots.size();
		if (!find_positive_root(neg_coeffs))
			break;
		// Последний добавленный корень — это положительный корень P(-x),
		// значит настоящий корень отрицателен
		if (all_roots.size() > old_size)
			all_roots.back() = -all_roots.back();
	}

	// Дефлируем исходный полином по найденным отрицательным корням
	for (size_t i = rational_roots.size(); i < all_roots.size(); ++i) {
		// Корни уже были дефлированы из рабочих копий
	}

	// Обработка оставшегося линейного множителя из обоих проходов
	if (coeffs.size() == 2 && abs_val(coeffs[1]) > eps) {
		T root = -coeffs[0] / coeffs[1];
		std::vector<T> orig_deriv = derivative(input_coeffs);
		root = newton_method(input_coeffs, orig_deriv, root, eps);
		if (abs_val(eval_poly(input_coeffs, root)) <
		    numeric_constants::scale_epsilon(eps,
		                                    numeric_constants::RELAXED_TOLERANCE))
			all_roots.push_back(root);
	}

	// Удаление дубликатов
	std::sort(all_roots.begin(), all_roots.end());
	auto last = std::unique(
	    all_roots.begin(), all_roots.end(), [eps](T a, T b) {
		    return abs_val(a - b) <
		           numeric_constants::scale_epsilon(
		               eps, numeric_constants::STRICT_TOLERANCE);
	    });
	all_roots.erase(last, all_roots.end());

	return all_roots;
}
