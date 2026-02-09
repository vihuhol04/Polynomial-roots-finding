// Заголовочный файл с реализацией дополнительных методов для основных методов
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "NumericConcepts.h"
#include "NumericConstants.h"
#include "mathUtils.h"
#include "polynomialUtils.h"

// МЕТОДЫ УТОЧНЕНИЯ КОРНЕЙ
/**
 * @brief Метод Ньютона-Рафсона для уточнения вещественных корней.
 *
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 * @param poly Коэффициенты многочлена.
 * @param deriv Коэффициенты производной.
 * @param guess Начальное приближение.
 * @param epsilon Требуемая точность.
 * @param max_iter Максимальное число итераций (по умолчанию 100).
 * @return Уточнённый корень.
 *
 * @note Использует контроль: проверяет isfinite после деления, а не пытается предсказать деление на ноль заранее.
 */
template <RealNumber T>
T newton_method(const std::vector<T>& poly, const std::vector<T>& deriv, T guess, T epsilon,
	size_t max_iter = numeric_constants::DEFAULT_MAX_ITERATIONS) {
	T f, f_prime, next, diff;
	for (size_t i = 0; i < max_iter; ++i) {
		f = eval_poly(poly, guess);
		f_prime = eval_poly(deriv, guess);

		next = guess - f / f_prime;

		// проверяем результат, а не делитель
		if (!isfinite_val(next))
			break;

		diff = abs_val(next - guess);
		if (diff < epsilon)
			return next;
		guess = next;
	}
	return guess;
}

/**
 * @brief Метод Ньютона для уточнения комплексных корней.
 *
 * @tparam Complex Комплексный тип (std::complex<T> или complex_precision<T>).
 * @tparam Real Вещественный тип для epsilon.
 *
 * Добавила: Семенидо Алина, КМБО-06-22 lina.semenido@yandex.ru
 */
template <typename Complex, RealNumber Real>
void newton_method(std::vector<Complex>& z, const std::vector<Complex>& c,
	Real eps) {
	for (Complex& root : z) {
		auto [f, df] = eval_poly_and_deriv(c, root);
		Complex next = root - f / df;
		// контроль результата
		if (isfinite_val(std::abs(next)))
			root = next;
	}
}

// АНАЛИЗ ЗНАКОВ И ГРАНИЦЫ КОРНЕЙ
/**
 * @brief Подсчёт количества перемен знаков в последовательности коэффициентов.
 *
 * Используется в правиле знаков Декарта.
 *
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 */
template <RealNumber T> int sign_changes(const std::vector<T>& coeffs) {
	if (coeffs.empty())
		return 0;

	int changes = 0;
	T prev = coeffs[0];

	for (size_t i = 1; i < coeffs.size(); ++i) {
		bool is_nonzero = (coeffs[i] != T(0));
		bool sign_changed = (prev > T(0) && coeffs[i] < T(0)) || (prev < T(0) && coeffs[i] > T(0));
		changes += is_nonzero * sign_changed;
		prev = is_nonzero ? coeffs[i] : prev;
	}
	return changes;
}

/**
 * @brief Вычисление верхней границы корней полинома по формуле Лагранжа.
 *
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 */
template <RealNumber T> T computeRootUpperBound(const std::vector<T>& poly) {
	if (poly.empty())
		return T(0);

	T max_neg = T(0);
	int n = static_cast<int>(poly.size()) - 1;

	for (int i = 0; i < n; ++i) {
		// Устранён if: используем max_val с условным выражением
		T ratio = (poly[i] < T(0)) ? abs_val(poly[i]) / abs_val(poly.back()) : T(0);
		max_neg = max_val(max_neg, ratio);
	}
	return T(1) + pow_val(max_neg, T(1) / T(n));
}

/**
 * @brief Граница Коши для модулей корней многочлена.
 *
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 */
template <RealNumber T> T computeRootBoundCauchy(const std::vector<T>& poly) {
	if (poly.empty())
		return T(0);

	// Игнорируем ведущие нули для правильной степени
	size_t actual_degree = poly.size() - 1;
	while (actual_degree > 0 && is_zero_val(poly[actual_degree]))
		--actual_degree;

	if (actual_degree == 0)
		return T(0);

	T max_coeff_abs = T(0);
	for (size_t i = 0; i < actual_degree; ++i) {
		max_coeff_abs = max_val(max_coeff_abs, abs_val(poly[i]));
	}

	T leading_coeff_abs = abs_val(poly[actual_degree]);
	if (is_zero_val(leading_coeff_abs)) {
		throw std::runtime_error("Leading coefficient is zero in computeRootBoundCauchy.");
	}
	return T(1) + max_coeff_abs / leading_coeff_abs;
}