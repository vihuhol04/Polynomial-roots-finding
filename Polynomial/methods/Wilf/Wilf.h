// Адаптер метода Вилфа для фреймворка Polynomial-roots-finding
// Метод Вилфа: локализация корней комплексного полинома через
// рекурсивное разбиение прямоугольников и последовательности Штурма.
// Оригинальная реализация: wilf_root_finder.hpp (C++17, header-only)
// Интеграция: автоматическая генерация

// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#include <complex>
#include <type_traits>
#include <vector>

#include "wilf_root_finder.hpp"

/**
 * @brief Поиск всех корней полинома методом Вилфа.
 *
 * Метод локализует корни в комплексной плоскости путём рекурсивного
 * разбиения прямоугольников с подсчётом числа корней внутри области
 * через последовательности Штурма (принцип аргумента).
 *
 * @tparam T Вещественный тип (float, double или long double).
 * @param coeffs Коэффициенты полинома (по возрастанию степеней, coeffs[0] = свободный член).
 * @param eps Требуемая точность локализации.
 * @return Вектор комплексных корней.
 *
 * @note Метод поддерживает только стандартные типы с плавающей точкой.
 */
template <typename T>
std::vector<std::complex<T>>
find_roots_by_Wilf(const std::vector<T> &coeffs, T eps) {
	static_assert(std::is_floating_point_v<T>,
		"find_roots_by_Wilf поддерживает только float, double или long double. "
		"Для повышенной точности используйте find_roots_by_Bairstow или find_roots_by_JenkinsTraub.");

	if (coeffs.size() < 2)
		return {};

	wilf::Options<T> opts;
	opts.coefficient_order = wilf::CoefficientOrder::ascending;

	wilf::Solver<T> solver;
	return solver.find_roots(coeffs, eps, opts);
}

/**
 * @brief Локализация корней полинома методом Вилфа (с кратностями).
 *
 * Возвращает прямоугольные области, содержащие корни, вместе с их
 * кратностями и оценкой погрешности.
 *
 * @tparam T Вещественный тип (float, double или long double).
 * @param coeffs Коэффициенты полинома (по возрастанию степеней).
 * @param eps Требуемая точность локализации.
 * @return Вектор RootBox (центр, размеры, кратность, погрешность).
 */
template <typename T>
std::vector<wilf::RootBox<T>>
localize_roots_by_Wilf(const std::vector<T> &coeffs, T eps) {
	static_assert(std::is_floating_point_v<T>,
		"localize_roots_by_Wilf поддерживает только float, double или long double. "
		"Для повышенной точности используйте find_roots_by_Bairstow или find_roots_by_JenkinsTraub.");

	if (coeffs.size() < 2)
		return {};

	wilf::Options<T> opts;
	opts.coefficient_order = wilf::CoefficientOrder::ascending;

	wilf::Solver<T> solver;
	return solver.localize(coeffs, eps, opts);
}
