// precisionSource/NumericConcepts.h
// C++20 Concepts для обеспечения корректности шаблонных параметров
// Применяются ко всем численным методам проекта
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#include <complex>
#include <concepts>
#include <cstddef>
#include <limits>

// БАЗОВЫЕ КОНЦЕПТЫ ДЛЯ ВЕЩЕСТВЕННЫХ ТИПОВ
/**
 * @brief Концепт для типов, поддерживающих базовые арифметические операции.
 *
 * Требует наличия операторов +, -, *, / и сравнений.
 * Подходит для float, double, long double, float_precision и др.
 */
template <typename T>
concept Arithmetic = requires(T a, T b) {
	{ a + b } -> std::convertible_to<T>;
	{ a - b } -> std::convertible_to<T>;
	{ a* b } -> std::convertible_to<T>;
	{ a / b } -> std::convertible_to<T>;
	{ -a } -> std::convertible_to<T>;
	{ a < b } -> std::convertible_to<bool>;
	{ a > b } -> std::convertible_to<bool>;
	{ a <= b } -> std::convertible_to<bool>;
	{ a >= b } -> std::convertible_to<bool>;
	{ a == b } -> std::convertible_to<bool>;
};

/**
 * @brief Концепт для вещественных числовых типов.
 *
 * Требует арифметику + возможность конструирования из целых чисел.
 * Охватывает: float, double, long double, float_precision, и др.
 */
template <typename T>
concept RealNumber = Arithmetic<T> && requires(T a) {
	// Конструктор из целого числа (для констант типа T(0), T(1))
	{ T(0) } -> std::convertible_to<T>;
	{ T(1) } -> std::convertible_to<T>;
	// Присваивание
	{ a = T(0) };
};

/**
 * @brief Концепт для типов с поддержкой std::numeric_limits.
 *
 * Необходим для кода, работающего с epsilon, infinity и т.д.
 */
template <typename T>
concept HasNumericLimits = requires {
	{ std::numeric_limits<T>::epsilon() } -> std::convertible_to<T>;
	{ std::numeric_limits<T>::infinity() } -> std::convertible_to<T>;
	{ std::numeric_limits<T>::min() } -> std::convertible_to<T>;
	{ std::numeric_limits<T>::max() } -> std::convertible_to<T>;
};

/**
 * @brief Концепт для типов, совместимых с arbitrary precision.
 *
 * Сочетает RealNumber + HasNumericLimits.
 * Это основной концепт для большинства алгоритмов.
 */
template <typename T>
concept PrecisionCompatible = RealNumber<T> && HasNumericLimits<T>;

// КОНЦЕПТЫ ДЛЯ КОМПЛЕКСНЫХ ТИПОВ
/**
 * @brief Проверяет, является ли тип комплексным числом.
 *
 * Работает с std::complex<T> и complex_precision<T>.
 */
template <typename T>
concept ComplexNumber = requires(T z) {
	{ z.real() };
	{ z.imag() };
	{ std::abs(z) };
}&& requires(T a, T b) {
	{ a + b } -> std::convertible_to<T>;
	{ a - b } -> std::convertible_to<T>;
	{ a * b } -> std::convertible_to<T>;
	{ a / b } -> std::convertible_to<T>;
};

/**
 * @brief Извлекает вещественный тип из комплексного.
 */
template <typename T> struct RealTypeOf {
	using type = T;
};

template <typename T> struct RealTypeOf<std::complex<T>> {
	using type = T;
};

// Для complex_precision — будет специализировано при необходимости
template <typename T> using RealType = typename RealTypeOf<T>::type;

// КОНЦЕПТЫ ДЛЯ ПОЛИНОМИАЛЬНЫХ ОПЕРАЦИЙ
/**
 * @brief Концепт для типа коэффициентов полинома (вещественные).
 */
template <typename T>
concept PolynomialCoefficient = RealNumber<T>;

/**
 * @brief Концепт для типа коэффициентов полинома (комплексные).
 */
template <typename T>
concept ComplexPolynomialCoefficient = ComplexNumber<T>;

/**
 * @brief Концепт для epsilon-параметра (точности).
 *
 * Должен быть положительным вещественным числом.
 */
template <typename T>
concept EpsilonType = RealNumber<T>;

// КОНЦЕПТЫ ДЛЯ КОНТЕЙНЕРОВ
/**
 * @brief Концепт для контейнера коэффициентов полинома.
 */
template <typename Container>
concept CoefficientContainer = requires(Container c) {
	{ c.size() } -> std::convertible_to<std::size_t>;
	{ c.begin() };
	{ c.end() };
	{ c.empty() } -> std::convertible_to<bool>;
	{ c[0] };
};

// ВСПОМОГАТЕЛЬНЫЕ КОНЦЕПТЫ
/**
 * @brief Концепт для типов, поддерживающих математические функции.
 *
 * Проверяет, что для типа определены abs_val, sqrt_val, log_val и т.д.
 */
template <typename T>
concept SupportsMathFunctions = requires(T x) {
	{ abs_val(x) } -> std::convertible_to<T>;
};

// Примечание: этот концепт проверяется отложенно, т.к. abs_val может быть определён позже через перегрузки
