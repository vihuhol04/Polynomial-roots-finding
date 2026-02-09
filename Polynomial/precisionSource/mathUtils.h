// Файл с математическими функциями для всех используемых типов
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>

#include "NumericConstants.h"
#include <type_traits>
#include <vector>

#include "complexprecision.h"
#include "fprecision.h"

// разрешение неоднозначности MSVC для сравнения float_precision и арифметических типов

// Поскольку fprecision.h (vendor) определяет operator<(float_precision, float_precision) и iprecision.h определяет 
// operator<(float_precision, int_precision), сравнения, подобные float_precision < int, становятся неоднозначными для MSVC std::complex
// Явные перегрузки для арифметических типов решают эту проблему

template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
operator==(const float_precision& a, const T& b) {
	return a == float_precision(b);
}

template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
operator==(const T& a, const float_precision& b) {
	return float_precision(a) == b;
}

template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
operator!=(const float_precision& a, const T& b) {
	return a != float_precision(b);
}

template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
operator!=(const T& a, const float_precision& b) {
	return float_precision(a) != b;
}

template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
operator<(const float_precision& a, const T& b) {
	return a < float_precision(b);
}

template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
operator<(const T& a, const float_precision& b) {
	return float_precision(a) < b;
}

template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
operator>(const float_precision& a, const T& b) {
	return a > float_precision(b);
}

template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
operator>(const T& a, const float_precision& b) {
	return float_precision(a) > b;
}

template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
operator<=(const float_precision& a, const T& b) {
	return a <= float_precision(b);
}

template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
operator<=(const T& a, const float_precision& b) {
	return float_precision(a) <= b;
}

template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
operator>=(const float_precision& a, const T& b) {
	return a >= float_precision(b);
}

template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
operator>=(const T& a, const float_precision& b) {
	return float_precision(a) >= b;
}

// ВСПОМОГАТЕЛЬНЫЕ МЕТАПРОГРАММИРУЮЩИЕ КОНСТРУКЦИИ

// Trait для извлечения базового вещественного типа из комплексного
template <typename T> struct underlying_real_type {
	using type = T;
};

template <typename T> struct underlying_real_type<std::complex<T>> {
	using type = T;
};

template <typename T> struct underlying_real_type<complex_precision<T>> {
	using type = T;
};

// МОДУЛЬ ЧИСЛА
// Для стандартных арифметических типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
abs_val(const T& x) {
	return std::abs(x);
}

// Для std::complex
template <typename T> inline T abs_val(const std::complex<T>& x) {
	return std::abs(x);
}

// Для float_precision
inline float_precision abs_val(const float_precision& x) { return abs(x); }

// Для complex_precision
template <typename T> inline T abs_val(const complex_precision<T>& z) {
	return abs(z);
}

// ЛОГАРИФМЫ И ЭКСПОНЕНТЫ
// Для стандартных арифметических типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
exp_val(const T& x) {
	return std::exp(x);
}

// Для float_precision
inline float_precision exp_val(const float_precision& x) { return exp(x); }

// Для complex_precision
template <typename T>
inline complex_precision<T> exp_val(const complex_precision<T>& z) {
	return exp(z);
}

// Для std::complex
template <typename T> inline std::complex<T> exp_val(const std::complex<T>& z) {
	return std::exp(z);
}

// exp(x)-1
template <typename T> inline T expm1_val(const T& x) {
	if constexpr (std::is_arithmetic_v<T>)
		return std::expm1(x);
	else
		return exp_val(x) - T(1);
}

// log_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
log_val(const T& x) {
	if (x <= T(0))
		return -std::numeric_limits<T>::infinity();
	return std::log(x);
}

// log_val для float_precision
inline float_precision log_val(const float_precision& x) {
	if (x <= float_precision(0))
		return -FP_INFINITY;
	return log(x);
}

// log_val для complex_precision
template <typename T>
inline complex_precision<T> log_val(const complex_precision<T>& z) {
	return log(z);
}

// log_val для std::complex
template <typename T> inline std::complex<T> log_val(const std::complex<T>& z) {
	return std::log(z);
}

// pow_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
pow_val(const T& base, const T& exp) {
	return std::pow(base, exp);
}

// pow_val для float_precision
inline float_precision pow_val(const float_precision& base,
	const float_precision& exp) {
	return pow(base, exp);
}

// pow_val для других типов
template <typename T, typename U>
inline typename std::enable_if_t<
	!std::is_arithmetic_v<T> && !std::is_same_v<T, float_precision>, T>
	pow_val(const T& base, const U& exp) {
	return exp_val(log_val(base) * T(exp));
}

// sqrt_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
sqrt_val(const T& x) {
	return std::sqrt(x);
}

// sqrt_val для float_precision
inline float_precision sqrt_val(const float_precision& x) { return sqrt(x); }

// sqrt_val для complex_precision
template <typename T>
inline complex_precision<T> sqrt_val(const complex_precision<T>& z) {
	return sqrt(z);
}

// sqrt_val для std::complex
template <typename T>
inline std::complex<T> sqrt_val(const std::complex<T>& z) {
	return std::sqrt(z);
}

// log10_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
log10_val(const T& x) {
	return std::log10(x);
}

// log10_val для float_precision
inline float_precision log10_val(const float_precision& x) { return log10(x); }

// log10_val для других типов
template <typename T>
inline typename std::enable_if_t<
	!std::is_arithmetic_v<T> && !std::is_same_v<T, float_precision>, T>
	log10_val(const T& x) {
	return log_val(x) / log_val(T(numeric_constants::FACTOR_10));
}

// ОКРУГЛЕНИЕ
// ceil_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
ceil_val(const T& x) {
	return std::ceil(x);
}

// ceil_val для float_precision
inline float_precision ceil_val(const float_precision& x) { return ceil(x); }

// floor_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
floor_val(const T& x) {
	return std::floor(x);
}

// floor_val для float_precision
inline float_precision floor_val(const float_precision& x) { return floor(x); }

// round_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
round_val(const T& x) {
	return std::round(x);
}

// round_val для float_precision
inline float_precision round_val(const float_precision& x) {
	// float_precision может не иметь round, используем floor(x + 0.5)
	return floor(x + float_precision("0.5"));
}

// ТРИГОНОМЕТРИЧЕСКИЕ ФУНКЦИИ
// sin_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
sin_val(const T& x) {
	return std::sin(x);
}

// sin_val для float_precision
inline float_precision sin_val(const float_precision& x) { return sin(x); }

// cos_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
cos_val(const T& x) {
	return std::cos(x);
}

// cos_val для float_precision
inline float_precision cos_val(const float_precision& x) { return cos(x); }

// atan2_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
atan2_val(const T& y, const T& x) {
	return std::atan2(y, x);
}

// atan2_val для float_precision
inline float_precision atan2_val(const float_precision& y,
	const float_precision& x) {
	return atan2(y, x);
}

// ПРОВЕРКИ
// isfinite_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
isfinite_val(const T& x) {
	return std::isfinite(x);
}

// isfinite_val для float_precision
inline bool isfinite_val(const float_precision& x) { return isfinite(x); }

template <typename T> inline bool is_finite_ld(const T& x) {
	return isfinite_val(x);
}

/**
 * @brief Проверка, является ли значение "нулевым" с учётом epsilon.
 *
 * По умолчанию использует epsilon, масштабированный относительно типа T.
 * Для arbitrary precision типов рекомендуется явно передавать eps.
 *
 * @tparam T Числовой тип.
 * @param x Проверяемое значение.
 * @param eps Порог для сравнения (по умолчанию: 1000 * numeric_limits<T>::epsilon()).
 */
template <typename T> inline bool is_zero_val(const T& x, T eps) {
	return abs_val(x) < eps;
}

// Перегрузка без явного eps: используем масштабированный epsilon типа
// Специализация для арифметических типов (double, float, etc.)
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
is_zero_val(const T& x) {
	// Используем 1000 * epsilon вместо хардкоженного 1e-12
	const T scale = T(numeric_constants::FACTOR_1K);
	T eps = scale * std::numeric_limits<T>::epsilon();
	return abs_val(x) < eps;
}

// Специализация для float_precision
inline bool is_zero_val(const float_precision& x) {
	float_precision scale(numeric_constants::FACTOR_1K);
	float_precision eps = scale * float_precision().epsilon();
	return abs_val(x) < eps;
}

// Специализация для других (non-arithmetic, non-float_precision) типов
template <typename T>
inline typename std::enable_if_t<
	!std::is_arithmetic_v<T> && !std::is_same_v<T, float_precision>, bool>
	is_zero_val(const T& x) {
	// Пробуем использовать numeric_limits, если доступен
	if constexpr (std::numeric_limits<T>::is_specialized) {
		T eps = numeric_constants::scale_epsilon(std::numeric_limits<T>::epsilon(),
			numeric_constants::FACTOR_1K);
		return abs_val(x) < eps;
	}
	else {
		// Последний резерв: масштабируем относительно единицы
		T eps = T(1) / T(numeric_constants::FACTOR_1B); // 1e-9 через целые числа
		return abs_val(x) < eps;
	}
}

/**
 * @brief Проверка, равно ли значение единице с учётом epsilon.
 */
template <typename T> inline bool is_one_val(const T& x, T eps) {
	return is_zero_val(x - T(1), eps);
}

template <typename T> inline bool is_one_val(const T& x) {
	return is_zero_val(x - T(1));
}

// is_positive_val для всех типов
template <typename T> inline bool is_positive_val(const T& x) {
	return x >= T(0);
}

// СРАВНЕНИЯ
// max_val для всех типов
template <typename T> inline T max_val(const T& a, const T& b) {
	return (a > b) ? a : b;
}

// min_val для всех типов
template <typename T> inline T min_val(const T& a, const T& b) {
	return (a < b) ? a : b;
}

// clamp_val для всех типов
template <typename T> inline T clamp_val(const T& v, const T& lo, const T& hi) {
	if (v < lo)
		return lo;
	if (v > hi)
		return hi;
	return v;
}

// LDExp и др.
// ldexp_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
ldexp_val(const T& x, int exp) {
	return std::ldexp(x, exp);
}

// ldexp_val для float_precision
inline float_precision ldexp_val(const float_precision& x, int exp) {
	return x * pow_val(float_precision(2), int(exp));
}

// ldexp_val для других типов
template <typename T>
inline typename std::enable_if_t<
	!std::is_arithmetic_v<T> && !std::is_same_v<T, float_precision>, T>
	ldexp_val(const T& x, int exp) {
	return x * pow_val(T(2), T(exp));
}

// EPS
// eps_for_degree для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
eps_for_degree(unsigned P) {
	// Минимальный порог: sqrt(epsilon) — разумный баланс между точностью и
	// устойчивостью
	T eps_small = std::sqrt(std::numeric_limits<T>::epsilon());
	T temp_1 = numeric_constants::scale_epsilon(std::numeric_limits<T>::epsilon(),
		numeric_constants::FACTOR_1K);
	T eps_machine = temp_1 * T(P + 1);
	return std::max(eps_small, eps_machine);
}

// eps_for_degree для float_precision
inline float_precision eps_for_degree(float_precision P) {
	// Для arbitrary precision: sqrt(epsilon)
	float_precision eps_small = sqrt(float_precision().epsilon());
	float_precision temp_1 = float_precision(numeric_constants::FACTOR_1K) *
		float_precision().epsilon();
	float_precision eps_machine = temp_1 * (P + float_precision(1));
	return max_val(eps_small, eps_machine);
}

// eps_for_degree для других типов
template <typename T>
inline typename std::enable_if_t<
	!std::is_arithmetic_v<T> && !std::is_same_v<T, float_precision>, T>
	eps_for_degree(unsigned P) {
	using real_t = typename underlying_real_type<T>::type;
	// sqrt(epsilon) как минимальный порог
	real_t eps_small = std::sqrt(std::numeric_limits<real_t>::epsilon());
	real_t temp_1 = numeric_constants::scale_epsilon(
		std::numeric_limits<real_t>::epsilon(), numeric_constants::FACTOR_1K);
	real_t eps_machine = temp_1 * real_t(P + 1);
	return T(std::max(eps_small, eps_machine));
}

// ПРИСВАИВАНИЕ БЕСКОНЕЧНОСТИ
// -бесконечность
template <typename T> inline T neg_inf_val() {
	if constexpr (std::is_arithmetic_v<T>)
		return -std::numeric_limits<T>::infinity();
	else if constexpr (std::is_same_v<T, float_precision>)
		return -FP_INFINITY;
	else if constexpr (std::numeric_limits<T>::has_infinity)
		return -std::numeric_limits<T>::infinity();
	else
		return std::numeric_limits<T>::lowest();
}

// +бесконечность
template <typename T> inline T pos_inf_val() {
	if constexpr (std::is_arithmetic_v<T>)
		return std::numeric_limits<T>::infinity();
	else if constexpr (std::is_same_v<T, float_precision>)
		return FP_INFINITY;
	else if constexpr (std::numeric_limits<T>::has_infinity)
		return std::numeric_limits<T>::infinity();
	else
		return std::numeric_limits<T>::max();
}

// ФАКТОРИАЛ ЧИСЛА
template <typename T> T factorial(int n) {
	T result = T(1);
	for (int i = 2; i <= n; ++i)
		result = result * T(i);
	return result;
}

/**
 * @brief Распределяет целое число M пропорционально весам weights.
 *
 * Использует метод наибольших остатков, чтобы сумма распределённых частей была в точности равна M.
 *
 * @param M Общее количество элементов для распределения.
 * @param weights Вектор весов (например, длины интервалов).
 * @return std::vector<unsigned> Вектор целочисленных количеств.
 */
inline std::vector<unsigned>
allocate_counts_proportional(unsigned M, const std::vector<double>& weights) {
	if (weights.empty())
		return {};

	double sum_w = 0.0;
	for (double w : weights)
		sum_w += w;

	std::vector<unsigned> result(weights.size(), 0);
	if (sum_w <= 0.0) {
		// Если сумма весов 0, просто кладём всё в первую корзину (или распределяем равномерно?) 
		// Выберем первую для определённости
		if (M > 0)
			result[0] = M;
		return result;
	}

	// 1. Считаем идеальные доли и их целые части
	std::vector<double> remainders(weights.size());
	unsigned allocated_sum = 0;

	for (size_t i = 0; i < weights.size(); ++i) {
		double ideal = static_cast<double>(M) * (weights[i] / sum_w);
		unsigned floor_val = static_cast<unsigned>(ideal);

		result[i] = floor_val;
		remainders[i] = ideal - floor_val;
		allocated_sum += floor_val;
	}

	// 2. Распределяем остаток (M - allocated_sum) методом наибольших остатков
	unsigned diff = M - allocated_sum;
	if (diff > 0) {
		// Индексы для сортировки по убыванию остатков
		std::vector<size_t> indices(weights.size());
		std::iota(indices.begin(), indices.end(), 0);

		// Сортировка индексов: у кого больше остаток -> тот раньше
		std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
			return remainders[a] > remainders[b];
			});

		// Раздаём по +1
		for (size_t k = 0; k < diff; ++k) {
			result[indices[k]]++;
		}
	}

	return result;
}