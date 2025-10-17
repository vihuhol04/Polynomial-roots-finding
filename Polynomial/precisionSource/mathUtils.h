#pragma once

#include <vector>
#include <limits>
#include <cmath>
#include <complex>
#include <type_traits>
#include <algorithm>

#include "fprecision.h"
#include "complexprecision.h"
#include "fractionprecision.h"
#include "intervalprecision.h"
#include "iprecision.h"
#include "mathprecision.h"
#include "polyprecision.h"

// -----------------------------
// Вспомогательные метапрограммы
// -----------------------------

template<typename U, typename = void>
struct has_to_double : std::false_type {
};

template<typename U>
struct has_to_double<U, std::void_t<decltype(std::declval<U>().to_double())>> : std::true_type {
};

template<typename U>
inline double to_double_fallback(const U&x) {
    if constexpr (has_to_double<U>::value) {
        return x.to_double();
    }
    else {
        return static_cast<double>(x);
    }
}

template<typename U, typename = void>
struct has_member_abs : std::false_type {
};

template<typename U>
struct has_member_abs<U, std::void_t<decltype(std::declval<U>().abs())>> : std::true_type {
};

template<typename T>
struct is_std_complex : std::false_type {
};

template<typename T>
struct is_std_complex<std::complex<T>> : std::true_type {
};

template<typename T>
struct is_complex_precision : std::false_type {
};

template<typename T>
struct is_complex_precision<complex_precision<T>> : std::true_type {
};

template<typename T>
struct is_complex_type : std::bool_constant<is_std_complex<T>::value || is_complex_precision<T>::value> {
};

template<typename T>
struct underlying_real_type {
    using type = T;
};

template<typename T>
struct underlying_real_type<std::complex<T>> {
    using type = T;
};

template<typename T>
struct underlying_real_type<complex_precision<T>> {
    using type = T;
};

// -----------------------------
// ABS / MATH utilities
// -----------------------------

// Для стандартных арифметических типов
template<typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
abs_val(const T&x) {
    return (x >= T(0)) ? x : -x;
}

// Для std::complex
template<typename T>
inline T abs_val(const std::complex<T>&x) {
    return std::abs(x);
}

// Для float_precision
inline float_precision abs_val(const float_precision&x) {
    return float_precision(std::fabs(to_double_fallback(x)));
}

// Для complex_precision
template<typename T>
inline T abs_val(const complex_precision<T>&z) {
    if constexpr (has_member_abs<complex_precision<T>>::value) {
        return z.abs();
    }
    else {
        double re = to_double_fallback(z.real());
        double im = to_double_fallback(z.imag());
        return static_cast<T>(std::sqrt(re * re + im * im));
    }
}

// -----------------------------
// Экспоненциальные и логарифмические функции
// -----------------------------

// exp_val - для стандартных типов
template<typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
exp_val(const T&x) {
    return std::exp(x);
}

// exp_val - для пользовательских типов
template<typename T>
inline typename std::enable_if_t<!std::is_arithmetic_v<T>, T>
exp_val(const T&x) {
    return static_cast<T>(std::exp(to_double_fallback(x)));
}

// exp_val - для std::complex
template<typename T>
inline std::complex<T> exp_val(const std::complex<T>&z) {
    return std::exp(z);
}

// exp_val - для complex_precision
template<typename T>
inline complex_precision<T> exp_val(const complex_precision<T>&z) {
    double re = to_double_fallback(z.real());
    double im = to_double_fallback(z.imag());
    std::complex<double> c = std::exp(std::complex<double>(re, im));
    return complex_precision<T>(static_cast<T>(c.real()), static_cast<T>(c.imag()));
}

// log_val - для стандартных типов
template<typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
log_val(const T&x) {
    if (x <= T(0)) return std::log(std::numeric_limits<T>::min());
    return std::log(x);
}

// log_val - для пользовательских типов
template<typename T>
inline typename std::enable_if_t<!std::is_arithmetic_v<T>, T>
log_val(const T&x) {
    double d = to_double_fallback(x);
    if (d <= 0.0) d = std::numeric_limits<double>::min();
    return static_cast<T>(std::log(d));
}

// log_val - для std::complex
template<typename T>
inline std::complex<T> log_val(const std::complex<T>&z) {
    return std::log(z);
}

// log_val - для complex_precision
template<typename T>
inline complex_precision<T> log_val(const complex_precision<T>&z) {
    double re = to_double_fallback(z.real());
    double im = to_double_fallback(z.imag());
    std::complex<double> c = std::log(std::complex<double>(re, im));
    return complex_precision<T>(static_cast<T>(c.real()), static_cast<T>(c.imag()));
}

// pow_val - для стандартных типов (одинаковые типы)
template<typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
pow_val(const T&base, const T&exponent) {
    return std::pow(base, exponent);
}

// pow_val - для стандартных типов (разные типы)
template<typename T, typename U>
inline typename std::enable_if_t<std::is_arithmetic_v<T> && std::is_arithmetic_v<U>, T>
pow_val(const T&base, const U&exponent) {
    return static_cast<T>(std::pow(base, static_cast<T>(exponent)));
}

// pow_val - для пользовательских типов
template<typename T, typename U>
inline typename std::enable_if_t<!std::is_arithmetic_v<T> || !std::is_arithmetic_v<U>, T>
pow_val(const T&base, const U&exponent) {
    return static_cast<T>(std::pow(to_double_fallback(base), to_double_fallback(exponent)));
}

// sqrt_val - для стандартных типов
template<typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
sqrt_val(const T&x) {
    return std::sqrt(x);
}

// sqrt_val - для пользовательских типов
template<typename T>
inline typename std::enable_if_t<!std::is_arithmetic_v<T>, T>
sqrt_val(const T&x) {
    return static_cast<T>(std::sqrt(to_double_fallback(x)));
}

// sqrt_val - для std::complex
template<typename T>
inline std::complex<T> sqrt_val(const std::complex<T>&z) {
    return std::sqrt(z);
}

// sqrt_val - для complex_precision
template<typename T>
inline complex_precision<T> sqrt_val(const complex_precision<T>&z) {
    double re = to_double_fallback(z.real());
    double im = to_double_fallback(z.imag());
    std::complex<double> c = std::sqrt(std::complex<double>(re, im));
    return complex_precision<T>(static_cast<T>(c.real()), static_cast<T>(c.imag()));
}

// log10_val - для стандартных типов
template<typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
log10_val(const T&x) {
    return std::log10(x);
}

// log10_val - для пользовательских типов
template<typename T>
inline typename std::enable_if_t<!std::is_arithmetic_v<T>, T>
log10_val(const T&x) {
    return log_val(x) / log_val(T(10));
}

// log10_val - для std::complex
template<typename T>
inline std::complex<T> log10_val(const std::complex<T>&z) {
    return std::log10(z);
}

// -----------------------------
// Проверки
// -----------------------------

// isfinite_val
inline bool isfinite_val(const float_precision&x) {
    return std::isfinite(to_double_fallback(x));
}

template<typename T>
inline bool isfinite_val(const complex_precision<T>&z) {
    return std::isfinite(to_double_fallback(z.real())) &&
           std::isfinite(to_double_fallback(z.imag()));
}

template<typename T>
inline bool isfinite_val(const T&x) {
    return std::isfinite(to_double_fallback(x));
}

// isnan_val
inline bool isnan_val(const float_precision&x) {
    return std::isnan(to_double_fallback(x));
}

template<typename T>
inline bool isnan_val(const complex_precision<T>&z) {
    return std::isnan(to_double_fallback(z.real())) ||
           std::isnan(to_double_fallback(z.imag()));
}

template<typename T>
inline bool isnan_val(const T&x) {
    return std::isnan(to_double_fallback(x));
}

// is_zero_val
template<typename T>
inline bool is_zero_val(const T&x, double eps = 1e-12) {
    return to_double_fallback(abs_val(x)) < eps;
}

template<typename T>
inline bool is_zero_val(const std::complex<T>&x, T eps = T(1e-12)) {
    return std::abs(x) < eps;
}

template<typename T>
inline bool is_zero_val(const complex_precision<T>&x, double eps = 1e-20) {
    return to_double_fallback(abs_val(x)) < eps;
}

// is_one_val
template<typename T>
inline bool is_one_val(const T&x, double eps = 1e-12) {
    return is_zero_val(x - T(1), eps);
}

// is_positive_val
template<typename T>
inline bool is_positive_val(const T&x) {
    return to_double_fallback(x) >= 0.0;
}

template<typename T>
inline bool is_positive_val(const std::complex<T>&x) {
    return std::real(x) >= T(0);
}

template<typename T>
inline bool is_positive_val(const complex_precision<T>&x) {
    return to_double_fallback(x.real()) >= 0.0;
}

// max_val
template<typename T>
inline T max_val(const T&a, const T&b) {
    return (a > b) ? a : b;
}

// eps_for_degree
template<typename T>
inline T eps_for_degree(unsigned P) {
    using real_t = typename underlying_real_type<T>::type;
    double eps_small = 1e-8;
    double temp_1 = 1000.0 * std::numeric_limits<real_t>::epsilon();
    double eps_machine = temp_1 * (P + 1);
    return static_cast<T>(std::max(eps_small, eps_machine));
}

// from_long_double
template<typename T>
T from_long_double(long double x) {
    return T(static_cast<double>(x));
}
