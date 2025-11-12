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

// =============================
//   ВСПОМОГАТЕЛЬНЫЕ ТИПЫ
// =============================

template <typename U, typename = void>
struct has_to_double : std::false_type {};

template <typename U>
struct has_to_double<U, std::void_t<decltype(std::declval<U>().to_double())>> : std::true_type {};

template <typename U, typename = void>
struct has_member_abs : std::false_type {};

template <typename U>
struct has_member_abs<U, std::void_t<decltype(std::declval<U>().abs())>> : std::true_type {};

template <typename U, typename = void>
struct has_member_exp : std::false_type {};

template <typename U>
struct has_member_exp<U, std::void_t<decltype(std::declval<U>().exp())>> : std::true_type {};

template <typename U, typename = void>
struct has_member_log : std::false_type {};

template <typename U>
struct has_member_log<U, std::void_t<decltype(std::declval<U>().log())>> : std::true_type {};

template <typename U, typename = void>
struct has_member_sqrt : std::false_type {};

template <typename U>
struct has_member_sqrt<U, std::void_t<decltype(std::declval<U>().sqrt())>> : std::true_type {};

template <typename T>
struct is_std_complex : std::false_type {};

template <typename T>
struct is_std_complex<std::complex<T>> : std::true_type {};

template <typename T>
struct is_complex_precision : std::false_type {};

template <typename T>
struct is_complex_precision<complex_precision<T>> : std::true_type {};

template <typename T>
struct is_complex_type : std::bool_constant<is_std_complex<T>::value || is_complex_precision<T>::value> {};

template <typename T>
struct underlying_real_type { using type = T; };

template <typename T>
struct underlying_real_type<std::complex<T>> { using type = T; };

template <typename T>
struct underlying_real_type<complex_precision<T>> { using type = T; };

// =============================
//       ABS / MATH UTILS
// =============================

// Для стандартных арифметических типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
abs_val(const T& x) { return std::abs(x); }

// Для std::complex
template <typename T>
inline T abs_val(const std::complex<T>& x) { return std::abs(x); }

// Для float_precision
inline float_precision abs_val(const float_precision& x) {
    return abs(x);
}

// Для complex_precision
template <typename T>
inline T abs_val(const complex_precision<T>& z) {
    return abs(z);
}

// =============================
//     ЛОГАРИФМЫ И ЭКСПОНЕНТЫ
// =============================

// Для стандартных арифметических типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
exp_val(const T& x) { return std::exp(x); }

// Для float_precision
inline float_precision exp_val(const float_precision& x) {
    return exp(x); 
}

// Для complex_precision
template <typename T>
inline complex_precision<T> exp_val(const complex_precision<T>& z) {
    return exp(z); 
}

// Для std::complex
template <typename T>
inline std::complex<T> exp_val(const std::complex<T>& z) { return std::exp(z); }

// expm1_val — безопасная версия exp(x)-1
template <typename T>
inline T expm1_val(const T& x) {
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
        return -std::numeric_limits<float_precision>::infinity();
    return log(x);
}

// log_val для complex_precision
template <typename T>
inline complex_precision<T> log_val(const complex_precision<T>& z) {
    return log(z);
}

// log_val для std::complex
template <typename T>
inline std::complex<T> log_val(const std::complex<T>& z) { return std::log(z); }

// pow_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
pow_val(const T& base, const T& exp) { return std::pow(base, exp); }

// pow_val для float_precision
inline float_precision pow_val(const float_precision& base, const float_precision& exp) {
    return pow(base, exp);
}

// pow_val для других типов
template <typename T, typename U>
inline typename std::enable_if_t<!std::is_arithmetic_v<T> && !std::is_same_v<T, float_precision>, T>
pow_val(const T& base, const U& exp) {
    return exp_val(log_val(base) * T(exp));
}

// sqrt_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
sqrt_val(const T& x) { return std::sqrt(x); }

// sqrt_val для float_precision
inline float_precision sqrt_val(const float_precision& x) {
    return sqrt(x);
}

// sqrt_val для complex_precision
template <typename T>
inline complex_precision<T> sqrt_val(const complex_precision<T>& z) {
    return sqrt(z);
}

// sqrt_val для std::complex
template <typename T>
inline std::complex<T> sqrt_val(const std::complex<T>& z) { return std::sqrt(z); }

// log10_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
log10_val(const T& x) { return std::log10(x); }

// log10_val для float_precision
inline float_precision log10_val(const float_precision& x) {
    return log10(x);
}

// log10_val для других типов
template <typename T>
inline typename std::enable_if_t<!std::is_arithmetic_v<T> && !std::is_same_v<T, float_precision>, T>
log10_val(const T& x) { return log_val(x) / log_val(T(10)); }

// =============================
//    ТРИГОНОМЕТРИЧЕСКИЕ ФУНКЦИИ
// =============================

// sin_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
sin_val(const T& x) { return std::sin(x); }

// sin_val для float_precision
inline float_precision sin_val(const float_precision& x) {
    return sin(x);
}


// cos_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
cos_val(const T& x) { return std::cos(x); }

// cos_val для float_precision
inline float_precision cos_val(const float_precision& x) {
    return cos(x);
}

// atan2_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
atan2_val(const T& y, const T& x) { return std::atan2(y, x); }

// atan2_val для float_precision
inline float_precision atan2_val(const float_precision& y, const float_precision& x) {
    return atan2(y, x);
}

// =============================
//          ПРОВЕРКИ
// =============================

// isfinite_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, bool>
isfinite_val(const T& x) { return std::isfinite(x); }

// isfinite_val для float_precision
inline bool isfinite_val(const float_precision& x) {
    return isfinite(x);
}

// isfinite_val для других типов
template <typename T>
inline typename std::enable_if_t<!std::is_arithmetic_v<T> && !std::is_same_v<T, float_precision>, bool>
isfinite_val(const T& x) {
    return abs_val(x) < std::numeric_limits<T>::max();
}

template <typename T>
inline bool is_finite_ld(const T& x) { return isfinite_val(x); }

// is_zero_val для всех типов
template <typename T>
inline bool is_zero_val(const T& x, T eps = T(1e-12)) {
    return abs_val(x) < eps;
}

// is_one_val для всех типов
template <typename T>
inline bool is_one_val(const T& x, T eps = T(1e-12)) {
    return is_zero_val(x - T(1), eps);
}

// is_positive_val для всех типов
template <typename T>
inline bool is_positive_val(const T& x) { return x >= T(0); }

// =============================
//       СРАВНЕНИЯ / КЛАМПЫ
// =============================

// max_val для всех типов
template <typename T>
inline T max_val(const T& a, const T& b) {
    return (a > b) ? a : b;
}

// min_val для всех типов
template <typename T>
inline T min_val(const T& a, const T& b) {
    return (a < b) ? a : b;
}

// clamp_val для всех типов
template <typename T>
inline T clamp_val(const T& v, const T& lo, const T& hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

// =============================
//        LDExp и др.
// =============================

// ldexp_val для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
ldexp_val(const T& x, int exp) { return std::ldexp(x, exp); }

// ldexp_val для float_precision
inline float_precision ldexp_val(const float_precision& x, int exp) {
    return x * pow_val(float_precision(2), int(exp));
}

// ldexp_val для других типов
template <typename T>
inline typename std::enable_if_t<!std::is_arithmetic_v<T> && !std::is_same_v<T, float_precision>, T>
ldexp_val(const T& x, int exp) {
    return x * pow_val(T(2), T(exp));
}

// =============================
//        EPS / HELPERS
// =============================

// eps_for_degree для стандартных типов
template <typename T>
inline typename std::enable_if_t<std::is_arithmetic_v<T>, T>
eps_for_degree(unsigned P) {
    T eps_small = T(1e-8);
    T temp_1 = T(1000) * std::numeric_limits<T>::epsilon();
    T eps_machine = temp_1 * T(P + 1);
    return std::max(eps_small, eps_machine);
}

// eps_for_degree для float_precision
inline float_precision eps_for_degree(float_precision P) {
    float_precision eps_small = float_precision(1e-8);
    float_precision temp_1 = float_precision(1000) * std::numeric_limits<float_precision>::epsilon();
    float_precision eps_machine = temp_1 * (P + float_precision(1));
    return max_val(eps_small, eps_machine);
}

// eps_for_degree для других типов
template <typename T>
inline typename std::enable_if_t<!std::is_arithmetic_v<T> && !std::is_same_v<T, float_precision>, T>
eps_for_degree(unsigned P) {
    using real_t = typename underlying_real_type<T>::type;
    real_t eps_small = real_t(1e-8);
    real_t temp_1 = real_t(1000) * std::numeric_limits<real_t>::epsilon();
    real_t eps_machine = temp_1 * real_t(P + 1);
    return T(std::max(eps_small, eps_machine));
}

template<typename T>
inline T neg_inf_val() {
    if constexpr (std::is_arithmetic_v<T>)
        return -std::numeric_limits<T>::infinity();
    else if constexpr (std::is_same_v<T, float_precision>)
        return -FP_INFINITY; // используем определённый константный FP_INFINITY для float_precision
    else
        return T(-1e300); // surrogate for -inf
}

template<typename T>
inline T pos_inf_val() {
    if constexpr (std::is_arithmetic_v<T>)
        return std::numeric_limits<T>::infinity();
    else if constexpr (std::is_same_v<T, float_precision>)
        return FP_INFINITY; // используем определённый константный FP_INFINITY для float_precision
    else
        return T(1e300); // surrogate for +inf
}