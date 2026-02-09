// Заголовочный файл для универсальных методов, связанных с полиномами
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <complex>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "NumericConcepts.h"
#include "complexprecision.h"
#include "fprecision.h"
#include "iprecision.h"

// ОПЕРАЦИИ С ПОЛИНОМАМИ
/**
 * @brief Вычисляет производную многочлена.
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 * @param poly Коэффициенты многочлена (от младших к старшим степеням).
 * @return Вектор коэффициентов производной.
 */
template <RealNumber T> std::vector<T> derivative(const std::vector<T> &poly) {
  std::vector<T> result;
  result.reserve(poly.size() > 1 ? poly.size() - 1 : 0);
  for (size_t i = 1; i < poly.size(); ++i)
    result.emplace_back(poly[i] * static_cast<T>(i));
  return result;
}

/**
 * @brief Деление многочлена на (x - r) — дефляция.
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 * @param poly Коэффициенты многочлена.
 * @param root Корень для дефляции.
 * @return Результирующий многочлен пониженной степени.
 * @throws std::invalid_argument если степень многочлена < 1.
 */
template <RealNumber T>
std::vector<T> deflate_poly(const std::vector<T> &poly, const T &root) {
  if (poly.size() < 2)
    throw std::invalid_argument("Polynomial degree must be at least 1");

  std::vector<T> result(poly.size() - 1);
  T rem = poly.back();
  result.back() = rem;

  for (int i = static_cast<int>(poly.size()) - 2; i > 0; --i) {
    rem = poly[i] + rem * root;
    result[i - 1] = rem;
  }
  return result;
}

/**
 * @brief Вычисляет значение многочлена в точке методом Горнера.
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 * @param coeffs Коэффициенты многочлена (от младших к старшим степеням).
 * @param x Точка вычисления.
 * @return Значение P(x).
 */
template <RealNumber T> T eval_poly(const std::vector<T> &coeffs, T x) {
  T res = T(0);
  for (int i = static_cast<int>(coeffs.size()) - 1; i >= 0; --i)
    res = res * x + coeffs[i];
  return res;
}

/**
 * @brief Вычисляет значение многочлена в комплексной точке.
 * @tparam T Базовый вещественный тип.
 */
template <typename T>
  requires std::is_arithmetic_v<T>
std::complex<T> eval_poly(const std::vector<T> &coeffs, std::complex<T> x) {
  std::complex<T> res(T(0));
  for (int i = static_cast<int>(coeffs.size()) - 1; i >= 0; --i)
    res = res * x + std::complex<T>(coeffs[i], T(0));
  return res;
}

inline float_precision eval_poly(const std::vector<float_precision> &coeffs,
                                 float_precision x) {
  float_precision res(0);
  for (int i = static_cast<int>(coeffs.size()) - 1; i >= 0; --i) {
    res = res * x + coeffs[i]; // обычные операторы * и +
  }
  return res;
}

/**
 * @brief Вычисляет значение вещественного полинома (float_precision) в
 * комплексной точке.
 *
 * Необходимо для методов типа Graeffe, которые работают с вещественными
 * коэффициентами, но вычисляют значения в комплексных точках.
 */
inline std::complex<float_precision>
eval_poly(const std::vector<float_precision> &coeffs,
          std::complex<float_precision> x) {
  std::complex<float_precision> res(float_precision(0), float_precision(0));
  for (int i = static_cast<int>(coeffs.size()) - 1; i >= 0; --i) {
    res =
        res * x + std::complex<float_precision>(coeffs[i], float_precision(0));
  }
  return res;
}

template <typename T>
inline complex_precision<T>
eval_poly(const std::vector<complex_precision<T>> &coeffs,
          complex_precision<T> x) {
  complex_precision<T> res(T(0), T(0));
  for (int i = static_cast<int>(coeffs.size()) - 1; i >= 0; --i)
    res = res * x + coeffs[i];
  return res;
}

/**
 * @brief Нормализация полинома (приведение старшего коэффициента к 1).
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 * @throws std::domain_error если старший коэффициент равен нулю.
 */
template <RealNumber T>
std::vector<T> normalizePolynomial(const std::vector<T> &poly) {
  if (poly.empty())
    return {};
  T leading = poly.back();
  if (leading == T(0))
    throw std::domain_error(
        "Cannot normalize polynomial with leading zero coefficient");

  std::vector<T> result;
  result.reserve(poly.size());
  for (const auto &coeff : poly)
    result.emplace_back(coeff / leading);
  return result;
}

// Функция проверки, что все коэффициенты полинома ненулевые
inline bool checkNonZeroCoefficients(const std::vector<int_precision> &coeffs) {
  for (const auto &coeff : coeffs) {
    if (coeff == int_precision(0))
      return false;
  }
  return true;
}

// Умножение двух полиномов (ascending порядок)
template <typename T>
std::vector<T> multiply_polynomials(const std::vector<T> &p1,
                                    const std::vector<T> &p2) {
  if (p1.empty() || p2.empty())
    return {};
  std::vector<T> result(p1.size() + p2.size() - 1, T(0));
  for (size_t i = 0; i < p1.size(); ++i)
    for (size_t j = 0; j < p2.size(); ++j)
      result[i + j] = result[i + j] + p1[i] * p2[j];
  return result;
}

// Сложение двух полиномов (ascending порядок)
template <typename T>
std::vector<T> add_polynomials(const std::vector<T> &p1,
                               const std::vector<T> &p2) {
  size_t max_size = std::max(p1.size(), p2.size());
  std::vector<T> result(max_size, T(0));
  for (size_t i = 0; i < p1.size(); ++i)
    result[i] = result[i] + p1[i];
  for (size_t i = 0; i < p2.size(); ++i)
    result[i] = result[i] + p2[i];
  return result;
}
