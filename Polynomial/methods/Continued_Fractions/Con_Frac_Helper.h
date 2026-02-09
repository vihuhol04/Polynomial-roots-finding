// Вспомогательные функции для метода цепных дробей
// Теория + основная реализация: Шаймарданов Санджар, КМБО-07-23
// Интеграция в фреймворк + исправление ошибок + настройка шаблонности для
// пользовательских типов: Павлова Анастасия, КМБО-01-22

#pragma once

#include "NumericConstants.h"
#include "polynomialUtils.h"
#include <cmath>
#include <vector>


// преобразование полинома для цепных дробей: P(x) -> P(a + 1/x)
template <typename T>
std::vector<T> transform(const std::vector<T> &coeffs, T a) {
  std::vector<T> res = coeffs;
  int n = int(res.size());
  for (int i = n - 2; i >= 0; --i) {
    for (int j = i; j >= 0; --j)
      res[j] += a * res[j + 1];
  }
  return res;
}

// вычисление частичного частного по правилу Коши
template <typename T> T compute_partial_quotient(const std::vector<T> &coeffs) {
  if (coeffs.size() < 2)
    return T(1);

  // ищем первый ненулевой коэффициент
  size_t k = coeffs.size() - 1;
  T threshold = numeric_constants::adaptive_epsilon<T>(
      numeric_constants::EPSILON_SCALE_PRECISE);
  while (k > 0 && abs(coeffs[k]) < threshold)
    k--;

  if (k == 0)
    return T(1);

  // берем отношение двух старших коэффициентов
  T a = coeffs[k] / coeffs[k - 1];
  return (a > 0) ? floor(a + 0.5)
                 : ceil(a - 0.5); // округление до ближайшего целого
}

// вычисление приближенного значения корня по цепной дроби
template <typename T>
T build_convergent(const std::vector<T> &a, T &approx, double eps) {
  T Pm2 = T(1), Pm1 = T(a[0]);
  T Qm2 = T(0), Qm1 = T(1);
  T prev_approx = Pm1 / Qm1;

  for (size_t k = 1; k < a.size(); ++k) {
    T P = T(a[k]) * Pm1 + Pm2;
    T Q = T(a[k]) * Qm1 + Qm2;
    approx = P / Q;

    if (abs(approx - prev_approx) < eps)
      break;
    prev_approx = approx;

    Pm2 = Pm1;
    Pm1 = P;
    Qm2 = Qm1;
    Qm1 = Q;
  }

  return approx;
}

// удаление рациональный корней
template <typename T>
std::vector<T>
remove_rational_roots(std::vector<T> &coeffs, std::vector<T> &rational_roots,
                      double eps = numeric_constants::adaptive_epsilon<double>(
                          numeric_constants::EPSILON_SCALE_STANDARD)) {
  std::vector<T> result = coeffs;
  for (int i = -100; i <= 100; ++i) {
    T r = T(i);
    int deflate_count = 0;
    while (abs(eval_poly(result, r)) < eps) {
      rational_roots.push_back(r);
      result = deflate_poly(result, r);
      if (++deflate_count > 20)
        break;
    }
  }
  return result;
}