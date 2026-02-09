// Метод цепных дробей (теория и реализация)
// Автор основной логики: Шаймарданов Санджар, КМБО-07-23
// Интеграция в фреймворк + исправление ошибок + настройка шаблонности для
// пользовательских типов: Павлова Анастасия, КМБО-01-22

#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <limits>
#include <vector>

#include "Con_Frac_Helper.h"
#include "Helper_for_all_methods.h"
#include "NumericConstants.h"

template <typename T>
std::vector<T> find_roots_by_Continued_Fractions(
    const std::vector<T> &input_coeffs,
    T eps = numeric_constants::adaptive_epsilon<T>(
        numeric_constants::EPSILON_SCALE_PRECISE)) {
  std::vector<T> coeffs = input_coeffs;
  std::vector<T> all_roots;

  // Удаление рациональных корней
  std::vector<T> rational_roots;
  coeffs = remove_rational_roots(coeffs, rational_roots, eps);
  all_roots.insert(all_roots.end(), rational_roots.begin(),
                   rational_roots.end());

  while (coeffs.size() > 2) {
    std::vector<T> current_coeffs = coeffs;
    std::vector<T> deriv = derivative(coeffs); // Символьное дифференцирование
    std::vector<T> partials;
    int iter = 0;
    const int max_iter = numeric_constants::DEFAULT_MAX_ITERATIONS;
    bool root_found = false;
    T best_approx = T(0);
    T min_error = std::numeric_limits<T>::max();

    while (iter < max_iter) {
      int sc = sign_changes(current_coeffs);
      if (sc == 0)
        break;

      T a = compute_partial_quotient(current_coeffs);
      partials.push_back(a);
      current_coeffs = transform(current_coeffs, -a);
      ++iter;

      // Расчет подходящей дроби
      T approx = T(0);
      build_convergent(partials, approx, eps);

      // Уточнение корня
      approx = newton_method(coeffs, deriv, approx, eps);

      T val = abs(eval_poly(coeffs, approx));
      if (val <
          numeric_constants::scale_epsilon(
              eps,
              numeric_constants::RELAXED_TOLERANCE)) { // Критерий сходимости
        all_roots.push_back(approx);
        coeffs = deflate_poly(coeffs, approx);
        root_found = true;
        break; // Перезапуск, так как полином изменился
      }

      if (val < min_error) {
        min_error = val;
        best_approx = approx;
      }
    }

    if (!root_found &&
        min_error < sqrt_val(eps)) { // Адаптивный порог вместо 1e-4
      best_approx = newton_method(coeffs, deriv, best_approx, eps);
      all_roots.push_back(best_approx);
      coeffs = deflate_poly(coeffs, best_approx);
      root_found = true;
    }

    if (!root_found)
      break;
  }

  // Обработка оставшегося линейного множителя
  if (coeffs.size() == 2 && abs(coeffs[1]) > eps) {
    T root = -coeffs[0] / coeffs[1];
    root = newton_method(input_coeffs, derivative(input_coeffs), root, eps);
    if (abs(eval_poly(input_coeffs, root)) <
        numeric_constants::scale_epsilon(eps,
                                         numeric_constants::RELAXED_TOLERANCE))
      all_roots.push_back(root);
  }

  // Удаление дубликатов
  std::sort(all_roots.begin(), all_roots.end());
  auto last = std::unique(all_roots.begin(), all_roots.end(), [eps](T a, T b) {
    return abs(a - b) <
           numeric_constants::scale_epsilon(
               eps, numeric_constants::STRICT_TOLERANCE); // Нечеткое сравнение
  });
  all_roots.erase(last, all_roots.end());

  return all_roots;
}