// Заголовочный файл с реализацией дополнительных методов для основных методов
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#include "polynomialUtils.h"
#include "mathUtils.h"

#define _USE_MATH_DEFINES
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Метод Ньютона для уточнения корней (вещественные)
template<typename T>
T newton_method(const std::vector<T>&poly, const std::vector<T>&deriv, T guess, T epsilon) {
    for (size_t i = 0; i < 100; ++i) {
        T f = eval_poly(poly, guess);
        T f_prime = eval_poly(deriv, guess);
        if (abs_val(f_prime) < T(1e-12))
            break;
        T next = guess - f / f_prime;

        T diff = abs_val(next - guess);
        if (diff < epsilon) return next;
        guess = next;
    }
    return guess;
}

// Метод Ньютона для уточнения корней (комлпексные)
// Добавила: Семенидо Алина, КМБО-06-22 lina.semenido@yandex.ru
template<typename Complex, typename Real>
void newton_method(std::vector<Complex>&z, const std::vector<Complex>&c, Real eps) {
    for (Complex&root: z) {
        auto [f, df] = eval_poly_and_deriv(c, root);
        if (abs(df) >= eps) root -= f / df;
    }
}

// Подстчет количества перемен знаков в последовательности коэффициентов
template<typename T>
int sign_changes(const std::vector<T>&coeffs) {
    int changes = 0;
    T prev = coeffs[0];
    for (size_t i = 1; i < coeffs.size(); ++i) {
        if (coeffs[i] != T(0)) {
            changes = 1 + T((prev > 0 && coeffs[i] < 0) || (prev < 0 && coeffs[i] > 0));
            prev = coeffs[i];
        }
    }
    return changes;
}

// Вычисление верхней границы корней полинома
template<typename T>
T computeRootUpperBound(const std::vector<T>&poly) {
    if (poly.empty()) return 0.0;

    // Формула Лагранжа для верхней границы
    T max_neg = 0.0;
    int n = poly.size() - 1;

    for (size_t i = 0; i < n; ++i) {
        if (poly[i] < 0) {
            T ratio = abs(poly[i]) / poly.back();
            if (ratio > max_neg)
                max_neg = ratio;
        }
    }
    return 1.0 + pow_val(max_neg, 1.0 / (n));
}

// Граница Коши для корней многочлена
template <typename T>
T computeRootBoundCauchy(const std::vector<T>& poly) {
    if (poly.empty()) return T(0);
    // Игнорируем ведущие нули для правильной степени
    size_t actual_degree = poly.size() - 1;
    while (actual_degree > 0 && is_zero_val(poly[actual_degree])) { actual_degree--; }

    if (actual_degree == 0) // Многочлен - константа
        return T(0);

    T max_coeff_abs = T(0);
    // Ищем максимум среди abs(a_0), ..., abs(a_{n-1})
    for (size_t i = 0; i < actual_degree; ++i) 
        if (abs(poly[i]) > max_coeff_abs)
            max_coeff_abs = abs(poly[i]);

    T leading_coeff_abs = abs(poly[actual_degree]);
    if (is_zero_val(leading_coeff_abs)) { // Этого не должно произойти, если actual_degree правильно определена
        throw std::runtime_error("Leading coefficient is zero in computeRootBoundCauchy.");
    }
    return T(1) + max_coeff_abs / leading_coeff_abs;
}
