// Заголовочный файл с реализацией дополнительных методов для основных методов
// Реализация: Павлова Анастасия, КМБО-01-22

#pragma once

#include "polynomialUtils.h"
#include "mathUtils.h"

#define _USE_MATH_DEFINES
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

template<typename CoefT, typename Real>
std::complex<Real> find_root_on_circle(const std::vector<CoefT>&poly,
                                       Real modulus, double eps, int N = 256);

// метод Ньютона для уточнения корней (вещественные)
template<typename T>
T newton_method(const std::vector<T>&poly, const std::vector<T>&deriv, T guess, double epsilon) {
    for (size_t i = 0; i < 100; ++i) {
        T f = eval_poly(poly, guess);
        T f_prime = eval_poly(deriv, guess);
        if (abs_val(f_prime) < T(1e-12))
            break;
        T next = guess - f / f_prime;

        // сравниваем в double, чтобы избежать operator<
        double diff = to_double_fallback(abs_val(next - guess));
        if (diff < epsilon)
            return next;

        guess = next;
    }
    return guess;
}

// метод Ньютона для уточнения корней (комлпексные)
// Добавила: Семенидо Алина, КМБО-06-22
template<typename Complex, typename Real>
void newton_method(std::vector<Complex>&z,
                   const std::vector<Complex>&c,
                   Real eps) {
    for (Complex&root: z) {
        auto [f, df] = eval_poly_and_deriv(c, root);

        if (abs(df) >= eps)
            root -= f / df;
    }
}

// автоматически приводит тип к вещественному (modulus -> вещественный тип)
template<typename T>
T find_root_on_circle_adapted(const std::vector<T>&poly, T modulus, double epsilon) {
    using value_type = typename underlying_real_type<T>::type;
    return find_root_on_circle<T, value_type>(poly, static_cast<value_type>(abs_val(modulus)), epsilon);
}

// поиск корней на круге (помогает при проверке и отборе корней)
template<typename CoefT, typename Real>
std::complex<Real> find_root_on_circle(const std::vector<CoefT>&poly,
                                       Real modulus, double eps,
                                       int N) {
    using C = std::complex<Real>;
    auto deriv = derivative(poly);

    C best = C(0);
    Real best_abs = std::numeric_limits<Real>::max();

    for (int i = 0; i < N; ++i) {
        double angle = 2.0 * M_PI * i / N;
        C guess(modulus * std::cos(angle), modulus * std::sin(angle));
        Real val_abs = std::abs(eval_poly(poly, guess));
        if (val_abs < best_abs) {
            best_abs = val_abs;
            best = guess;
        }
    }
    return newton_method(poly, deriv, best, eps);
}

// подстчет количества перемен знаков в последовательности коэффициентов
template<typename T>
int sign_changes(const std::vector<T>&coeffs) {
    int changes = 0;
    T prev = coeffs[0];
    for (size_t i = 1; i < coeffs.size(); ++i) {
        if (coeffs[i] != T(0)) {
            if ((prev > 0 && coeffs[i] < 0) || (prev < 0 && coeffs[i] > 0))
                ++changes;
            prev = coeffs[i];
        }
    }
    return changes;
}

// вычисление верхней границы корней полинома
template<typename T>
T computeRootUpperBound(const std::vector<T>&poly) {
    if (poly.empty())
        return 0.0;

    // метод Лагранжа для верхней границы
    T max_neg = 0.0;
    int n = poly.size() - 1;

    for (size_t i = 0; i < n; ++i) {
        if (poly[i] < 0) {
            T ratio = abs(poly[i]) / poly.back();
            if (ratio > max_neg)
                max_neg = ratio;
        }
    }

    return 1.0 + pow(max_neg, 1.0 / (n));
}
