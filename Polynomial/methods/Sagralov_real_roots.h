// Метод Сагралова для изоляции вещественных корней
// Теория: Павлова Анастасия, КМБО-01-22
// Реализация: Корастелин Никита, КМБО-01-23

#pragma once

#include "mathUtils.h"
#include "polynomialUtils.h"
#include "Helper_for_all_methods.h"

// Функция для изоляции вещественных корней многочлена методом Сагралова
template <typename T>
std::vector<std::pair<T, T>> find_real_roots_by_Sagralov(const std::vector<T>& input_coeffs, double epsilon = 1e-8) {
    if (input_coeffs.empty()) {
        throw std::invalid_argument("Polynomial must not be empty");
    }

    // Нормализуем многочлен (делаем старший коэффициент равным 1)
    std::vector<T> coeffs = normalizePolynomial(input_coeffs);
    std::vector<std::pair<T, T>> isolated_intervals;

    // Проверяем, есть ли корень в нуле
    if (is_zero(eval_poly(coeffs, T(0)), epsilon)) {
        // Если есть корень в нуле, обрабатываем отдельно
        isolated_intervals.emplace_back(-epsilon, epsilon);

        // Дефлируем многочлен (делим на x)
        coeffs = deflate_poly(coeffs, T(0));
    }

    // Определяем границы поиска корней
    T upperBound = computeRootBoundCauchy(coeffs);
    T lowerBound = -upperBound;

    // Строим последовательность Штурма
    std::vector<std::vector<T>> sturm_sequence = { coeffs, derivative(coeffs) };

    while (sturm_sequence.back().size() > 1) {
        const auto& p1 = sturm_sequence[sturm_sequence.size() - 2];
        const auto& p2 = sturm_sequence[sturm_sequence.size() - 1];

        // Вычисляем остаток от деления p1 на p2 с обратным знаком
        auto remainder = polyRemainder(p1, p2);
        for (auto& coeff : remainder) {
            coeff = -coeff;
        }

        sturm_sequence.push_back(remainder);
    }

    // Рекурсивно изолируем корни
    isolateRootsRecursive(coeffs, sturm_sequence, lowerBound, upperBound, isolated_intervals, epsilon);

    // Уточняем границы интервалов
    for (auto& interval : isolated_intervals) {
        refineInterval(interval, coeffs, epsilon);
    }

    return isolated_intervals;
}

// Вспомогательная функция для вычисления остатка от деления многочленов
template <typename T>
std::vector<T> polyRemainder(const std::vector<T>& a, const std::vector<T>& b) {
    if (b.empty()) {
        throw std::invalid_argument("Division by zero polynomial");
    }

    std::vector<T> remainder = a;
    std::vector<T> divisor = b;

    while (remainder.size() >= divisor.size()) {
        T scale = remainder.back() / divisor.back();
        size_t shift = remainder.size() - divisor.size();

        for (size_t i = 0; i < divisor.size(); ++i) {
            remainder[i + shift] -= scale * divisor[i];
        }

        // Удаляем ведущие нули
        while (!remainder.empty() && is_zero(remainder.back())) {
            remainder.pop_back();
        }
    }

    return remainder.empty() ? std::vector<T>{T(0)} : remainder;
}

// Вспомогательная функция для подсчета перемен знаков в последовательности Штурма
template <typename T>
int countSignChanges(const std::vector<std::vector<T>>& sturm_sequence, T x) {
    int changes = 0;
    T prev = eval_poly(sturm_sequence[0], x);

    for (size_t i = 1; i < sturm_sequence.size(); ++i) {
        T current = eval_poly(sturm_sequence[i], x);
        changes += (!is_zero(current) && ((prev > T(0) && current < T(0)) || (prev < T(0) && current > T(0)))); // тут убран if
        prev = current;
    }

    return changes;
}

// Рекурсивная функция для изоляции корней
template <typename T>
void isolateRootsRecursive(
    const std::vector<T>& poly,
    const std::vector<std::vector<T>>& sturm_sequence,
    T a, T b,
    std::vector<std::pair<T, T>>& result,
    double epsilon) {

    int rootCount = countSignChanges(sturm_sequence, a) - countSignChanges(sturm_sequence, b);

    if (rootCount == 0) {
        return; // Нет корней на интервале
    }

    if (rootCount == 1) {
        // Нашли интервал с ровно одним корнем
        result.emplace_back(a, b);
        return;
    }

    // Делим интервал пополам и рекурсивно обрабатываем каждую половину
    T mid = (a + b) / T(2);

    isolateRootsRecursive(poly, sturm_sequence, a, mid, result, epsilon);
    isolateRootsRecursive(poly, sturm_sequence, mid, b, result, epsilon);
}

// Функция для уточнения границ интервала
template <typename T>
void refineInterval(std::pair<T, T>& interval, const std::vector<T>& poly, double epsilon) {
    T a = interval.first;
    T b = interval.second;

    // Проверяем знаки на концах интервала
    T fa = eval_poly(poly, a);
    T fb = eval_poly(poly, b);

    if (is_zero(fa, epsilon)) {
        interval.second = a + epsilon;
        return;
    }

    if (is_zero(fb, epsilon)) {
        interval.first = b - epsilon;
        return;
    }

    // Используем метод Ньютона для уточнения
    auto deriv = derivative(poly);
    T guess = (a + b) / T(2);
    T root = newton_method(poly, deriv, guess, epsilon);

    // Обновляем интервал с небольшим запасом
    interval.first = root - epsilon;
    interval.second = root + epsilon;
}

template <typename T>
T computeRootBoundCauchy(const std::vector<T>& poly) {
    if (poly.empty()) return 0.0;
    T max_coeff = 0.0;
    for (size_t i = 0; i < poly.size() - 1; ++i) {
        if (abs(poly[i]) > max_coeff) {
            max_coeff = abs(poly[i]);
        }
    }
    return 1.0 + max_coeff / abs(poly.back());
}