// Алгоритм CAD
// Теория + реализация: Трудолюбов Никита, КМБО-01-22

#pragma once

#include <map>
#include <set>

// Подключение заголовочных файлов
#include "multipolyprecision.h"
#include "Sagralov_real_roots.h"

using Float = float_precision;

// Шаблонный класс для матрицы
template<typename T>
class Matrix {
public:
    std::vector<std::vector<T>> mat;
    size_t rows = 0, cols = 0;

    Matrix(size_t r = 0, size_t c = 0) : rows(r), cols(c), mat(r, std::vector<T>(c, T())) {
    }

    Matrix getMinor(size_t r, size_t c) const {
        Matrix minor(rows - 1, cols - 1);
        size_t m_r = 0;
        for (size_t i = 0; i < rows; ++i) {
            if (i == r)
                continue;
            size_t m_c = 0;
            for (size_t j = 0; j < cols; ++j) {
                if (j == c)
                    continue;
                minor.mat[m_r][m_c] = mat[i][j];
                m_c++;
            }
            m_r++;
        }
        return minor;
    }
};

// Функции для алгоритма CAD
// Вспомогательные функции
template<typename CoeffType>
bool operator<(const mpolynomial<CoeffType>&a, const mpolynomial<CoeffType>&b) {
    return a.toString() < b.toString();
}

// Вычисление детерминанта (предварительное объявление)
template<typename CoeffType>
mpolynomial<CoeffType> determinant(const Matrix<mpolynomial<CoeffType>>&m);

// Извлечение коэффициентов
template<typename CoeffType>
std::map<size_t, mpolynomial<CoeffType>> get_coefficients(const mpolynomial<CoeffType>&p, size_t var_to_elim) {
    std::map<size_t, mpolynomial<CoeffType>> coeffs;
    if (p.empty())
        return coeffs;
    for (const auto&term: p.get_terms()) {
        const auto&exponents = term.first;
        const CoeffType&coeff_val = term.second;
        size_t power = (var_to_elim < exponents.size()) ? exponents[var_to_elim] : 0;
        mpolynomial<CoeffType> coeff_poly_term;
        std::vector<unsigned int> new_exponents = exponents;
        if (var_to_elim < new_exponents.size())
            new_exponents[var_to_elim] = 0;
        coeff_poly_term.add_term(new_exponents, coeff_val);
        coeffs[power] += coeff_poly_term;
    }
    return coeffs;
}

// Общая функция для вычисления результанта
template<typename CoeffType>
mpolynomial<CoeffType> resultant(const mpolynomial<CoeffType>&p, const mpolynomial<CoeffType>&q, size_t var) {
    size_t n = p.degree(var);
    size_t m = q.degree(var);

    if (p.empty() || q.empty())
        return mpolynomial<CoeffType>();
    if (n == 0) {
        return pow(p, m);
    }
    if (m == 0) {
        return pow(q, n);
    }

    auto p_coeffs = get_coefficients(p, var);
    auto q_coeffs = get_coefficients(q, var);

    Matrix<mpolynomial<CoeffType>> sylvester(n + m, n + m);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j <= n; ++j) {
            sylvester.mat[i][i + j] = p_coeffs.count(n - j) ? p_coeffs.at(n - j) : mpolynomial<CoeffType>();
        }
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= m; ++j) {
            sylvester.mat[i + m][i + j] = q_coeffs.count(m - j) ? q_coeffs.at(m - j) : mpolynomial<CoeffType>();
        }
    }

    return determinant(sylvester);
}

// Реализация детерминанта
template<typename CoeffType>
mpolynomial<CoeffType> determinant(const Matrix<mpolynomial<CoeffType>>&m) {
    if (m.rows != m.cols)
        throw std::runtime_error("Matrix must be square");
    size_t n = m.rows;
    if (n == 0)
        return mpolynomial<CoeffType>(CoeffType(1));
    if (n == 1)
        return m.mat[0][0];
    if (n == 2)
        return (m.mat[0][0] * m.mat[1][1]) - (m.mat[0][1] * m.mat[1][0]);

    mpolynomial<CoeffType> det;
    mpolynomial<CoeffType> sign(CoeffType(1));
    for (size_t j = 0; j < n; ++j) {
        det += sign * m.mat[0][j] * determinant(m.getMinor(0, j));
        sign *= mpolynomial<CoeffType>(CoeffType(-1));
    }
    return det;
}

// Дискриминант по общей формуле
template<typename CoeffType>
mpolynomial<CoeffType> discriminant(const mpolynomial<CoeffType>&p, size_t var_to_elim) {
    size_t n = p.degree(var_to_elim);
    if (n < 1)
        return mpolynomial<CoeffType>();

    auto coeffs = get_coefficients(p, var_to_elim);
    mpolynomial<CoeffType> an = coeffs.count(n) ? coeffs.at(n) : mpolynomial<CoeffType>();
    if (an.empty())
        return mpolynomial<CoeffType>();

    mpolynomial<CoeffType> p_deriv = p.derivative(var_to_elim);

    mpolynomial<CoeffType> sign(CoeffType(1));
    if ((n * (n - 1) / 2) % 2 != 0) {
        sign = mpolynomial<CoeffType>(CoeffType(-1));
    }

    mpolynomial<CoeffType> res = resultant(p, p_deriv, var_to_elim);

    // Формула: sign * (1/an) * res
    // Для CAD важны корни, которые содержатся в res, поэтому деление на an (полином) опускаем.
    return sign * res;
}

// Главная функция проекции
template<typename CoeffType>
std::vector<mpolynomial<CoeffType>> project_one_polynomial(const mpolynomial<CoeffType>&p, size_t var_to_elim) {
    std::set<mpolynomial<CoeffType>> projection_set;
    auto add_to_set = [](std::set<mpolynomial<CoeffType>>&s, const mpolynomial<CoeffType>&poly) {
        if (!poly.empty())
            s.insert(poly);
    };

    std::cout << "  1. Проецирование коэффициентов..." << std::endl;
    auto coeffs = get_coefficients(p, var_to_elim);
    for (const auto&pair: coeffs) {
        add_to_set(projection_set, pair.second);
    }

    std::cout << "  2. Вычисление и проецирование дискриминанта..." << std::endl;
    mpolynomial<CoeffType> disc = discriminant(p, var_to_elim);
    add_to_set(projection_set, disc);
    std::cout << "     (Результат(P, P'): " << disc.toString() << ")" << std::endl;

    std::vector<mpolynomial<CoeffType>> result(projection_set.begin(), projection_set.end());
    return result;
}

// Функция на основе метода Сагралова
template<typename CoeffType>
std::vector<Float> find_real_roots(const mpolynomial<CoeffType>&p, size_t var_index) {
    if (p.empty()) {
        return {};
    }

    // Шаг 1: Преобразуем полином в вектор коэффициентов
    std::map<size_t, CoeffType> coeffs_map;
    size_t max_deg = 0;
    for (const auto&term: p.get_terms()) {
        bool is_univariate_term = true;
        size_t current_deg = 0;
        for (size_t i = 0; i < term.first.size(); ++i) {
            if (i == var_index) {
                current_deg = term.first[i];
            }
            else if (term.first[i] != 0) {
                is_univariate_term = false;
                break;
            }
        }
        if (!is_univariate_term) {
            std::cerr << "Ошибка: find_real_roots вызван для полинома, не являющегося унивариантным: " << p.toString()
                    << std::endl;
            return {};
        }
        if (current_deg > max_deg)
            max_deg = current_deg;
        coeffs_map[current_deg] += term.second;
    }

    if (max_deg == 0)
        return {};

    // Шаг 2: Конвертируем коэффициенты в std::vector<Float>
    std::vector<Float> float_coeffs(max_deg + 1, Float(0));
    for (const auto&pair: coeffs_map) {
        float_coeffs[pair.first] = Float(pair.second) * Float(pair.first <= max_deg); // Замена условия if
    }

    while (float_coeffs.size() > 1 && is_zero_val(float_coeffs.back())) {
        float_coeffs.pop_back();
    }

    try {
        // Шаг 3: Вызываем метод изоляции корней из Sagralov.h
        std::vector<std::pair<Float, Float>> isolated_intervals = find_real_roots_by_Sagralov<Float>(float_coeffs);

        // Шаг 4: Преобразуем интервалы в вектор корней (берем середины)
        std::vector<Float> roots;
        for (const auto&interval: isolated_intervals) {
            roots.push_back((interval.first + interval.second) / Float(2));
        }
        return roots;
    }
    catch (const std::exception&e) {
        std::cerr << "Произошла ошибка при поиске корней: " << e.what() << std::endl;
        return {};
    }
}

template<typename CoeffType>
void lift_recursive(
    // Тип карты для текущей точки изменен на Float
    std::map<size_t, Float>&current_sample,
    const std::vector<std::vector<mpolynomial<CoeffType>>>&lifting_polys,
    const mpolynomial<CoeffType>&original_poly,
    size_t target_level,
    std::vector<std::map<size_t, Float>>&solutions) {
    // Базовый случай: достигли нужной глубины рекурсии
    if (target_level == lifting_polys.size()) {
        mpolynomial<CoeffType> final_eval = original_poly.evaluate(current_sample);
        if (final_eval.is_zero()) {
            solutions.push_back(current_sample);
        }
        return;
    }

    // Шаг рекурсии
    const std::vector<mpolynomial<CoeffType>>&polys_at_level = lifting_polys[target_level];

    std::set<Float> unique_roots;
    for (const auto&poly: polys_at_level) {
        mpolynomial<CoeffType> p_eval = poly.evaluate(current_sample);

        std::vector<Float> roots = find_real_roots(p_eval, target_level);
        unique_roots.insert(roots.begin(), roots.end());
    }

    std::vector<Float> samples;
    if (!unique_roots.empty()) {
        std::vector<Float> sorted_roots(unique_roots.begin(), unique_roots.end());

        samples.push_back(sorted_roots.front() - Float(1));
        for (size_t i = 0; i < sorted_roots.size(); ++i) {
            samples.push_back(sorted_roots[i]);
            if (i + 1 < sorted_roots.size()) {
                samples.push_back((sorted_roots[i] + sorted_roots[i + 1]) / Float(2));
            }
        }
        samples.push_back(sorted_roots.back() + Float(1));
    }
    else {
        samples.push_back(Float(0));
    }

    // Рекурсивный вызов для каждой пробной точки
    for (const Float&s: samples) {
        current_sample[target_level] = s;
        lift_recursive(current_sample, lifting_polys, original_poly, target_level + 1, solutions);
        current_sample.erase(target_level);
    }
}

template<typename CoeffType>
void find_roots_by_CAD(const mpolynomial<CoeffType>&initial_poly, int num_vars) {
    std::cout << "\n========================================================" << std::endl;
    std::cout << "ЗАПУСК МЕТОДА CAD для P = " << initial_poly.toString() << std::endl;
    std::cout << "========================================================" << std::endl;

    if (num_vars == 1) {
        std::cout << "\n--- РЕШЕНИЕ УРАВНЕНИЯ С ОДНОЙ ПЕРЕМЕННОЙ ---" << std::endl;

        std::vector<Float> roots = find_real_roots(initial_poly, 0);

        if (roots.empty()) {
            std::cout << "\nРешений не найдено." << std::endl;
        }
        else {
            std::cout << "\n>>> НАЙДЕНЫ СЛЕДУЮЩИЕ РЕШЕНИЯ (" << roots.size() << " шт.):" << std::endl;
            int count = 1;
            for (const auto&r: roots) {
                std::cout << "  " << count++ << ". x_0 = " << r << ";" << std::endl;
            }
            std::cout << "<<<" << std::endl;
        }
        return;
    }

    // ПОЛНЫЙ АЛГОРИТМ ДЛЯ НЕСКОЛЬКИХ ПЕРЕМЕННЫХ
    std::cout << "\n--- ФАЗА 1: ПРОЕКЦИЯ ---" << std::endl;
    std::vector<std::vector<mpolynomial<CoeffType>>> lifting_polys(num_vars);
    if (num_vars > 0) {
        lifting_polys[num_vars - 1] = {initial_poly};
    }
    std::vector<mpolynomial<CoeffType>> current_level_polys = {initial_poly};
    for (int i = num_vars - 1; i > 0; --i) {
        std::set<mpolynomial<CoeffType>> next_level_set;
        for (const auto&p: current_level_polys) {
            if (p.degree(i) == 0) {
                next_level_set.insert(p);
                continue;
            }
            std::vector<mpolynomial<CoeffType>> projected = project_one_polynomial(p, i);
            for (const auto&proj_p: projected) {
                if (!proj_p.empty() && !proj_p.is_zero()) {
                    next_level_set.insert(proj_p);
                }
            }
        }
        if (next_level_set.empty())
            break;
        current_level_polys.assign(next_level_set.begin(), next_level_set.end());
        lifting_polys[i - 1] = current_level_polys;
    }
    std::cout << "\n--- ФАЗА 2: ПОДНЯТИЕ ---" << std::endl;

    std::vector<std::map<size_t, Float>> all_solutions;
    std::map<size_t, Float> sample_point;

    lift_recursive(sample_point, lifting_polys, initial_poly, 0, all_solutions);

    if (all_solutions.empty()) {
        std::cout << "\nРешений в рамках выполненного поиска не найдено." << std::endl;
    }
    else {
        std::cout << "\n\n>>> НАЙДЕНЫ СЛЕДУЮЩИЕ РЕШЕНИЯ (" << all_solutions.size() << " шт.):" << std::endl;
        int count = 1;
        for (const auto&sol: all_solutions) {
            std::cout << "  " << count++ << ". ";

            std::map<size_t, CoeffType> sorted_sol(sol.begin(), sol.end());
            for (const auto&[var, val]: sorted_sol) {
                std::cout << "x_" << var << " = " << val << "; ";
            }
            std::cout << std::endl;
        }
        std::cout << "<<<" << std::endl;
    }
}
