// Модификация Хосмана для комплексных корней на основе метода Грефе
// hosman_modification_graeffe - поиск корней без учета кратности (порождает кучу мусора)
// hosman_with_multiplicities - верся с учетом кратности 

#pragma once

#include "mathUtils.h"
#include "polynomialUtils.h"
#include "Helper_for_all_methods.h"
#include "Graeffe_real_roots.h"
#include <algorithm>
#include <map>
#include <set>

// Вспомогательные функции для метода Хосмана
// Вычисление полинома в комплексной точке
template <typename T>
std::complex<T> eval_complex_poly(const std::vector<T>& coeffs_desc, const std::complex<T>& z)
{
    if (coeffs_desc.empty())
        return std::complex<T>(T(0), T(0));
    std::complex<T> result = std::complex<T>(coeffs_desc[0], T(0));
    for (size_t i = 1; i < coeffs_desc.size(); ++i)
        result = result * z + std::complex<T>(coeffs_desc[i], T(0));
    return result;
}

// Факториал
template <typename T>
T factorial(int n)
{
    T result = T(1);
    for (int i = 2; i <= n; ++i)
        result = result * T(i);
    return result;
}

// Умножение двух полиномов (ascending порядок)
template <typename T>
std::vector<T> multiply_polynomials(const std::vector<T>& p1, const std::vector<T>& p2)
{
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
std::vector<T> add_polynomials(const std::vector<T>& p1, const std::vector<T>& p2)
{
    size_t max_size = std::max(p1.size(), p2.size());
    std::vector<T> result(max_size, T(0));
    for (size_t i = 0; i < p1.size(); ++i)
        result[i] = result[i] + p1[i];
    for (size_t i = 0; i < p2.size(); ++i)
        result[i] = result[i] + p2[i];
    return result;
}

// Вычисление (r^2 - x^2)^k как полинома от x (ascending порядок)
template <typename T>
std::vector<T> power_of_binomial(T r_squared, int k)
{
    if (k == 0)
        return { T(1) };

    // Начинаем с (r^2 - x^2)
    std::vector<T> result = { r_squared, T(0), T(-1) }; // (r^2 - x^2)

    // Возводим в степень последовательным умножением
    for (int i = 1; i < k; ++i)
    {
        std::vector<T> temp(result.size() + 2, T(0));
        for (size_t j = 0; j < result.size(); ++j)
        {
            temp[j] = temp[j] + result[j] * r_squared;
            if (j + 2 < temp.size())
                temp[j + 2] = temp[j + 2] - result[j];
        }
        result = temp;
    }
    return result;
}

// Основная функция: метод Хосмана
template <typename T>
std::vector<std::complex<T>> hosman_modification_graeffe (const std::vector<T>& coeffs_desc, T epsilon = T(1e-10L), int maxIter = 100) {
    // Находим модули корней методом Греффе
    std::vector<T> graeffe_moduli = find_moduli_roots_by_graeffe(coeffs_desc, epsilon, maxIter);
    if (graeffe_moduli.empty()) return {}; // Graeffe не нашёл модулей

    // Подсчитываем кратности модулей
    std::vector<T> moduli = graeffe_moduli;

    std::sort(moduli.begin(), moduli.end(), [](const T& a, const T& b) { return a < b; });

    // Группируем модули: считаем количество для каждого уникального значения
    std::vector<std::pair<T, int>> unique_moduli; // (модуль, count)

    size_t i = 0;
    while (i < moduli.size()) {
        T current_mod = moduli[i];
        int count = 0;
        size_t start_i = i;

        // Считаем все модули близкие к current_mod (порог больше для накопленных ошибок Греффе)
        while (i < moduli.size() && abs_val(moduli[i] - current_mod) < T(1e-7)) {
            count++;
            i++;
        }

        unique_moduli.push_back({ current_mod, count });
    }

    // Нормализация полинома
    size_t first_nz = 0;
    while (first_nz < coeffs_desc.size() && is_zero_val(coeffs_desc[first_nz]))
        ++first_nz;
    if (first_nz >= coeffs_desc.size())
        return {};

    std::vector<T> coeffs(coeffs_desc.begin() + first_nz, coeffs_desc.end());
    size_t deg = coeffs.size() - 1;
    if (deg == 0)
        return {};

    T lead = coeffs.front();
    for (auto& c : coeffs)
        c = c / lead;

    // Преобразуем в ascending порядок
    std::vector<T> f_asc(coeffs.rbegin(), coeffs.rend());

    // Вычисляем все производные
    std::vector<std::vector<T>> derivatives;
    derivatives.push_back(f_asc);
    for (size_t ord = 1; ord <= deg; ++ord) {
        derivatives.push_back(derivative(derivatives.back()));
        if (derivatives.back().empty())
            break;
    }

    // Вычисляем масштаб полинома для адаптивных порогов
    T max_coeff = T(0);
    for (const auto& c : coeffs_desc)
        max_coeff = std::max(max_coeff, abs_val(c));

    std::vector<std::complex<T>> roots;
    roots.reserve(graeffe_moduli.size() * 2);

    // Для каждого уникального модуля применяем метод Хосмана
    for (const auto& [r, multiplicity] : unique_moduli) {
        if (r <= epsilon) continue;

        T r_squared = r * r;

        // Вычисляем ЛОКАЛЬНЫЕ адаптивные пороги для текущего модуля r
        T r_pow2 = r * r;
        T local_scale_factor = max_coeff * r_pow2;
        if (local_scale_factor < T(1))
            local_scale_factor = T(1);

        // Базовые множители
        T base_eval_mult = T(1000000000); // 1e-3 базовый
        T base_resid_mult = T(100000);    // 1e-7 базовый

        // Для больших local_scale_factor добавляем дополнительный множитель
        if (local_scale_factor > T(1e6))
            base_resid_mult = base_resid_mult * T(1000000); // дополнительно 1e6
        else if (local_scale_factor > T(1e3))
            base_resid_mult = base_resid_mult * T(1000); // дополнительно 1e3

        const T evaluation_threshold = epsilon * base_eval_mult * local_scale_factor;
        const T residual_threshold = epsilon * base_resid_mult * local_scale_factor;
        const T relaxed_residual_threshold = epsilon * base_eval_mult * local_scale_factor;

        // Построение уравнения Хосмана
        // Разложение Тейлора
        std::vector<T> hosman_eq = f_asc;

        for (int k = 2; k < static_cast<int>(derivatives.size()); k += 2)
        {
            if (derivatives[k].empty())
                break;

            std::vector<T> binomial_power = power_of_binomial(r_squared, k / 2);
            std::vector<T> term = multiply_polynomials(derivatives[k], binomial_power);

            T fact_k = factorial<T>(k);
            for (auto& coef : term)
                coef = coef / fact_k;

            // Знак: (-1)^(k/2)
            if ((k / 2) % 2 == 1)
                for (auto& coef : term)
                    coef = -coef;

            hosman_eq = add_polynomials(hosman_eq, term);
        }

        // Переводим в descending для find_moduli_roots_by_graeffe
        std::vector<T> hosman_eq_desc(hosman_eq.rbegin(), hosman_eq.rend());

        // Находим кандидаты для x
        std::vector<std::pair<T, bool>> candidate_reals; // (value, from_grid_sampling)

        // Проверяем x=0 явно
        T val_at_zero = eval_poly(hosman_eq, T(0));
        if (abs_val(val_at_zero) < evaluation_threshold)
            candidate_reals.push_back({ T(0), false });

        // Находим модули корней уравнения Хосмана
        std::vector<T> candidate_moduli = find_moduli_roots_by_graeffe(hosman_eq_desc, epsilon, maxIter);

        // Тестируем +- x_mod
        for (const T& x_mod : candidate_moduli)
        {
            T val_pos = eval_poly(hosman_eq, x_mod);
            T val_neg = eval_poly(hosman_eq, -x_mod);

            // Тестируем кандидатов из Hosman equation
            if (abs_val(val_pos) < evaluation_threshold)
                candidate_reals.push_back({ x_mod, false });
            if (abs_val(val_neg) < evaluation_threshold && abs_val(x_mod) > epsilon)
                candidate_reals.push_back({ -x_mod, false });
        }

        // Если ни один кандидат не прошёл, делаем grid sampling с мягким порогом
        if (candidate_reals.empty())
        {
            // Grid sampling активирован
            int num_samples = 100;
            T grid_threshold = evaluation_threshold * T(10); // В 10 раз мягче для grid

            for (int i = 1; i <= num_samples; ++i)
            {
                T x_test = r * T(i) / T(num_samples + 1);
                T val_pos = eval_poly(hosman_eq, x_test);
                T val_neg = eval_poly(hosman_eq, -x_test);

                if (abs_val(val_pos) < grid_threshold)
                {
                    candidate_reals.push_back({ x_test, true }); // от grid
                }
                if (abs_val(val_neg) < grid_threshold && abs_val(x_test) > epsilon)
                {
                    candidate_reals.push_back({ -x_test, true }); // от grid
                }
            }
            // Grid sampling завершён
        }

        if (candidate_reals.empty())
        {
            continue; // Нет кандидатов для данного r
        }

        // Кандидаты найдены

        // ШАГ 2.3: Валидация кандидатов уравнением (1.6)
        std::vector<std::pair<T, T>> valid_roots; // (a, residual)

        for (const auto& [a, from_grid] : candidate_reals)
        {
            T b_squared = r_squared - a * a;
            if (b_squared < -epsilon)
                continue;
            if (b_squared < T(0))
                b_squared = T(0);

            // Для комплексных корней: проверяем f(a±bi) на исходном полиноме
            if (b_squared >= epsilon)
            {
                T b = sqrt_val(b_squared);
                std::complex<T> z_plus(a, b);
                std::complex<T> z_minus(a, -b);

                T f_plus_abs = abs_val(eval_complex_poly(coeffs_desc, z_plus));
                T f_minus_abs = abs_val(eval_complex_poly(coeffs_desc, z_minus));

                if (f_plus_abs > evaluation_threshold || f_minus_abs > evaluation_threshold)
                {
                    continue; // Кандидат не проходит pre-check
                }
            }

            T abs_residual;
            T current_threshold = from_grid ? relaxed_residual_threshold : residual_threshold;

            // Для вещественных корней (b²≈0) проверяем f(a)=0
            if (b_squared < epsilon)
            {
                // Вещественный корень
                T f_a = eval_poly(derivatives[0], a);
                abs_residual = abs_val(f_a);
            }
            else
            {
                // Комплексный корень: уравнение (1.6)
                // f'(a) - f'''(a)·b²/3! + f⁵(a)·b⁴/5! - ...
                T residual = eval_poly(derivatives[1], a);
                T b_power_k = b_squared;

                for (int k = 3; k < static_cast<int>(derivatives.size()); k += 2)
                {
                    if (derivatives[k].empty())
                        break;

                    T term = eval_poly(derivatives[k], a) * b_power_k / factorial<T>(k);
                    if (((k - 1) / 2) % 2 == 1)
                        term = -term;

                    residual = residual + term;
                    b_power_k = b_power_k * b_squared;
                }

                abs_residual = abs_val(residual);
            }

            if (abs_residual < current_threshold)
                valid_roots.push_back({ a, abs_residual });
        }

        // Если после валидации не нашли корней, пробуем grid sampling
        if (valid_roots.empty())
        {
            // Grid rescue: пробуем найти корни вручную

            // Сначала проверяем вещественные корни ±r
            T f_plus_r = eval_poly(derivatives[0], r);
            T f_minus_r = eval_poly(derivatives[0], -r);

            if (abs_val(f_plus_r) < evaluation_threshold)
            {
                valid_roots.push_back({ r, abs_val(f_plus_r) });
            }
            if (abs_val(f_minus_r) < evaluation_threshold)
            {
                valid_roots.push_back({ -r, abs_val(f_minus_r) });
            }

            // Если вещественных не нашли, продолжаем grid sampling для комплексных
            if (valid_roots.empty())
            {
                int num_samples = 100;
                T grid_eval_threshold = evaluation_threshold * T(1000); // ОЧЕНЬ мягко для grid rescue
                T grid_resid_threshold = relaxed_residual_threshold;    // используем relaxed

                for (int i = 1; i <= num_samples; ++i)
                {
                    T a_test = r * T(i) / T(num_samples + 1);
                    T b_squared = r_squared - a_test * a_test;

                    if (b_squared < T(0))
                        continue;

                    // Проверяем Hosman equation
                    T val_pos = eval_poly(hosman_eq, a_test);
                    T val_neg = eval_poly(hosman_eq, -a_test);

                    // Проверяем Hosman equation
                    if (abs_val(val_pos) < grid_eval_threshold)
                    {
                        // Если b² мало, возможен вещественный корень ≈±r (уже проверили выше)
                        if (b_squared >= r * epsilon)
                        {
                            T b = sqrt_val(b_squared);
                            std::complex<T> z_plus(a_test, b);
                            std::complex<T> z_minus(a_test, -b);

                            T f_plus_abs = abs_val(eval_complex_poly(coeffs_desc, z_plus));
                            T f_minus_abs = abs_val(eval_complex_poly(coeffs_desc, z_minus));

                            if (f_plus_abs <= evaluation_threshold && f_minus_abs <= evaluation_threshold)
                            {
                                valid_roots.push_back({ a_test, abs_val(val_pos) });
                            }
                        }
                    }

                    // Для отрицательного -a_test
                    if (abs_val(val_neg) < grid_eval_threshold && abs_val(a_test) > epsilon)
                    {
                        // Если b² мало, возможен вещественный корень ≈±r (уже проверили выше)
                        if (b_squared >= r * epsilon)
                        {
                            T b = sqrt_val(b_squared);
                            std::complex<T> z_plus(-a_test, b);
                            std::complex<T> z_minus(-a_test, -b);

                            T f_plus_abs = abs_val(eval_complex_poly(coeffs_desc, z_plus));
                            T f_minus_abs = abs_val(eval_complex_poly(coeffs_desc, z_minus));

                            if (f_plus_abs <= evaluation_threshold && f_minus_abs <= evaluation_threshold)
                            {
                                valid_roots.push_back({ -a_test, abs_val(val_neg) });
                            }
                        }
                    }
                }
            } // Конец grid sampling
        }

        if (valid_roots.empty())
            continue;

        // ШАГ 2.4: Сортируем по residual и добавляем все валидные корни
        std::sort(valid_roots.begin(), valid_roots.end(),
            [](const std::pair<T, T>& x, const std::pair<T, T>& y)
            {
                return x.second < y.second;
            });

        // Для каждого валидного 'a', добавляем соответствующие корни
        // Избегаем дубликатов по 'a'
        std::set<T> processed_a;
        for (const auto& [a, res] : valid_roots)
        {
            // Проверяем, не обработали ли уже это 'a'
            bool duplicate = false;
            for (const T& prev_a : processed_a)
            {
                if (abs_val(prev_a - a) < epsilon * T(5))
                {
                    duplicate = true;
                    break;
                }
            }
            if (duplicate)
                continue;

            processed_a.insert(a);

            T b_squared = r_squared - a * a;
            if (b_squared < T(0))
                b_squared = T(0);
            T b = sqrt_val(b_squared);

            if (abs_val(b) > epsilon * T(5))
            {
                // Комплексная пара a±bi
                roots.push_back(std::complex<T>(a, b));
                roots.push_back(std::complex<T>(a, -b));
            }
            else
            {
                // Реальный корень a
                roots.push_back(std::complex<T>(a, T(0)));
            }
        }
    }

    // ШАГ 3: Улучшенная дедупликация с выбором по минимальной невязке
    struct RootWithResidual
    {
        std::complex<T> root;
        T residual;
    };

    std::vector<RootWithResidual> roots_with_residuals;
    roots_with_residuals.reserve(roots.size());

    // Вычисляем невязку для каждого корня
    for (const auto& root : roots)
    {
        std::complex<T> root_t(root.real(), root.imag());
        auto eval = eval_complex_poly(coeffs_desc, root_t);
        T residual = abs_val(eval);
        roots_with_residuals.push_back({ root, residual });
    }

    // Дедупликация с улучшенным порогом (адаптивный к модулю корня)
    auto nearly_equal = [&](const std::complex<T>& a, const std::complex<T>& b)
        {
            T modulus_a = sqrt_val(a.real() * a.real() + a.imag() * a.imag());
            T modulus_b = sqrt_val(b.real() * b.real() + b.imag() * b.imag());
            T avg_modulus = (modulus_a + modulus_b) / T(2);

            // Используем умеренный относительный порог: 0.5% от модуля + абсолютный порог
            T relative_threshold = avg_modulus * T(0.005);
            T absolute_threshold = epsilon * T(100);
            T threshold = absolute_threshold + relative_threshold;

            T diff_real = abs_val(a.real() - b.real());
            T diff_imag = abs_val(a.imag() - b.imag());

            return diff_real < threshold && diff_imag < threshold;
        };

    std::vector<std::complex<T>> unique_roots;
    for (size_t i = 0; i < roots_with_residuals.size(); ++i)
    {
        bool is_duplicate = false;
        int best_existing_idx = -1;

        // Проверяем, есть ли похожий корень в unique_roots
        for (size_t j = 0; j < unique_roots.size(); ++j)
        {
            if (nearly_equal(roots_with_residuals[i].root, unique_roots[j]))
            {
                is_duplicate = true;
                best_existing_idx = static_cast<int>(j);
                break;
            }
        }

        if (!is_duplicate)
        {
            // Новый уникальный корень
            unique_roots.push_back(roots_with_residuals[i].root);
        }
        else
        {
            // Есть дубликат - выбираем корень с меньшей невязкой
            // Найдем невязку текущего корня в unique_roots
            std::complex<T> existing_root_t(
                unique_roots[best_existing_idx].real(),
                unique_roots[best_existing_idx].imag());
            auto current_eval = eval_complex_poly(coeffs_desc, existing_root_t);
            T current_residual = abs_val(current_eval);

            if (roots_with_residuals[i].residual < current_residual)
            {
                // Новый корень лучше - заменяем
                unique_roots[best_existing_idx] = roots_with_residuals[i].root;
            }
        }
    }

    return unique_roots;
}

// ============ Версия метода Хосмана с учётом кратностей модулей ============

/**
 * @brief Модифицированный метод Хосмана, использующий информацию о кратностях модулей
 *
 * Ключевое отличие: для каждого модуля с кратностью k находит всех кандидатов,
 * сортирует их по невязке |f(z)| и возвращает ровно k лучших.
 *
 * Решает проблему избыточного числа корней: вместо 24+ корней для полинома степени 5
 * возвращает ровно 5 корней (соответствует сумме кратностей).
 *
 * @param coeffs_desc Коэффициенты полинома в порядке убывания степеней
 * @param moduli_with_multiplicities Модули с кратностями от find_moduli_with_multiplicities_by_graeffe()
 * @param epsilon Порог для определения нуля
 * @return Вектор комплексных корней, количество = sum(multiplicities)
 */
template <typename T>
std::vector<std::complex<T>> hosman_with_multiplicities(
    const std::vector<T>& coeffs_desc,
    const std::vector<ModulusWithMultiplicity<T>>& moduli_with_multiplicities,
    T epsilon = T(1e-10L))
{

    if (moduli_with_multiplicities.empty())
    {
        return {};
    }

    // Нормализация полинома
    size_t first_nz = 0;
    while (first_nz < coeffs_desc.size() && is_zero_val(coeffs_desc[first_nz]))
        ++first_nz;
    if (first_nz >= coeffs_desc.size())
        return {};

    std::vector<T> coeffs(coeffs_desc.begin() + first_nz, coeffs_desc.end());
    size_t deg = coeffs.size() - 1;
    if (deg == 0)
        return {};

    T lead = coeffs.front();
    for (auto& c : coeffs)
        c = c / lead;

    std::vector<T> f_asc(coeffs.rbegin(), coeffs.rend());

    // Вычисление производных полинома
    std::vector<std::vector<T>> derivatives;
    derivatives.push_back(f_asc);

    for (int deriv_order = 1; deriv_order <= static_cast<int>(deg); ++deriv_order)
    {
        const auto& prev = derivatives.back();
        if (prev.size() <= 1)
        {
            derivatives.push_back({});
            continue;
        }

        std::vector<T> current_deriv(prev.size() - 1);
        for (size_t i = 1; i < prev.size(); ++i)
            current_deriv[i - 1] = prev[i] * T(i);

        derivatives.push_back(current_deriv);
    }

    struct CandidateRoot
    {
        std::complex<T> root;
        T residual;
    };

    std::vector<CandidateRoot> all_final_roots;

    T evaluation_threshold = epsilon * T(1000);
    T residual_threshold = epsilon * T(10);
    T relaxed_residual_threshold = epsilon * T(100);

    // Обработка каждого модуля с учётом его кратности
    for (const auto& mod_with_mult : moduli_with_multiplicities)
    {
        T r = mod_with_mult.value;
        int multiplicity = mod_with_mult.multiplicity;

        if (r < epsilon)
        {
            for (int i = 0; i < multiplicity; ++i)
            {
                all_final_roots.push_back({ std::complex<T>(T(0), T(0)), T(0) });
            }
            continue;
        }

        T r_squared = r * r;

        // Построение уравнения Хосмана H(x)
        std::vector<T> hosman_eq = derivatives[0];
        T r_power_k = r_squared;

        for (int k = 1; k <= static_cast<int>(deg) / 2; ++k)
        {
            int deriv_order = 2 * k;
            if (deriv_order >= static_cast<int>(derivatives.size()))
                break;
            if (derivatives[deriv_order].empty())
                break;

            std::vector<T> binomial_term = power_of_binomial(r_squared, k);
            std::vector<T> derivative_term = derivatives[deriv_order];

            T factorial_val = factorial<T>(deriv_order);
            for (auto& c : derivative_term)
                c = c / factorial_val;

            std::vector<T> combined = multiply_polynomials(derivative_term, binomial_term);

            T sign_val = (k % 2 == 0) ? T(1) : T(-1);
            for (auto& c : combined)
                c = c * sign_val;

            hosman_eq = add_polynomials(hosman_eq, combined);
            r_power_k = r_power_k * r_squared;
        }

        // Поиск корней H(x) = 0 в диапазоне [-r, r]
        std::vector<std::pair<T, bool>> candidate_reals;

        // Проверка специальных точек
        {
            T f_zero = eval_poly(hosman_eq, T(0));
            T f_plus_r = eval_poly(hosman_eq, r);
            T f_minus_r = eval_poly(hosman_eq, -r);

            if (abs_val(f_zero) < evaluation_threshold)
                candidate_reals.push_back({ T(0), false });
            if (abs_val(f_plus_r) < evaluation_threshold)
                candidate_reals.push_back({ r, false });
            if (abs_val(f_minus_r) < evaluation_threshold)
                candidate_reals.push_back({ -r, false });
        }

        // Grid sampling с увеличенной плотностью
        {
            int num_samples = 500;
            T grid_threshold = evaluation_threshold * T(50);

            for (int i = 1; i <= num_samples; ++i)
            {
                T x_test = r * T(i) / T(num_samples + 1);
                T val_pos = eval_poly(hosman_eq, x_test);
                T val_neg = eval_poly(hosman_eq, -x_test);

                if (abs_val(val_pos) < grid_threshold)
                    candidate_reals.push_back({ x_test, true });
                if (abs_val(val_neg) < grid_threshold && abs_val(x_test) > epsilon)
                    candidate_reals.push_back({ -x_test, true });
            }
        }

        // Валидация кандидатов
        std::vector<CandidateRoot> valid_candidates;

        for (const auto& [a, from_grid] : candidate_reals)
        {
            T b_squared = r_squared - a * a;
            if (b_squared < -epsilon)
                continue;
            if (b_squared < T(0))
                b_squared = T(0);

            if (b_squared >= epsilon)
            {
                // Комплексная пара a ± bi
                T b = sqrt_val(b_squared);

                std::complex<T> z_plus(a, b);
                std::complex<T> z_minus(a, -b);

                T f_plus = abs_val(eval_complex_poly(coeffs_desc, z_plus));
                T f_minus = abs_val(eval_complex_poly(coeffs_desc, z_minus));
                T eval_thr = evaluation_threshold;

                if (f_plus > eval_thr && f_minus > eval_thr)
                    continue;

                if (f_plus <= eval_thr)
                    valid_candidates.push_back({ std::complex<T>(a, b), T(f_plus) });
                if (f_minus <= eval_thr)
                    valid_candidates.push_back({ std::complex<T>(a, -b), T(f_minus) });
            }
            else
            {
                // Вещественный корень
                T f_a = eval_poly(derivatives[0], a);
                T abs_residual = abs_val(f_a);
                T current_threshold = from_grid ? relaxed_residual_threshold : residual_threshold;

                if (abs_residual < current_threshold)
                    valid_candidates.push_back({ std::complex<T>(a, T(0)), abs_residual });
            }
        }

        // Отбор лучших кандидатов по кратности
        if (valid_candidates.empty())
        {
            for (int i = 0; i < multiplicity; ++i)
            {
                all_final_roots.push_back({ std::complex<T>(r, T(0)), T(1e100) });
            }
            continue;
        }

        std::sort(valid_candidates.begin(), valid_candidates.end(),
            [](const CandidateRoot& a, const CandidateRoot& b)
            {
                return a.residual < b.residual;
            });

        int num_to_take = std::min(multiplicity, static_cast<int>(valid_candidates.size()));
        for (int i = 0; i < num_to_take; ++i)
        {
            all_final_roots.push_back(valid_candidates[i]);
        }

        for (int i = num_to_take; i < multiplicity; ++i)
        {
            all_final_roots.push_back({ std::complex<T>(r, T(0)), T(1e100) });
        }
    }

    std::vector<std::complex<T>> result;
    result.reserve(all_final_roots.size());
    for (const auto& candidate : all_final_roots)
    {
        result.push_back(candidate.root);
    }

    return result;
}