// Метод Сагралова CIsolate для изоляции комплексных корней
// Реализация согласно теории (раздел 2.2)

#pragma once

#include <cmath>
#include <vector>
#include <complex>
#include <optional>
#include <algorithm>

#include "Sagralov_detail.h"
#include "mathUtils.h"

template <typename T>
using cplx = std::complex<T>;

template <typename T>
struct Disk
{
    cplx<T> center;
    T radius;
    int num_roots;

    Disk() : center(0, 0), radius(0), num_roots(0) {}
    Disk(cplx<T> c, T r, int k = 0) : center(c), radius(r), num_roots(k) {}
};

template <typename T>
struct Component
{
    std::vector<cplx<T>> squares;
    T halfw;
    Disk<T> disk;
    int N_C; // 2^(2*n_C)
    int n_C;

    Component() : halfw(0), N_C(4), n_C(1) {}
};

// Тест Пелле после N итераций Грефа
// Согласно теории: проверяем условие НА КАЖДОЙ итерации (п.2 алгоритма TG)
template <typename T>
int T_G(const std::vector<cplx<T>>& a, const cplx<T>& m, T r)
{
    int n = (int)a.size() - 1;
    if (n < 1)
        return 0;

    // Вычисление N = ⌈log(1+log n+5)⌉
    T n_val = T(max_val(2, n));
    T log_n = log_val(n_val);
    T inner = T(1) + log_n + T(5);
    T log_inner = log_val(inner);
    T log2_inner = log_inner / log_val(T(2));
    int N = (int)ceil_val(log2_inner);
    N = max_val(1, min_val(N, 12));

    auto c = detail_sagralov::build_Fdelta(a, m, r);

    // Выполняем N итераций Грефа, проверяя условие Пелле после каждой
    for (int iter = 0; iter <= N; ++iter)
    {
        std::vector<T> mag(c.size());
        T sum = T(0);
        for (size_t i = 0; i < c.size(); ++i)
        {
            mag[i] = abs_val(c[i]);
            sum += mag[i];
        }

        // Проверка условия Пелле: |c_k| > sum_{i != k} |c_i|
        // Тест Пелле: |c_k| > K * sum_{i≠k}|c_i|, где K≥1
        // Используем K=1.01 для надёжности при численных ошибках
        const T K = T(1.01);
        for (size_t k = 0; k < mag.size(); ++k)
            if (mag[k] > K * (sum - mag[k]))
                return (int)k;

        // Следующая итерация Грефа
        if (iter < N)
        {
            c = detail_sagralov::graeffe_iteration(c);

            // Масштабирование степенями двойки
            // Находим максимальный коэффициент и нормализуем все коэффициенты
            T max_mag = T(0);
            for (const auto& z : c)
            {
                max_mag = max_val(max_mag, abs_val(z));
            }

            if (max_mag > T(0) && isfinite_val(max_mag))
            {
                // Делим все коэффициенты на max_mag для нормализации
                // Приводим максимальный к диапазону [0.5, 1) делением на степени 2
                T scale = T(1);
                T normalized = max_mag;

                // Приводим к диапазону [0.5, 1) умножением/делением на 2
                while (normalized > T(1))
                {
                    normalized /= T(2);
                    scale *= T(2);
                }
                while (normalized <= T(0.5) && normalized > T(0))
                {
                    normalized *= T(2);
                    scale /= T(2);
                }

                // Масштабируем все коэффициенты
                if (scale != T(1))
                {
                    for (auto& z : c)
                    {
                        z /= scale;
                    }
                }
            }
        }
    }

    return -1;
}

// Адаптивный тест T* с переменным числом итераций Греффе
// При неудаче базового T_G постепенно увеличиваем N
template <typename T>
int T_star(const std::vector<cplx<T>>& a, const cplx<T>& m, T r)
{
    int n = (int)a.size() - 1;
    if (n < 1)
        return 0;

    // Базовое значение N согласно теории
    int N_base = (int)std::ceil(std::log2(1.0L + std::log((long double)max_val(2, n)) + 5.0L));
    N_base = max_val(1, N_base);

    // Пытаемся с возрастающим N (имитация переменной точности)
    for (int N = N_base; N <= min_val(N_base + 6, 18); ++N)
    {
        auto c = detail_sagralov::build_Fdelta(a, m, r);

        // Выполняем N итераций Греффе с проверкой Пелле после каждой
        for (int iter = 0; iter <= N; ++iter)
        {
            std::vector<T> mag(c.size());
            T sum = T(0);
            for (size_t i = 0; i < c.size(); ++i)
            {
                mag[i] = abs_val(c[i]);
                sum += mag[i];
            }

            // Тест Пелле: |c_k| > K * sum_{i≠k}|c_i|
            const T K = T(1.01);
            for (size_t k = 0; k < mag.size(); ++k)
                if (mag[k] > K * (sum - mag[k]))
                    return (int)k;

            // Следующая итерация Греффе с нормализацией
            if (iter < N)
            {
                c = detail_sagralov::graeffe_iteration(c);

                // Численно стабильная нормализация степенями двойки
                T max_mag = T(0);
                for (const auto& z : c)
                    max_mag = max_val(max_mag, abs_val(z));

                if (max_mag > T(0) && isfinite_val(max_mag))
                {
                    T scale_factor = T(1) / max_mag;
                    for (auto& z : c)
                        z *= scale_factor;
                }
            }
        }
    }

    // Все попытки исчерпаны - возвращаем -1 (неопределённость)
    return -1;
}

// Метод Ньютона-Шрёдера для кластеров
template <typename T>
std::optional<Disk<T>> NewtonTest(const std::vector<cplx<T>>& a, const Disk<T>& D, int N_C)
{
    if (D.num_roots <= 0)
        return std::nullopt;

    int k = D.num_roots;
    cplx<T> x_C = D.center;
    T r = D.radius;

    cplx<T> F_val = detail_sagralov::eval_poly(a, x_C);
    auto deriv = detail_sagralov::derivative(a);
    cplx<T> Fp_val = detail_sagralov::eval_poly(deriv, x_C);

    if (abs_val(Fp_val) < T(1e-18))
        return std::nullopt;

    cplx<T> x_new = x_C - cplx<T>(T(k), T(0)) * (F_val / Fp_val);
    T r_new = r / T(2 * N_C);

    int k_new = T_star(a, x_new, r_new);
    if (k_new == k)
        return Disk<T>(x_new, r_new, k);

    return std::nullopt;
}

// Бисекция квадрата на 4 подквадрата
template <typename T>
std::vector<std::pair<cplx<T>, T>> BisectionSquare(const cplx<T>& center, T halfw)
{
    T h2 = halfw / static_cast<T>(2);
    return {
        {cplx<T>(center.real() - h2, center.imag() - h2), h2},
        {cplx<T>(center.real() + h2, center.imag() - h2), h2},
        {cplx<T>(center.real() - h2, center.imag() + h2), h2},
        {cplx<T>(center.real() + h2, center.imag() + h2), h2} };
}

// CIsolate - изоляция комплексных корней с ненулевыми Re и Im частями
// ВАЖНО: метод НЕ предназначен для корней на вещественной или мнимой осях!
// Для вещественных корней используйте Sagralov_real_roots
// Для чисто мнимых корней рекомендуется преобразование P(ix) и применение метода вещественных корней
template <typename T_in, typename T = long double>
std::vector<Disk<T>> CIsolate(
    const std::vector<T_in>& coeffs_in,
    cplx<T> center = cplx<T>(0, 0),
    T halfw = static_cast<T>(0.0),
    T min_radius = static_cast<T>(1e-3),
    int max_iterations = 500)
{
    // Преобразуем в комплексные коэффициенты
    std::vector<cplx<T>> a;
    a.reserve(coeffs_in.size());
    for (const auto& v : coeffs_in)
    {
        a.emplace_back(static_cast<T>(v), static_cast<T>(0));
    }

    int n = (int)a.size() - 1;
    if (n < 1)
        return {};

    // Вычисляем границу Коши ДО масштабирования
    T max_coef = T(0);
    for (size_t i = 0; i < a.size() - 1; ++i)
        max_coef = max_val(max_coef, abs_val(a[i]));
    T leading = abs_val(a.back());
    if (leading < T(1e-10))
        return {}; // Вырожденный полином

    T cauchy_bound = T(1) + max_coef / leading;

    // Проблема: CIsolate плохо работает для |z| < 1 (underflow в итерациях Грефе)
    //
    // Решение: Трансформация P(x) → Q(y) = P(y/k)
    // - Если z корень P, то w = k·z корень Q
    // - Модули умножаются: |w| = k·|z|
    // - После нахождения дисков делим на k
    //
    // Выбор k: целевая граница Коши ≈ 2.5-3.0 (стабильная область для итераций Грефе)
    T scale_k = T(1);

    if (cauchy_bound < T(2.0))
    {
        // Малые/средние корни: увеличиваем для стабильности
        scale_k = T(2.5) / cauchy_bound;

        // Q(y) = P(y/k): a_i → a_i / k^i
        T k_pow = T(1);
        T inv_k = T(1) / scale_k;
        for (size_t i = 0; i < a.size(); ++i)
        {
            a[i] *= k_pow;
            k_pow *= inv_k;
        }
    }
    else if (cauchy_bound > T(10.0))
    {
        // Очень большие корни: уменьшаем для предотвращения overflow
        scale_k = T(3.0) / cauchy_bound;

        T k_pow = T(1);
        T inv_k = T(1) / scale_k;
        for (size_t i = 0; i < a.size(); ++i)
        {
            a[i] *= k_pow;
            k_pow *= inv_k;
        }
    }

    // Нормализация: 1/4 < |a_n| <= 1 (требование теории)
    a = detail_sagralov::normalize_polynomial(a);

    // Пересчитываем границу Коши ПОСЛЕ масштабирования
    if (halfw <= T(0))
    {
        max_coef = T(0);
        for (size_t i = 0; i < a.size() - 1; ++i)
            max_coef = max_val(max_coef, abs_val(a[i]));
        leading = abs_val(a.back());

        T new_cauchy_bound = T(1) + max_coef / leading;

        // Увеличиваем область для надёжности
        halfw = new_cauchy_bound * T(1.8);
    }

    std::vector<Component<T>> work_queue;

    auto process_square = [&](const cplx<T>& sq_center, T sq_half) -> std::optional<Component<T>>
        {
            T radius = sq_half * sqrt_val(static_cast<T>(2)) * static_cast<T>(1.15);
            int k = T_star(a, sq_center, radius);

            if (k > 0)
            {
                Component<T> comp;
                comp.squares.push_back(sq_center);
                comp.halfw = sq_half;
                comp.disk = Disk<T>(sq_center, radius, k);
                comp.N_C = 4;
                comp.n_C = 1;
                return comp;
            }
            return std::nullopt;
        };

    // Начальное деление: делим область на 4 квадрата
    auto initial_subs = BisectionSquare(center, halfw);
    for (const auto& [sq_center, sq_half] : initial_subs)
    {
        auto comp_opt = process_square(sq_center, sq_half);
        if (comp_opt)
            work_queue.push_back(*comp_opt);
    }

    // Если ни один начальный квадрат не содержит корней, делим мельче
    if (work_queue.empty() && n > 0)
    {
        // Попробуем деление на 16 подквадратов (двойная бисекция)
        T half2 = halfw / T(2);
        for (const auto& [c1, h1] : initial_subs)
        {
            auto sub2 = BisectionSquare(c1, h1);
            for (const auto& [sq_center, sq_half] : sub2)
            {
                auto comp_opt = process_square(sq_center, sq_half);
                if (comp_opt)
                    work_queue.push_back(*comp_opt);
            }
        }
    }

    std::vector<Disk<T>> isolated_roots;
    int iteration = 0;

    while (!work_queue.empty() && iteration < max_iterations)
    {
        ++iteration;

        Component<T> current = work_queue.back();
        work_queue.pop_back();

        if (current.disk.radius <= min_radius)
        {
            isolated_roots.push_back(current.disk);
            continue;
        }

        // Согласно теории: проверяем T*(Δ), если k≥0 добавляем в isolated_roots
        int k_current = T_star(a, current.disk.center, current.disk.radius);

        if (k_current <= 0)
        {
            // Если T_star не подтверждает корни, пропускаем
            continue;
        }

        if (k_current == 1)
        {
            // Один корень - сразу изолирован
            isolated_roots.push_back(current.disk);
            continue;
        }

        // k_current > 1: несколько корней, попробуем NewtonTest
        auto newton_disk = NewtonTest(a, current.disk, current.N_C);
        if (newton_disk && newton_disk->radius < current.disk.radius * T(0.8))
        {
            // Успешное уточнение
            Component<T> refined;
            refined.squares.push_back(newton_disk->center);
            refined.halfw = newton_disk->radius / sqrt_val(T(2));
            refined.disk = *newton_disk;
            refined.n_C = current.n_C + 1;
            refined.N_C = (1 << (2 * refined.n_C));
            work_queue.push_back(refined);
            continue;
        }

        // NewtonTest не помог, применяем бисекцию по теории

        std::vector<std::pair<cplx<T>, T>> subsquares;
        for (const auto& sq : current.squares)
        {
            auto subs = BisectionSquare(sq, current.halfw);
            subsquares.insert(subsquares.end(), subs.begin(), subs.end());
        }

        struct SquareData
        {
            cplx<T> center;
            T halfw;
            T radius;
            int roots;
        };
        std::vector<SquareData> nonempty_squares;

        for (const auto& [sq_center, sq_half] : subsquares)
        {
            T radius = sq_half * sqrt_val(T(2)) * T(1.15);
            int k = T_star(a, sq_center, radius);
            if (k > 0)
                nonempty_squares.push_back({ sq_center, sq_half, radius, k });
        }

        if (nonempty_squares.empty())
            continue;

        std::vector<bool> processed(nonempty_squares.size(), false);

        for (size_t i = 0; i < nonempty_squares.size(); ++i)
        {
            if (processed[i])
                continue;

            std::vector<cplx<T>> component_centers;
            component_centers.push_back(nonempty_squares[i].center);
            int total_roots = nonempty_squares[i].roots;
            T comp_halfw = nonempty_squares[i].halfw;
            processed[i] = true;

            for (size_t j = i + 1; j < nonempty_squares.size(); ++j)
            {
                if (processed[j])
                    continue;

                T dist = abs_val(nonempty_squares[i].center - nonempty_squares[j].center);
                T threshold = comp_halfw * T(2.2);

                if (dist <= threshold)
                {
                    component_centers.push_back(nonempty_squares[j].center);
                    total_roots += nonempty_squares[j].roots;
                    processed[j] = true;
                }
            }

            cplx<T> centroid(0, 0);
            for (const auto& c : component_centers)
                centroid += c;
            centroid /= cplx<T>(T(component_centers.size()), T(0));

            T max_distance = T(0);
            for (const auto& c : component_centers)
                max_distance = max_val(max_distance, abs_val(c - centroid));

            T component_radius = max_distance + comp_halfw * sqrt_val(T(2)) * T(1.15);

            Component<T> new_component;
            new_component.squares = component_centers;
            new_component.halfw = comp_halfw;
            new_component.disk = Disk<T>(centroid, component_radius, total_roots);
            new_component.n_C = max_val(1, current.n_C - 1);
            new_component.N_C = (1 << (2 * new_component.n_C));

            work_queue.push_back(new_component);
        }
    }

    // Если достигнут лимит итераций, добавляем оставшиеся компоненты как результаты
    // с увеличенным радиусом для покрытия возможных корней
    if (iteration >= max_iterations && !work_queue.empty())
    {
        for (const auto& comp : work_queue)
        {
            if (comp.disk.num_roots > 0)
            {
                // Добавляем диск с увеличенным радиусом для надежности
                T safe_radius = comp.disk.radius * T(2);
                isolated_roots.emplace_back(comp.disk.center, safe_radius, comp.disk.num_roots);
            }
        }
    }

    // Обратная трансформация: диски для Q(y) → диски для P(x)
    // Если w = k·z, то z = w/k
    // Делим центры и радиусы на k
    if (scale_k != T(1))
    {
        for (auto& disk : isolated_roots)
        {
            disk.center /= scale_k;
            disk.radius /= scale_k;
        }
    }

    return isolated_roots;
}