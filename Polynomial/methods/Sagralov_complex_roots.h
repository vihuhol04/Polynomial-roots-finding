// Метод Сагралова для изоляции вещественных корней
// Теория: Шишигина Алина, КМБО-01-22
// Реализация: Павлова Анастасия, КМБО-01-22

#pragma once

#include "mathUtils.h"
#include "polynomialUtils.h"
#include "Helper_for_all_methods.h"

#include <queue>
#include <optional>
#include <cassert>

// Приближённая реализация CIsolate (Сагралова) — тест Пелле + Graeffe + Newton/Bisection
// Работает на long double + std::complex<long double>
// Ограничение: предполагает простые корни внутри B и 2B (как в теории).

// Вынесли, чтоб код был более читаемым
template<typename T>
using cplx = std::complex<T>;

template<typename T>
struct Disk {
    cplx<T> center;
    T radius;
    Disk() : center(0, 0), radius(0) {}
    Disk(cplx<T> c, T r) : center(c), radius(r) {}
};

namespace detail {

    // Вычисление значения многочлена с комплексными коэффициентами в точке z (метод Горнера)
    template<typename T>
    cplx<T> eval_poly(const std::vector<cplx<T>>& a, const cplx<T>& z) {
        cplx<T> res = cplx<T>(0, 0);
        for (size_t i = 0; i < a.size(); ++i) {
            res = res * z + a[i];
        }
        return res;
    }

    // Вычисление многочлена с вещественными коэффициентами (переведёнными в комплексные) в точке z
    template<typename T>
    cplx<T> eval_poly_real(const std::vector<T>& a_real, const cplx<T>& z) {
        std::vector<cplx<T>> a;
        a.reserve(a_real.size());
        for (auto v : a_real) a.emplace_back(v, 0);
        return eval_poly(a, z);
    }

    // Вычисление коэффициентов k-й производной многочлена.
    // Вход: коэффициенты a_0 ... a_n, представляющие многочлен a_0 x^n + a_1 x^{n-1} + ... + a_n
    // Возвращает коэффициенты F^{(k)} в том же формате (степень уменьшается на k).
    template<typename T>
    std::vector<cplx<T>> derivative_coeffs(const std::vector<cplx<T>>& a, int k) {
        int n = (int)a.size() - 1;
        if (k == 0) return a;
        if (k > n) return std::vector<cplx<T>>(1, cplx<T>(0, 0));
        std::vector<cplx<T>> cur = a;
        for (int iter = 0; iter < k; ++iter) {
            int m = (int)cur.size() - 1;
            std::vector<cplx<T>> next(m);
            for (int i = 0; i < m; ++i) {
                // cur[i] соответствует коэффициенту при x^{m-i}
                next[i] = cur[i] * cplx<T>((T)(m - i), 0);
            }
            cur.swap(next);
        }
        return cur;
    }

    // Факториал для небольших k
    inline long double fact_ld(int k) {
        static std::vector<long double> f = { 1.0L };
        while ((int)f.size() <= k) f.push_back(f.back() * (long double)f.size());
        return f[k];
    }

    // Построение коэффициентов c_i для F_delta(x) = F(m + r x):
    // c_i = F^{(i)}(m) / i! * r^i
    // Вход: коэффициенты многочлена a (a[0]..a[n], где a[0] при x^n).
    template<typename T>
    std::vector<cplx<T>> build_Fdelta_coeffs(const std::vector<cplx<T>>& a, const cplx<T>& m, T r) {
        int n = (int)a.size() - 1;
        std::vector<cplx<T>> c(n + 1);
        for (int i = 0; i <= n; ++i) {
            // вычисляем i-ю производную в точке m
            std::vector<cplx<T>> der = derivative_coeffs(a, i);
            cplx<T> val = eval_poly(der, m);
            long double denom = fact_ld(i);
            c[i] = val / (denom)*std::pow(r, (long double)i);
        }
        return c;
    }

    // Разбиение квадрата на 4 подквадрата.
    // Квадрат задаётся центром и полуразмером h.
    // Возвращает центры подквадратов и новый полуразмер h/2.
    template<typename T>
    std::vector<std::pair<cplx<T>, T>> subdivide_square(const cplx<T>& center, T halfw) {
        std::vector<std::pair<cplx<T>, T>> out;
        T h2 = halfw / 2;
        // смещения: (-h2, -h2), (h2, -h2), (-h2, h2), (h2, h2)
        out.emplace_back(cplx<T>(center.real() - h2, center.imag() - h2), h2);
        out.emplace_back(cplx<T>(center.real() + h2, center.imag() - h2), h2);
        out.emplace_back(cplx<T>(center.real() - h2, center.imag() + h2), h2);
        out.emplace_back(cplx<T>(center.real() + h2, center.imag() + h2), h2);
        return out;
    }

    // Итерация Греффе для комплексных коэффициентов:
    // Реализует F^{[1]}(x) = (-1)^n (F_e(x)^2 - x F_o(x)^2)
    // Вход: коэффициенты a[0..n] (a[0] при x^n), возвращает преобразованный многочлен
    template<typename T>
    std::vector<cplx<T>> graeffe_one_iter(const std::vector<cplx<T>>& a) {
        int n = (int)a.size() - 1;

        // Строим чётную и нечётную части относительно убывающих степеней
        std::vector<cplx<T>> a_even, a_odd;
        for (int i = 0; i <= n; ++i) {
            if ((i % 2) == 0) a_even.push_back(a[i]);
            else a_odd.push_back(a[i]);
        }

        // Возведение в квадрат (свёртка)
        auto square_poly = [&](const std::vector<cplx<T>>& p) {
            int deg = (int)p.size() - 1;
            int deg_result = 2 * deg;
            std::vector<cplx<T>> res(deg_result + 1, cplx<T>(0, 0));
            for (int i = 0; i <= deg; ++i)
                for (int j = 0; j <= deg; ++j)
                    res[i + j] += p[i] * p[j];
            return res;
            };

        auto sq_even = square_poly(a_even);
        auto sq_odd = square_poly(a_odd);

        // Строим многочлены Fe и Fo в полном формате
        int deg_full = n;
        std::vector<cplx<T>> Fe_full(deg_full + 1, cplx<T>(0, 0));
        std::vector<cplx<T>> Fo_full(deg_full + 1, cplx<T>(0, 0));
        for (int i = 0; i <= n; ++i) {
            int power = n - i;
            if ((power % 2) == 0) Fe_full[i] = a[i];
            else Fo_full[i] = a[i];
        }

        auto poly_square_full = [&](const std::vector<cplx<T>>& p_full) {
            int degp = (int)p_full.size() - 1;
            std::vector<cplx<T>> res(2 * degp + 1, cplx<T>(0, 0));
            for (int i = 0; i <= degp; ++i) {
                for (int j = 0; j <= degp; ++j) {
                    res[i + j] += p_full[i] * p_full[j];
                }
            }
            return res;
            };

        auto Fe2 = poly_square_full(Fe_full);
        auto Fo2 = poly_square_full(Fo_full);

        // Вычисляем Fe2 - x*Fo2
        int deg_res = 2 * deg_full;
        std::vector<cplx<T>> res(deg_res + 1, cplx<T>(0, 0));
        for (int k = 0; k <= deg_res; ++k) {
            cplx<T> fe = (k < (int)Fe2.size() ? Fe2[k] : cplx<T>(0, 0));
            cplx<T> fo_shift = (k >= 1 && (k - 1) < (int)Fo2.size() ? Fo2[k - 1] : cplx<T>(0, 0));
            res[k] = fe - fo_shift;
        }

        // Умножаем на (-1)^n
        if ((n % 2) != 0) {
            for (auto& x : res) x = -x;
        }
        return res;
    }

    // После итерации Греффе степень удваивается. 
    // Для теста Пелле нужны только модули коэффициентов → возвращаем как есть.
} // namespace detail

// Тест T_G: возвращает k ∈ {0..n}, если найдено, иначе -1
// Параметры:
//  - a : коэффициенты многочлена как комплексные числа (a[0] при x^n ... a[n])
//  - m : центр диска
//  - r : радиус
//  - maxGraeffeIter : число итераций (по умолчанию вычисляется из n)
//  - K : пороговый коэффициент (в Пелле K>=1, здесь K=1 по умолчанию)
template<typename T>
int T_G(const std::vector<cplx<T>>& a, const cplx<T>& m, T r, int maxGraeffeIter = -1, T K = 1.0L) {
    int n = (int)a.size() - 1;
    if (n < 1) return 0;
    if (maxGraeffeIter < 0) {
        // N = floor(log2(1 + log n)) + 5  приблизительный размер с натуральными логарифмами
        long double ln = std::log(std::max(1.0, (double)n));
        int N = (int)std::floor(std::log2(1.0L + std::log((long double)n))) + 5;
        if (N < 1) N = 1;
        maxGraeffeIter = N;
    }

    // Строим коэффициенты F_delta c_i = F^{(i)}(m)/i! * r^i
    auto c = detail::build_Fdelta_coeffs(a, m, r); // длина вектора n+1 комплекс

    // Работаем с vector<cplx<T>>
    std::vector<cplx<T>> current = c;

    for (int iter = 0; iter <= maxGraeffeIter; ++iter) {
        // Проверяем состояние Пелле: для каждого k проверить |c_k| > k * sum_{i != k} |c_i|
        std::vector<long double> mags(current.size());
        long double total = 0;
        for (size_t i = 0; i < current.size(); ++i) {
            mags[i] = std::abs(current[i]);
            total += mags[i];
        }
        for (size_t k = 0; k < mags.size(); ++k) {
            long double others = total - mags[k];
            if (mags[k] > (long double)K * others) {
                // Пелле указывает ровно на k корней внутри единичного диска (преобразованного многочлена)
                return (int)k;
            }
        }

        if (iter == maxGraeffeIter) break;
        // выполняем одну итерацию Грефа для текущего (комплексного полинома)
        current = detail::graeffe_one_iter(current);
        // Необязательно нормализовать, чтобы избежать overflow/underflow: масштабировать по максимальной величине
        long double maxm = 0;
        for (auto& z : current) maxm = std::max(maxm, (long double)std::abs(z));
        if (maxm > 0) {
            long double scale = std::pow(10.0L, std::floor(std::log10(maxm)));
            if (std::isfinite(scale) && scale != 0.0L) {
                for (auto& z : current) z /= scale;
            }
        }
    }
    return -1;
}

// Тест T_star: обёртка над T_G с подбором параметров для повышения надёжности
// Возвращает k ≥ 0, если диск содержит ровно k корней, иначе -1
template<typename T>
int T_star(const std::vector<cplx<T>>& a, const cplx<T>& m, T r, long double base_eps = 1e-12L) {
    // Попробуем несколько комбинаций: увеличиваем количество итераций Грефа, немного ослабляем порог K и т.д.
    int try_iters[] = { 1, 3, 5, 8, 12 };
    long double tryKs[] = { 1.0L, 1.0L, 1.2L, 1.5L };
    for (int Ni = 0; Ni < (int)(sizeof(try_iters) / sizeof(int)); ++Ni) {
        for (int Ki = 0; Ki < (int)(sizeof(tryKs) / sizeof(long double)); ++Ki) {
            int k = T_G(a, m, r, try_iters[Ni], (T)tryKs[Ki]);
            if (k >= 0) return k;
        }
    }
    // сдаемся и возвращаем -1
    return -1;
}

// NewtonTest: попытка сузить компоненту с центром centerC и полуразмером h, содержащую k корней
// Возвращает уменьшенный диск, если успешно, иначе nullopt
template<typename T>
std::optional<Disk<T>> NewtonTest(const std::vector<cplx<T>>& a, const Disk<T>& D, int k) {
    // Базовый Newton-like шаг: x' = x - k * F(x)/F'(x)
    cplx<T> xC = D.center;
    cplx<T> Fx = detail::eval_poly(a, xC);
    // вычисляем коэффициенты первой производной:
    std::vector<cplx<T>> deriv = detail::derivative_coeffs(a, 1);
    cplx<T> Fpx = detail::eval_poly(deriv, xC);

    // избегаем маленьких значений производных
    long double denom = std::abs(Fpx);
    if (denom < 1e-18L) return std::nullopt;

    cplx<T> x_new = xC - (T)k * (Fx / Fpx);

    // выбираем новый радиус, уменьшенный в несколько раз
    T r_new = D.radius * (T)0.5L; // сокращение вдвое
    // проверяем T_star на диске с центром в x_new с радиусом r_new
    int k2 = T_star(a, x_new, r_new);
    if (k2 == k) return Disk<T>(x_new, r_new);
    return std::nullopt;
}

// Bisection: делит квадратную область на 4 части
// запускает T_star для каждого описанного диска и возвращает непустые (k > 0)
template<typename T>
std::vector<std::pair<cplx<T>, T>> BisectionComponent(const std::vector<cplx<T>>& a, const cplx<T>& center, T halfw) {
    std::vector<std::pair<cplx<T>, T>> out;
    auto subs = detail::subdivide_square(center, halfw);
    for (auto& p : subs) {
        cplx<T> c = p.first;
        T h = p.second;
        // ограниченный радиус диска r = sqrt(2) * h
        T r = h * (T)std::sqrt(2.0L);
        int k = T_star(a, c, r);
        if (k > 0) {
            out.emplace_back(c, h);
        }
    }
    return out;
}

// Основная функция CIsolate:
// Вход:
//   coeffs : коэффициенты многочлена (вещественные или комплексные).
//   center : центр начального квадрата.
//   halfw : полуразмер квадрата.
//   min_radius : остановка при радиусе ≤ min_radius.
// Возвращает вектор дисков, каждый из которых содержит ровно один простой корень.
template<typename T_in, typename T = long double>
std::vector<Disk<T>> CIsolate(const std::vector<T_in>& coeffs_in,
    cplx<T> center,
    T halfw,
    T min_radius = (T)1e-6,
    int max_iterations = 200) {
    // преобразуем входные коэффициенты в комплексные коэффициенты типа cplx<T>
    std::vector<cplx<T>> a;
    a.reserve(coeffs_in.size());
    for (auto& v : coeffs_in) a.emplace_back((T)v, (T)0);

    int deg = (int)a.size() - 1;
    if (deg < 1) return {};

    // необязательная нормализация: делаем так, чтобы модуль старшего коэффициента был в диапазоне (1/4, 1]
    long double leading = std::abs(a[0]);
    if (leading != 0.0L) {
        // масштабируем многочлен так, чтобы |a0 * s| ∈ (1/4, 1], где s ≈ степень числа 2
        long double scale = 1.0L;
        // при необходимости умножаем на подходящую степень 2
        while (leading * scale > 1.0L) scale /= 2.0L;
        while (leading * scale <= 0.25L) scale *= 2.0L;
        if (scale != 1.0L) {
            for (auto& c : a) c *= (T)scale;
        }
    }

    // начальная компонента — весь квадрат
    struct Component { cplx<T> center; T halfw; int estimated_k; int nC; };
    std::vector<Component> components;
    // вычисляем начальное T_star для описанной окружности квадрата
    T r_init = halfw * (T)std::sqrt(2.0L);
    int k_init = T_star(a, center, r_init);
    if (k_init <= 0) {
        // возможно, в квадрате нет корней — возвращаем пусто
        if (k_init == 0) return {};
        // если результат неопределён: продолжаем деление квадрата на части
    }
    // начинаем с разбиения квадрата на 4 части
    auto first_sub = detail::subdivide_square(center, halfw);
    for (auto& p : first_sub) {
        T r = p.second * (T)std::sqrt(2.0L);
        int k = T_star(a, p.first, r);
        if (k > 0) components.push_back({ p.first, p.second, k, 1 });
    }

    std::vector<Disk<T>> isolating_disks;
    int iter = 0;
    // основной цикл: обрабатываем каждую компоненту, пока не найдены все изолирующие диски или не превышено max_iterations
    while (!components.empty() && iter < max_iterations) {
        ++iter;
        Component C = components.back();
        components.pop_back();

        // вычисляем описанную окружность для текущей компоненты
        cplx<T> mC = C.center;
        T rC = C.halfw * (T)std::sqrt(2.0L);
        int kC = C.estimated_k;
        // быстрый тест: если диск уже достаточно мал и проверен с помощью T_star
        if (C.halfw <= min_radius) {
            // перепроверяем T_star (возможно, на удвоенном диске)
            int k_now = T_star(a, mC, rC);
            if (k_now == kC && kC == 1) {
                isolating_disks.emplace_back(Disk<T>(mC, rC));
                continue;
            }
            // если не прошло проверку — продолжаем уточнение
        }

        if (kC == 1) {
            // пробуем применить NewtonTest для быстрого нахождения изолирующего диска
            auto opt = NewtonTest(a, Disk<T>(mC, rC), kC);
            if (opt) {
                // если радиус полученного диска достаточно мал — добавляем его
                if (opt->radius <= min_radius) {
                    isolating_disks.push_back(*opt);
                }
                else {
                    // иначе добавляем как новую компоненту для дальнейшего уменьшения
                    components.push_back({ opt->center, opt->radius / (T)std::sqrt(2.0L), kC, C.nC + 1 });
                }
                continue;
            }
        }
        else if (kC > 1) {
            // пробуем NewtonTest для кластера корней
            auto optc = NewtonTest(a, Disk<T>(mC, rC), kC);
            if (optc) {
                // если успешно, добавляем уменьшенную компоненту
                components.push_back({ optc->center, optc->radius / (T)std::sqrt(2.0L), kC, C.nC + 1 });
                continue;
            }
        }

        // если NewtonTest не сработал или не применим → делим квадрат (Bisection)
        auto children = BisectionComponent(a, C.center, C.halfw);
        for (auto& ch : children) {
            // вычисляем k для каждой дочерней области
            T rch = ch.second * (T)std::sqrt(2.0L);
            int kk = T_star(a, ch.first, rch);
            if (kk <= 0) continue;
            // если kk == 1 и размер маленький → сразу добавляем изолирующий диск
            if (kk == 1 && ch.second <= min_radius) {
                isolating_disks.emplace_back(Disk<T>(ch.first, rch));
            }
            else {
                components.push_back({ ch.first, ch.second, kk, std::max(1, C.nC - 1) });
            }
        }
    }

    // Постобработка: пробуем уточнить центры изолирующих дисков с помощью итераций Ньютона (для простых корней)
    for (auto& D : isolating_disks) {
        cplx<T> x = D.center;
        for (int it = 0; it < 10; ++it) {
            cplx<T> Fx = detail::eval_poly(a, x);
            cplx<T> Fpx = detail::eval_poly(detail::derivative_coeffs(a, 1), x);
            if (std::abs(Fpx) < 1e-18L) break;
            cplx<T> xn = x - Fx / Fpx;
            if (std::abs(xn - x) < 1e-12L) { x = xn; break; }
            x = xn;
        }
        D.center = x;
    }

    return isolating_disks;
}