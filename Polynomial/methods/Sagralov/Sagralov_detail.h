// Вспомогательные функции для метода Сагралова CIsolate
// Разделение для уменьшения размера основного файла

#pragma once

#include <cmath>
#include <vector>
#include <complex>
#include "mathUtils.h"

template <typename T>
using cplx = std::complex<T>;

namespace detail_sagralov
{
    // Факториал
    inline long double factorial(int k)
    {
        static std::vector<long double> f = { 1.0L };
        while ((int)f.size() <= k)
            f.push_back(f.back() * (long double)f.size());
        return (k >= 0 && k < (int)f.size()) ? f[k] : 1.0L;
    }

    // Вычисление многочлена методом Горнера
    template <typename T>
    cplx<T> eval_poly(const std::vector<cplx<T>>& a, const cplx<T>& z)
    {
        cplx<T> result(0, 0);
        for (const auto& coef : a)
            result = result * z + coef;
        return result;
    }

    // Производная многочлена
    template <typename T>
    std::vector<cplx<T>> derivative(const std::vector<cplx<T>>& a)
    {
        int n = (int)a.size() - 1;
        if (n < 1)
            return { cplx<T>(0, 0) };
        std::vector<cplx<T>> d(n);
        for (int i = 0; i < n; ++i)
            d[i] = a[i] * cplx<T>((T)(n - i), 0);
        return d;
    }

    // k-я производная
    template <typename T>
    std::vector<cplx<T>> derivative_k(const std::vector<cplx<T>>& a, int k)
    {
        std::vector<cplx<T>> cur = a;
        for (int iter = 0; iter < k; ++iter)
            cur = derivative(cur);
        return cur;
    }

    template <typename T>
    std::vector<cplx<T>> build_Fdelta(const std::vector<cplx<T>>& a, const cplx<T>& m, T r)
    {
        int n = (int)a.size() - 1;
        if (n < 0)
            return {};

        std::vector<cplx<T>> c(n + 1);
        std::vector<cplx<T>> curr_poly = a;

        T r_pow = T(1);
        for (int i = 0; i <= n; ++i)
        {
            cplx<T> val = eval_poly(curr_poly, m);
            long double fac_ld = factorial(i);
            T fac = T(fac_ld);

            c[i] = val * (r_pow / fac);

            if (i < n)
            {
                curr_poly = derivative(curr_poly);
                r_pow *= r;
            }
        }

        return c;
    }

    // Итерация Грефа: F^[1](x) = (-1)^n (F_e^2 - x*F_o^2), сжатие до степени n
    template <typename T>
    std::vector<cplx<T>> graeffe_iteration(const std::vector<cplx<T>>& a)
    {
        // Преобразование Греффе: F^[1](x) = (-1)^n * (F_e(x)^2 - x*F_o(x)^2)
        // где F(x) = F_e(x^2) + x*F_o(x^2)
        // Корни нового полинома: квадраты корней исходного
        // Степень сохраняется: n

        int n = (int)a.size() - 1;
        if (n == 0)
            return a;

        // Разделяем на четную и нечетную части
        // a[i] соответствует коэфф при x^(n-i)
        // Fe и Fo имеют размер n+1 (с нулями на местах противоположной четности)
        std::vector<cplx<T>> Fe(n + 1, cplx<T>(0, 0));
        std::vector<cplx<T>> Fo(n + 1, cplx<T>(0, 0));

        for (int i = 0; i <= n; ++i)
        {
            int power = n - i;
            if (power % 2 == 0)
                Fe[i] = a[i];
            else
                Fo[i] = a[i];
        }

        // Возводим в квадрат
        auto sq = [](const std::vector<cplx<T>>& p) -> std::vector<cplx<T>>
            {
                int d = (int)p.size() - 1;
                std::vector<cplx<T>> r(2 * d + 1, cplx<T>(0, 0));
                for (int i = 0; i <= d; ++i)
                    for (int j = 0; j <= d; ++j)
                        r[i + j] += p[i] * p[j];
                return r;
            };

        auto Fe2 = sq(Fe);
        auto Fo2 = sq(Fo);

        // Fe2 - x*Fo2 в полной степени 2n
        // Fo2 после сдвига может дать степень 2n+1, поэтому резервируем 2n+2
        std::vector<cplx<T>> full(2 * n + 2, cplx<T>(0, 0));
        for (size_t i = 0; i < Fe2.size(); ++i)
            full[i] += Fe2[i];
        for (size_t i = 0; i < Fo2.size(); ++i)
            full[i + 1] -= Fo2[i];

        if (n % 2 == 1)
            for (auto& c : full)
                c = -c;

        // Сжатие: берём только чётные степени (т.к. подстановка y=x^2)
        // Результат степени n
        std::vector<cplx<T>> result(n + 1);
        for (int j = 0; j <= n; ++j)
            result[j] = full[2 * j];

        return result;
    }

    // Нормализация старшего коэффициента: 1/4 < |a_n| <= 1
    template <typename T>
    std::vector<cplx<T>> normalize_polynomial(const std::vector<cplx<T>>& a)
    {
        std::vector<cplx<T>> result = a;
        if (result.empty())
            return result;

        T leading = abs_val(result[0]);
        if (is_zero_val(leading))
            return result;

        double leading_dbl = T(leading);
        long double scale = 1.0L;
        while (leading_dbl * scale > 1.0L)
            scale /= 2.0L;
        while (leading_dbl * scale <= 0.25L)
            scale *= 2.0L;

        if (scale != 1.0L)
        {
            T scale_t = T(scale);
            for (auto& c : result)
                c *= scale_t;
        }

        return result;
    }

} // namespace detail_sagralov