// Алгоритм Дюрана-Кернера для поиска комплексных корней
// Теория + полная реализация: Семенидо Алина, КМБО-06-22

#ifndef ROOT_FINDER_H
#define ROOT_FINDER_H

#include <vector>
#include <complex>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <cmath>
#include <algorithm>
#include <stdexcept>


namespace detail {
    constexpr std::size_t DEFAULT_ITERATIONS_MULTIPLIER = 64;
    constexpr std::size_t DEFAULT_NEWTON_PASSES = 2;
    constexpr std::size_t MAX_SKIP_THRESHOLD = 10;
}

template <typename Complex>
using real_of_t = std::decay_t<decltype(std::abs(Complex{})) > ;

namespace detail {

    // ---------- вспомогательные решения низких степеней ----------
    template <typename Complex>
    std::vector<Complex> solve_linear(const std::vector<Complex>& c) {
        // c = {a0, 1}, т.е. x + a0 = 0
        return { -c[0] };
    }

    template <typename Complex>
    std::vector<Complex> solve_quadratic(const std::vector<Complex>& c) {
        // c = {a0, a1, 1}  =>  x² + a1 x + a0 = 0
        const Complex& a0 = c[0];
        const Complex& a1 = c[1];

        Complex d = sqrt(a1 * a1 - Complex(4) * a0);                // дискриминант
        Complex q = (std::real(a1) >= 0) ? (-a1 - d) : (-a1 + d);    // численно устойчиво
        Complex r1 = q / Complex(2);
        Complex r2 = a0 / r1;

        return { r1, r2 };
    }

    // ---------- утилиты алгоритма Дюрана-Кернера ----------
    template <typename Complex, typename Real>
    Real bounding_radius(const std::vector<Complex>& c) {
        Real R = Real(1);
        for (std::size_t i = 0, n = c.size() - 1; i < n; ++i)
            R = std::max(R, Real(1) + std::abs(c[i]));
        if (!std::isfinite(R))
            throw std::overflow_error("Переполнение при оценке радиуса");
        return R;
    }

    template <typename Complex, typename Real>
    std::vector<Complex> initial_roots(std::size_t n, Real R) {
        const Real PI = std::acos(Real(-1));
        std::vector<Complex> z(n);
        for (std::size_t k = 0; k < n; ++k) {
            Real ang = Real(2) * PI * Real(k) / Real(n);
            z[k] = Complex(std::cos(ang), std::sin(ang)) * R;
        }
        return z;
    }


    template <typename Complex, typename Real>
    bool durand_kerner_step(std::vector<Complex>& z,
        const std::vector<Complex>& c,
        Real eps) {
        bool ok = true;
        std::size_t skips = 0;
        const std::size_t n = z.size();

        for (std::size_t i = 0; i < n; ++i) {
            Complex f = eval_poly(c, z[i]);

            Complex denom = 1;
            for (std::size_t j = 0; j < n; ++j)
                if (j != i) denom *= (z[i] - z[j]);

            if (std::abs(denom) < eps) {
                if (++skips > detail::MAX_SKIP_THRESHOLD)
                    throw std::runtime_error("Нестабильность Дюрана-Кернера");
                continue;
            }

            Complex corr = f / denom;
            z[i] -= corr;
            if (std::abs(corr) / std::max<Real>(1, std::abs(z[i])) > eps)
                ok = false;
        }
        return ok;
    }

    template <typename Complex, typename Real>
    void newton_pass(std::vector<Complex>& z,
        const std::vector<Complex>& c,
        Real eps) {
        for (Complex& root : z) {
            auto [f, df] = eval_poly_and_deriv(c, root);
            if (std::abs(df) >= eps)
                root -= f / df;
        }
    }

} // namespace detail

// ------------------ ПУБЛИЧНАЯ ФУНКЦИЯ ------------------
template <typename Complex, typename Eps = real_of_t<Complex>>
std::vector<Complex> find_roots_by_Durand_Kerner(const std::vector<Complex>& coeffs,
    Eps eps = std::numeric_limits<Eps>::epsilon(),
    std::size_t max_iter = 0,
    std::size_t newton = detail::DEFAULT_NEWTON_PASSES)
{
    using Real = Eps;

    const std::size_t m = coeffs.size();
    if (m < 2) return {};

    const Complex leading = coeffs.front();           // aₙ
    if (std::abs(leading) < eps)
        throw std::invalid_argument("Старший коэффициент равен нулю");

    // Переворачиваем {aₙ…a₀} → {a₀…aₙ} и нормируем
    std::vector<Complex> c(coeffs.rbegin(), coeffs.rend());
    for (Complex& a : c) a /= leading;                // моническая форма

    const std::size_t n = m - 1;

    if (n == 1) return detail::solve_linear(c);
    if (n == 2) return detail::solve_quadratic(c);

    if (max_iter == 0) max_iter = n * detail::DEFAULT_ITERATIONS_MULTIPLIER;

    Real R = detail::bounding_radius<Complex, Real>(c);
    std::vector<Complex> z = detail::initial_roots<Complex, Real>(n, R);

    for (std::size_t k = 0; k < max_iter; ++k)
        if (detail::durand_kerner_step(z, c, eps)) break;

    for (std::size_t p = 0; p < newton; ++p)
        detail::newton_pass(z, c, eps);

    return z;
}


#endif // ROOT_FINDER_H
