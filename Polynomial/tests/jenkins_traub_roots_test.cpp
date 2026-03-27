// Тест метода Дженкинса-Трауба (CPOLY) для поиска корней полиномов
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "JenkinsTraub/JenkinsTraub.h"

using Complex = std::complex<double>;

bool approx_contains(const std::vector<Complex> &roots,
                     Complex expected, double tol) {
    for (const auto &r : roots) {
        if (std::abs(r - expected) < tol)
            return true;
    }
    return false;
}

int main() {
    setlocale(LC_ALL, "ru_RU.UTF-8");

    int passed = 0, failed = 0;

    // Тест 1: Вещественные коэффициенты: x^3 - 6x^2 + 11x - 6 = 0, корни: 1, 2, 3
    {
        std::vector<double> coeffs = {-6.0, 11.0, -6.0, 1.0};
        auto roots = find_roots_by_JenkinsTraub(coeffs, 1e-10);
        bool ok = roots.size() == 3 &&
                  approx_contains(roots, {1.0, 0.0}, 1e-4) &&
                  approx_contains(roots, {2.0, 0.0}, 1e-4) &&
                  approx_contains(roots, {3.0, 0.0}, 1e-4);
        std::cout << "Тест 1 (вещ. x^3-6x^2+11x-6): " << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
        if (!ok) {
            std::cout << "  Найдено " << roots.size() << " корней:";
            for (auto &r : roots) std::cout << " (" << r.real() << "," << r.imag() << ")";
            std::cout << std::endl;
        }
        ok ? ++passed : ++failed;
    }

    // Тест 2: x^2 + 1 = 0, корни: +i, -i
    {
        std::vector<double> coeffs = {1.0, 0.0, 1.0};
        auto roots = find_roots_by_JenkinsTraub(coeffs, 1e-10);
        bool ok = roots.size() == 2 &&
                  approx_contains(roots, {0.0, 1.0}, 1e-4) &&
                  approx_contains(roots, {0.0, -1.0}, 1e-4);
        std::cout << "Тест 2 (x^2+1): " << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
        if (!ok) {
            std::cout << "  Найдено " << roots.size() << " корней:";
            for (auto &r : roots) std::cout << " (" << r.real() << "," << r.imag() << ")";
            std::cout << std::endl;
        }
        ok ? ++passed : ++failed;
    }

    // Тест 3: Комплексные коэффициенты из C-Poly.pdf:
    // P(z) = z^4 + (2-3i)z^3 + (-1+4i)z^2 + (3+2i)z + (5-i)
    // Эталонные корни:
    //   r1 = -0.779014 - 0.159185i
    //   r2 =  0.352931 + 0.870215i
    //   r3 =  1.198171 - 1.029129i
    //   r4 = -2.772087 + 3.318099i
    {
        std::vector<Complex> coeffs = {
            {5.0, -1.0},   // a_0
            {3.0, 2.0},    // a_1
            {-1.0, 4.0},   // a_2
            {2.0, -3.0},   // a_3
            {1.0, 0.0}     // a_4
        };
        auto roots = find_roots_by_JenkinsTraub(coeffs, 1e-10);
        bool ok = roots.size() == 4 &&
                  approx_contains(roots, {-0.779014, -0.159185}, 1e-2) &&
                  approx_contains(roots, {0.352931, 0.870215}, 1e-2) &&
                  approx_contains(roots, {1.198171, -1.029129}, 1e-2) &&
                  approx_contains(roots, {-2.772087, 3.318099}, 1e-2);
        std::cout << "Тест 3 (компл. полином из C-Poly.pdf): " << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
        if (!ok) {
            std::cout << "  Найдено " << roots.size() << " корней:";
            for (auto &r : roots) std::cout << " (" << r.real() << "," << r.imag() << ")";
            std::cout << std::endl;
        }
        ok ? ++passed : ++failed;
    }

    // Тест 4: x^4 - 1 = 0, корни: 1, -1, i, -i
    {
        std::vector<double> coeffs = {-1.0, 0.0, 0.0, 0.0, 1.0};
        auto roots = find_roots_by_JenkinsTraub(coeffs, 1e-10);
        bool ok = roots.size() == 4 &&
                  approx_contains(roots, {1.0, 0.0}, 1e-4) &&
                  approx_contains(roots, {-1.0, 0.0}, 1e-4) &&
                  approx_contains(roots, {0.0, 1.0}, 1e-4) &&
                  approx_contains(roots, {0.0, -1.0}, 1e-4);
        std::cout << "Тест 4 (x^4-1): " << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
        if (!ok) {
            std::cout << "  Найдено " << roots.size() << " корней:";
            for (auto &r : roots) std::cout << " (" << r.real() << "," << r.imag() << ")";
            std::cout << std::endl;
        }
        ok ? ++passed : ++failed;
    }

    // Тест 5: Линейный полином x - 5 = 0, корень: 5
    {
        std::vector<double> coeffs = {-5.0, 1.0};
        auto roots = find_roots_by_JenkinsTraub(coeffs, 1e-10);
        bool ok = roots.size() == 1 &&
                  approx_contains(roots, {5.0, 0.0}, 1e-8);
        std::cout << "Тест 5 (x-5): " << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
        ok ? ++passed : ++failed;
    }

    std::cout << "\nИтого: " << passed << " пройдено, " << failed << " провалено." << std::endl;
    return failed > 0 ? 1 : 0;
}
