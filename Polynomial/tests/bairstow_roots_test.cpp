// Тест метода Берстоу для поиска корней полиномов
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <algorithm>

#include "Bairstow/Bairstow.h"

bool approx_contains(const std::vector<std::complex<double>> &roots,
                     std::complex<double> expected, double tol) {
    for (const auto &r : roots) {
        if (std::abs(r - expected) < tol)
            return true;
    }
    return false;
}

int main() {
    int passed = 0, failed = 0;

    // Тест 1: x^2 - 1 = 0, корни: -1, 1
    {
        std::vector<double> coeffs = {-1.0, 0.0, 1.0};
        auto roots = find_roots_by_Bairstow(coeffs, 1e-10);
        bool ok = roots.size() == 2 &&
                  approx_contains(roots, {1.0, 0.0}, 1e-6) &&
                  approx_contains(roots, {-1.0, 0.0}, 1e-6);
        std::cout << "Тест 1 (x^2-1): " << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
        if (!ok) {
            std::cout << "  Найдено " << roots.size() << " корней:";
            for (auto &r : roots) std::cout << " (" << r.real() << "," << r.imag() << ")";
            std::cout << std::endl;
        }
        ok ? ++passed : ++failed;
    }

    // Тест 2: x^3 - 6x^2 + 11x - 6 = 0, корни: 1, 2, 3
    {
        std::vector<double> coeffs = {-6.0, 11.0, -6.0, 1.0};
        auto roots = find_roots_by_Bairstow(coeffs, 1e-10);
        bool ok = roots.size() == 3 &&
                  approx_contains(roots, {1.0, 0.0}, 1e-4) &&
                  approx_contains(roots, {2.0, 0.0}, 1e-4) &&
                  approx_contains(roots, {3.0, 0.0}, 1e-4);
        std::cout << "Тест 2 (x^3-6x^2+11x-6): " << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
        if (!ok) {
            std::cout << "  Найдено " << roots.size() << " корней:";
            for (auto &r : roots) std::cout << " (" << r.real() << "," << r.imag() << ")";
            std::cout << std::endl;
        }
        ok ? ++passed : ++failed;
    }

    // Тест 3: x^2 + 1 = 0, корни: +i, -i
    {
        std::vector<double> coeffs = {1.0, 0.0, 1.0};
        auto roots = find_roots_by_Bairstow(coeffs, 1e-10);
        bool ok = roots.size() == 2 &&
                  approx_contains(roots, {0.0, 1.0}, 1e-6) &&
                  approx_contains(roots, {0.0, -1.0}, 1e-6);
        std::cout << "Тест 3 (x^2+1): " << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
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
        auto roots = find_roots_by_Bairstow(coeffs, 1e-10);
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

    // Тест 5: x^4 + (2-3i) ... => полином со всеми вещественными коэффициентами
    // x^4 - 10x^3 + 35x^2 - 50x + 24, корни: 1, 2, 3, 4
    {
        std::vector<double> coeffs = {24.0, -50.0, 35.0, -10.0, 1.0};
        auto roots = find_roots_by_Bairstow(coeffs, 1e-10);
        bool ok = roots.size() == 4 &&
                  approx_contains(roots, {1.0, 0.0}, 1e-4) &&
                  approx_contains(roots, {2.0, 0.0}, 1e-4) &&
                  approx_contains(roots, {3.0, 0.0}, 1e-4) &&
                  approx_contains(roots, {4.0, 0.0}, 1e-4);
        std::cout << "Тест 5 (x^4-10x^3+35x^2-50x+24): " << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
        if (!ok) {
            std::cout << "  Найдено " << roots.size() << " корней:";
            for (auto &r : roots) std::cout << " (" << r.real() << "," << r.imag() << ")";
            std::cout << std::endl;
        }
        ok ? ++passed : ++failed;
    }

    std::cout << "\nИтого: " << passed << " пройдено, " << failed << " провалено." << std::endl;
    return failed > 0 ? 1 : 0;
}
