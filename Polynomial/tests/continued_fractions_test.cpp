// Тест метода цепных дробей для поиска вещественных корней полиномов
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

#include "Continued_Fractions/Con_Frac.h"

bool approx_contains(const std::vector<double> &roots, double expected, double tol) {
    for (const auto &r : roots) {
        if (std::abs(r - expected) < tol)
            return true;
    }
    return false;
}

int main() {
    setlocale(LC_ALL, "ru_RU.UTF-8");

    int passed = 0, failed = 0;

    // Тест 1: x^3 - 6x^2 + 11x - 6 = 0, корни: 1, 2, 3
    {
        std::vector<double> coeffs = {-6.0, 11.0, -6.0, 1.0};
        auto roots = find_roots_by_Continued_Fractions(coeffs, 1e-6);
        bool ok = roots.size() == 3 &&
                  approx_contains(roots, 1.0, 1e-3) &&
                  approx_contains(roots, 2.0, 1e-3) &&
                  approx_contains(roots, 3.0, 1e-3);
        std::cout << "Тест 1 (x^3-6x^2+11x-6, корни 1,2,3): "
                  << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
        if (!ok) {
            std::cout << "  Найдено " << roots.size() << " корней:";
            for (auto r : roots) std::cout << " " << r;
            std::cout << std::endl;
        }
        ok ? ++passed : ++failed;
    }

    // Тест 2: x^2 - 5x + 6 = 0, корни: 2, 3
    {
        std::vector<double> coeffs = {6.0, -5.0, 1.0};
        auto roots = find_roots_by_Continued_Fractions(coeffs, 1e-6);
        bool ok = roots.size() == 2 &&
                  approx_contains(roots, 2.0, 1e-3) &&
                  approx_contains(roots, 3.0, 1e-3);
        std::cout << "Тест 2 (x^2-5x+6, корни 2,3): "
                  << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
        if (!ok) {
            std::cout << "  Найдено " << roots.size() << " корней:";
            for (auto r : roots) std::cout << " " << r;
            std::cout << std::endl;
        }
        ok ? ++passed : ++failed;
    }

    // Тест 3: x^2 - 1 = 0, корни: -1, 1
    {
        std::vector<double> coeffs = {-1.0, 0.0, 1.0};
        auto roots = find_roots_by_Continued_Fractions(coeffs, 1e-6);
        bool ok = roots.size() == 2 &&
                  approx_contains(roots, -1.0, 1e-3) &&
                  approx_contains(roots, 1.0, 1e-3);
        std::cout << "Тест 3 (x^2-1, корни -1,1): "
                  << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
        if (!ok) {
            std::cout << "  Найдено " << roots.size() << " корней:";
            for (auto r : roots) std::cout << " " << r;
            std::cout << std::endl;
        }
        ok ? ++passed : ++failed;
    }

    // Тест 4: x^4 - 10x^3 + 35x^2 - 50x + 24 = 0, корни: 1, 2, 3, 4
    {
        std::vector<double> coeffs = {24.0, -50.0, 35.0, -10.0, 1.0};
        auto roots = find_roots_by_Continued_Fractions(coeffs, 1e-6);
        bool ok = roots.size() == 4 &&
                  approx_contains(roots, 1.0, 1e-3) &&
                  approx_contains(roots, 2.0, 1e-3) &&
                  approx_contains(roots, 3.0, 1e-3) &&
                  approx_contains(roots, 4.0, 1e-3);
        std::cout << "Тест 4 (x^4-10x^3+35x^2-50x+24, корни 1,2,3,4): "
                  << (ok ? "ПРОЙДЕН" : "ПРОВАЛЕН") << std::endl;
        if (!ok) {
            std::cout << "  Найдено " << roots.size() << " корней:";
            for (auto r : roots) std::cout << " " << r;
            std::cout << std::endl;
        }
        ok ? ++passed : ++failed;
    }

    std::cout << "\nИтого: " << passed << " пройдено, " << failed << " провалено." << std::endl;
    return failed > 0 ? 1 : 0;
}
