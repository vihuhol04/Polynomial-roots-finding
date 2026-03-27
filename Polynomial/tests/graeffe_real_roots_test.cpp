// Файл для тестирования метода Грефе (вещественные корни). Тестируем улучшенный
// алгоритм с отбросом дублирующихся модулей Реализация: Павлова Анастасия,
// КМБО-01-22

#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>

#include "Graeffe_real_roots.h"
#include "NumericConstants.h"
#include "generate_high_degree_polynomial.h"

using namespace std;

// Основная функция тестирования - генерирует полиномы и тестирует метод Грефе
void run_graeffe_tests() {
    cout << "\n=========================================\n";
    cout << "РАСШИРЕННОЕ ТЕСТИРОВАНИЕ МЕТОДА ГРЕФЕ\n";
    cout << "=========================================\n";

    int test_id = 0;

    // степени
    vector<unsigned> degrees = { 5, 10, 15, 20, 30, 40, 50, 75, 100 };

    // типы сценариев
    enum ScenarioType {
        SIMPLE,
        MULTIPLE,
        CLUSTERED,
        MIXED
    };

    vector<ScenarioType> scenarios = {
        SIMPLE, MULTIPLE, CLUSTERED, MIXED };

    // =========================
    // ТЕСТЫ ДЛЯ double
    // =========================
    cout << "\n===== ТЕСТЫ (double) =====\n";

    for (auto P : degrees) {
        for (auto scenario : scenarios) {
            cout << "\n-----------------------------------------\n";
            cout << "ТЕСТ #" << test_id++
                << " | degree=" << P << " | type=double\n";

            unsigned num_complex_pairs = 0;
            unsigned num_clusters = 0;
            vector<unsigned> cluster_counts;
            vector<double> cluster_radii;
            vector<pair<unsigned, unsigned>> multiplicity_groups;

            switch (scenario) {
            case SIMPLE:
                cout << "Сценарий: простые корни\n";
                break;

            case MULTIPLE:
                cout << "Сценарий: кратные корни\n";
                multiplicity_groups = { {2, P / 4}, {3, P / 6} };
                break;

            case CLUSTERED:
                cout << "Сценарий: кластеризованные корни\n";
                num_clusters = 2;
                cluster_counts = { P / 4, P / 5 };
                cluster_radii = { 1e-3, 1e-2 };
                break;

            case MIXED:
                cout << "Сценарий: смешанный\n";
                num_clusters = 2;
                cluster_counts = { P / 5, P / 6 };
                cluster_radii = { 1e-3, 5e-3 };
                multiplicity_groups = { {2, P / 6} };
                break;
            }

            double default_cluster_radius = 0.1;
            bool normalize_coeffs = true;
            uint64_t seed = 12345 + P * 10 + scenario;

            vector<double> coefficients;
            vector<double> real_roots_repeated;
            vector<double> unique_real_roots;
            vector<unsigned> real_root_multiplicities;
            vector<complex<double>> complex_roots;

            generate_high_degree_polynomial(
                P, num_complex_pairs, num_clusters,
                cluster_counts, cluster_radii,
                multiplicity_groups, default_cluster_radius,
                normalize_coeffs, seed,
                coefficients, real_roots_repeated,
                unique_real_roots, real_root_multiplicities,
                complex_roots);

            double epsilon =
                numeric_constants::adaptive_epsilon<double>(
                    numeric_constants::EPSILON_SCALE_STANDARD);

            int maxIter = (P > 50) ? 10 : 20;

            vector<double> coeffs = coefficients;
            reverse(coeffs.begin(), coeffs.end());

            auto moduli = find_moduli_roots_by_graeffe(
                coeffs, epsilon, maxIter);

            auto moduli_mult =
                find_moduli_with_multiplicities_by_graeffe(
                    coeffs, epsilon, maxIter);

            cout << "Найдено модулей: " << moduli.size() << endl;

            if (!moduli.empty()) {
                cout << "Первые 5: ";
                for (size_t i = 0; i < moduli.size() && i < 5; ++i)
                    cout << moduli[i] << " ";
                cout << endl;
            }

            if (!moduli_mult.empty()) {
                cout << "С кратностями (до 3):\n";
                for (size_t i = 0; i < moduli_mult.size() && i < 3; ++i)
                    cout << "  |r|=" << moduli_mult[i].value
                    << " mult=" << moduli_mult[i].multiplicity << endl;
            }

            cout << "РЕЗУЛЬТАТ: "
                << (!moduli.empty() ? "PASSED" : "FAILED")
                << endl;
        }
    }

    // =========================
    // ТЕСТЫ ДЛЯ float_precision
    // =========================
    cout << "\n===== ТЕСТЫ (float_precision) =====\n";

    for (auto P : degrees) {
        for (auto scenario : scenarios) {
            cout << "\n-----------------------------------------\n";
            cout << "ТЕСТ #" << test_id++
                << " | degree=" << P << " | type=high_precision\n";

            unsigned num_complex_pairs = 0;
            unsigned num_clusters = 0;
            vector<unsigned> cluster_counts;
            vector<float_precision> cluster_radii;
            vector<pair<unsigned, unsigned>> multiplicity_groups;

            switch (scenario) {
            case SIMPLE:
                cout << "Сценарий: простые корни\n";
                break;

            case MULTIPLE:
                cout << "Сценарий: кратные корни\n";
                multiplicity_groups = { {2, P / 4}, {3, P / 6} };
                break;

            case CLUSTERED:
                cout << "Сценарий: кластеризованные корни\n";
                num_clusters = 2;
                cluster_counts = { P / 4, P / 5 };
                cluster_radii = { 1e-4, 1e-3 };
                break;

            case MIXED:
                cout << "Сценарий: смешанный\n";
                num_clusters = 2;
                cluster_counts = { P / 5, P / 6 };
                cluster_radii = { 1e-4, 5e-4 };
                multiplicity_groups = { {2, P / 6} };
                break;
            }

            float_precision default_cluster_radius = 0.1;
            bool normalize_coeffs = true;
            uint64_t seed = 9999 + P * 7 + scenario;

            vector<float_precision> coefficients;
            vector<float_precision> real_roots_repeated;
            vector<float_precision> unique_real_roots;
            vector<unsigned> real_root_multiplicities;
            vector<complex<float_precision>> complex_roots;

            generate_high_degree_polynomial(
                P, num_complex_pairs, num_clusters,
                cluster_counts, cluster_radii,
                multiplicity_groups, default_cluster_radius,
                normalize_coeffs, seed,
                coefficients, real_roots_repeated,
                unique_real_roots, real_root_multiplicities,
                complex_roots);

            float_precision epsilon =
                numeric_constants::adaptive_epsilon<float_precision>(
                    numeric_constants::EPSILON_SCALE_PRECISE);

            int maxIter = (P > 50) ? 8 : 15;

            vector<float_precision> coeffs = coefficients;
            reverse(coeffs.begin(), coeffs.end());

            auto moduli = find_moduli_roots_by_graeffe(
                coeffs, epsilon, maxIter);

            auto moduli_mult =
                find_moduli_with_multiplicities_by_graeffe(
                    coeffs, epsilon, maxIter);

            cout << "Найдено модулей: " << moduli.size() << endl;

            if (!moduli.empty()) {
                cout << "Первые 5: ";
                for (size_t i = 0; i < moduli.size() && i < 5; ++i)
                    cout << moduli[i] << " ";
                cout << endl;
            }

            if (!moduli_mult.empty()) {
                cout << "С кратностями (до 3):\n";
                for (size_t i = 0; i < moduli_mult.size() && i < 3; ++i)
                    cout << "  |r|=" << moduli_mult[i].value
                    << " mult=" << moduli_mult[i].multiplicity << endl;
            }

            cout << "РЕЗУЛЬТАТ: "
                << (!moduli.empty() ? "PASSED" : "FAILED")
                << endl;
        }
    }

    cout << "\n=========================================\n";
    cout << "ВСЕ ТЕСТЫ ЗАВЕРШЕНЫ\n";
    cout << "=========================================\n";
}

int main() {
  // Установка русской локали для корректного вывода
  setlocale(LC_ALL, "ru_RU.UTF-8");

  try {
    cout << "ТЕСТИРОВАНИЕ МЕТОДА ГРЕФЕ" << endl;
    cout << "Генерация полиномов и вычисление модулей корней" << endl << endl;

    run_graeffe_tests();

    cout << "ТЕСТИРОВАНИЕ ЗАВЕРШЕНО " << endl;
    return EXIT_SUCCESS;
  } catch (const exception &e) {
    cerr << "ОШИБКА: " << e.what() << endl;
    return EXIT_FAILURE;
  }
}