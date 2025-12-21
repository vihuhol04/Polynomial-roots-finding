// Файл для тестирования модификации Хосмана метода Грефе
// Реализация: Павлова Анастасия, КМБО-01-22

#include "mathUtils.h"
#include "polynomialUtils.h"
#include "generate_high_degree_polynomial.h"
#include "Graeffe_complex_roots.h"

using namespace std;

// Основная функция тестирования - генерирует полиномы и тестирует метод Хосмана
void run_hosman_tests() {
    // ТЕСТ 1
    {
        cout << "ТЕСТ 1: Полином с 2 сопряженными парами" << endl;

        unsigned P = 4;
        unsigned num_complex_pairs = 2;
        unsigned num_clusters = 0;
        vector<unsigned> cluster_counts = {};
        vector<double> cluster_radii = {};
        vector<pair<unsigned, unsigned>> multiplicity_groups = {};
        double default_cluster_radius = 0.1;
        bool normalize_coeffs = true;
        uint64_t seed = 12345;

        vector<double> coefficients;
        vector<double> real_roots_repeated;
        vector<double> unique_real_roots;
        vector<unsigned> real_root_multiplicities;
        vector<complex<double>> complex_roots;

        generate_high_degree_polynomial(
            P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
            multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
            coefficients, real_roots_repeated, unique_real_roots,
            real_root_multiplicities, complex_roots
        );

        cout << "СГЕНЕРИРОВАННЫЙ ПОЛИНОМ:" << endl;
        cout << "Степень: " << P << endl;
        cout << "Полином: ";
        print_polynomial(coefficients);
        cout << endl;

        reverse(coefficients.begin(), coefficients.end());
        vector<double> coefficients_1 = coefficients;

        cout << "ИСХОДНЫЕ КОРНИ:" << endl;
        cout << "Вещественные корни: ";
        for (auto root : real_roots_repeated) {
            cout << root << " ";
        }
        cout << endl;

        cout << "Комплексные корни (сопряженные пары):" << endl;
        if (complex_roots.empty()) {
            cout << "  нет" << endl;
        } else {
            for (size_t i = 0; i < complex_roots.size(); i += 2) {
                cout << "  " << complex_roots[i] << " и " << complex_roots[i + 1] << endl;
            }
        }
        cout << endl;

        cout << "ВЫЧИСЛЕНИЕ МЕТОДОМ ХОСМАНА (базовая версия):" << endl;
        double epsilon = 1e-6;

        // Вычисление корней методом Хосмана (базовая версия)
        vector<complex<double>> hosman_roots = hosman_modification_graeffe(coefficients, epsilon);
        cout << "Найдено корней: " << hosman_roots.size() << endl;
        cout << "Вычисленные корни:" << endl;
        for (size_t i = 0; i < hosman_roots.size(); ++i) {
            cout << "  " << i + 1 << ": " << hosman_roots[i] << endl;
        }
        cout << endl;

        // ТЕСТИРОВАНИЕ МЕТОДА HOSMAN_WITH_MULTIPLICITIES
        cout << "ВЫЧИСЛЕНИЕ МЕТОДОМ ХОСМАНА С КРАТНОСТЯМИ:" << endl;

        auto found_graeffe = find_moduli_with_multiplicities_by_graeffe(coefficients_1, epsilon);
        std::cout << "Корни вещественные (Греф):\n";
        for (const auto& root : found_graeffe) {
            std::cout << "root: " << root.value << "    mult: " << root.multiplicity << "\n";
        }
        std::cout << "\n";

        auto found_graeffe_complex = hosman_with_multiplicities(coefficients_1, found_graeffe, epsilon);
        std::cout << "Корни комплексные (Хосман):\n";
        for (auto& r : found_graeffe_complex) {
            std::cout << r.real() << " " << r.imag() << "i\n";
        }
    }

    // ТЕСТ 2
    {
        cout << "ТЕСТ 2: Полином с 10 сопряженными парами (precision source)" << endl;

        unsigned P = 20;
        unsigned num_complex_pairs = 10;
        unsigned num_clusters = 0;
        vector<unsigned> cluster_counts = {};
        vector<float_precision> cluster_radii = {};
        vector<pair<unsigned, unsigned>> multiplicity_groups = {};
        float_precision default_cluster_radius = 0.1;
        bool normalize_coeffs = true;
        uint64_t seed = 951;

        vector<float_precision> coefficients;
        vector<float_precision> real_roots_repeated;
        vector<float_precision> unique_real_roots;
        vector<unsigned> real_root_multiplicities;
        vector<complex<float_precision>> complex_roots;

        generate_high_degree_polynomial(
            P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
            multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
            coefficients, real_roots_repeated, unique_real_roots,
            real_root_multiplicities, complex_roots
        );

        cout << "СГЕНЕРИРОВАННЫЙ ПОЛИНОМ:" << endl;
        cout << "Степень: " << P << endl;
        cout << "Полином: ";
        print_polynomial(coefficients);
        cout << endl;

        reverse(coefficients.begin(), coefficients.end());
        vector<float_precision> coefficients_1 = coefficients;

        cout << "ИСХОДНЫЕ КОРНИ:" << endl;
        cout << "Вещественные корни: ";
        for (auto root : real_roots_repeated) {
            cout << root << " ";
        }
        cout << endl;

        cout << "Комплексные корни (сопряженные пары):" << endl;
        if (complex_roots.empty()) {
            cout << "  нет" << endl;
        }
        else {
            for (size_t i = 0; i < complex_roots.size(); i += 2) {
                cout << "  " << complex_roots[i] << " и " << complex_roots[i + 1] << endl;
            }
        }
        cout << endl;

        cout << "ВЫЧИСЛЕНИЕ МЕТОДОМ ХОСМАНА (базовая версия):" << endl;
        float_precision epsilon = 1e-12;

        // Вычисление корней методом Хосмана (базовая версия)
        vector<complex<float_precision>> hosman_roots = hosman_modification_graeffe(coefficients, epsilon);
        cout << "Найдено корней: " << hosman_roots.size() << endl;
        cout << "Вычисленные корни:" << endl;
        for (size_t i = 0; i < hosman_roots.size(); ++i) {
            cout << "  " << i + 1 << ": " << hosman_roots[i] << endl;
        }
        cout << endl;

        // ТЕСТИРОВАНИЕ МЕТОДА HOSMAN_WITH_MULTIPLICITIES
        cout << "ВЫЧИСЛЕНИЕ МЕТОДОМ ХОСМАНА С КРАТНОСТЯМИ:" << endl;

        auto found_graeffe = find_moduli_with_multiplicities_by_graeffe(coefficients_1, epsilon);
        std::cout << "Корни вещественные (Греф):\n";
        for (const auto& root : found_graeffe) {
            std::cout << "root: " << root.value << "    mult: " << root.multiplicity << "\n";
        }
        std::cout << "\n";

        auto found_graeffe_complex = hosman_with_multiplicities(coefficients_1, found_graeffe, epsilon);
        std::cout << "Корни комплексные (Хосман):\n";
        for (auto& r : found_graeffe_complex) {
            std::cout << r.real() << " " << r.imag() << "i\n";
        }
    }

    // ТЕСТ 3
    {
        cout << "ТЕСТ 3: Полином со смешанными корнями" << endl;

        unsigned P = 10;
        unsigned num_complex_pairs = 3;
        unsigned num_clusters = 0;
        vector<unsigned> cluster_counts = {};
        vector<double> cluster_radii = {};
        vector<pair<unsigned, unsigned>> multiplicity_groups = {};
        double default_cluster_radius = 0.1;
        bool normalize_coeffs = true;
        uint64_t seed = 5291;

        vector<double> coefficients;
        vector<double> real_roots_repeated;
        vector<double> unique_real_roots;
        vector<unsigned> real_root_multiplicities;
        vector<complex<double>> complex_roots;

        generate_high_degree_polynomial(
            P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
            multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
            coefficients, real_roots_repeated, unique_real_roots,
            real_root_multiplicities, complex_roots
        );

        cout << "СГЕНЕРИРОВАННЫЙ ПОЛИНОМ:" << endl;
        cout << "Степень: " << P << endl;
        cout << "Полином: ";
        print_polynomial(coefficients);
        cout << endl;

        reverse(coefficients.begin(), coefficients.end());
        vector<double> coefficients_1 = coefficients;

        cout << "ИСХОДНЫЕ КОРНИ:" << endl;
        cout << "Вещественные корни: ";
        for (auto root : real_roots_repeated) {
            cout << root << " ";
        }
        cout << endl;

        cout << "Комплексные корни (сопряженные пары):" << endl;
        if (complex_roots.empty()) {
            cout << "  нет" << endl;
        }
        else {
            for (size_t i = 0; i < complex_roots.size(); i += 2) {
                cout << "  " << complex_roots[i] << " и " << complex_roots[i + 1] << endl;
            }
        }
        cout << endl;

        cout << "ВЫЧИСЛЕНИЕ МЕТОДОМ ХОСМАНА (базовая версия):" << endl;
        double epsilon = 1e-6;

        // Вычисление корней методом Хосмана (базовая версия)
        vector<complex<double>> hosman_roots = hosman_modification_graeffe(coefficients, epsilon);
        cout << "Найдено корней: " << hosman_roots.size() << endl;
        cout << "Вычисленные корни:" << endl;
        for (size_t i = 0; i < hosman_roots.size(); ++i) {
            cout << "  " << i + 1 << ": " << hosman_roots[i] << endl;
        }
        cout << endl;

        // ТЕСТИРОВАНИЕ МЕТОДА HOSMAN_WITH_MULTIPLICITIES
        cout << "ВЫЧИСЛЕНИЕ МЕТОДОМ ХОСМАНА С КРАТНОСТЯМИ:" << endl;

        auto found_graeffe = find_moduli_with_multiplicities_by_graeffe(coefficients_1, epsilon);
        std::cout << "Корни вещественные (Греф):\n";
        for (const auto& root : found_graeffe) {
            std::cout << "root: " << root.value << "    mult: " << root.multiplicity << "\n";
        }
        std::cout << "\n";

        auto found_graeffe_complex = hosman_with_multiplicities(coefficients_1, found_graeffe, epsilon);
        std::cout << "Корни комплексные (Хосман):\n";
        for (auto& r : found_graeffe_complex) {
            std::cout << r.real() << " " << r.imag() << "i\n";
        }
    }

    // ТЕСТ 4
    {
        cout << "ТЕСТ 4: Полином со смешанными корнями (precision source)" << endl;

        unsigned P = 20;
        unsigned num_complex_pairs = 7;
        unsigned num_clusters = 0;
        vector<unsigned> cluster_counts = {};
        vector<float_precision> cluster_radii = {};
        vector<pair<unsigned, unsigned>> multiplicity_groups = {};
        float_precision default_cluster_radius = 0.1;
        bool normalize_coeffs = true;
        uint64_t seed = 7618;

        vector<float_precision> coefficients;
        vector<float_precision> real_roots_repeated;
        vector<float_precision> unique_real_roots;
        vector<unsigned> real_root_multiplicities;
        vector<complex<float_precision>> complex_roots;

        generate_high_degree_polynomial(
            P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
            multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
            coefficients, real_roots_repeated, unique_real_roots,
            real_root_multiplicities, complex_roots
        );

        cout << "СГЕНЕРИРОВАННЫЙ ПОЛИНОМ:" << endl;
        cout << "Степень: " << P << endl;
        cout << "Полином: ";
        print_polynomial(coefficients);
        cout << endl;

        reverse(coefficients.begin(), coefficients.end());
        vector<float_precision> coefficients_1 = coefficients;

        cout << "ИСХОДНЫЕ КОРНИ:" << endl;
        cout << "Вещественные корни: ";
        for (auto root : real_roots_repeated) {
            cout << root << " ";
        }
        cout << endl;

        cout << "Комплексные корни (сопряженные пары):" << endl;
        if (complex_roots.empty()) {
            cout << "  нет" << endl;
        }
        else {
            for (size_t i = 0; i < complex_roots.size(); i += 2) {
                cout << "  " << complex_roots[i] << " и " << complex_roots[i + 1] << endl;
            }
        }
        cout << endl;

        cout << "ВЫЧИСЛЕНИЕ МЕТОДОМ ХОСМАНА (базовая версия):" << endl;
        float_precision epsilon = 1e-6;

        // Вычисление корней методом Хосмана (базовая версия)
        vector<complex<float_precision>> hosman_roots = hosman_modification_graeffe(coefficients, epsilon);
        cout << "Найдено корней: " << hosman_roots.size() << endl;
        cout << "Вычисленные корни:" << endl;
        for (size_t i = 0; i < hosman_roots.size(); ++i) {
            cout << "  " << i + 1 << ": " << hosman_roots[i] << endl;
        }
        cout << endl;

        // ТЕСТИРОВАНИЕ МЕТОДА HOSMAN_WITH_MULTIPLICITIES
        cout << "ВЫЧИСЛЕНИЕ МЕТОДОМ ХОСМАНА С КРАТНОСТЯМИ:" << endl;

        auto found_graeffe = find_moduli_with_multiplicities_by_graeffe(coefficients_1, epsilon);
        std::cout << "Корни вещественные (Греф):\n";
        for (const auto& root : found_graeffe) {
            std::cout << "root: " << root.value << "    mult: " << root.multiplicity << "\n";
        }
        std::cout << "\n";

        auto found_graeffe_complex = hosman_with_multiplicities(coefficients_1, found_graeffe, epsilon);
        std::cout << "Корни комплексные (Хосман):\n";
        for (auto& r : found_graeffe_complex) {
            std::cout << r.real() << " " << r.imag() << "i\n";
        }
    }
}

int main() {
    // Установка русской локали для корректного вывода
    setlocale(LC_ALL, "Russian");

    try {
        cout << "ТЕСТИРОВАНИЕ МОДИФИКАЦИИ ХОСМАНА МЕТОДА ГРЕФЕ" << endl;
        cout << "Генерация полиномов и вычисление комплексных корней" << endl << endl;

        run_hosman_tests();

        cout << "ТЕСТИРОВАНИЕ ЗАВЕРШЕНО" << endl;
        return EXIT_SUCCESS;
    }
    catch (const exception& e) {
        cerr << "ОШИБКА: " << e.what() << endl;
        return EXIT_FAILURE;
    }
}