// Файл для тестирования метода Грефе (вещественные корни). Тестируем улучшенный алгоритм с отбросом дублирующихся модулей
// Реализация: Павлова Анастасия, КМБО-01-22

#include "mathUtils.h"
#include "polynomialUtils.h"
#include "Graeffe_real_roots.h"
#include "generate_high_degree_polynomial.h"

using namespace std;

// Основная функция тестирования - генерирует полиномы и тестирует метод Грефе
void run_graeffe_tests() {
    cout << endl << endl;

    // ТЕСТ 1
    {
        cout << "ТЕСТ 1: Полином со степенью 5" << endl;

        unsigned P = 5;                                                // степень полинома
        unsigned num_complex_pairs = 0;                                 // кол-во комплексных пар
        unsigned num_clusters = 0;                                      // кол-во кластеров
        vector<unsigned> cluster_counts = {};                           // сколько корней в каждом кластере
        vector<double> cluster_radii = {};                              // радиусы кластеров
        vector<pair<unsigned, unsigned>> multiplicity_groups = {};      // сколько корней каждой кратности
        double default_cluster_radius = 0.1;                            // радиус кластера по умолчанию
        bool normalize_coeffs = true;                                   // приводить старший коэффициент к 1?
        uint64_t seed = 12345;                                          // тип рандома

        vector<double> coefficients;                                    // вектор коэффициентов на выходе
        vector<double> real_roots_repeated;                             // вектор всех вещественных корней
        vector<double> unique_real_roots;                               // вектор уникальных вещественных корней
        vector<unsigned> real_root_multiplicities;                      // вектор кратных вещественных корней
        vector<complex<double>> complex_roots;                          // вектор комплексных корней

        generate_high_degree_polynomial(
            P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
            multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
            coefficients, real_roots_repeated, unique_real_roots,
            real_root_multiplicities, complex_roots
        );

        cout << "СГЕНЕРИРОВАННЫЙ ПОЛИНОМ:" << endl;
        cout << "Степень: " << P << endl;

        cout << endl << "Полином: ";
        print_polynomial(coefficients);
        cout << endl;

        cout << "ИСХОДНЫЕ КОРНИ:" << endl;
        cout << "Вещественные корни (с кратностями): ";
        for (auto root : real_roots_repeated) {
            cout << root << " ";
        }
        cout << endl;

        cout << "Уникальные вещественные корни и их кратности:" << endl;
        for (size_t i = 0; i < unique_real_roots.size(); ++i) {
            cout << "  " << unique_real_roots[i] << " : кратность " << real_root_multiplicities[i] << endl;
        }

        cout << "Комплексные корни (сопряженные пары):" << endl;
        if (complex_roots.empty()) {
            cout << "  нет" << endl;
        }
        else {
            for (size_t i = 0; i < complex_roots.size(); i += 2) {
                cout << "  " << complex_roots[i] << " и " << complex_roots[i + 1] << endl;
            }
        }

        cout << endl << "ВЫЧИСЛЕНИЕ МЕТОДОМ ГРЕФЕ" << endl;
        double epsilon = 1e-12;

        // меняем порядок коэффициентов
        vector<double> coefficients_for_graeffe = coefficients;
        reverse(coefficients_for_graeffe.begin(), coefficients_for_graeffe.end());

        // ищем корни методом Грефе
        vector<double> computed_moduli = find_moduli_roots_by_graeffe(coefficients_for_graeffe, epsilon);
        cout << "Вычисленные модули корней: ";
        if (computed_moduli.empty()) {
            cout << "метод не нашел корней";
        } else {
            for (auto modulus : computed_moduli)
                cout << modulus << " ";
        }
        cout << endl;

        // ищем корни методом Грефе с учетом их кратности
        auto moduli_with_mult = find_moduli_with_multiplicities_by_graeffe(coefficients_for_graeffe, epsilon);
        cout << "Модули с кратностями:" << endl;
        if (moduli_with_mult.empty()) {
            cout << "  не найдены" << endl;
        } else {
            for (const auto& m : moduli_with_mult) 
                cout << "  |r| = " << m.value << ", кратность = " << m.multiplicity << endl;
        }
        cout << endl << endl;
    }

    // ТЕСТ 2
    {
        cout << "ТЕСТ 2: Полином со степенью 20" << endl;

        unsigned P = 20;                                                // степень полинома
        unsigned num_complex_pairs = 0;                                 // кол-во комплексных пар
        unsigned num_clusters = 0;                                      // кол-во кластеров
        vector<unsigned> cluster_counts = {};                           // сколько корней в каждом кластере
        vector<double> cluster_radii = {};                              // радиусы кластеров
        vector<pair<unsigned, unsigned>> multiplicity_groups = {};      // сколько корней каждой кратности
        double default_cluster_radius = 0.1;                            // радиус кластера по умолчанию
        bool normalize_coeffs = true;                                   // приводить старший коэффициент к 1?
        uint64_t seed = 12345;                                          // тип рандома

        vector<double> coefficients;                                    // вектор коэффициентов на выходе
        vector<double> real_roots_repeated;                             // вектор всех вещественных корней
        vector<double> unique_real_roots;                               // вектор уникальных вещественных корней
        vector<unsigned> real_root_multiplicities;                      // вектор кратных вещественных корней
        vector<complex<double>> complex_roots;                          // вектор комплексных корней

        generate_high_degree_polynomial(
            P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
            multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
            coefficients, real_roots_repeated, unique_real_roots,
            real_root_multiplicities, complex_roots
        );

        cout << "СГЕНЕРИРОВАННЫЙ ПОЛИНОМ:" << endl;
        cout << "Степень: " << P << endl;
        cout << "Размер вектора коэффициентов" << coefficients.size() << endl;

        cout << endl << "Полином: ";
        print_polynomial(coefficients);
        cout << endl;

        cout << "ИСХОДНЫЕ КОРНИ:" << endl;
        cout << "Вещественные корни (с кратностями): ";
        for (auto root : real_roots_repeated) {
            cout << root << " ";
        }
        cout << endl;

        cout << "Уникальные вещественные корни и их кратности:" << endl;
        for (size_t i = 0; i < unique_real_roots.size(); ++i) {
            cout << "  " << unique_real_roots[i] << " : кратность " << real_root_multiplicities[i] << endl;
        }

        cout << "Комплексные корни (сопряженные пары):" << endl;
        if (complex_roots.empty()) {
            cout << "  нет" << endl;
        }
        else {
            for (size_t i = 0; i < complex_roots.size(); i += 2) {
                cout << "  " << complex_roots[i] << " и " << complex_roots[i + 1] << endl;
            }
        }

        cout << endl << "ВЫЧИСЛЕНИЕ МЕТОДОМ ГРЕФЕ" << endl;
        double epsilon = 1e-12;

        // меняем порядок коэффициентов
        vector<double> coefficients_for_graeffe = coefficients;
        reverse(coefficients_for_graeffe.begin(), coefficients_for_graeffe.end());

        // ищем корни методом Грефе
        vector<double> computed_moduli = find_moduli_roots_by_graeffe(coefficients_for_graeffe, epsilon);
        cout << "Вычисленные модули корней: ";
        if (computed_moduli.empty()) {
            cout << "метод не нашел корней";
        }
        else {
            for (auto modulus : computed_moduli)
                cout << modulus << " ";
        }
        cout << endl;

        // ищем корни методом Грефе с учетом их кратности
        auto moduli_with_mult = find_moduli_with_multiplicities_by_graeffe(coefficients_for_graeffe, epsilon);
        cout << "Модули с кратностями:" << endl;
        if (moduli_with_mult.empty()) {
            cout << "  не найдены" << endl;
        }
        else {
            for (const auto& m : moduli_with_mult)
                cout << "  |r| = " << m.value << ", кратность = " << m.multiplicity << endl;
        }
        cout << endl << endl;
    }

    // ТЕСТ 3
    {
        cout << "ТЕСТ 3: Полином с кратными корнями и степень 5 (precision source)" << endl;

        unsigned P = 5;                                                // степень полинома
        unsigned num_complex_pairs = 0;                                 // кол-во комплексных пар
        unsigned num_clusters = 0;                                      // кол-во кластеров
        vector<unsigned> cluster_counts = {};                           // сколько корней в каждом кластере
        vector<float_precision> cluster_radii = {};                              // радиусы кластеров
        vector<pair<unsigned, unsigned>> multiplicity_groups = { {2, 3}, {1, 3} };      // сколько корней каждой кратности
        float_precision default_cluster_radius = 0.1;                            // радиус кластера по умолчанию
        bool normalize_coeffs = true;                                   // приводить старший коэффициент к 1?
        uint64_t seed = 12345;                                          // тип рандома

        vector<float_precision> coefficients;                                    // вектор коэффициентов на выходе
        vector<float_precision> real_roots_repeated;                             // вектор всех вещественных корней
        vector<float_precision> unique_real_roots;                               // вектор уникальных вещественных корней
        vector<unsigned> real_root_multiplicities;                      // вектор кратных вещественных корней
        vector<complex<float_precision>> complex_roots;                          // вектор комплексных корней

        generate_high_degree_polynomial(
            P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
            multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
            coefficients, real_roots_repeated, unique_real_roots,
            real_root_multiplicities, complex_roots
        );

        cout << "СГЕНЕРИРОВАННЫЙ ПОЛИНОМ:" << endl;
        cout << "Степень: " << P << endl;

        cout << endl << "Полином: ";
        print_polynomial(coefficients);
        cout << endl;

        cout << "ИСХОДНЫЕ КОРНИ:" << endl;
        cout << "Вещественные корни (с кратностями): ";
        for (auto root : real_roots_repeated) {
            cout << root << " ";
        }
        cout << endl;

        cout << "Уникальные вещественные корни и их кратности:" << endl;
        for (size_t i = 0; i < unique_real_roots.size(); ++i) {
            cout << "  " << unique_real_roots[i] << " : кратность " << real_root_multiplicities[i] << endl;
        }

        cout << "Комплексные корни (сопряженные пары):" << endl;
        if (complex_roots.empty()) {
            cout << "  нет" << endl;
        }
        else {
            for (size_t i = 0; i < complex_roots.size(); i += 2) {
                cout << "  " << complex_roots[i] << " и " << complex_roots[i + 1] << endl;
            }
        }

        cout << endl << "ВЫЧИСЛЕНИЕ МЕТОДОМ ГРЕФЕ" << endl;
        float_precision epsilon = 1e-6;

        // меняем порядок коэффициентов
        vector<float_precision> coefficients_for_graeffe = coefficients;
        reverse(coefficients_for_graeffe.begin(), coefficients_for_graeffe.end());

        // ищем корни методом Грефе
        vector<float_precision> computed_moduli = find_moduli_roots_by_graeffe(coefficients_for_graeffe, epsilon);
        cout << "Вычисленные модули корней: ";
        if (computed_moduli.empty()) {
            cout << "метод не нашел корней";
        }
        else {
            for (auto modulus : computed_moduli)
                cout << modulus << " ";
        }
        cout << endl;

        // ищем корни методом Грефе с учетом их кратности
        auto moduli_with_mult = find_moduli_with_multiplicities_by_graeffe(coefficients_for_graeffe, epsilon);
        cout << "Модули с кратностями:" << endl;
        if (moduli_with_mult.empty()) {
            cout << "  не найдены" << endl;
        }
        else {
            for (const auto& m : moduli_with_mult)
                cout << "  |r| = " << m.value << ", кратность = " << m.multiplicity << endl;
        }
        cout << endl << endl;
    }

    // ТЕСТ 4
    {
        cout << "ТЕСТ 4: Полином с кластеризованными корнями (precision source)" << endl;

        unsigned P = 10;                                                // степень полинома
        unsigned num_complex_pairs = 0;                                 // кол-во комплексных пар
        unsigned num_clusters = 0;                                      // кол-во кластеров
        vector<unsigned> cluster_counts = {  };                           // сколько корней в каждом кластере
        vector<float_precision> cluster_radii = { };                              // радиусы кластеров
        vector<pair<unsigned, unsigned>> multiplicity_groups = {};      // сколько корней каждой кратности
        float_precision default_cluster_radius = 0.1;                            // радиус кластера по умолчанию
        bool normalize_coeffs = true;                                   // приводить старший коэффициент к 1?
        uint64_t seed = 8351;                                          // тип рандома

        vector<float_precision> coefficients;                                    // вектор коэффициентов на выходе
        vector<float_precision> real_roots_repeated;                             // вектор всех вещественных корней
        vector<float_precision> unique_real_roots;                               // вектор уникальных вещественных корней
        vector<unsigned> real_root_multiplicities;                      // вектор кратных вещественных корней
        vector<complex<float_precision>> complex_roots;                          // вектор комплексных корней

        generate_high_degree_polynomial(
            P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
            multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
            coefficients, real_roots_repeated, unique_real_roots,
            real_root_multiplicities, complex_roots
        );

        cout << "СГЕНЕРИРОВАННЫЙ ПОЛИНОМ:" << endl;
        cout << "Степень: " << P << endl;

        cout << endl << "Полином: ";
        print_polynomial(coefficients);
        cout << endl;

        cout << "ИСХОДНЫЕ КОРНИ:" << endl;
        cout << "Вещественные корни (с кратностями): ";
        for (auto root : real_roots_repeated) {
            cout << root << " ";
        }
        cout << endl;

        cout << "Уникальные вещественные корни и их кратности:" << endl;
        for (size_t i = 0; i < unique_real_roots.size(); ++i) {
            cout << "  " << unique_real_roots[i] << " : кратность " << real_root_multiplicities[i] << endl;
        }

        cout << "Комплексные корни (сопряженные пары):" << endl;
        if (complex_roots.empty()) {
            cout << "  нет" << endl;
        }
        else {
            for (size_t i = 0; i < complex_roots.size(); i += 2) {
                cout << "  " << complex_roots[i] << " и " << complex_roots[i + 1] << endl;
            }
        }

        cout << endl << "ВЫЧИСЛЕНИЕ МЕТОДОМ ГРЕФЕ" << endl;
        float_precision epsilon = 1e-12;

        // меняем порядок коэффициентов
        vector<float_precision> coefficients_for_graeffe = coefficients;
        reverse(coefficients_for_graeffe.begin(), coefficients_for_graeffe.end());

        // ищем корни методом Грефе
        vector<float_precision> computed_moduli = find_moduli_roots_by_graeffe(coefficients_for_graeffe, epsilon);
        cout << "Вычисленные модули корней: ";
        if (computed_moduli.empty()) {
            cout << "метод не нашел корней";
        }
        else {
            for (auto modulus : computed_moduli)
                cout << modulus << " ";
        }
        cout << endl;

        // ищем корни методом Грефе с учетом их кратности
        auto moduli_with_mult = find_moduli_with_multiplicities_by_graeffe(coefficients_for_graeffe, epsilon);
        cout << "Модули с кратностями:" << endl;
        if (moduli_with_mult.empty()) {
            cout << "  не найдены" << endl;
        }
        else {
            for (const auto& m : moduli_with_mult)
                cout << "  |r| = " << m.value << ", кратность = " << m.multiplicity << endl;
        }
        cout << endl << endl;
    }

    // ТЕСТ 5
    {
        cout << "ТЕСТ 5: Полином степени 30 с кластерами (3 шт.) и precision source" << endl;

        unsigned P = 30;                                                // степень полинома
        unsigned num_complex_pairs = 0;                                 // кол-во комплексных пар
        unsigned num_clusters = 3;                                      // кол-во кластеров
        vector<unsigned> cluster_counts = { 2, 3, 6 };                           // сколько корней в каждом кластере
        vector<float_precision> cluster_radii = {0.001, 0.01, 0.0001};                              // радиусы кластеров
        vector<pair<unsigned, unsigned>> multiplicity_groups = {};      // сколько корней каждой кратности
        float_precision default_cluster_radius = 0.1;                            // радиус кластера по умолчанию
        bool normalize_coeffs = true;                                   // приводить старший коэффициент к 1?
        uint64_t seed = 12345;                                          // тип рандома

        vector<float_precision> coefficients;                                    // вектор коэффициентов на выходе
        vector<float_precision> real_roots_repeated;                             // вектор всех вещественных корней
        vector<float_precision> unique_real_roots;                               // вектор уникальных вещественных корней
        vector<unsigned> real_root_multiplicities;                      // вектор кратных вещественных корней
        vector<complex<float_precision>> complex_roots;                          // вектор комплексных корней

        generate_high_degree_polynomial(
            P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
            multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
            coefficients, real_roots_repeated, unique_real_roots,
            real_root_multiplicities, complex_roots
        );

        cout << "СГЕНЕРИРОВАННЫЙ ПОЛИНОМ:" << endl;
        cout << "Степень: " << P << endl;

        cout << endl << "Полином: ";
        print_polynomial(coefficients);
        cout << endl;

        cout << "ИСХОДНЫЕ КОРНИ:" << endl;
        cout << "Вещественные корни (с кратностями): ";
        for (auto root : real_roots_repeated) {
            cout << root << " ";
        }
        cout << endl;

        cout << "Уникальные вещественные корни и их кратности:" << endl;
        for (size_t i = 0; i < unique_real_roots.size(); ++i) {
            cout << "  " << unique_real_roots[i] << " : кратность " << real_root_multiplicities[i] << endl;
        }

        cout << "Комплексные корни (сопряженные пары):" << endl;
        if (complex_roots.empty()) {
            cout << "  нет" << endl;
        }
        else {
            for (size_t i = 0; i < complex_roots.size(); i += 2) {
                cout << "  " << complex_roots[i] << " и " << complex_roots[i + 1] << endl;
            }
        }

        cout << endl << "ВЫЧИСЛЕНИЕ МЕТОДОМ ГРЕФЕ" << endl;
        float_precision epsilon = 1e-12;

        // меняем порядок коэффициентов
        vector<float_precision> coefficients_for_graeffe = coefficients;
        reverse(coefficients_for_graeffe.begin(), coefficients_for_graeffe.end());

        // ищем корни методом Грефе
        vector<float_precision> computed_moduli = find_moduli_roots_by_graeffe(coefficients_for_graeffe, epsilon);
        cout << "Вычисленные модули корней: ";
        if (computed_moduli.empty()) {
            cout << "метод не нашел корней";
        }
        else {
            for (auto modulus : computed_moduli)
                cout << modulus << " ";
        }
        cout << endl;

        // ищем корни методом Грефе с учетом их кратности
        auto moduli_with_mult = find_moduli_with_multiplicities_by_graeffe(coefficients_for_graeffe, epsilon);
        cout << "Модули с кратностями:" << endl;
        if (moduli_with_mult.empty()) {
            cout << "  не найдены" << endl;
        }
        else {
            for (const auto& m : moduli_with_mult)
                cout << "  |r| = " << m.value << ", кратность = " << m.multiplicity << endl;
        }
        cout << endl << endl;
    }
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
    }
    catch (const exception& e) {
        cerr << "ОШИБКА: " << e.what() << endl;
        return EXIT_FAILURE;
    }
}