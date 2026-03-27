// Файл для обширного тестирования метода Сагралова (ANewDsc) – поиск вещественных корней
// Реализация: Павлова Анастасия, КМБО-01-22
// МАСШТАБНОЕ ТЕСТИРОВАНИЕ для дипломной работы

#include "NumericConstants.h"
#include "Sagralov_ANewDsc.h"
#include "generate_high_degree_polynomial.h"
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <map>
#include <algorithm>
#include <numeric>

using namespace std;
using namespace chrono;

// Структура для хранения результатов теста
struct TestResult {
    string test_name;
    unsigned degree;
    unsigned expected_real_roots;
    unsigned found_intervals;
    double execution_time_ms;
    bool all_roots_found;
    double average_interval_width;
    double max_interval_width;
    double min_interval_width;
    vector<double> expected_roots;
    vector<pair<double, double>> found_intervals_list;
};

// Глобальная статистика
map<string, vector<TestResult>> all_test_results;

// Вспомогательные функции для форматированного вывода
static void print_header(const string& text) {
    cout << "\n" << string(100, '=') << endl;
    cout << ">>> " << text << " <<<" << endl;
    cout << string(100, '=') << endl;
}

static void print_subheader(const string& text) {
    cout << "\n" << string(80, '-') << endl;
    cout << "--- " << text << " ---" << endl;
    cout << string(80, '-') << endl;
}

// Функция для оценки качества найденных интервалов
template<typename fp_t>
static void evaluate_root_detection(
    const vector<fp_t>& expected_roots,
    const vector<pair<fp_t, fp_t>>& found_intervals,
    fp_t tolerance,
    TestResult& result) {

    result.expected_roots.clear();
    for (const auto& r : expected_roots) {
        result.expected_roots.push_back(static_cast<double>(r));
    }

    result.found_intervals_list.clear();
    result.expected_real_roots = static_cast<unsigned>(expected_roots.size());
    result.found_intervals = static_cast<unsigned>(found_intervals.size());
    result.all_roots_found = true;

    vector<double> widths;
    for (const auto& interval : found_intervals) {
        double width = static_cast<double>(interval.second - interval.first);
        widths.push_back(width);
        result.found_intervals_list.push_back({
            static_cast<double>(interval.first),
            static_cast<double>(interval.second)
            });
    }

    if (!widths.empty()) {
        double sum = 0.0;
        for (double w : widths) {
            sum += w;
        }
        result.average_interval_width = sum / static_cast<double>(widths.size());
        result.max_interval_width = *max_element(widths.begin(), widths.end());
        result.min_interval_width = *min_element(widths.begin(), widths.end());
    }
    else {
        result.average_interval_width = 0.0;
        result.max_interval_width = 0.0;
        result.min_interval_width = 0.0;
    }

    // Проверяем, покрыты ли все ожидаемые корни интервалами
    vector<bool> root_found(expected_roots.size(), false);

    for (size_t i = 0; i < expected_roots.size(); ++i) {
        for (const auto& interval : found_intervals) {
            double root = static_cast<double>(expected_roots[i]);
            double left = static_cast<double>(interval.first);
            double right = static_cast<double>(interval.second);
            double tol = static_cast<double>(tolerance);

            if (root >= left - tol && root <= right + tol) {
                root_found[i] = true;
                break;
            }
        }
    }

    for (bool found : root_found) {
        if (!found) {
            result.all_roots_found = false;
            break;
        }
    }
}

// Функция для печати подробных результатов теста
template<typename fp_t>
static void print_detailed_results(
    const string& test_name,
    const vector<fp_t>& coefficients,
    const vector<fp_t>& expected_roots,
    const vector<pair<fp_t, fp_t>>& found_intervals,
    double execution_time,
    const TestResult& result) {

    cout << "\n=== ДЕТАЛИ ТЕСТА: " << test_name << " ===\n";

    cout << "Время: " << fixed << setprecision(3)
        << execution_time << " мс\n";

    cout << "Ожидалось корней: " << expected_roots.size()
        << " | Найдено интервалов: " << found_intervals.size() << "\n";

    cout << (result.all_roots_found ? "✔ ВСЕ НАЙДЕНЫ\n"
        : "✘ ПРОПУЩЕНЫ КОРНИ\n");

    // ================= ROOTS =================
    cout << "\nОжидаемые корни:\n";
    for (size_t i = 0; i < expected_roots.size(); ++i) {
        cout << "  r" << i << " = " << expected_roots[i] << "\n";
    }

    // ================= INTERVALS =================
    cout << "\nИнтервалы:\n";

    vector<pair<double, double>> intervals;
    for (auto& it : found_intervals) {
        intervals.push_back({
            (double)it.first,
            (double)it.second
            });
    }

    sort(intervals.begin(), intervals.end());

    for (size_t i = 0; i < intervals.size(); ++i) {
        double w = intervals[i].second - intervals[i].first;
        cout << "  I" << i
            << " = [" << intervals[i].first
            << ", " << intervals[i].second << "]"
            << " width=" << scientific << w << "\n";
    }

    // ================= ПЕРЕСЕЧЕНИЯ =================
    cout << "\nАнализ интервалов:\n";

    int overlaps = 0;
    int duplicates = 0;
    int nested = 0;

    for (size_t i = 0; i < intervals.size(); ++i) {
        for (size_t j = i + 1; j < intervals.size(); ++j) {

            auto A = intervals[i];
            auto B = intervals[j];

            // пересечение
            if (A.second >= B.first && B.second >= A.first) {
                overlaps++;
                cout << "  ⚠ Пересечение: I" << i << " & I" << j << "\n";
            }

            // почти одинаковые
            if (fabs(A.first - B.first) < 1e-12 &&
                fabs(A.second - B.second) < 1e-12) {
                duplicates++;
                cout << "  ⚠ Дубликаты: I" << i << " & I" << j << "\n";
            }

            // вложенность
            if (A.first >= B.first && A.second <= B.second) {
                nested++;
                cout << "  ⚠ Вложен: I" << i << " в I" << j << "\n";
            }
        }
    }

    cout << "\nСтатистика:\n";
    cout << "  Пересечений: " << overlaps << "\n";
    cout << "  Дубликатов: " << duplicates << "\n";
    cout << "  Вложенных: " << nested << "\n";

    cout << "  Средняя ширина: " << scientific << result.average_interval_width << "\n";
    cout << "  Мин ширина: " << result.min_interval_width << "\n";
    cout << "  Макс ширина: " << result.max_interval_width << "\n";

    cout << string(80, '=') << "\n";
}

// ==================== КАТЕГОРИЯ 1: ПРОСТЫЕ КОРНИ ====================
static void test_simple_roots() {
    print_subheader("КАТЕГОРИЯ 1: ПРОСТЫЕ ВЕЩЕСТВЕННЫЕ КОРНИ");

    vector<unsigned> degrees = { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100 };
    vector<uint64_t> seeds = { 12345, 23456, 34567, 45678, 56789 };

    for (unsigned degree : degrees) {
        for (uint64_t seed : seeds) {
            string test_name = "Простые_корни_степень" + to_string(degree) + "_seed" + to_string(seed);

            unsigned num_complex_pairs = 0;
            unsigned num_clusters = 0;
            vector<unsigned> cluster_counts = {};
            vector<double> cluster_radii = {};
            vector<pair<unsigned, unsigned>> multiplicity_groups = {};
            double default_cluster_radius = 0.1;
            bool normalize_coeffs = true;

            vector<double> coefficients;
            vector<double> real_roots_repeated;
            vector<double> unique_real_roots;
            vector<unsigned> real_root_multiplicities;
            vector<complex<double>> complex_roots;

            generate_high_degree_polynomial(
                degree, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
                multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
                coefficients, real_roots_repeated, unique_real_roots,
                real_root_multiplicities, complex_roots);

            double TOLERANCE = numeric_constants::adaptive_epsilon<double>(
                numeric_constants::EPSILON_SCALE_PRECISE);

            auto start_time = high_resolution_clock::now();
            auto found = find_real_roots_ANewDsc(coefficients, TOLERANCE);
            auto end_time = high_resolution_clock::now();

            double execution_time = duration<double, milli>(end_time - start_time).count();

            TestResult result;
            result.test_name = test_name;
            result.degree = degree;
            result.execution_time_ms = execution_time;

            evaluate_root_detection(real_roots_repeated, found, TOLERANCE, result);

            all_test_results["1. Простые корни"].push_back(result);

            if (degree <= 30 || seed == seeds[0]) {
                print_detailed_results(test_name, coefficients, real_roots_repeated,
                    found, execution_time, result);
            }
        }
    }
}

// ==================== КАТЕГОРИЯ 2: КРАТНЫЕ КОРНИ ====================
static void test_multiple_roots() {
    print_subheader("КАТЕГОРИЯ 2: КРАТНЫЕ ВЕЩЕСТВЕННЫЕ КОРНИ");

    struct MultipleRootsConfig {
        unsigned degree;
        vector<pair<unsigned, unsigned>> multiplicities;
        uint64_t seed;
    };

    // Создаем вектор конфигураций через push_back
    vector<MultipleRootsConfig> configs;

    MultipleRootsConfig cfg1;
    cfg1.degree = 8;
    cfg1.multiplicities = { {2,2}, {1,3} };
    cfg1.seed = 11111;
    configs.push_back(cfg1);

    MultipleRootsConfig cfg2;
    cfg2.degree = 12;
    cfg2.multiplicities = { {3,2}, {2,2}, {1,2} };
    cfg2.seed = 22222;
    configs.push_back(cfg2);

    MultipleRootsConfig cfg3;
    cfg3.degree = 15;
    cfg3.multiplicities = { {2,3}, {3,2}, {1,2} };
    cfg3.seed = 33333;
    configs.push_back(cfg3);

    MultipleRootsConfig cfg4;
    cfg4.degree = 20;
    cfg4.multiplicities = { {4,2}, {3,2}, {2,2}, {1,2} };
    cfg4.seed = 44444;
    configs.push_back(cfg4);

    MultipleRootsConfig cfg5;
    cfg5.degree = 24;
    cfg5.multiplicities = { {2,4}, {3,3}, {4,2} };
    cfg5.seed = 55555;
    configs.push_back(cfg5);

    MultipleRootsConfig cfg6;
    cfg6.degree = 30;
    cfg6.multiplicities = { {5,2}, {4,2}, {3,2}, {2,2}, {1,2} };
    cfg6.seed = 66666;
    configs.push_back(cfg6);

    MultipleRootsConfig cfg7;
    cfg7.degree = 36;
    cfg7.multiplicities = { {6,2}, {5,2}, {4,2}, {3,2} };
    cfg7.seed = 77777;
    configs.push_back(cfg7);

    MultipleRootsConfig cfg8;
    cfg8.degree = 42;
    cfg8.multiplicities = { {7,2}, {6,2}, {5,2}, {4,2} };
    cfg8.seed = 88888;
    configs.push_back(cfg8);

    MultipleRootsConfig cfg9;
    cfg9.degree = 50;
    cfg9.multiplicities = { {10,2}, {8,2}, {6,2}, {5,2} };
    cfg9.seed = 99999;
    configs.push_back(cfg9);

    MultipleRootsConfig cfg10;
    cfg10.degree = 60;
    cfg10.multiplicities = { {12,2}, {10,2}, {8,2}, {6,2} };
    cfg10.seed = 10101;
    configs.push_back(cfg10);

    MultipleRootsConfig cfg11;
    cfg11.degree = 75;
    cfg11.multiplicities = { {15,2}, {12,2}, {10,2}, {8,2} };
    cfg11.seed = 20202;
    configs.push_back(cfg11);

    MultipleRootsConfig cfg12;
    cfg12.degree = 100;
    cfg12.multiplicities = { {20,2}, {16,2}, {12,2}, {10,2} };
    cfg12.seed = 30303;
    configs.push_back(cfg12);

    for (const auto& config : configs) {
        string test_name = "Кратные_корни_степень" + to_string(config.degree) +
            "_мультипл" + to_string(config.multiplicities.size());

        unsigned num_complex_pairs = 0;
        unsigned num_clusters = 0;
        vector<unsigned> cluster_counts = {};
        vector<double> cluster_radii = {};
        double default_cluster_radius = 0.1;
        bool normalize_coeffs = true;

        vector<double> coefficients;
        vector<double> real_roots_repeated;
        vector<double> unique_real_roots;
        vector<unsigned> real_root_multiplicities;
        vector<complex<double>> complex_roots;

        generate_high_degree_polynomial(
            config.degree, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
            config.multiplicities, default_cluster_radius, normalize_coeffs, config.seed,
            coefficients, real_roots_repeated, unique_real_roots,
            real_root_multiplicities, complex_roots);

        double TOLERANCE = numeric_constants::adaptive_epsilon<double>(
            numeric_constants::EPSILON_SCALE_PRECISE * 2);

        auto start_time = high_resolution_clock::now();
        auto found = find_real_roots_ANewDsc(coefficients, TOLERANCE);
        auto end_time = high_resolution_clock::now();

        double execution_time = duration<double, milli>(end_time - start_time).count();

        TestResult result;
        result.test_name = test_name;
        result.degree = config.degree;
        result.execution_time_ms = execution_time;

        evaluate_root_detection(real_roots_repeated, found, TOLERANCE, result);

        all_test_results["2. Кратные корни"].push_back(result);

        if (config.degree <= 30) {
            print_detailed_results(test_name, coefficients, real_roots_repeated,
                found, execution_time, result);
        }
    }
}

// ==================== КАТЕГОРИЯ 3: КЛАСТЕРЫ КОРНЕЙ ====================
static void test_clustered_roots() {
    print_subheader("КАТЕГОРИЯ 3: КЛАСТЕРЫ КОРНЕЙ");

    struct ClusterConfig {
        unsigned degree;
        vector<unsigned> cluster_sizes;
        vector<double> cluster_radii;
        uint64_t seed;
    };

    // Создаем вектор конфигураций через push_back
    vector<ClusterConfig> configs;

    ClusterConfig cfg1; cfg1.degree = 10; cfg1.cluster_sizes = { 2, 3 }; cfg1.cluster_radii = { 1e-2, 1e-3 }; cfg1.seed = 12121; configs.push_back(cfg1);
    ClusterConfig cfg2; cfg2.degree = 15; cfg2.cluster_sizes = { 3, 4 }; cfg2.cluster_radii = { 1e-3, 1e-4 }; cfg2.seed = 23232; configs.push_back(cfg2);
    ClusterConfig cfg3; cfg3.degree = 20; cfg3.cluster_sizes = { 4, 5 }; cfg3.cluster_radii = { 1e-4, 1e-5 }; cfg3.seed = 34343; configs.push_back(cfg3);
    ClusterConfig cfg4; cfg4.degree = 25; cfg4.cluster_sizes = { 5, 5 }; cfg4.cluster_radii = { 1e-5, 1e-6 }; cfg4.seed = 45454; configs.push_back(cfg4);
    ClusterConfig cfg5; cfg5.degree = 30; cfg5.cluster_sizes = { 6, 6 }; cfg5.cluster_radii = { 1e-6, 1e-7 }; cfg5.seed = 56565; configs.push_back(cfg5);
    ClusterConfig cfg6; cfg6.degree = 35; cfg6.cluster_sizes = { 7, 7 }; cfg6.cluster_radii = { 1e-7, 1e-8 }; cfg6.seed = 67676; configs.push_back(cfg6);
    ClusterConfig cfg7; cfg7.degree = 40; cfg7.cluster_sizes = { 8, 8 }; cfg7.cluster_radii = { 1e-8, 1e-9 }; cfg7.seed = 78787; configs.push_back(cfg7);
    ClusterConfig cfg8; cfg8.degree = 45; cfg8.cluster_sizes = { 9, 9 }; cfg8.cluster_radii = { 1e-9, 1e-10 }; cfg8.seed = 89898; configs.push_back(cfg8);
    ClusterConfig cfg9; cfg9.degree = 50; cfg9.cluster_sizes = { 10, 10 }; cfg9.cluster_radii = { 1e-10, 1e-11 }; cfg9.seed = 90909; configs.push_back(cfg9);
    ClusterConfig cfg10; cfg10.degree = 60; cfg10.cluster_sizes = { 12, 12 }; cfg10.cluster_radii = { 1e-11, 1e-12 }; cfg10.seed = 12121; configs.push_back(cfg10);
    ClusterConfig cfg11; cfg11.degree = 75; cfg11.cluster_sizes = { 15, 15 }; cfg11.cluster_radii = { 1e-12, 1e-13 }; cfg11.seed = 23232; configs.push_back(cfg11);
    ClusterConfig cfg12; cfg12.degree = 100; cfg12.cluster_sizes = { 20, 20 }; cfg12.cluster_radii = { 1e-13, 1e-14 }; cfg12.seed = 34343; configs.push_back(cfg12);

    for (const auto& config : configs) {
        string test_name = "Кластеры_степень" + to_string(config.degree) +
            "_кластеров" + to_string(config.cluster_sizes.size());

        unsigned num_complex_pairs = 0;
        unsigned num_clusters = static_cast<unsigned>(config.cluster_sizes.size());
        vector<pair<unsigned, unsigned>> multiplicity_groups = {};
        double default_cluster_radius = 0.1;
        bool normalize_coeffs = true;

        vector<double> coefficients;
        vector<double> real_roots_repeated;
        vector<double> unique_real_roots;
        vector<unsigned> real_root_multiplicities;
        vector<complex<double>> complex_roots;

        generate_high_degree_polynomial(
            config.degree, num_complex_pairs, num_clusters,
            config.cluster_sizes, config.cluster_radii,
            multiplicity_groups, default_cluster_radius, normalize_coeffs, config.seed,
            coefficients, real_roots_repeated, unique_real_roots,
            real_root_multiplicities, complex_roots);

        double TOLERANCE = numeric_constants::adaptive_epsilon<double>(
            numeric_constants::EPSILON_SCALE_COARSE * 0.1);

        auto start_time = high_resolution_clock::now();
        auto found = find_real_roots_ANewDsc(coefficients, TOLERANCE);
        auto end_time = high_resolution_clock::now();

        double execution_time = duration<double, milli>(end_time - start_time).count();

        TestResult result;
        result.test_name = test_name;
        result.degree = config.degree;
        result.execution_time_ms = execution_time;

        evaluate_root_detection(real_roots_repeated, found, TOLERANCE, result);

        all_test_results["3. Кластеры корней"].push_back(result);

        if (config.degree <= 30) {
            print_detailed_results(test_name, coefficients, real_roots_repeated,
                found, execution_time, result);
        }
    }
}

// ==================== КАТЕГОРИЯ 4: СМЕШАННЫЕ СЛУЧАИ ====================
static void test_mixed_cases() {
    print_subheader("КАТЕГОРИЯ 4: СМЕШАННЫЕ СЛУЧАИ (КРАТНЫЕ + КЛАСТЕРЫ)");

    struct MixedConfig {
        unsigned degree;
        unsigned num_complex_pairs;
        vector<unsigned> cluster_sizes;
        vector<double> cluster_radii;
        vector<pair<unsigned, unsigned>> multiplicities;
        uint64_t seed;
    };

    // Создаем вектор конфигураций через push_back
    vector<MixedConfig> configs;

    MixedConfig cfg1;
    cfg1.degree = 15;
    cfg1.num_complex_pairs = 2;
    cfg1.cluster_sizes = { 2, 3 };
    cfg1.cluster_radii = { 1e-3, 1e-4 };
    cfg1.multiplicities = { {2,2}, {1,2} };
    cfg1.seed = 13579;
    configs.push_back(cfg1);

    MixedConfig cfg2;
    cfg2.degree = 20;
    cfg2.num_complex_pairs = 3;
    cfg2.cluster_sizes = { 3, 4 };
    cfg2.cluster_radii = { 1e-4, 1e-5 };
    cfg2.multiplicities = { {3,2}, {2,2} };
    cfg2.seed = 24680;
    configs.push_back(cfg2);

    MixedConfig cfg3;
    cfg3.degree = 25;
    cfg3.num_complex_pairs = 4;
    cfg3.cluster_sizes = { 4, 4 };
    cfg3.cluster_radii = { 1e-5, 1e-6 };
    cfg3.multiplicities = { {4,2}, {3,2} };
    cfg3.seed = 35791;
    configs.push_back(cfg3);

    MixedConfig cfg4;
    cfg4.degree = 30;
    cfg4.num_complex_pairs = 5;
    cfg4.cluster_sizes = { 5, 5 };
    cfg4.cluster_radii = { 1e-6, 1e-7 };
    cfg4.multiplicities = { {5,2}, {4,2} };
    cfg4.seed = 46802;
    configs.push_back(cfg4);

    MixedConfig cfg5;
    cfg5.degree = 35;
    cfg5.num_complex_pairs = 5;
    cfg5.cluster_sizes = { 6, 6 };
    cfg5.cluster_radii = { 1e-7, 1e-8 };
    cfg5.multiplicities = { {6,2}, {5,2} };
    cfg5.seed = 57913;
    configs.push_back(cfg5);

    MixedConfig cfg6;
    cfg6.degree = 40;
    cfg6.num_complex_pairs = 6;
    cfg6.cluster_sizes = { 7, 7 };
    cfg6.cluster_radii = { 1e-8, 1e-9 };
    cfg6.multiplicities = { {7,2}, {6,2} };
    cfg6.seed = 68024;
    configs.push_back(cfg6);

    MixedConfig cfg7;
    cfg7.degree = 45;
    cfg7.num_complex_pairs = 6;
    cfg7.cluster_sizes = { 8, 8 };
    cfg7.cluster_radii = { 1e-9, 1e-10 };
    cfg7.multiplicities = { {8,2}, {7,2} };
    cfg7.seed = 79135;
    configs.push_back(cfg7);

    MixedConfig cfg8;
    cfg8.degree = 50;
    cfg8.num_complex_pairs = 7;
    cfg8.cluster_sizes = { 9, 9 };
    cfg8.cluster_radii = { 1e-10, 1e-11 };
    cfg8.multiplicities = { {9,2}, {8,2} };
    cfg8.seed = 80246;
    configs.push_back(cfg8);

    MixedConfig cfg9;
    cfg9.degree = 60;
    cfg9.num_complex_pairs = 8;
    cfg9.cluster_sizes = { 10, 10 };
    cfg9.cluster_radii = { 1e-11, 1e-12 };
    cfg9.multiplicities = { {10,2}, {9,2} };
    cfg9.seed = 91357;
    configs.push_back(cfg9);

    MixedConfig cfg10;
    cfg10.degree = 75;
    cfg10.num_complex_pairs = 10;
    cfg10.cluster_sizes = { 12, 12 };
    cfg10.cluster_radii = { 1e-12, 1e-13 };
    cfg10.multiplicities = { {12,2}, {10,2} };
    cfg10.seed = 12468;
    configs.push_back(cfg10);

    MixedConfig cfg11;
    cfg11.degree = 100;
    cfg11.num_complex_pairs = 15;
    cfg11.cluster_sizes = { 15, 15 };
    cfg11.cluster_radii = { 1e-13, 1e-14 };
    cfg11.multiplicities = { {15,2}, {12,2} };
    cfg11.seed = 23579;
    configs.push_back(cfg11);

    for (const auto& config : configs) {
        string test_name = "Смешанный_степень" + to_string(config.degree) +
            "_компл" + to_string(config.num_complex_pairs);

        unsigned num_clusters = static_cast<unsigned>(config.cluster_sizes.size());
        double default_cluster_radius = 0.1;
        bool normalize_coeffs = true;

        vector<double> coefficients;
        vector<double> real_roots_repeated;
        vector<double> unique_real_roots;
        vector<unsigned> real_root_multiplicities;
        vector<complex<double>> complex_roots;

        generate_high_degree_polynomial(
            config.degree, config.num_complex_pairs, num_clusters,
            config.cluster_sizes, config.cluster_radii,
            config.multiplicities, default_cluster_radius, normalize_coeffs, config.seed,
            coefficients, real_roots_repeated, unique_real_roots,
            real_root_multiplicities, complex_roots);

        double TOLERANCE = numeric_constants::adaptive_epsilon<double>(
            numeric_constants::EPSILON_SCALE_COARSE * 0.5);

        auto start_time = high_resolution_clock::now();
        auto found = find_real_roots_ANewDsc(coefficients, TOLERANCE);
        auto end_time = high_resolution_clock::now();

        double execution_time = duration<double, milli>(end_time - start_time).count();

        TestResult result;
        result.test_name = test_name;
        result.degree = config.degree;
        result.execution_time_ms = execution_time;

        evaluate_root_detection(real_roots_repeated, found, TOLERANCE, result);

        all_test_results["4. Смешанные случаи"].push_back(result);

        if (config.degree <= 30) {
            print_detailed_results(test_name, coefficients, real_roots_repeated,
                found, execution_time, result);
        }
    }
}

// ==================== КАТЕГОРИЯ 5: ЭКСТРЕМАЛЬНЫЕ СЛУЧАИ ====================
static void test_extreme_cases() {
    print_subheader("КАТЕГОРИЯ 5: ЭКСТРЕМАЛЬНЫЕ СЛУЧАИ");

    // Тест 1: Полиномы очень высокой степени
    {
        vector<unsigned> high_degrees = { 80, 90, 100, 120, 150, 180, 200 };

        for (unsigned degree : high_degrees) {
            string test_name = "Экстремально_высокая_степень_" + to_string(degree);

            unsigned num_complex_pairs = degree / 4; // 25% комплексных корней
            unsigned num_clusters = 0;
            vector<unsigned> cluster_counts = {};
            vector<double> cluster_radii = {};
            vector<pair<unsigned, unsigned>> multiplicity_groups = {};
            double default_cluster_radius = 0.1;
            bool normalize_coeffs = true;
            uint64_t seed = 98765 + degree;

            vector<double> coefficients;
            vector<double> real_roots_repeated;
            vector<double> unique_real_roots;
            vector<unsigned> real_root_multiplicities;
            vector<complex<double>> complex_roots;

            generate_high_degree_polynomial(
                degree, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
                multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
                coefficients, real_roots_repeated, unique_real_roots,
                real_root_multiplicities, complex_roots);

            double TOLERANCE = numeric_constants::adaptive_epsilon<double>(
                numeric_constants::EPSILON_SCALE_COARSE * 10);

            auto start_time = high_resolution_clock::now();
            auto found = find_real_roots_ANewDsc(coefficients, TOLERANCE);
            auto end_time = high_resolution_clock::now();

            double execution_time = duration<double, milli>(end_time - start_time).count();

            TestResult result;
            result.test_name = test_name;
            result.degree = degree;
            result.execution_time_ms = execution_time;

            evaluate_root_detection(real_roots_repeated, found, TOLERANCE, result);

            all_test_results["5. Экстремальные случаи"].push_back(result);

            if (degree <= 100) {
                print_detailed_results(test_name, coefficients, real_roots_repeated,
                    found, execution_time, result);
            }
        }
    }

    // Тест 2: Полиномы с очень близкими корнями
    {
        vector<pair<unsigned, double>> extreme_clusters = {
            {20, 1e-10}, {30, 1e-12}, {40, 1e-14}, {50, 1e-15}
        };

        for (const auto& cluster : extreme_clusters) {
            string test_name = "Сверхблизкие_корни_степень" + to_string(cluster.first) +
                "_радиус" + to_string(cluster.second);

            unsigned num_complex_pairs = 0;
            unsigned num_clusters = 2;
            vector<unsigned> cluster_counts = { cluster.first / 2, cluster.first / 2 };
            vector<double> cluster_radii = { cluster.second, cluster.second * 0.1 };
            vector<pair<unsigned, unsigned>> multiplicity_groups = {};
            double default_cluster_radius = 0.1;
            bool normalize_coeffs = true;
            uint64_t seed = 55555;

            vector<double> coefficients;
            vector<double> real_roots_repeated;
            vector<double> unique_real_roots;
            vector<unsigned> real_root_multiplicities;
            vector<complex<double>> complex_roots;

            generate_high_degree_polynomial(
                cluster.first, num_complex_pairs, num_clusters,
                cluster_counts, cluster_radii,
                multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
                coefficients, real_roots_repeated, unique_real_roots,
                real_root_multiplicities, complex_roots);

            double TOLERANCE = numeric_constants::adaptive_epsilon<double>(
                numeric_constants::EPSILON_SCALE_COARSE * cluster.second);

            auto start_time = high_resolution_clock::now();
            auto found = find_real_roots_ANewDsc(coefficients, TOLERANCE);
            auto end_time = high_resolution_clock::now();

            double execution_time = duration<double, milli>(end_time - start_time).count();

            TestResult result;
            result.test_name = test_name;
            result.degree = cluster.first;
            result.execution_time_ms = execution_time;

            evaluate_root_detection(real_roots_repeated, found, TOLERANCE, result);

            all_test_results["5. Экстремальные случаи"].push_back(result);

            print_detailed_results(test_name, coefficients, real_roots_repeated,
                found, execution_time, result);
        }
    }
}

// ==================== КАТЕГОРИЯ 6: РАЗНЫЕ ТИПЫ ДАННЫХ ====================
static void test_different_precisions() {
    print_subheader("КАТЕГОРИЯ 6: ТЕСТИРОВАНИЕ РАЗНЫХ ТИПОВ ДАННЫХ");

    struct PrecisionConfig {
        unsigned degree;
        unsigned num_complex_pairs;
        uint64_t seed;
    };

    // Создаем вектор конфигураций через push_back
    vector<PrecisionConfig> configs;

    PrecisionConfig cfg1; cfg1.degree = 10; cfg1.num_complex_pairs = 2; cfg1.seed = 111; configs.push_back(cfg1);
    PrecisionConfig cfg2; cfg2.degree = 15; cfg2.num_complex_pairs = 3; cfg2.seed = 222; configs.push_back(cfg2);
    PrecisionConfig cfg3; cfg3.degree = 20; cfg3.num_complex_pairs = 4; cfg3.seed = 333; configs.push_back(cfg3);
    PrecisionConfig cfg4; cfg4.degree = 25; cfg4.num_complex_pairs = 5; cfg4.seed = 444; configs.push_back(cfg4);
    PrecisionConfig cfg5; cfg5.degree = 30; cfg5.num_complex_pairs = 6; cfg5.seed = 555; configs.push_back(cfg5);

    for (const auto& config : configs) {
        // Тест с float
        {
            string test_name = "Float_точность_степень" + to_string(config.degree);

            unsigned num_clusters = 0;
            vector<unsigned> cluster_counts = {};
            vector<float> cluster_radii = {};
            vector<pair<unsigned, unsigned>> multiplicity_groups = {};
            float default_cluster_radius = 0.1f;
            bool normalize_coeffs = true;

            vector<float> coefficients;
            vector<float> real_roots_repeated;
            vector<float> unique_real_roots;
            vector<unsigned> real_root_multiplicities;
            vector<complex<float>> complex_roots;

            generate_high_degree_polynomial(
                config.degree, config.num_complex_pairs, num_clusters,
                cluster_counts, cluster_radii,
                multiplicity_groups, default_cluster_radius, normalize_coeffs, config.seed,
                coefficients, real_roots_repeated, unique_real_roots,
                real_root_multiplicities, complex_roots);

            float TOLERANCE = numeric_constants::adaptive_epsilon<float>(
                numeric_constants::EPSILON_SCALE_PRECISE);

            auto start_time = high_resolution_clock::now();
            auto found = find_real_roots_ANewDsc(coefficients, TOLERANCE);
            auto end_time = high_resolution_clock::now();

            double execution_time = duration<double, milli>(end_time - start_time).count();

            TestResult result;
            result.test_name = test_name;
            result.degree = config.degree;
            result.execution_time_ms = execution_time;

            vector<double> expected_roots_double(real_roots_repeated.begin(),
                real_roots_repeated.end());
            vector<pair<double, double>> found_intervals_double;
            for (const auto& interval : found) {
                found_intervals_double.push_back({ static_cast<double>(interval.first),
                                                  static_cast<double>(interval.second) });
            }

            result.expected_roots = expected_roots_double;
            result.expected_real_roots = static_cast<unsigned>(expected_roots_double.size());
            result.found_intervals = static_cast<unsigned>(found.size());
            result.all_roots_found = true;

            // Вычисляем ширины интервалов
            vector<double> widths;
            for (const auto& interval : found_intervals_double) {
                widths.push_back(interval.second - interval.first);
            }

            if (!widths.empty()) {
                double sum = 0.0;
                for (double w : widths) {
                    sum += w;
                }
                result.average_interval_width = sum / static_cast<double>(widths.size());
                result.max_interval_width = *max_element(widths.begin(), widths.end());
                result.min_interval_width = *min_element(widths.begin(), widths.end());
            }
            else {
                result.average_interval_width = 0.0;
                result.max_interval_width = 0.0;
                result.min_interval_width = 0.0;
            }

            all_test_results["6. Разные точности"].push_back(result);
        }

        // Тест с double
        {
            string test_name = "Double_точность_степень" + to_string(config.degree);

            unsigned num_clusters = 0;
            vector<unsigned> cluster_counts = {};
            vector<double> cluster_radii = {};
            vector<pair<unsigned, unsigned>> multiplicity_groups = {};
            double default_cluster_radius = 0.1;
            bool normalize_coeffs = true;

            vector<double> coefficients;
            vector<double> real_roots_repeated;
            vector<double> unique_real_roots;
            vector<unsigned> real_root_multiplicities;
            vector<complex<double>> complex_roots;

            generate_high_degree_polynomial(
                config.degree, config.num_complex_pairs, num_clusters,
                cluster_counts, cluster_radii,
                multiplicity_groups, default_cluster_radius, normalize_coeffs, config.seed,
                coefficients, real_roots_repeated, unique_real_roots,
                real_root_multiplicities, complex_roots);

            double TOLERANCE = numeric_constants::adaptive_epsilon<double>(
                numeric_constants::EPSILON_SCALE_PRECISE);

            auto start_time = high_resolution_clock::now();
            auto found = find_real_roots_ANewDsc(coefficients, TOLERANCE);
            auto end_time = high_resolution_clock::now();

            double execution_time = duration<double, milli>(end_time - start_time).count();

            TestResult result;
            result.test_name = test_name;
            result.degree = config.degree;
            result.execution_time_ms = execution_time;

            evaluate_root_detection(real_roots_repeated, found, TOLERANCE, result);

            all_test_results["6. Разные точности"].push_back(result);
        }
    }
}

// ==================== ОСНОВНАЯ ФУНКЦИЯ ТЕСТИРОВАНИЯ ====================
static void run_comprehensive_sagralov_tests() {
    print_header("МАСШТАБНОЕ ТЕСТИРОВАНИЕ МЕТОДА САГРАЛОВА (ANewDsc)");
    cout << "Тестирование полиномов от 5 до 100 степени\n";
    cout << "Всего будет выполнено несколько сотен тестов...\n\n";

    auto global_start_time = high_resolution_clock::now();

    // Запускаем все категории тестов
    test_simple_roots();
    test_multiple_roots();
    test_clustered_roots();
    test_mixed_cases();
    test_extreme_cases();
    test_different_precisions();

    auto global_end_time = high_resolution_clock::now();
    double total_execution_time = duration<double, milli>(global_end_time - global_start_time).count();

    // Выводим сводную статистику
    print_header("СВОДНАЯ СТАТИСТИКА");

    cout << "\nОбщее время выполнения всех тестов: " << fixed << setprecision(3)
        << total_execution_time << " мс (" << total_execution_time / 1000 << " с)\n\n";

    int total_tests = 0;
    int successful_tests = 0;

    for (const auto& category : all_test_results) {
        cout << category.first << ":\n";
        cout << "  Выполнено тестов: " << category.second.size() << "\n";

        int cat_success = 0;
        double cat_time = 0.0;
        for (const auto& result : category.second) {
            if (result.all_roots_found) cat_success++;
            cat_time += result.execution_time_ms;
        }

        cout << "  Успешных: " << cat_success << " ("
            << fixed << setprecision(1) << (100.0 * cat_success / static_cast<double>(category.second.size())) << "%)\n";
        cout << "  Среднее время: " << fixed << setprecision(3)
            << cat_time / static_cast<double>(category.second.size()) << " мс\n\n";

        total_tests += static_cast<int>(category.second.size());
        successful_tests += cat_success;
    }

    cout << "\nИТОГО:\n";
    cout << "  Всего тестов: " << total_tests << "\n";
    cout << "  Успешных тестов: " << successful_tests << "\n";
    cout << "  Процент успеха: " << fixed << setprecision(2)
        << (total_tests > 0 ? (100.0 * successful_tests / total_tests) : 0.0) << "%\n";

}

static void test_noisy_polynomials() {
    print_subheader("КАТЕГОРИЯ 7: УСТОЙЧИВОСТЬ К ШУМУ");

    vector<unsigned> degrees = { 10, 20, 30, 50 };
    vector<double> noise_levels = { 1e-10, 1e-8, 1e-6, 1e-4 };

    for (auto P : degrees) {
        for (auto noise : noise_levels) {
            string test_name = "Шум_deg" + to_string(P) + "_noise" + to_string(noise);

            vector<double> coefficients;
            vector<double> real_roots_repeated;
            vector<double> unique_real_roots;
            vector<unsigned> real_root_multiplicities;
            vector<complex<double>> complex_roots;

            generate_high_degree_polynomial(
                P, 0, 0, {}, {}, {},
                0.1, true, 12345,
                coefficients, real_roots_repeated,
                unique_real_roots, real_root_multiplicities, complex_roots);

            // Добавляем шум
            for (auto& c : coefficients) {
                c += noise * ((rand() % 2000 - 1000) / 1000.0);
            }

            double tol = numeric_constants::adaptive_epsilon<double>(
                numeric_constants::EPSILON_SCALE_COARSE);

            auto start = high_resolution_clock::now();
            auto found = find_real_roots_ANewDsc(coefficients, tol);
            auto end = high_resolution_clock::now();

            TestResult result;
            result.test_name = test_name;
            result.degree = P;
            result.execution_time_ms =
                duration<double, milli>(end - start).count();

            evaluate_root_detection(real_roots_repeated, found, tol, result);

            all_test_results["7. Шум в коэффициентах"].push_back(result);
        }
    }
}

static void test_ill_conditioned() {
    print_subheader("КАТЕГОРИЯ 8: ПЛОХО ОБУСЛОВЛЕННЫЕ ПОЛИНОМЫ");

    vector<unsigned> degrees = { 10, 20, 30, 40 };

    for (auto P : degrees) {
        string test_name = "Ill_conditioned_deg" + to_string(P);

        vector<double> roots;
        for (unsigned i = 0; i < P; ++i) {
            roots.push_back(1.0 + i * 1e-6); // почти совпадающие
        }

        // генерим через обычный генератор
        vector<double> coefficients;
        vector<double> real_roots_repeated = roots;
        vector<double> unique_real_roots = roots;
        vector<unsigned> mult(roots.size(), 1);
        vector<complex<double>> complex_roots;

        generate_high_degree_polynomial(
            P, 0, 0, {}, {}, {},
            0.1, true, 9999,
            coefficients, real_roots_repeated,
            unique_real_roots, mult, complex_roots);

        double tol = numeric_constants::adaptive_epsilon<double>(
            numeric_constants::EPSILON_SCALE_COARSE);

        auto start = high_resolution_clock::now();
        auto found = find_real_roots_ANewDsc(coefficients, tol);
        auto end = high_resolution_clock::now();

        TestResult result;
        result.test_name = test_name;
        result.degree = P;
        result.execution_time_ms =
            duration<double, milli>(end - start).count();

        evaluate_root_detection(real_roots_repeated, found, tol, result);

        all_test_results["8. Плохая обусловленность"].push_back(result);

        print_detailed_results(test_name, coefficients,
            real_roots_repeated, found,
            result.execution_time_ms, result);
    }
}

static void test_random_stress() {
    print_subheader("КАТЕГОРИЯ 9: СТРЕСС-ТЕСТ (100+ ПОЛИНОМОВ)");

    int N = 120;

    for (int i = 0; i < N; ++i) {
        unsigned P = 5 + rand() % 95;

        string test_name = "Stress_" + to_string(i) + "_deg" + to_string(P);

        vector<double> coefficients;
        vector<double> real_roots_repeated;
        vector<double> unique_real_roots;
        vector<unsigned> real_root_multiplicities;
        vector<complex<double>> complex_roots;

        generate_high_degree_polynomial(
            P,
            rand() % (P / 3),
            rand() % 3,
            {},
            {},
            {},
            0.1,
            true,
            1000 + i,
            coefficients,
            real_roots_repeated,
            unique_real_roots,
            real_root_multiplicities,
            complex_roots);

        double tol = numeric_constants::adaptive_epsilon<double>(
            numeric_constants::EPSILON_SCALE_STANDARD);

        auto start = high_resolution_clock::now();
        auto found = find_real_roots_ANewDsc(coefficients, tol);
        auto end = high_resolution_clock::now();

        TestResult result;
        result.test_name = test_name;
        result.degree = P;
        result.execution_time_ms =
            duration<double, milli>(end - start).count();

        evaluate_root_detection(real_roots_repeated, found, tol, result);

        all_test_results["9. Стресс-тест"].push_back(result);
    }
}

static void test_tolerance_sensitivity() {
    print_subheader("КАТЕГОРИЯ 10: ЧУВСТВИТЕЛЬНОСТЬ К EPSILON");

    vector<unsigned> degrees = { 10, 20, 30 };
    vector<double> scales = { 0.1, 1.0, 10.0, 100.0 };

    for (auto P : degrees) {
        for (auto scale : scales) {
            string test_name = "Tol_deg" + to_string(P) + "_scale" + to_string(scale);

            vector<double> coefficients;
            vector<double> real_roots_repeated;
            vector<double> unique_real_roots;
            vector<unsigned> real_root_multiplicities;
            vector<complex<double>> complex_roots;

            generate_high_degree_polynomial(
                P, 0, 0, {}, {}, {},
                0.1, true, 777,
                coefficients, real_roots_repeated,
                unique_real_roots, real_root_multiplicities, complex_roots);

            double tol = numeric_constants::adaptive_epsilon<double>(
                numeric_constants::EPSILON_SCALE_STANDARD * scale);

            auto start = high_resolution_clock::now();
            auto found = find_real_roots_ANewDsc(coefficients, tol);
            auto end = high_resolution_clock::now();

            TestResult result;
            result.test_name = test_name;
            result.degree = P;
            result.execution_time_ms =
                duration<double, milli>(end - start).count();

            evaluate_root_detection(real_roots_repeated, found, tol, result);

            all_test_results["10. Sensitivity epsilon"].push_back(result);
        }
    }
}

int main() {
    setlocale(LC_ALL, "ru_RU.UTF-8");

    try {
        run_comprehensive_sagralov_tests();
        test_noisy_polynomials();
        test_ill_conditioned();
        test_random_stress();
        test_tolerance_sensitivity();

        cout << "\n\nТестирование успешно завершено!\n";
        return EXIT_SUCCESS;
    }
    catch (const exception& e) {
        cerr << "\nКРИТИЧЕСКАЯ ОШИБКА: " << e.what() << endl;
        return EXIT_FAILURE;
    }
    catch (...) {
        cerr << "\nНЕИЗВЕСТНАЯ ОШИБКА!" << endl;
        return EXIT_FAILURE;
    }
}