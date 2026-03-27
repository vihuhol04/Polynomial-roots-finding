// Файл для тестирования метода Сагралова – поиск комплексных корней
// Реализация: Павлова Анастасия, КМБО-01-22 (расширенная версия)

#include "Sagralov_complex_roots.h"
#include "generate_high_degree_polynomial.h"

#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;

// =========================================
// ВЫВОД ДИСКОВ
// =========================================
template<typename fp_t>
void print_disks(const vector<Disk<fp_t>>& disks) {

    cout << "Найдено дисков: " << disks.size() << endl;

    fp_t total_radius = fp_t(0);
    fp_t max_radius = fp_t(0);
    int total_roots = 0;

    for (const auto& d : disks) {
        cout << "Центр: " << d.center
            << " | Радиус: " << d.radius
            << " | Корней: " << d.num_roots << endl;

        total_radius += d.radius;

        if (d.radius > max_radius)
            max_radius = d.radius;

        total_roots += d.num_roots;
    }

    cout << "\nСУММАРНО:\n";
    cout << "  Всего корней: " << total_roots << endl;

    if (!disks.empty()) {
        cout << "  Средний радиус: "
            << total_radius / static_cast<fp_t>(disks.size()) << endl;
    }
    else {
        cout << "  Средний радиус: 0\n";
    }

    cout << "  Макс радиус: " << max_radius << endl;
}

// =========================================
// УНИВЕРСАЛЬНЫЙ ТЕСТ
// =========================================
template<typename fp_t>
static void run_single_test(
    const string& name,
    unsigned P,
    unsigned num_complex_pairs,
    unsigned num_clusters,
    const vector<unsigned>& cluster_counts,
    const vector<fp_t>& cluster_radii,
    const vector<pair<unsigned, unsigned>>& multiplicity_groups,
    fp_t default_cluster_radius,
    uint64_t seed) {

    cout << "\n=========================================\n";
    cout << name << endl;
    cout << "Степень: " << P << endl;
    cout << "=========================================\n";

    vector<fp_t> coefficients;
    vector<fp_t> real_roots_repeated;
    vector<fp_t> unique_real_roots;
    vector<unsigned> real_root_multiplicities;
    vector<complex<fp_t>> complex_roots;

    generate_high_degree_polynomial(
        P,
        num_complex_pairs,
        num_clusters,
        cluster_counts,
        cluster_radii,
        multiplicity_groups,
        default_cluster_radius,
        true,
        seed,
        coefficients,
        real_roots_repeated,
        unique_real_roots,
        real_root_multiplicities,
        complex_roots);

    cout << "\nПолином:\n";
    print_polynomial(coefficients);

    auto disks = CIsolate<fp_t, fp_t>(coefficients);

    cout << "\nРЕЗУЛЬТАТ:\n";
    print_disks(disks);
}

// =========================================
// ОСНОВНОЙ РАННЕР
// =========================================
static void run_sagralov_tests() {

    cout << "\n=========================================\n";
    cout << "ТЕСТИРОВАНИЕ МЕТОДА САГРАЛОВА (КОМПЛЕКСНЫЕ КОРНИ)\n";
    cout << "=========================================\n";

    vector<unsigned> degrees = { 5, 10, 20, 30, 50 };

    // =========================================
    // DOUBLE
    // =========================================
    cout << "\n===== DOUBLE =====\n";

    for (unsigned P : degrees) {

        uint64_t base_seed = 12345ULL + static_cast<uint64_t>(P) * 10ULL;

        // Простые
        run_single_test<double>(
            "Простые комплексные корни",
            P,
            P / 2,
            0,
            {},
            {},
            {},
            0.1,
            base_seed + 1ULL);

        // Кластеры
        run_single_test<double>(
            "Кластеризованные корни",
            P,
            P / 3,
            2,
            { P / 5, P / 6 },
            { 1e-3, 1e-2 },
            {},
            0.1,
            base_seed + 2ULL);

        // Кратные
        run_single_test<double>(
            "Кратные комплексные корни",
            P,
            P / 3,
            0,
            {},
            {},
            { {2, P / 6}, {3, P / 8} },
            0.1,
            base_seed + 3ULL);

        // Смешанные
        run_single_test<double>(
            "Смешанный случай",
            P,
            P / 4,
            2,
            { P / 6, P / 7 },
            { 1e-3, 5e-3 },
            { {2, P / 7} },
            0.1,
            base_seed + 4ULL);
    }

    // =========================================
    // HIGH PRECISION
    // =========================================
    cout << "\n===== FLOAT_PRECISION =====\n";

    vector<unsigned> hp_degrees = { 10, 20, 30 };

    for (unsigned P : hp_degrees) {

        uint64_t base_seed = 54321ULL + static_cast<uint64_t>(P) * 10ULL;

        run_single_test<float_precision>(
            "HP: Простые",
            P,
            P / 2,
            0,
            {},
            {},
            {},
            static_cast<float_precision>(0.1),
            base_seed + 1ULL);

        run_single_test<float_precision>(
            "HP: Кластеры",
            P,
            P / 3,
            2,
            { P / 5, P / 6 },
            { static_cast<float_precision>(1e-4),
             static_cast<float_precision>(1e-3) },
            {},
            static_cast<float_precision>(0.1),
            base_seed + 2ULL);

        run_single_test<float_precision>(
            "HP: Смешанный",
            P,
            P / 4,
            2,
            { P / 6, P / 7 },
            { static_cast<float_precision>(1e-4),
             static_cast<float_precision>(5e-4) },
            { {2, P / 6} },
            static_cast<float_precision>(0.1),
            base_seed + 3ULL);
    }

    cout << "\n=========================================\n";
    cout << "ВСЕ ТЕСТЫ ЗАВЕРШЕНЫ\n";
    cout << "=========================================\n";
}

// =========================================
// MAIN
// =========================================
int main() {
    setlocale(LC_ALL, "ru_RU.UTF-8");

    try {
        run_sagralov_tests();
        return EXIT_SUCCESS;
    }
    catch (const exception& e) {
        cerr << "ОШИБКА: " << e.what() << endl;
        return EXIT_FAILURE;
    }
}