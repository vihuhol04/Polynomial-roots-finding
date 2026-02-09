// Файл для тестирования метода Сагралова (ANewDsc) – поиск вещественных корней
// Реализация: Павлова Анастасия, КМБО-01-22

#include "NumericConstants.h"
#include "Sagralov_ANewDsc.h"
#include "generate_high_degree_polynomial.h"
#include <complex>
#include <iostream>
#include <vector>

using namespace std;

template <typename fp_t>
void print_real_roots_ground_truth(
    const vector<fp_t> &real_roots_repeated,
    const vector<fp_t> &unique_real_roots,
    const vector<unsigned> &real_root_multiplicities) {
  cout << "Исходные вещественные корни (с кратностями): ";
  for (auto r : real_roots_repeated)
    cout << r << " ";
  cout << endl;

  cout << "Уникальные вещественные корни: ";
  for (size_t i = 0; i < unique_real_roots.size(); ++i)
    cout << unique_real_roots[i] << " (мультипл. "
         << real_root_multiplicities[i] << ")  ";
  cout << endl;
}

void run_sagralov_tests() {
  cout << endl << "=== ТЕСТИРОВАНИЕ МЕТОДА САГРАЛОВА (ANewDsc) ===" << endl;

  // ТЕСТ 1
  {
    cout << "\nТЕСТ 1: Простые вещественные корни (степень 10)" << endl;

    unsigned P = 10;
    unsigned num_complex_pairs = 0;
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
        real_root_multiplicities, complex_roots);

    double TOLERANCE = numeric_constants::adaptive_epsilon<double>(
        numeric_constants::EPSILON_SCALE_PRECISE);

    cout << "\nСГЕНИРИРОВАННЫЙ ПОЛМНОМ: \n";
    print_polynomial(coefficients);
    print_real_roots_ground_truth(real_roots_repeated, unique_real_roots,
                                  real_root_multiplicities);

    auto found = find_real_roots_ANewDsc(coefficients, TOLERANCE);
    cout << "Интервалы (Сагралов):\n";
    for (auto &i : found)
      cout << "  [" << i.first << ", " << i.second << "]\n";
  }

  // ТЕСТ 2
  {
    cout << "\nТЕСТ 2: Кратные корни float_precision" << endl;

    unsigned P = 12;
    unsigned num_complex_pairs = 0;
    unsigned num_clusters = 0;
    vector<unsigned> cluster_counts = {};
    vector<float_precision> cluster_radii = {};
    vector<pair<unsigned, unsigned>> multiplicity_groups = {{2, 3}, {1, 2}};
    float_precision default_cluster_radius = 0.1;
    bool normalize_coeffs = true;
    uint64_t seed = 12345;

    vector<float_precision> coefficients;
    vector<float_precision> real_roots_repeated;
    vector<float_precision> unique_real_roots;
    vector<unsigned> real_root_multiplicities;
    vector<complex<float_precision>> complex_roots;

    generate_high_degree_polynomial(
        P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
        multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
        coefficients, real_roots_repeated, unique_real_roots,
        real_root_multiplicities, complex_roots);

    float_precision TOLERANCE =
        numeric_constants::adaptive_epsilon<float_precision>(
            numeric_constants::EPSILON_SCALE_PRECISE);

    cout << "\nСГЕНИРИРОВАННЫЙ ПОЛМНОМ: \n";
    print_polynomial(coefficients);
    print_real_roots_ground_truth(real_roots_repeated, unique_real_roots,
                                  real_root_multiplicities);

    auto found = find_real_roots_ANewDsc(coefficients, TOLERANCE);
    cout << "Интервалы (Сагралов):\n";
    for (auto &i : found)
      cout << "  [" << i.first << ", " << i.second << "]\n";
  }

  // ТЕСТ 3
  {
    cout << "\nТЕСТ 3: Кластеры корней float_precision" << endl;

    unsigned P = 15;
    unsigned num_complex_pairs = 0;
    unsigned num_clusters = 0;
    vector<unsigned> cluster_counts = {3, 4};
    vector<float_precision> cluster_radii = {0.0001, 0.001};
    vector<pair<unsigned, unsigned>> multiplicity_groups = {};
    float_precision default_cluster_radius = 0.1;
    bool normalize_coeffs = true;
    uint64_t seed = 9876;

    vector<float_precision> coefficients;
    vector<float_precision> real_roots_repeated;
    vector<float_precision> unique_real_roots;
    vector<unsigned> real_root_multiplicities;
    vector<complex<float_precision>> complex_roots;

    generate_high_degree_polynomial(
        P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
        multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
        coefficients, real_roots_repeated, unique_real_roots,
        real_root_multiplicities, complex_roots);

    float_precision TOLERANCE =
        numeric_constants::adaptive_epsilon<float_precision>(
            numeric_constants::EPSILON_SCALE_COARSE * 5);

    cout << "\nСГЕНИРИРОВАННЫЙ ПОЛМНОМ: \n";
    print_polynomial(coefficients);
    print_real_roots_ground_truth(real_roots_repeated, unique_real_roots,
                                  real_root_multiplicities);

    auto found = find_real_roots_ANewDsc(coefficients, TOLERANCE);
    cout << "Интервалы (Сагралов):\n";
    for (auto &i : found)
      cout << "  [" << i.first << ", " << i.second << "]\n";
  }

  // ТЕСТ 4
  {
    cout << "\nТЕСТ 4: Высокая степень (30)" << endl;

    unsigned P = 30;
    unsigned num_complex_pairs = 0;
    unsigned num_clusters = 0;
    vector<unsigned> cluster_counts = {};
    vector<float_precision> cluster_radii = {};
    vector<pair<unsigned, unsigned>> multiplicity_groups = {};
    float_precision default_cluster_radius = 0.1;
    bool normalize_coeffs = true;
    uint64_t seed = 4567;

    vector<float_precision> coefficients;
    vector<float_precision> real_roots_repeated;
    vector<float_precision> unique_real_roots;
    vector<unsigned> real_root_multiplicities;
    vector<complex<float_precision>> complex_roots;

    generate_high_degree_polynomial(
        P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
        multiplicity_groups, default_cluster_radius, normalize_coeffs, seed,
        coefficients, real_roots_repeated, unique_real_roots,
        real_root_multiplicities, complex_roots);

    float_precision TOLERANCE =
        numeric_constants::adaptive_epsilon<float_precision>(
            numeric_constants::EPSILON_SCALE_COARSE * 5);

    cout << "\nСГЕНИРИРОВАННЫЙ ПОЛМНОМ: \n";
    print_polynomial(coefficients);
    print_real_roots_ground_truth(real_roots_repeated, unique_real_roots,
                                  real_root_multiplicities);

    auto found = find_real_roots_ANewDsc(coefficients, TOLERANCE);
    cout << "Интервалы (Сагралов):\n";
    for (auto &i : found)
      cout << "  [" << i.first << ", " << i.second << "]\n";
  }
}

int main() {
    setlocale(LC_ALL, "ru_RU.UTF-8");
  try {
    run_sagralov_tests();
    return EXIT_SUCCESS;
  } catch (const exception &e) {
    cerr << "ОШИБКА: " << e.what() << endl;
    return EXIT_FAILURE;
  }
}