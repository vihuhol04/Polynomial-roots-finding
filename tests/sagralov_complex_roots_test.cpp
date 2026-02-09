// Файл для тестирования метода Сагралова – поиск комплексных корней
// Реализация: Павлова Анастасия, КМБО-01-22

#include "Sagralov_complex_roots.h"
#include "generate_high_degree_polynomial.h"
#include <complex>
#include <iostream>
#include <vector>

using namespace std;

void run_sagralov_tests() {
  cout << endl
       << "ТЕСТИРОВАНИЕ МЕТОДА САГРАЛОВА ДЛЯ КОМПЛЕКСНЫХ КОРНЕЙ" << endl;

  // ТЕСТ 1
  {
    cout << "\nТЕСТ 1: Простые комплексные корни (степень 10)" << endl;

    unsigned P = 10;
    unsigned num_complex_pairs = 5;
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

    cout << "\nСГЕНИРИРОВАННЫЙ ПОЛМНОМ: \n";
    print_polynomial(coefficients);

    auto disks = CIsolate<double, long double>(coefficients);
    cout << "Интервалы (Сагралов):\n";
    for (auto &i : disks)
      cout << "Центр: " << i.center
           << "      Количество корней в диске: " << i.num_roots
           << "    Радиус диска: " << i.radius << "\n";
  }

  // ТЕСТ 2
  // ТЕСТ 2
  {
    cout << "\nТЕСТ 2: Высокая степень (30) float_precision" << endl;

    unsigned P = 30;
    unsigned num_complex_pairs = 15;
    unsigned num_clusters = 0;
    vector<unsigned> cluster_counts = {};
    vector<float_precision> cluster_radii = {};
    vector<pair<unsigned, unsigned>> multiplicity_groups = {};
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

    float_precision TOLERANCE = 1e-6;

    cout << "\nСГЕНИРИРОВАННЫЙ ПОЛМНОМ: \n";
    print_polynomial(coefficients);

    auto disks = CIsolate<float_precision, float_precision>(coefficients);
    cout << "Интервалы (Сагралов):\n";
    for (auto &i : disks)
      cout << "Центр: " << i.center
           << "      Количество корней в диске: " << i.num_roots
           << "    Радиус диска: " << i.radius << "\n";
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