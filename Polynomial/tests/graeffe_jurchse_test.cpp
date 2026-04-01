// Тестирование метода Юркша (нормализованный Грефе)
// Комплексные тесты: степени 5-100, вещественные/комплексные/кластеризованные/кратные корни
// Типы: double, long double, float_precision
// Результаты сохраняются в .txt файлы

// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#ifdef _WIN32
#include <windows.h>
#endif

#include <algorithm>
#include <chrono>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Graeffe_Jurchse.h"
#include "Graeffe_real_roots.h"
#include "NumericConstants.h"
#include "generate_high_degree_polynomial.h"

using namespace std;

// =====================================================================
// Утилиты для отчёта
// =====================================================================

template <typename T>
double to_double(const T &val) {
	return static_cast<double>(val);
}

template <typename T>
string type_name() {
	if constexpr (is_same_v<T, double>)
		return "double";
	else if constexpr (is_same_v<T, long double>)
		return "long double";
	else
		return "float_precision";
}

template <typename T>
void sort_descending(vector<T> &v) {
	sort(v.begin(), v.end(), [](const T &a, const T &b) { return a > b; });
}

// Вычисление модулей из списка корней (вещественных + комплексных)
template <typename T>
vector<T> expected_moduli(const vector<T> &real_roots,
                          const vector<complex<T>> &complex_roots) {
	vector<T> moduli;
	for (auto &r : real_roots)
		moduli.push_back(abs_val(r));
	for (auto &c : complex_roots)
		moduli.push_back(
		    sqrt_val(c.real() * c.real() + c.imag() * c.imag()));
	sort_descending(moduli);
	return moduli;
}

// Оценка ошибки: для каждого ожидаемого модуля ищем ближайший найденный
// Считаем min(относительная, абсолютная) ошибку
template <typename T>
double compute_max_relative_error(const vector<T> &computed,
                                  const vector<T> &expected, T eps) {
	int ne = int(expected.size());
	int nc = int(computed.size());
	if (ne == 0)
		return 0.0;
	double max_err = 0.0;
	vector<bool> used(nc, false);
	for (int i = 0; i < ne; ++i) {
		double e = to_double(expected[i]);
		double best = 1e30;
		int best_j = -1;
		for (int j = 0; j < nc; ++j) {
			if (used[j])
				continue;
			double c = to_double(computed[j]);
			double diff = abs(c - e);
			if (diff < best) {
				best = diff;
				best_j = j;
			}
		}
		if (best_j >= 0) {
			used[best_j] = true;
			double denom = max(abs(e), 0.01);
			double rel = best / denom;
			max_err = max(max_err, rel);
		} else {
			max_err = max(max_err, abs(e) > 0.01 ? 1.0 : 0.01);
		}
	}
	return max_err;
}

// Шаблонный тест для одного полинома (Юркша)

template <typename T>
struct TestResult {
	string test_name;
	string data_type;
	int degree;
	int roots_found;
	int roots_expected;
	double max_relative_error;
	double elapsed_ms;
	bool passed;
};

template <typename T>
TestResult<T> run_single_jurchse_test(
    const string &test_name, unsigned P, unsigned num_complex_pairs,
    unsigned num_clusters, const vector<unsigned> &cluster_counts,
    const vector<T> &cluster_radii,
    const vector<pair<unsigned, unsigned>> &multiplicity_groups,
    T default_cluster_radius, uint64_t seed, T tolerance, int maxIter = 40,
    ostream &log = cout) {

	TestResult<T> result;
	result.test_name = test_name;
	result.data_type = type_name<T>();
	result.degree = P;

	vector<T> coefficients;
	vector<T> real_roots_repeated;
	vector<T> unique_real_roots;
	vector<unsigned> real_root_multiplicities;
	vector<complex<T>> complex_roots;

	generate_high_degree_polynomial(
	    P, num_complex_pairs, num_clusters, cluster_counts, cluster_radii,
	    multiplicity_groups, default_cluster_radius, true, seed, coefficients,
	    real_roots_repeated, unique_real_roots, real_root_multiplicities,
	    complex_roots);

	// ascending -> descending для Грефе
	vector<T> coeffs_desc(coefficients.rbegin(), coefficients.rend());

	T epsilon =
	    numeric_constants::adaptive_epsilon<T>(numeric_constants::EPSILON_SCALE_PRECISE);

	auto start = chrono::high_resolution_clock::now();
	vector<T> moduli = find_moduli_by_Jurchse(coeffs_desc, epsilon, maxIter);
	auto end = chrono::high_resolution_clock::now();

	result.elapsed_ms =
	    chrono::duration<double, milli>(end - start).count();

	vector<T> exp_mod = expected_moduli(real_roots_repeated, complex_roots);
	result.roots_expected = int(exp_mod.size());
	result.roots_found = int(moduli.size());

	sort_descending(moduli);

	result.max_relative_error =
	    compute_max_relative_error(moduli, exp_mod, tolerance);
	result.passed = (result.max_relative_error < to_double(tolerance)) &&
	                (result.roots_found > 0);

	log << "  " << test_name << " [" << type_name<T>() << "] deg=" << P
	    << "  корней: " << result.roots_found << "/" << result.roots_expected
	    << "  макс.отн.ошибка=" << scientific << setprecision(3)
	    << result.max_relative_error << "  время=" << fixed << setprecision(1)
	    << result.elapsed_ms << "мс  "
	    << (result.passed ? "PASSED" : "FAILED") << endl;

	return result;
}

// Набор тестов для конкретного типа данных

template <typename T>
vector<TestResult<T>> run_all_jurchse_tests(ostream &log) {
	vector<TestResult<T>> results;
	T tol;
	if constexpr (is_same_v<T, float_precision>)
		tol = T("0.5");
	else
		tol = T(0.5);

	int maxIter = 30;
	if constexpr (is_same_v<T, float_precision>)
		maxIter = 15;

	log << "\n========== МЕТОД ЮРКША: тесты для " << type_name<T>()
	    << " ==========\n";

	// --- 1. Только вещественные корни ---
	for (unsigned deg : {5u, 10u, 20u, 30u, 50u}) {
		if constexpr (is_same_v<T, float_precision>) {
			if (deg > 20) continue;
		}

		T dcr;
		if constexpr (is_same_v<T, float_precision>)
			dcr = T("0.1");
		else
			dcr = T(0.1);

		vector<T> cr;
		results.push_back(run_single_jurchse_test<T>(
		    "Веществ.корни", deg, 0, 0, {}, cr, {}, dcr, 12345 + deg, tol,
		    maxIter, log));
	}

	// --- 2. С комплексными парами ---
	for (unsigned deg : {10u, 20u, 30u}) {
		if constexpr (is_same_v<T, float_precision>) {
			if (deg > 20) continue;
		}
		unsigned cp = deg / 4;

		T dcr;
		if constexpr (is_same_v<T, float_precision>)
			dcr = T("0.1");
		else
			dcr = T(0.1);

		vector<T> cr;
		results.push_back(run_single_jurchse_test<T>(
		    "Компл.пары", deg, cp, 0, {}, cr, {}, dcr, 54321 + deg, tol,
		    maxIter, log));
	}

	// --- 3. Кластеризованные корни ---
	for (unsigned deg : {10u, 20u, 30u}) {
		if constexpr (is_same_v<T, float_precision>) {
			if (deg > 20) continue;
		}
		unsigned nc = 2;
		vector<unsigned> cc = {3, 3};

		T r1, r2, dcr;
		if constexpr (is_same_v<T, float_precision>) {
			r1 = T("0.01");
			r2 = T("0.01");
			dcr = T("0.1");
		} else {
			r1 = T(0.01);
			r2 = T(0.01);
			dcr = T(0.1);
		}
		vector<T> cr = {r1, r2};

		results.push_back(run_single_jurchse_test<T>(
		    "Кластериз.", deg, 0, nc, cc, cr, {}, dcr, 99999 + deg, tol,
		    maxIter, log));
	}

	// --- 4. Кратные корни ---
	for (unsigned deg : {10u, 20u}) {
		if constexpr (is_same_v<T, float_precision>) {
			if (deg > 15) continue;
		}
		vector<pair<unsigned, unsigned>> mg = {{2, 2}, {3, 1}};

		T dcr;
		if constexpr (is_same_v<T, float_precision>)
			dcr = T("0.1");
		else
			dcr = T(0.1);

		vector<T> cr;
		results.push_back(run_single_jurchse_test<T>(
		    "Кратные", deg, 0, 0, {}, cr, mg, dcr, 77777 + deg, tol, maxIter,
		    log));
	}

	// --- 5. Смешанный: компл. + кластеры + кратные ---
	for (unsigned deg : {20u, 30u, 50u}) {
		if constexpr (is_same_v<T, float_precision>) {
			if (deg > 20) continue;
		}
		unsigned cp = 2;
		unsigned nc = 1;
		vector<unsigned> cc = {2};
		vector<pair<unsigned, unsigned>> mg = {{2, 1}};

		T r1, dcr;
		if constexpr (is_same_v<T, float_precision>) {
			r1 = T("0.01");
			dcr = T("0.1");
		} else {
			r1 = T(0.01);
			dcr = T(0.1);
		}
		vector<T> cr = {r1};

		results.push_back(run_single_jurchse_test<T>(
		    "Смешанный", deg, cp, nc, cc, cr, mg, dcr, 11111 + deg, tol,
		    maxIter, log));
	}

	// --- 6. Высокие степени (только double/long double) ---
	if constexpr (!is_same_v<T, float_precision>) {
		for (unsigned deg : {60u, 80u, 100u}) {
			T dcr = T(0.1);
			vector<T> cr;
			unsigned cp = deg / 5;
			results.push_back(run_single_jurchse_test<T>(
			    "Высок.степ.", deg, cp, 0, {}, cr, {}, dcr, 33333 + deg, tol,
			    maxIter, log));
		}
	}

	return results;
}

// Сохранение результатов в .txt

template <typename T>
void save_results_to_file(const vector<TestResult<T>> &results,
                          const string &filename) {
	ofstream out(filename);
	if (!out.is_open()) {
		cerr << "Ошибка записи файла: " << filename << endl;
		return;
	}

	out << "======================================================================"
	       "====\n";
	out << "РЕЗУЛЬТАТЫ ТЕСТИРОВАНИЯ: МЕТОД ЮРКША (Нормализованный Грефе)\n";
	out << "Тип данных: " << type_name<T>() << "\n";
	out << "======================================================================"
	       "====\n\n";

	out << left << setw(22) << "Тест" << setw(8) << "Степ." << setw(14)
	    << "Корней(найд)" << setw(14) << "Корней(ожид)" << setw(18)
	    << "Макс.отн.ошиб." << setw(12) << "Время(мс)" << setw(10)
	    << "Статус"
	    << "\n";
	out << string(98, '-') << "\n";

	int passed = 0, total = 0;
	for (auto &r : results) {
		out << left << setw(22) << r.test_name << setw(8) << r.degree
		    << setw(14) << r.roots_found << setw(14) << r.roots_expected
		    << setw(18) << scientific << setprecision(3)
		    << r.max_relative_error << setw(12) << fixed << setprecision(1)
		    << r.elapsed_ms << setw(10) << (r.passed ? "PASSED" : "FAILED")
		    << "\n";
		++total;
		if (r.passed)
			++passed;
	}

	out << "\n" << string(98, '=') << "\n";
	out << "Итого: " << passed << " из " << total << " тестов пройдено.\n";
	out.close();

	cout << "Результаты сохранены в: " << filename << endl;
}

// main

int main() {
#ifdef _WIN32
	SetConsoleOutputCP(65001);
	SetConsoleCP(65001);
#endif

	cout << "========================================\n";
	cout << " ТЕСТИРОВАНИЕ МЕТОДА ЮРКША\n";
	cout << " (Нормализованная модификация Грефе)\n";
	cout << "========================================\n\n";

	// double
	{
		stringstream log;
		auto results = run_all_jurchse_tests<double>(log);
		cout << log.str();
		save_results_to_file(results, "jurchse_results_double.txt");
	}

	// long double
	{
		stringstream log;
		auto results = run_all_jurchse_tests<long double>(log);
		cout << log.str();
		save_results_to_file(results, "jurchse_results_long_double.txt");
	}

	// float_precision
	{
		stringstream log;
		auto results = run_all_jurchse_tests<float_precision>(log);
		cout << log.str();
		save_results_to_file(results, "jurchse_results_float_precision.txt");
	}

	cout << "\nВсе тесты метода Юркша завершены.\n";
	return 0;
}
