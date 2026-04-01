// Тестирование модификации Сан-Хуана метода Греффе
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

#include "Graeffe_SanJuan.h"
#include "Graeffe_real_roots.h"
#include "NumericConstants.h"
#include "generate_high_degree_polynomial.h"

using namespace std;

// =====================================================================
// Утилиты
// =====================================================================

template <typename T>
double to_double_val(const T &val) {
	return static_cast<double>(val);
}

template <typename T>
string type_name_sj() {
	if constexpr (is_same_v<T, double>)
		return "double";
	else if constexpr (is_same_v<T, long double>)
		return "long double";
	else
		return "float_precision";
}

template <typename T>
void sort_desc(vector<T> &v) {
	sort(v.begin(), v.end(), [](const T &a, const T &b) { return a > b; });
}

template <typename T>
vector<T> get_expected_moduli(const vector<T> &real_roots,
                              const vector<complex<T>> &complex_roots) {
	vector<T> moduli;
	for (auto &r : real_roots)
		moduli.push_back(abs_val(r));
	for (auto &c : complex_roots)
		moduli.push_back(
		    sqrt_val(c.real() * c.real() + c.imag() * c.imag()));
	sort_desc(moduli);
	return moduli;
}

template <typename T>
double max_relative_error(const vector<T> &computed,
                          const vector<T> &expected) {
	int ne = int(expected.size());
	int nc = int(computed.size());
	if (ne == 0)
		return 0.0;
	double max_err = 0.0;
	vector<bool> used(nc, false);
	for (int i = 0; i < ne; ++i) {
		double e = to_double_val(expected[i]);
		double best = 1e30;
		int best_j = -1;
		for (int j = 0; j < nc; ++j) {
			if (used[j])
				continue;
			double c = to_double_val(computed[j]);
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

// =====================================================================
// Структура результата теста
// =====================================================================

struct SJTestResult {
	string test_name;
	string data_type;
	int degree;
	int m_dominant;
	int roots_found;
	int roots_expected;
	double max_rel_err;
	double elapsed_ms;
	bool passed;
};

// =====================================================================
// Запуск одного теста Сан-Хуана
// =====================================================================

template <typename T>
SJTestResult run_single_sanjuan_test(
    const string &test_name, unsigned P, unsigned num_complex_pairs,
    unsigned num_clusters, const vector<unsigned> &cluster_counts,
    const vector<T> &cluster_radii,
    const vector<pair<unsigned, unsigned>> &multiplicity_groups,
    T default_cluster_radius, uint64_t seed, T tolerance,
    int graeffe_iters = 10, ostream &log = cout) {

	SJTestResult result;
	result.test_name = test_name;
	result.data_type = type_name_sj<T>();
	result.degree = P;
	result.m_dominant = 0;

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

	// ascending -> descending
	vector<T> coeffs_desc(coefficients.rbegin(), coefficients.rend());

	T epsilon =
	    numeric_constants::adaptive_epsilon<T>(numeric_constants::EPSILON_SCALE_PRECISE);

	// Тест find_all_moduli_by_SanJuan (все корни)
	auto start = chrono::high_resolution_clock::now();
	vector<T> moduli =
	    find_all_moduli_by_SanJuan(coeffs_desc, graeffe_iters, epsilon);
	auto end = chrono::high_resolution_clock::now();

	result.elapsed_ms =
	    chrono::duration<double, milli>(end - start).count();

	vector<T> exp_mod = get_expected_moduli(real_roots_repeated, complex_roots);
	result.roots_expected = int(exp_mod.size());
	result.roots_found = int(moduli.size());
	result.m_dominant = result.roots_found;

	sort_desc(moduli);

	result.max_rel_err = max_relative_error(moduli, exp_mod);
	result.passed =
	    (result.max_rel_err < to_double_val(tolerance)) && (result.roots_found > 0);

	log << "  " << test_name << " [" << type_name_sj<T>() << "] deg=" << P
	    << "  корней: " << result.roots_found << "/" << result.roots_expected
	    << "  макс.отн.ошибка=" << scientific << setprecision(3)
	    << result.max_rel_err << "  время=" << fixed << setprecision(1)
	    << result.elapsed_ms << "мс  "
	    << (result.passed ? "PASSED" : "FAILED") << endl;

	// Дополнительный тест: доминирующие корни (m = n/3)
	if (P > 5) {
		int m_dom = max(1, int(P) / 3);
		auto start2 = chrono::high_resolution_clock::now();
		vector<T> dom_moduli = find_dominant_moduli_by_SanJuan(
		    coeffs_desc, m_dom, graeffe_iters, epsilon);
		auto end2 = chrono::high_resolution_clock::now();

		double dom_ms = chrono::duration<double, milli>(end2 - start2).count();

		// Первые m_dom модулей из ожидаемых
		vector<T> exp_dom(exp_mod.begin(),
		                  exp_mod.begin() + min(m_dom, int(exp_mod.size())));
		sort_desc(dom_moduli);

		double dom_err = max_relative_error(dom_moduli, exp_dom);

		log << "    -> доминир.(m=" << m_dom << "): "
		    << dom_moduli.size() << " корней, ошибка=" << scientific
		    << setprecision(3) << dom_err << "  время=" << fixed
		    << setprecision(1) << dom_ms << "мс" << endl;
	}

	return result;
}

// =====================================================================
// Полный набор тестов для типа T
// =====================================================================

template <typename T>
vector<SJTestResult> run_all_sanjuan_tests(ostream &log) {
	vector<SJTestResult> results;
	T tol;
	if constexpr (is_same_v<T, float_precision>)
		tol = T("0.5");
	else
		tol = T(0.5);

	int giter = 10;
	if constexpr (is_same_v<T, float_precision>)
		giter = 8;

	log << "\n========== МЕТОД САН-ХУАНА: тесты для " << type_name_sj<T>()
	    << " ==========\n";

	// --- 1. Только вещественные корни ---
	for (unsigned deg : {5u, 10u, 20u, 30u, 50u}) {
		if constexpr (is_same_v<T, float_precision>) {
			if (deg > 20)
				continue;
		}
		T dcr;
		if constexpr (is_same_v<T, float_precision>)
			dcr = T("0.1");
		else
			dcr = T(0.1);

		vector<T> cr;
		results.push_back(run_single_sanjuan_test<T>(
		    "Веществ.корни", deg, 0, 0, {}, cr, {}, dcr, 12345 + deg, tol,
		    giter, log));
	}

	// --- 2. С комплексными парами ---
	for (unsigned deg : {10u, 20u, 30u}) {
		if constexpr (is_same_v<T, float_precision>) {
			if (deg > 20)
				continue;
		}
		unsigned cp = deg / 4;
		T dcr;
		if constexpr (is_same_v<T, float_precision>)
			dcr = T("0.1");
		else
			dcr = T(0.1);

		vector<T> cr;
		results.push_back(run_single_sanjuan_test<T>(
		    "Компл.пары", deg, cp, 0, {}, cr, {}, dcr, 54321 + deg, tol,
		    giter, log));
	}

	// --- 3. Кластеризованные корни ---
	for (unsigned deg : {10u, 20u, 30u}) {
		if constexpr (is_same_v<T, float_precision>) {
			if (deg > 20)
				continue;
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

		results.push_back(run_single_sanjuan_test<T>(
		    "Кластериз.", deg, 0, nc, cc, cr, {}, dcr, 99999 + deg, tol,
		    giter, log));
	}

	// --- 4. Кратные корни ---
	for (unsigned deg : {10u, 20u}) {
		if constexpr (is_same_v<T, float_precision>) {
			if (deg > 15)
				continue;
		}
		vector<pair<unsigned, unsigned>> mg = {{2, 2}, {3, 1}};
		T dcr;
		if constexpr (is_same_v<T, float_precision>)
			dcr = T("0.1");
		else
			dcr = T(0.1);

		vector<T> cr;
		results.push_back(run_single_sanjuan_test<T>(
		    "Кратные", deg, 0, 0, {}, cr, mg, dcr, 77777 + deg, tol, giter,
		    log));
	}

	// --- 5. Смешанный ---
	for (unsigned deg : {20u, 30u, 50u}) {
		if constexpr (is_same_v<T, float_precision>) {
			if (deg > 20)
				continue;
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

		results.push_back(run_single_sanjuan_test<T>(
		    "Смешанный", deg, cp, nc, cc, cr, mg, dcr, 11111 + deg, tol,
		    giter, log));
	}

	// --- 6. Высокие степени ---
	if constexpr (!is_same_v<T, float_precision>) {
		for (unsigned deg : {60u, 80u, 100u}) {
			T dcr = T(0.1);
			vector<T> cr;
			unsigned cp = deg / 5;
			results.push_back(run_single_sanjuan_test<T>(
			    "Высок.степ.", deg, cp, 0, {}, cr, {}, dcr, 33333 + deg,
			    tol, giter, log));
		}
	}

	return results;
}

// =====================================================================
// Сохранение в .txt
// =====================================================================

void save_sj_results(const vector<SJTestResult> &results,
                     const string &filename) {
	ofstream out(filename);
	if (!out.is_open()) {
		cerr << "Ошибка записи: " << filename << endl;
		return;
	}

	out << "======================================================================"
	       "====\n";
	out << "РЕЗУЛЬТАТЫ ТЕСТИРОВАНИЯ: МОДИФИКАЦИЯ САН-ХУАНА метода Грефе\n";
	if (!results.empty())
		out << "Тип данных: " << results[0].data_type << "\n";
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
		    << setw(18) << scientific << setprecision(3) << r.max_rel_err
		    << setw(12) << fixed << setprecision(1) << r.elapsed_ms
		    << setw(10) << (r.passed ? "PASSED" : "FAILED") << "\n";
		++total;
		if (r.passed)
			++passed;
	}

	out << "\n" << string(98, '=') << "\n";
	out << "Итого: " << passed << " из " << total << " тестов пройдено.\n";
	out.close();

	cout << "Результаты сохранены в: " << filename << endl;
}

// =====================================================================
// main
// =====================================================================

int main() {
#ifdef _WIN32
	SetConsoleOutputCP(65001);
	SetConsoleCP(65001);
#endif

	cout << "==============================================\n";
	cout << " ТЕСТИРОВАНИЕ МОДИФИКАЦИИ САН-ХУАНА\n";
	cout << " (Усечённый полином после итераций Грефе)\n";
	cout << "==============================================\n\n";

	// double
	{
		stringstream log;
		auto results = run_all_sanjuan_tests<double>(log);
		cout << log.str();
		save_sj_results(results, "sanjuan_results_double.txt");
	}

	// long double
	{
		stringstream log;
		auto results = run_all_sanjuan_tests<long double>(log);
		cout << log.str();
		save_sj_results(results, "sanjuan_results_long_double.txt");
	}

	// float_precision
	{
		stringstream log;
		auto results = run_all_sanjuan_tests<float_precision>(log);
		cout << log.str();
		save_sj_results(results, "sanjuan_results_float_precision.txt");
	}

	cout << "\nВсе тесты Сан-Хуана завершены.\n";
	return 0;
}
