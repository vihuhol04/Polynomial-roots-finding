// Метод Греффе для вычисления модулей вещественных корней
// find_moduli_roots_by_graeffe - реализация метода без отсеивания повторяющихся корней (порождает кучу мусора) 
// find_moduli_with_multiplicities_by_graeffe - отсеивает дублирующиеся корни 
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <algorithm>
#include <chrono>
#include <iostream>
#include <utility>
#include <vector>

#include "NumericConcepts.h"
#include "NumericConstants.h"
#include "mathUtils.h"

// ОСНОВНАЯ ФУНКЦИЯ МЕТОДА ГРЕФФЕ

/**
 * @brief Вычисление модулей корней методом Греффе (логарифмическая модификация).
 *
 * Идея: хранить коэффициенты в логарифмическом виде для повышенной точности.
 * Векторы L и S хранят логарифмы модулей и знаки коэффициентов соответственно.
 *
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 * @param coeffs_desc Коэффициенты многочлена (от старших к младшим степеням).
 * @param epsilon Требуемая точность сходимости.
 * @param maxIter Максимальное число итераций Греффе.
 * @return Вектор модулей корней.
 */
template <RealNumber T>
std::vector<T> find_moduli_roots_by_graeffe(const std::vector<T>& coeffs_desc,
	T epsilon = numeric_constants::adaptive_epsilon<T>(numeric_constants::EPSILON_SCALE_PRECISE),
	int maxIter = numeric_constants::DEFAULT_MAX_ITERATIONS, bool debug = false)
{
	std::cerr << "[TRACE] Entering find_moduli_roots_by_graeffe. Deg=" << coeffs_desc.size() - 1 << " MaxIter=" << maxIter << std::endl;

	const T NEG_INF = neg_inf_val<T>();

	// Безопасный логарифм: использует isfinite для проверки результата
	auto safe_log = [&](const T& x) -> T {
		T v = abs_val(x);
		T result = log_val(v);
		return isfinite_val(result) ? result : NEG_INF;
		};

	// Вычисление log(sum(exp(vals))) — устойчиво при больших разностях
	auto log_sum_exp = [&](const std::vector<T>& vals) -> T {
		if (vals.empty())
			return NEG_INF;
		T m = *::std::max_element(vals.begin(), vals.end());
		if (m == NEG_INF)
			return NEG_INF;

		T s = T(0);
		for (auto v : vals)
			s += exp_val(v - m) * T(v != NEG_INF);

		T log_s = log_val(s);
		return isfinite_val(log_s) ? m + log_s : NEG_INF;
		};

	// Пропуск нулевых коэффициентов при старших степенях многочлена
	// Коэффициент считается "нулевым", если его логарифм не является конечным числом, что автоматически определяется через safe_log и isfinite_val
	size_t first_nz = 0;
	while (first_nz < coeffs_desc.size() && !isfinite_val(safe_log(coeffs_desc[first_nz])))
		++first_nz;
	if (first_nz >= coeffs_desc.size())
		return {};

	// Берем только значимую часть
	std::vector<T> coeffs(coeffs_desc.begin() + first_nz, coeffs_desc.end());
	size_t deg = coeffs.size() - 1;
	if (deg == 0)
		return {};

	// Нормализация (ведущий коэффициент = 1)
	T lead = coeffs.front();
	if (!isfinite_val(safe_log(lead)))
		return {};

	std::vector<T> p;
	p.reserve(coeffs.size());
	for (auto& c : coeffs)
		p.push_back(c / lead);

	// Меняем порядок коэффициентов: было от старших к младшим → станет от младших к старшим
	// Это необходимо для корректной работы логарифмической модификации Греффе
	std::reverse(p.begin(), p.end());

	// Инициализация логарифмов и знаков коэффициентов
	std::vector<T> L(p.size());
	std::vector<int> S(p.size());
	for (size_t i = 0; i < p.size(); ++i) {
		L[i] = safe_log(p[i]);
		S[i] = (p[i] >= T(0)) * 2 - 1; // арифметико-логическое выражение вместо if
	}

	int iters_done = 0;
	bool converged = false;

	if (debug) {
		std::cerr << "[DEBUG] Graeffe: начало, степень полинома = " << (p.size() - 1) << ", maxIter = " << maxIter << std::endl;
	}

	// Векторы для расщепления (аллоцируем один раз вовне цикла)
	std::vector<T> Lpe, Lpo;
	std::vector<int> Spe, Spo;

	// Основной цикл метода Греффе
	for (int iter = 0; iter < maxIter; ++iter) {
		auto iter_start = std::chrono::high_resolution_clock::now();

		// Разделение на чётные/нечётные индексы
		// Используем шаг 2 вместо проверки i % 2 == 0 для устранения ветвления
		Lpe.clear();
		Spe.clear();
		Lpo.clear();
		Spo.clear();

		// Предварительно резервируем память для избежания реаллокаций
		size_t half_size = (L.size() + 1) / 2;
		Lpe.reserve(half_size);
		Spe.reserve(half_size);
		Lpo.reserve(half_size);
		Spo.reserve(half_size);

		for (size_t i = 0; i < L.size(); i += 2) {
			Lpe.push_back(L[i]);
			Spe.push_back(S[i]);
		}
		for (size_t i = 1; i < L.size(); i += 2) {
			Lpo.push_back(L[i]);
			Spo.push_back(S[i]);
		}

		// Свёртка: возведение в квадрат (в логарифмической форме)
		auto conv_square_log = [&](const std::vector<T>& Lvec, const std::vector<int>& Svec) {
			size_t msize = (Lvec.empty()) ? 0 : (2 * Lvec.size() - 1);
			std::vector<T> Lout(msize, NEG_INF);
			std::vector<int> Sout(msize, 0);

			for (size_t m = 0; m < msize; ++m) {
				std::vector<std::pair<T, int>> terms;

				// Свёртка: b[m] = сумма (a[i]*a[m-i])
				for (size_t i = 0; i < Lvec.size(); ++i) {
					size_t j = m - i;
					if (j < Lvec.size()) {
						T lv = Lvec[i] + Lvec[j];
						terms.emplace_back(lv, Svec[i] * Svec[j]);
					}
				}

				if (terms.empty())
					continue;

				// Если все слагаемые одного знака, используем log_sum_exp
				bool all_same_sign = true;
				int first_sign = terms[0].second;
				for (const auto& t : terms) {
					if (t.second != first_sign) {
						all_same_sign = false;
						break;
					}
				}

				if (all_same_sign) {
					std::vector<T> vals;
					for (const auto& t : terms)
						vals.push_back(t.first);
					Lout[m] = log_sum_exp(vals);
					Sout[m] = first_sign;
				}
				else {
					T sum = T(0);
					for (const auto& t : terms)
						sum += t.second * exp_val(t.first);

					// Используем isfinite вместо сравнения с магической константой
					T log_abs = log_val(abs_val(sum));
					if (!isfinite_val(log_abs)) {
						Lout[m] = NEG_INF;
						Sout[m] = 0;
					}
					else {
						Lout[m] = log_abs;
						Sout[m] = (sum > T(0)) ? 1 : -1;
					}
				}
			}
			return ::std::make_pair(Lout, Sout);
			};

		// Получаем квадраты для чётной и нечётной частей
		auto [Lpe2, Spe2] = conv_square_log(Lpe, Spe);
		auto [Lpo2, Spo2] = conv_square_log(Lpo, Spo);

		// Формируем новый многочлен: Q(y) = Pe(y)^2 - y * Po2(y)
		std::vector<T> Lq(L.size(), NEG_INF);
		std::vector<int> Sq(L.size(), 0);

		for (size_t m = 0; m < L.size(); ++m) {
			bool hasPe = (m < Lpe2.size() && is_finite_ld(Lpe2[m]));
			bool hasPo = (m >= 1 && (m - 1) < Lpo2.size() && is_finite_ld(Lpo2[m - 1]));

			if (!hasPe && !hasPo)
				continue;

			if (hasPe && !hasPo) {
				Lq[m] = Lpe2[m];
				Sq[m] = Spe2[m];
			}
			else if (!hasPe && hasPo) {
				Lq[m] = Lpo2[m - 1];
				Sq[m] = -Spo2[m - 1];
			}
			else {
				// Оба есть: вычисляем комбинацию Pe^2 - y*Po^2
				T A = Lpe2[m];
				T B = Lpo2[m - 1];
				int sA = Spe2[m];
				int sB = Spo2[m - 1];

				if (A == NEG_INF && B == NEG_INF)
					continue;

				// Если A примерно равен B и знаки противоположные
				T log_diff = abs_val(A - B);
				if (sA != sB && log_diff < numeric_constants::log_diff_threshold<T>()) {
					// Почти одинаковые величины с противоположными знаками
					T val_A = sA * exp_val(A);
					T val_B = -sB * exp_val(B);
					T comb = val_A + val_B;

					T log_comb = log_val(abs_val(comb));
					// Используем isfinite вместо сравнения с магической константой
					Lq[m] = isfinite_val(log_comb) ? log_comb : NEG_INF;
					Sq[m] = isfinite_val(log_comb) ? ((comb > T(0)) * 2 - 1) : 0;
				}
				else {
					// Обычное логарифмическое сложение
					T M = max_val(A, B);
					T u = sA * exp_val(A - M);
					T v = -sB * exp_val(B - M);
					T comb = u + v;

					T log_comb = log_val(abs_val(comb));
					Lq[m] = isfinite_val(log_comb) ? (M + log_comb) : NEG_INF;
					Sq[m] = isfinite_val(log_comb) ? ((comb > T(0)) * 2 - 1) : 0;
				}
			}
		}

		// нормализация по старшему коэффициенту
		T lead_log = Lq.back();
		if (!is_finite_ld(lead_log)) {
			break;
		}

		for (size_t i = 0; i < Lq.size(); ++i)
			Lq[i] -= lead_log * T(is_finite_ld(Lq[i]));

		// нормализация старого L
		T lead_log_old = L.back();
		if (is_finite_ld(lead_log_old)) {
			for (size_t i = 0; i < L.size(); ++i)
				L[i] -= lead_log_old * T(is_finite_ld(L[i]));
		}

		// проверка сходимости
		T max_rel = T(0);
		bool local_nonfinite = false;

		// вычисляем log-разности для текущего и предыдущего полиномов
		std::vector<T> log_ratios_old, log_ratios_new;
		for (size_t i = 1; i < L.size(); ++i) {
			if (is_finite_ld(L[i]) && is_finite_ld(L[i - 1]))
				log_ratios_old.push_back(L[i] - L[i - 1]);
			if (is_finite_ld(Lq[i]) && is_finite_ld(Lq[i - 1]))
				log_ratios_new.push_back(Lq[i] - Lq[i - 1]);
		}

		// Сравниваем изменение log-разностей
		if (log_ratios_old.size() == log_ratios_new.size() &&
			!log_ratios_old.empty()) {
			for (size_t i = 0; i < log_ratios_old.size(); ++i) {
				T ratio_old = log_ratios_old[i];
				T ratio_new = log_ratios_new[i];

				if (!is_finite_ld(ratio_old) || !is_finite_ld(ratio_new)) {
					local_nonfinite = true;
					max_rel = pos_inf_val<T>();
					break;
				}

				T delta = abs_val(ratio_new - T(2) * ratio_old);
				T scale = max_val(abs_val(ratio_new), abs_val(ratio_old));

				T rel_change = (scale > T(0)) * (delta / scale);
				max_rel = max_val(max_rel, rel_change);
			}
		}
		else {
			// Сравниваем просто изменение коэффициентов
			size_t min_size = min_val(L.size(), Lq.size());
			for (size_t i = 0; i < min_size; ++i) {
				if (!is_finite_ld(L[i]) || !is_finite_ld(Lq[i]))
					continue;

				T delta = Lq[i] - L[i];
				if (!is_finite_ld(delta)) {
					local_nonfinite = true;
					max_rel = pos_inf_val<T>();
					break;
				}

				// Ограничиваем delta для безопасного exp()
				// Используем константу из NumericConstants.h
				constexpr double LOG_MAX_SAFE_D = numeric_constants::LOG_MAX_SAFE;
				T LOG_MAX_SAFE = T(LOG_MAX_SAFE_D);
				T clamped = clamp_val(delta, -LOG_MAX_SAFE, LOG_MAX_SAFE);
				T diff = abs_val(exp_val(clamped) - T(1));

				// Проверка после вычисления вместо предварительной проверки
				max_rel = isfinite_val(diff) ? max_val(max_rel, diff) : pos_inf_val<T>();
				local_nonfinite = local_nonfinite || !isfinite_val(diff);
			}
		}

		// Обновляем L и S
		L = std::move(Lq);
		S = std::move(Sq);
		iters_done = iter + 1;

		// DEBUG: вывод времени итерации
		if (debug) {
			auto iter_end = std::chrono::high_resolution_clock::now();
			auto iter_ms = std::chrono::duration_cast<std::chrono::milliseconds>(iter_end - iter_start).count();
			std::cerr << "[DEBUG] Graeffe iter " << iter << ": max_rel = " << static_cast<double>(max_rel) << ", time = " << iter_ms << " ms" << std::endl;
		}

		// Проверка критериев остановки
		if (local_nonfinite || !isfinite_val(max_rel))
			return {};

		if (max_rel < epsilon) {
			if (debug) {
				std::cerr << "[DEBUG] Graeffe: сходимость достигнута на итерации " << iter << std::endl;
			}
			converged = true;
			break;
		}
	}

	// обратно переворачиваем
	std::vector<T> Ldesc(L.rbegin(), L.rend());
	std::vector<int> Sdesc(S.rbegin(), S.rend());
	size_t n = Ldesc.size() - 1;

	if (!converged || n == 0)
		return {};

	// вычисляем степень двойки для итераций
	const int effective_iters = max_val(iters_done, 1);
	T pow_two_k = pow_val(T(2), static_cast<T>(effective_iters));

	// Проверяем наличие нулевых средних коэффициентов
	// Если некоторые средние коэффициенты равны нулю (или -беск. в лог-форме), это означает, что все корни имеют одинаковый модуль
	bool has_zero_middle = false;
	for (size_t i = 1; i + 1 < Ldesc.size(); ++i) {
		if (!is_finite_ld(Ldesc[i]) || Ldesc[i] == NEG_INF) {
			has_zero_middle = true;
			break;
		}
	}

	// если есть нулевые средние коэффициенты, все корни имеют одинаковый модуль
	if (has_zero_middle) {
		T L_first = Ldesc[0];
		T L_last = Ldesc.back();

		if (is_finite_ld(L_first) && is_finite_ld(L_last)) {
			T log_ratio = (L_last - L_first) / static_cast<T>(n);
			T common_modulus = exp_val(log_ratio / pow_two_k);
			// если модуль неконечный или отрицательный, заменяем на 0
			common_modulus = (!is_finite_ld(common_modulus) || common_modulus < T(0)) ? T(0) : common_modulus;
			return std::vector<T>(n, common_modulus);
		}
	}

	// извлечение модулей по соседним коэффициентам
	std::vector<T> moduli;
	moduli.reserve(n);

	if (debug) {
		std::cerr << "[DEBUG] Graeffe: извлечение модулей" << std::endl;
		std::cerr << "[DEBUG]   n = " << n << ", iters_done = " << iters_done << ", pow_two_k = " << static_cast<double>(pow_two_k) << std::endl;
		std::cerr << "[DEBUG]   Ldesc (коэффициенты в log): ";
		for (size_t i = 0; i < Ldesc.size(); ++i)
			std::cerr << static_cast<double>(Ldesc[i]) << " ";
		std::cerr << std::endl;
	}

	// Переменные для логарифмов соседних коэффициентов — аллоцируем вне цикла
	T Lprev, Lcurr;
	for (size_t j = 1; j <= n; ++j) {
		Lprev = Ldesc[j - 1];
		Lcurr = Ldesc[j];

		// если оба нулевые, модуль = 0
		if (!is_finite_ld(Lprev) && !is_finite_ld(Lcurr)) {
			if (debug) {
				std::cerr << "[DEBUG]   j=" << j << ": оба -inf, модуль = 0" << std::endl;
			}
			moduli.push_back(T(0));
			continue;
		}

		// ВАЖНО: если Lcurr = -inf (нулевой коэффициент при младшей степени), это означает что полином делится на x, т.е. один из корней = 0
		// В этом случае модуль этого корня = 0
		if (!is_finite_ld(Lcurr)) {
			if (debug) {
				std::cerr << "[DEBUG]   j=" << j << ": Lcurr=-inf (нулевой коэфф.), модуль = 0" << std::endl;
			}
			moduli.push_back(T(0));
			continue;
		}

		// Подстановка для -inf только для Lprev: если предыдущий коэффициент был нулевой, но текущий ненулевой, используем LOG_TINY для стабильности
		constexpr double LOG_TINY_D = numeric_constants::LOG_TINY;
		const T LOG_TINY = T(LOG_TINY_D);
		if (!is_finite_ld(Lprev))
			Lprev = LOG_TINY;

		// log|x_j| = (log|A_j| - log|A_{j-1}|) / 2^k
		T log_modulus = (Lcurr - Lprev) / pow_two_k;

		if (debug) {
			std::cerr << "[DEBUG]   j=" << j << ": Lprev=" << static_cast<double>(Ldesc[j - 1]) << ", Lcurr=" << static_cast<double>(Ldesc[j])
				<< ", log_modulus=" << static_cast<double>(log_modulus) << std::endl;
		}

		if (!is_finite_ld(log_modulus)) {
			moduli.push_back(T(0));
			continue;
		}

		T modulus = exp_val(log_modulus);
		modulus = (isfinite_val(modulus) && modulus >= T(0)) ? modulus : T(0);

		if (debug) {
			std::cerr << "[DEBUG]   j=" << j << ": modulus=" << static_cast<double>(modulus) << std::endl;
		}

		moduli.push_back(modulus);
	}
	return moduli;
}

// СТРУКТУРЫ И ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
/**
 * @brief Структура для хранения модуля корня с его кратностью.
 * @tparam T Вещественный тип, должен удовлетворять концепту RealNumber.
 */
template <RealNumber T> struct ModulusWithMultiplicity {
	T value;          // значение модуля
	int multiplicity; // сколько раз этот модуль встретился

	ModulusWithMultiplicity(T val, int mult = 1) : value(val), multiplicity(mult) {}
};

/**
 * @brief Версия метода Греффе, возвращающая модули с кратностями.
 *
 * Группирует близкие по значению модули и подсчитывает их кратность.
 *
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 * @param coeffs_desc Коэффициенты многочлена.
 * @param epsilon Требуемая точность.
 * @param maxIter Максимальное число итераций.
 * @param debug Выводить отладочную информацию.
 * @return Вектор модулей с кратностями.
 */
template <RealNumber T>
std::vector<ModulusWithMultiplicity<T>>
find_moduli_with_multiplicities_by_graeffe(const std::vector<T>& coeffs_desc, T epsilon = numeric_constants::adaptive_epsilon<T>(
	numeric_constants::EPSILON_SCALE_PRECISE), int maxIter = numeric_constants::DEFAULT_MAX_ITERATIONS, bool debug = false) {
	std::cerr << "[TRACE] Entering find_moduli_with_multiplicities_by_graeffe" << std::endl;
	// Сначала находим модули стандартным методом
	std::vector<T> raw_moduli = find_moduli_roots_by_graeffe(coeffs_desc, epsilon, maxIter, debug);

	if (raw_moduli.empty())
		return {};

	// Группируем модули с учётом погрешности
	std::vector<ModulusWithMultiplicity<T>> grouped;
	// Порог группировки: используем функцию из NumericConstants.h
	const T tolerance = numeric_constants::grouping_tolerance(epsilon);

	for (const auto& modulus : raw_moduli) {
		bool found = false;

		for (auto& group : grouped) {
			T diff = abs_val(group.value - modulus);

			if (diff < tolerance) {
				// Найден похожий модуль — увеличиваем кратность
				group.multiplicity++;
				// Обновляем значение модуля (скользящее среднее)
				group.value = (group.value * T(group.multiplicity - 1) + modulus) /
					T(group.multiplicity);
				found = true;
				break;
			}
		}

		if (!found)
			grouped.emplace_back(modulus, 1);
	}
	return grouped;
}