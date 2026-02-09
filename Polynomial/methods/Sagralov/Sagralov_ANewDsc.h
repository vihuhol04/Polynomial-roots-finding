// Алгоритм ANewDsc
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#include "NumericConcepts.h"
#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

#include "NumericConstants.h"
#include "mathUtils.h"
#include "polynomialUtils.h"

// СТРУКТУРА ИНТЕРВАЛА
/**
 * @brief Интервал для метода ANewDsc с информацией о вариабельности.
 *
 * Хранит число изменений знака, вычисленное по правилу Декарта
 * для полиномиального преобразования.
 *
 * @tparam T Вещественный тип, должен удовлетворять концепту RealNumber.
 */
template <RealNumber T> struct Interval01 {
	T left;          // левая граница интервала
	T right;         // правая граница интервала
	int variability; // количество возможных корней внутри интервала

	// Конструктор: по умолчанию количество возможных корней = -1 (ещё не вычислено)
	Interval01(T l, T r, int v = -1) : left(l), right(r), variability(v) {}

	T width() const { return abs_val(right - left); }
	T midpoint() const { return (left + right) / T(2); }
};

// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
namespace anewdsc_detail {
	/**
	 * @brief Декартово преобразование интервала [a, b].
	 *
	 * Преобразует исходный полином P(x) на интервале [a,b] в новый Q(y),
	 * где корни в [a,b] отображаются в [0, +inf).
	 */
	template <RealNumber T>
	std::vector<T> descartes_transform(const std::vector<T>& poly, T a, T b) {
		int n = static_cast<int>(poly.size()) - 1;
		if (n < 0)
			return {};

		std::vector<T> result(n + 1, T(0));

		std::vector<T> numerator(n + 1);
		std::vector<T> denominator(n + 1);

		for (int i = 0; i <= n; ++i) {
			// (a*y + b)^i
			std::fill(numerator.begin(), numerator.begin() + i + 1, T(0));
			numerator[0] = T(1);
			for (int k = 0; k < i; ++k) {
				for (int j = k + 1; j >= 1; --j)
					numerator[j] = numerator[j] * b + numerator[j - 1] * a;
				numerator[0] *= b;
			}

			// (y + 1)^(n - i)
			std::fill(denominator.begin(), denominator.begin() + (n - i) + 1, T(0));
			denominator[0] = T(1);
			for (int k = 0; k < n - i; ++k) {
				for (int j = k + 1; j >= 1; --j)
					denominator[j] = denominator[j] + denominator[j - 1];
			}

			// Перемножаем обе части и добавляем в итоговый полином
			for (int p = 0; p <= i; ++p)
				for (int q = 0; q <= n - i; ++q)
					result[p + q] += poly[i] * numerator[p] * denominator[q];
		}
		return result;
	}

	/**
	 * @brief Подсчёт изменений знака в последовательности коэффициентов.
	 *
	 * Используется правило Декарта: количество положительных корней
	 * не больше числа изменений знака.
	 */
	template <RealNumber T>
	int count_sign_changes(const std::vector<T>& coeffs, T eps) {
		int changes = 0;
		bool have_sign = false;
		bool last_sign = false;

		for (const auto& c : coeffs) {
			if (is_zero_val(c, eps))
				continue;

			bool current_sign = (c > T(0));
			changes += have_sign * (current_sign != last_sign);
			last_sign = current_sign;
			have_sign = true;
		}
		return changes;
	}

	/**
	 * @brief 01-тест из метода ANewDsc.
	 *
	 * Проверяет количество корней полинома в интервале [a,b]:
	 * - 0: нет корней
	 * - 1: ровно один корень
	 * - -1: неопределённо (несколько корней или кратные)
	 */
	template <RealNumber T>
	int test_01(const std::vector<T>& poly, T a, T b, T eps) {
		if (abs_val(b - a) < eps)
			return -1;

		T fa = eval_poly(poly, a);
		T fb = eval_poly(poly, b);

		// Точный корень на границе
		if (is_zero_val(fa, eps) || is_zero_val(fb, eps))
			return 1;

		auto transformed = descartes_transform(poly, a, b);
		int v = count_sign_changes(transformed, eps);

		return (v <= 1) ? v : -1;
	}

	/**
	 * @brief Поиск допустимой точки внутри интервала.
	 *
	 * Ищет точку, где значение полинома максимально далеко от нуля.
	 *
	 * @param grid_size Размер сетки для поиска (параметр алгоритма).
	 */
	template <RealNumber T>
	std::optional<T> find_admissible_point(const std::vector<T>& poly, T a, T b, int grid_size, T eps) {
		T mid = (a + b) / T(2);
		T step = (b - a) / T(grid_size);

		T max_val_found = T(0);
		T best_point = mid;

		for (int i = 0; i <= grid_size; ++i) {
			T x = a + step * T(i);
			T val = abs_val(eval_poly(poly, x));
			best_point = (val > max_val_found) ? x : best_point;
			max_val_found = max_val(max_val_found, val);
		}

		T threshold = max_val_found / T(4);

		for (int i = grid_size / 2 - 2; i <= grid_size / 2 + 2; ++i) {
			if (i < 0 || i > grid_size)
				continue;

			T x = a + step * T(i);
			T val = abs_val(eval_poly(poly, x));

			if (val >= threshold)
				return x;
		}
		return mid;
	}

	/**
	 * @brief Тест Ньютона для ускорения сходимости.
	 *
	 * @param N_param Параметр, задающий насколько узким должен стать интервал.
	 */
	template <RealNumber T>
	std::optional<Interval01<T>> newton_test(const std::vector<T>& poly, const Interval01<T>& interval, int N_param, T eps) {
		if (interval.variability <= 1)
			return std::nullopt;

		T width = interval.width();
		T target_width = width / T(N_param);

		auto deriv = derivative(poly);

		std::vector<T> sample_points = { interval.left + width * T(0.25), interval.midpoint(), interval.left + width * T(0.75) };

		T avg_center = T(0);
		int valid_samples = 0;

		for (const auto& xi : sample_points) {
			T f_val = eval_poly(poly, xi);
			T fp_val = eval_poly(deriv, xi);

			// проверяем результат деления
			T correction = f_val / fp_val;
			if (!isfinite_val(correction))
				continue;

			T candidate = xi - correction;

			// Корректная точка Ньютона должна лежать в интервале
			bool in_interval = (candidate >= interval.left && candidate <= interval.right);
			avg_center += candidate * T(in_interval);
			valid_samples += in_interval;
		}

		if (valid_samples == 0)
			return std::nullopt;

		T lambda = avg_center / T(valid_samples);
		T new_left = lambda - target_width / T(2);
		T new_right = lambda + target_width / T(2);

		new_left = max_val(interval.left, new_left);
		new_right = min_val(interval.right, new_right);

		int test_left = test_01(poly, interval.left, new_left, eps);
		int test_right = test_01(poly, new_right, interval.right, eps);

		if (test_left == 0 && test_right == 0) {
			int test_center = test_01(poly, new_left, new_right, eps);
			if (test_center >= 0)
				return Interval01<T>(new_left, new_right, test_center);
		}

		T delta = target_width / T(2);
		if (test_01(poly, interval.left + delta, interval.right, eps) == 0)
			return Interval01<T>(interval.left, interval.left + delta);
		if (test_01(poly, interval.left, interval.right - delta, eps) == 0)
			return Interval01<T>(interval.right - delta, interval.right);

		return std::nullopt;
	}
} // namespace anewdsc_detail

// ОСНОВНАЯ ФУНКЦИЯ АЛГОРИТМА
/**
 * @brief Алгоритм ANewDsc для изоляции вещественных корней.
 *
 * Реализует улучшенный алгоритм Декарта с тестом Ньютона для ускорения.
 *
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 * @param input Коэффициенты многочлена.
 * @param epsilon Требуемая точность изоляции.
 * @param max_iterations Максимальное число итераций.
 * @return Вектор изолирующих интервалов (left, right).
 */
template <RealNumber T>
std::vector<std::pair<T, T>>
find_real_roots_ANewDsc(const std::vector<T>& input, T epsilon = numeric_constants::adaptive_epsilon<T>(numeric_constants::EPSILON_SCALE_STANDARD), int max_iterations = 1000) {
	if (input.empty() || input.size() == 1)
		return {};

	std::vector<T> poly = normalizePolynomial(input);

	// Оценка границы корней по правилу Коши
	T maxc = T(0);
	for (size_t i = 0; i + 1 < poly.size(); ++i)
		maxc = max_val(maxc, abs_val(poly[i]));

	T lead = abs_val(poly.back());
	if (is_zero_val(lead))
		throw std::runtime_error("zero leading coefficient");

	T bound = T(1) + maxc / lead;

	std::vector<Interval01<T>> active_intervals;
	active_intervals.emplace_back(-bound, bound, -1);

	std::vector<std::pair<T, T>> isolated_roots;
	int iteration = 0;

	// Параметр размера сетки для find_admissible_point
	constexpr int GRID_SIZE = numeric_constants::ADMISSIBLE_POINT_MAX_TRIES;
	// Параметр для теста Ньютона
	constexpr int NEWTON_PARAM = 4;

	while (!active_intervals.empty() && iteration < max_iterations) {
		++iteration;

		Interval01<T> current = active_intervals.back();
		active_intervals.pop_back();

		if (current.width() < epsilon) {
			isolated_roots.emplace_back(current.left, current.right);
			continue;
		}

		int test_result = anewdsc_detail::test_01(poly, current.left, current.right, epsilon);

		if (test_result == 0)
			continue;

		if (test_result == 1) {
			if (current.width() < epsilon) {
				isolated_roots.emplace_back(current.left, current.right);
				continue;
			}
			test_result = -1;
		}

		// Ускоряем поиск с помощью теста Ньютона
		if (current.variability > 1) {
			auto newton_result = anewdsc_detail::newton_test(poly, current, NEWTON_PARAM, epsilon);

			if (newton_result) {
				active_intervals.push_back(*newton_result);
				continue;
			}
		}

		// Бисекция
		auto m_opt = anewdsc_detail::find_admissible_point(poly, current.left, current.right, GRID_SIZE, epsilon);
		T m = m_opt.value_or(current.midpoint());

		int left_test = anewdsc_detail::test_01(poly, current.left, m, epsilon);
		int right_test = anewdsc_detail::test_01(poly, m, current.right, epsilon);

		if (left_test != 0)
			active_intervals.emplace_back(current.left, m, -1);
		if (right_test != 0)
			active_intervals.emplace_back(m, current.right, -1);
	}

	if (isolated_roots.size() <= 1)
		return isolated_roots;

	// Сортировка и слияние перекрывающихся интервалов
	std::sort(isolated_roots.begin(), isolated_roots.end(),
		[](const auto& a, const auto& b) { return a.first < b.first; });

	std::vector<std::pair<T, T>> merged;
	merged.push_back(isolated_roots[0]);

	for (size_t i = 1; i < isolated_roots.size(); ++i) {
		auto& last = merged.back();
		const auto& current = isolated_roots[i];

		if (current.first <= last.second + epsilon)
			last.second = max_val(last.second, current.second);
		else
			merged.push_back(current);
	}
	return merged;
}