// Константы для численных методов поиска корней полиномов
// Вынесены из кода для устранения "магических чисел"

#pragma once

#include <limits>
#include <type_traits>

namespace numeric_constants {
	// КОНСТАНТЫ ИТЕРАЦИЙ

	/// Максимальное число итераций по умолчанию для итерационных методов
	constexpr int DEFAULT_MAX_ITERATIONS = 100;

	/// Максимальное число итераций для метода Ньютона
	constexpr int NEWTON_MAX_ITERATIONS = 100;

	/// Максимальное число итераций для метода Сагралова (уточнение корней)
	constexpr int SAGRALOV_MAX_ITERATIONS = 100;

	/// Максимальное число итераций для бисекции
	constexpr int BISECTION_MAX_ITERATIONS = 60;

	/// Максимальное число попыток поиска допустимой точки
	constexpr int ADMISSIBLE_POINT_MAX_TRIES = 20;

	// МАСШТАБНЫЕ МНОЖИТЕЛИ (FACTOR_*)
	// Используются для замены магических чисел (10, 100, 1000 и т.д.) при масштабировании epsilon и задания порогов

	constexpr long long FACTOR_5 = 5;
	constexpr long long FACTOR_10 = 10;
	constexpr long long FACTOR_20 = 20;
	constexpr long long FACTOR_50 = 50;
	constexpr long long FACTOR_100 = 100;
	constexpr long long FACTOR_1K = 1000;
	constexpr long long FACTOR_10K = 10000;
	constexpr long long FACTOR_100K = 100000;
	constexpr long long FACTOR_1M = 1000000;
	constexpr long long FACTOR_100M = 100000000;
	constexpr long long FACTOR_1B = 1000000000;

	// Семантические алиасы для часто используемых порогов
	constexpr long long RELAXED_TOLERANCE = FACTOR_100;
	constexpr long long STRICT_TOLERANCE = FACTOR_10;
	constexpr long long GRID_SEARCH_RELAXATION = FACTOR_10;
	constexpr long long GRID_SEARCH_RELAXATION_HUGE = FACTOR_1K;

	// ПОРОГИ ТОЧНОСТИ (адаптивные, зависят от типа T)
	// МАСШТАБЫ ДЛЯ АДАПТИВНОГО EPSILON
	//
	// Эти константы определяют, насколько "грубым" или "точным" будет порог
	// Результат: epsilon() * scale (чем больше scale, тем больше результат)
	//
	// Для типа double: epsilon() =(примерно) 2.2e-16
	// - COARSE (1e3)  -> eps =(примерно) 2.2e-13 (грубый порог для проверки "почти ноль")
	// - STANDARD (1e8) -> eps =(примерно) 2.2e-8  (стандартный порог для методов Сагралова)
	// - PRECISE (1e10) -> eps =(примерно) 2.2e-6  (точный порог для методов Греффе)
	//
	// Для arbitrary precision типов: epsilon() зависит от заданной точности, и эти масштабы автоматически адаптируются

	/// Масштаб для грубых проверок (is_zero_val без параметра)
	constexpr long long EPSILON_SCALE_COARSE = 1000LL;

	/// Масштаб для стандартных методов (Сагралов, ANewDsc) — обеспечивает 1e-8 для double
	constexpr long long EPSILON_SCALE_STANDARD = 100000000LL;

	/// Масштаб для точных методов (Грефе, Хосман) — обеспечивает 1e-6 для double
	constexpr long long EPSILON_SCALE_PRECISE = 10000000000LL;

	/// Масштаб для сверхточных проверок (шаг Ньютона) — обеспечивает 1e-2 для double
	constexpr long long EPSILON_SCALE_ULTRA = 1000000000000000000LL;

	/**
	 * @brief Вычисляет адаптивный epsilon для типа T.
	 *
	 * Для стандартных типов использует std::numeric_limits<T>::epsilon().
	 * Для arbitrary precision типов может использовать их собственные методы.
	 *
	 * @tparam T Числовой тип (double, float, float_precision, etc.)
	 * @param scale Множитель — используйте EPSILON_SCALE_* константы.
	 */
	template <typename T>
	inline T adaptive_epsilon(long long scale = EPSILON_SCALE_COARSE) {
		if constexpr (std::numeric_limits<T>::is_specialized) {
			return T(scale) * std::numeric_limits<T>::epsilon();
		}
		else {
			// Для типов без numeric_limits: используем 1/scale
			return T(1) / T(scale);
		}
	}

	/**
	 * @brief Универсальная функция для масштабирования epsilon.
	 *
	 * @param eps Базовый epsilon.
	 * @param factor Множитель (используйте FACTOR_* константы).
	 */
	template <typename T> inline T scale_epsilon(T eps, long long factor) {
		return eps * T(factor);
	}

	/**
	 * @brief Порог для сравнения близких значений.
	 *
	 * Используется для группировки модулей корней, сравнения коэффициентов и т.д.
	 */
	template <typename T> inline T grouping_tolerance(T epsilon) {
		return epsilon * T(100);
	}

	// КОНСТАНТЫ ДЛЯ БЕЗОПАСНЫХ ВЫЧИСЛЕНИЙ
	/**
	 * @brief Безопасный максимум для exp() — чтобы избежать overflow.
	 *
	 * ln(DBL_MAX) =(примерно) 709.8, используем 700 для безопасности.
	 * Для arbitrary precision это все равно безопасное значение.
	 */
	constexpr double LOG_MAX_SAFE = 700.0;

	/**
	 * @brief Безопасный минимум для log() — замена для -inf.
	 *
	 * ln(DBL_MIN) =(примерно) -708, используем -700 для консистентности.
	 */
	constexpr double LOG_TINY = -700.0;

	// КОНСТАНТЫ ДЛЯ МЕТОДА ГРЕФЕ
	/**
	 * @brief Порог для сравнения логарифмических разностей в методе Грефе.
	 *
	 * Если |A - B| < LOG_DIFF_THRESHOLD, считаем величины примерно равными.
	 */
	template <typename T> inline T log_diff_threshold() { return T(2); }

	// КОНСТАНТЫ ДЛЯ МЕТОДА САГРАЛОВА
	/**
	 * @brief Коэффициент надёжности для радиуса компоненты.
	 *
	 * Значение 1.15 обеспечивает запас при численных ошибках.
	 */
	template <typename T> inline T component_radius_safety_factor() {
		return T(115) / T(100); // 1.15 через целые числа
	}

	/**
	 * @brief Коэффициент K для надёжности при численных ошибках (тест Пелле).
	 *
	 * Значение 1.01 используется как множитель строгости проверки.
	 */
	template <typename T> inline T sagralov_safety_factor_k() {
		return T(101) / T(100); // 1.01 через целые числа
	}
} // namespace numeric_constants
