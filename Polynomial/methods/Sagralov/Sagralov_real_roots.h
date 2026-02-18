// Детерминированный алгоритм Сагралова для изоляции вещественных корней
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#include <algorithm>
#include <limits>
#include <map>
#include <utility>
#include <vector>

#include "NumericConstants.h"
#include "mathUtils.h"
#include "polynomialUtils.h"
#include "Helper_for_all_methods.h"
#include "NumericConcepts.h"

// РАЗРЕЖЁННЫЙ ПОЛИНОМ
/**
 * @brief Разрежённый полином, хранящий только ненулевые коэффициенты.
 *
 * Эффективен для полиномов высокой степени с малым числом членов.
 *
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 */
template <RealNumber T> struct SparsePolynomial {
	std::map<int, T> terms;

	SparsePolynomial() = default;

	/**
	 * @brief Конструктор из плотного представления полинома.
	 * @param dense Вектор коэффициентов.
	 * @param eps Порог для определения нулевого коэффициента.
	 */
	explicit SparsePolynomial(const std::vector<T>& dense, T eps) {
		for (size_t i = 0; i < dense.size(); ++i)
			if (!is_zero_val(dense[i], eps))
				terms[static_cast<int>(i)] = dense[i];
	}

	/// Конструктор с eps по умолчанию на основе numeric_limits
	explicit SparsePolynomial(const std::vector<T>& dense) : SparsePolynomial(dense, numeric_constants::scale_epsilon(
		std::numeric_limits<T>::epsilon(), numeric_constants::FACTOR_1K)) {
	}
	std::vector<T> toDense() const {
		if (terms.empty())
			return { T(0) };

		int md = terms.rbegin()->first;
		std::vector<T> out(md + 1, T(0));
		for (auto& kv : terms)
			out[kv.first] = kv.second;
		return out;
	}
	int degree() const { return terms.empty() ? -1 : terms.rbegin()->first; }
	size_t num_terms() const { return terms.size(); }
	T eval(T x) const {
		if (terms.empty())
			return T(0);
		T res(0);
		T x_pow(1);
		int prev_deg = 0;
		for (auto& kv : terms) {
			for (int d = prev_deg; d < kv.first; ++d)
				x_pow *= x;
			res += kv.second * x_pow;
			prev_deg = kv.first;
		}
		return res;
	}

	SparsePolynomial<T> operator+(const SparsePolynomial<T>& other) const {
		SparsePolynomial<T> res;
		for (auto& kv : terms)
			res.terms[kv.first] = kv.second;
		for (auto& kv : other.terms)
			res.terms[kv.first] += kv.second;
		return res;
	}
	SparsePolynomial<T> operator*(const SparsePolynomial<T>& other) const {
		SparsePolynomial<T> res;
		for (auto& kv1 : terms)
			for (auto& kv2 : other.terms)
				res.terms[kv1.first + kv2.first] += kv1.second * kv2.second;
		return res;
	}
	SparsePolynomial<T> shift(T a) const {
		if (terms.empty())
			return *this;

		SparsePolynomial<T> result;

		for (const auto& [deg, coef] : terms) {
			std::vector<T> binom_row(deg + 1);
			binom_row[0] = T(1);
			for (int k = 1; k <= deg; ++k)
				binom_row[k] = binom_row[k - 1] * T(deg - k + 1) / T(k);

			T a_pow = T(1);
			for (int j = 0; j <= deg; ++j) {
				result.terms[j] += coef * binom_row[j] * a_pow;
				a_pow *= (j < deg) ? a : T(1);
			}
		}
		return result;
	}
	SparsePolynomial<T> reciprocal() const {
		SparsePolynomial<T> res;
		if (terms.empty())
			return res;
		int d = degree();
		for (auto& kv : terms)
			res.terms[d - kv.first] = kv.second;
		return res;
	}
};

// СТРУКТУРА РЕЗУЛЬТАТА
/**
 * @brief Корень с кратностью, представленный изолирующим интервалом.
 * @tparam T Вещественный тип, должен удовлетворять концепту RealNumber.
 */
template <RealNumber T> struct RootWithMultiplicity {
	std::pair<T, T> interval;
	int multiplicity;
	RootWithMultiplicity(T a, T b, int m = 1) : interval(a, b), multiplicity(m) {}
};

// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ (ДЕТАЛИ РЕАЛИЗАЦИИ)
namespace detail {
	template <RealNumber T>
	std::vector<T> polyRemainder(const std::vector<T>& a, const std::vector<T>& b, T eps);
	template <RealNumber T>
	int countSignChanges(const std::vector<std::vector<T>>& seq, T x);
	template <RealNumber T>
	void isolateRootsRecursive(const std::vector<T>& poly,
		const std::vector<std::vector<T>>& sturm, T a, T b,
		std::vector<std::pair<T, T>>& out, T eps);
	template <RealNumber T>
	void refineInterval(std::pair<T, T>& interval, const std::vector<T>& poly,
		T epsilon);
	template <RealNumber T>
	int fastMultiplicityDeflation(const std::vector<T>& poly, T root, T eps);
	template <RealNumber T>
	int checkRootMultiplicity(const std::vector<T>& poly, T root, T eps);
	template <RealNumber T> struct RecursiveDerivativeChain {
		std::vector<SparsePolynomial<T>> polynomials;
		std::vector<int> min_degrees;
		std::vector<T> free_terms;
	};
	template <RealNumber T>
	RecursiveDerivativeChain<T>
		buildRecursiveDerivativeChain(const SparsePolynomial<T>& poly);
	template <RealNumber T>
	void isolateRootsByBackwardPass(const RecursiveDerivativeChain<T>& chain, T lo,
		T hi, std::vector<std::pair<T, T>>& out, T eps);
	template <RealNumber T>
	SparsePolynomial<T> sparseDerivative(const SparsePolynomial<T>& poly);
	template <RealNumber T>
	std::vector<T> findRealRoots(const SparsePolynomial<T>& poly, T eps);
} // namespace detail

// ГЛАВНАЯ ФУНКЦИЯ
/**
 * @brief Алгоритм Сагралова для изоляции вещественных корней с кратностями.
 *
 * @tparam T Тип коэффициентов, должен удовлетворять концепту RealNumber.
 * @param input Коэффициенты многочлена (от младших к старшим степеням).
 * @param epsilon Требуемая точность изоляции.
 * @return Вектор корней с кратностями.
 */
template <RealNumber T>
std::vector<RootWithMultiplicity<T>>
find_real_roots_by_Sagralov(const std::vector<T>& input, T epsilon = numeric_constants::adaptive_epsilon<T>(numeric_constants::EPSILON_SCALE_STANDARD)) {
	if (input.empty())
		throw std::invalid_argument("empty poly");
	if (input.size() == 1)
		return {};

	std::vector<T> coeffs = input;
	std::vector<RootWithMultiplicity<T>> results;

	// Корень в нуле с кратностью
	int mult0 = 0;
	while (coeffs.size() > 1 && is_zero_val(eval_poly(coeffs, T(0)), epsilon)) {
		mult0++;
		auto df = deflate_poly(coeffs, T(0));
		if (df.empty())
			break;
		coeffs = normalizePolynomial(df);
	}
	if (mult0 > 0) {
		// Используем static_cast вместо создания переменной e
		results.emplace_back(-static_cast<T>(epsilon), static_cast<T>(epsilon), mult0);
	}
	if (coeffs.size() <= 1)
		return results;

	int deg = (int)coeffs.size() - 1;
	if (deg <= 3) {
		// нормализуем по ведущему коэффициенту
		if (!is_zero_val(coeffs.back())) {
			T lc = coeffs.back();
			for (auto& c : coeffs)
				c /= lc;
		}

		std::vector<T> roots;
		if (deg == 1) {
			if (!is_zero_val(coeffs[1]))
				roots.push_back(-coeffs[0] / coeffs[1]);
		}
		else if (deg == 2) {
			T a = coeffs[2], b = coeffs[1], c = coeffs[0];
			if (!is_zero_val(a)) {
				T disc = b * b - T(4) * a * c;
				if (disc > T(0)) {
					T sd = sqrt_val(disc);
					// Численно устойчивая формула для квадратного уравнения:
					// Сначала находим корень с большим модулем (избегаем потери точности при вычитании близких чисел), затем второй по формуле Виета
					T root1 =
						(b >= T(0)) ? (-b - sd) / (T(2) * a) : (-b + sd) / (T(2) * a);
					// По формуле Виета: x1 * x2 = c/a
					T root2 = c / (a * root1);
					roots.emplace_back(root1);
					roots.emplace_back(root2);
				}
				else if (is_zero_val(disc, epsilon))
					roots.push_back(-b / (T(2) * a));
			}
		}
		else if (deg == 3) {
			// пробуем рациональные корни
			std::vector<T> candidates{ T(-3), T(-2), T(-1), T(0), T(1), T(2), T(3) };
			std::vector<T> work = coeffs;
			for (auto cand : candidates) {
				if (is_zero_val(eval_poly(work, cand), epsilon)) {
					roots.push_back(cand);
					work = deflate_poly(work, cand);
					if (work.size() == 2) { // линейный
						if (!is_zero_val(work[1]))
							roots.push_back(-work[0] / work[1]);
						break;
					}
					// ищем второй рациональный
					for (auto cand2 : candidates) {
						if (is_zero_val(eval_poly(work, cand2), epsilon)) {
							roots.push_back(cand2);
							work = deflate_poly(work, cand2);
							if (work.size() == 2 && !is_zero_val(work[1]))
								roots.push_back(-work[0] / work[1]);
							cand2 = candidates.back();
							break;
						}
					}
					break;
				}
			}
			// Бисекция для поиска оставшихся корней (если рациональные не нашлись)
			if (static_cast<int>(roots.size()) < deg) {
				T boundLocal = T(5);
				T step = T(0.25);
				T prev = -boundLocal;
				T prevVal = eval_poly(coeffs, prev);
				for (T x = prev + step;
					x <= boundLocal && static_cast<int>(roots.size()) < deg;
					x += step) {
					T val = eval_poly(coeffs, x);
					// Проверка смены знака: арифметико-логическое выражение
					bool sign_change =
						(prevVal < T(0) && val > T(0)) || (prevVal > T(0) && val < T(0));
					if (sign_change) {
						T L = prev, R = x;
						T fL = prevVal;
						// Максимум итераций бисекции: используем константу
						constexpr int MAX_BISECTION_ITERS =
							numeric_constants::BISECTION_MAX_ITERATIONS;
						for (int it = 0; it < MAX_BISECTION_ITERS; ++it) {
							T m = (L + R) / T(2);
							T fm = eval_poly(coeffs, m);
							// Условие остановки: значение близко к нулю или интервал
							// достаточно мал
							if (is_zero_val(fm, epsilon) || abs_val(R - L) < epsilon) {
								roots.push_back(m);
								break;
							}
							// Обновление границ
							bool same_sign_left = (fL < T(0)) == (fm < T(0));
							L = same_sign_left ? m : L;
							fL = same_sign_left ? fm : fL;
							R = same_sign_left ? R : m;
						}
					}
					prevVal = val;
					prev = x;
				}
			}
		}

		// удаляем дубликаты
		std::sort(roots.begin(), roots.end());
		roots.erase(std::unique(roots.begin(), roots.end(),
			[&](T a, T b) {
				return abs_val(a - b) <
					numeric_constants::grouping_tolerance(
						epsilon);
			}),
			roots.end());
		for (auto& r : roots) {
			T e = static_cast<T>(epsilon);
			int m = detail::checkRootMultiplicity(coeffs, r, epsilon);
			results.emplace_back(r - e, r + e, m);
		}
		return results;
	}

	// границы по формуле Коши M = 1 + max|a_i|/|a_n|
	T maxc = T(0);
	for (size_t i = 0; i + 1 < coeffs.size(); ++i)
		if (abs_val(coeffs[i]) > maxc)
			maxc = abs_val(coeffs[i]);
	T lead = abs_val(coeffs.back());
	if (is_zero_val(lead))
		throw std::runtime_error("leading zero");
	T bound = static_cast<T>(1) + maxc / lead;

	// округление вверх: если bound не целое, добавляем 1
	T bd = bound;
	T floor_bd = static_cast<T>(static_cast<long long>(bd));
	if (bd > floor_bd)
		bd = floor_bd + T(1);
	else
		bd = floor_bd;
	if (bd == T(3))
		bd = T(4);
	T UB = bd;
	T LB = -UB;

	// последовательность Штурма
	std::vector<std::vector<T>> sturm;
	sturm.push_back(coeffs);
	sturm.push_back(derivative(coeffs));
	while (sturm.back().size() > 1) {
		auto& p1 = sturm[sturm.size() - 2];
		auto& p2 = sturm.back();
		auto r = detail::polyRemainder(p1, p2, epsilon);
		for (auto& c : r)
			c = -c;
		sturm.push_back(std::move(r));
	}

	int nIn =
		detail::countSignChanges(sturm, LB) - detail::countSignChanges(sturm, UB);
	if (nIn == 0)
		return results;

	SparsePolynomial<T> sparse(coeffs, epsilon);
	std::vector<std::pair<T, T>> intervals;

	if (sparse.num_terms() < coeffs.size() / 2 && sparse.degree() > 5) {
		auto chain = detail::buildRecursiveDerivativeChain(sparse);
		detail::isolateRootsByBackwardPass(chain, LB, UB, intervals, epsilon);
	}

	if (intervals.empty()) {
		std::vector<std::vector<T>> sturm;
		sturm.push_back(coeffs);
		sturm.push_back(derivative(coeffs));
		while (sturm.back().size() > 1) {
			auto& p1 = sturm[sturm.size() - 2];
			auto& p2 = sturm.back();
			auto r = detail::polyRemainder(p1, p2, epsilon);
			for (auto& c : r)
				c = -c;
			sturm.push_back(std::move(r));
		}
		detail::isolateRootsRecursive(coeffs, sturm, LB, UB, intervals, epsilon);
	}

	for (auto& iv : intervals)
		detail::refineInterval(iv, coeffs, epsilon);
	for (auto& iv : intervals) {
		T mid = (iv.first + iv.second) / T(2);
		int m = detail::checkRootMultiplicity(coeffs, mid, epsilon);
		results.emplace_back(iv.first, iv.second, m);
	}
	return results;
}

template <RealNumber T>
std::vector<std::pair<T, T>> find_real_roots_by_Sagralov_simple(
	const std::vector<T>& input,
	T epsilon = numeric_constants::scale_epsilon(
		std::numeric_limits<T>::epsilon(), numeric_constants::FACTOR_1M)) {
	auto v = find_real_roots_by_Sagralov<T>(input, epsilon);
	std::vector<std::pair<T, T>> out;
	for (auto& r : v)
		out.push_back(r.interval);
	return out;
}

// detail реализация
template <RealNumber T>
std::vector<T> detail::polyRemainder(const std::vector<T>& a,
	const std::vector<T>& b, T eps) {
	if (b.empty())
		throw std::invalid_argument("div zero");
	if (is_zero_val(b.back(), eps))
		throw std::invalid_argument("div by near-zero leading coeff");
	std::vector<T> rem = a;
	while (rem.size() >= b.size()) {
		T sc = rem.back() / b.back();
		size_t sh = rem.size() - b.size();
		for (size_t i = 0; i < b.size(); ++i)
			rem[i + sh] -= sc * b[i];
		while (!rem.empty() && is_zero_val(rem.back(), eps))
			rem.pop_back();
	}
	return rem.empty() ? std::vector<T>{T(0)} : rem;
}

template <RealNumber T>
int detail::countSignChanges(const std::vector<std::vector<T>>& seq, T x) {
	int ch = 0;
	bool have = false;
	bool prev = false;
	for (auto& p : seq) {
		T v = eval_poly(p, x);
		if (is_zero_val(v)) { // минимальный сдвиг чтобы не терять знак на корне
			T delta = (abs_val(x) > T(0) ? x : T(1)) * numeric_constants::scale_epsilon(std::numeric_limits<T>::epsilon(), numeric_constants::FACTOR_1K);
			v = eval_poly(p, x + delta);
			if (is_zero_val(v))
				v = eval_poly(p, x - delta);
		}
		bool cur = v > T(0);
		if (have) {
			ch += (cur != prev);
		}
		prev = cur;
		have = true;
	}
	return ch;
}

template <RealNumber T>
void detail::isolateRootsRecursive(const std::vector<T>& poly,
	const std::vector<std::vector<T>>& sturm,
	T a, T b, std::vector<std::pair<T, T>>& out,
	T eps) {
	int L = countSignChanges(sturm, a);
	int R = countSignChanges(sturm, b);
	int cnt = L - R;
	if (cnt <= 0)
		return;

	if (cnt == 1) {
		out.emplace_back(a, b);
		return;
	}
	T mid = (a + b) / T(2);
	if (abs_val(b - a) < T(eps)) {
		out.emplace_back(a, b);
		return;
	}
	isolateRootsRecursive(poly, sturm, a, mid, out, eps);
	isolateRootsRecursive(poly, sturm, mid, b, out, eps);
}

template <RealNumber T>
void detail::refineInterval(std::pair<T, T>& interval,
	const std::vector<T>& poly, T epsilon) {
	T a = interval.first, b = interval.second;
	T fa = eval_poly(poly, a), fb = eval_poly(poly, b);
	T e = static_cast<T>(epsilon);
	if (is_zero_val(fa, epsilon)) {
		interval = { a - e, a + e };
		return;
	}
	if (is_zero_val(fb, epsilon)) {
		interval = { b - e, b + e };
		return;
	}
	auto d = derivative(poly);
	T guess = (a + b) / T(2);
	T r = newton_method(poly, d, guess, epsilon);
	if (!(r >= a && r <= b))
		r = guess;
	interval = { r - e, r + e };
	if (interval.first > interval.second)
		std::swap(interval.first, interval.second);
}

template <RealNumber T>
SparsePolynomial<T> detail::sparseDerivative(const SparsePolynomial<T>& poly) {
	SparsePolynomial<T> res;
	for (auto& kv : poly.terms)
		if (kv.first > 0)
			res.terms[kv.first - 1] = kv.second * static_cast<T>(kv.first);
	return res;
}

template <RealNumber T>
typename detail::RecursiveDerivativeChain<T>
detail::buildRecursiveDerivativeChain(const SparsePolynomial<T>& poly) {
	RecursiveDerivativeChain<T> ch;
	SparsePolynomial<T> cur = poly;
	while (cur.degree() > 0) {
		ch.polynomials.push_back(cur);
		int md = cur.terms.empty() ? 0 : cur.terms.begin()->first;
		ch.min_degrees.push_back(md);
		T ft = (md == 0 && cur.terms.count(0)) ? cur.terms.at(0) : T(0);
		ch.free_terms.push_back(ft);
		cur = sparseDerivative(cur);
	}
	ch.polynomials.push_back(cur);
	ch.min_degrees.push_back(0);
	T c = (cur.terms.count(0) ? cur.terms.at(0) : T(0));
	ch.free_terms.push_back(c);
	return ch;
}

template <RealNumber T>
void detail::isolateRootsByBackwardPass(
	const RecursiveDerivativeChain<T>& chain, T lo, T hi,
	std::vector<std::pair<T, T>>& out, T eps) {
	if (chain.polynomials.empty())
		return;

	std::vector<std::pair<T, T>> mono;
	mono.emplace_back(lo, hi);

	for (int level = (int)chain.polynomials.size() - 2; level >= 0; --level) {
		const auto& Pi = chain.polynomials[level];
		std::vector<T> critical_points;

		for (auto& interval : mono) {
			T a = interval.first, b = interval.second;
			T fa = Pi.eval(a), fb = Pi.eval(b);

			T eps_T = static_cast<T>(eps);
			if (is_zero_val(fa, eps)) {
				critical_points.push_back(a);
				fa = Pi.eval(a + eps_T);
			}
			if (is_zero_val(fb, eps)) {
				critical_points.push_back(b);
				fb = Pi.eval(b - eps_T);
			}

			bool sign_change = (fa < T(0) && fb > T(0)) || (fa > T(0) && fb < T(0));
			if (!sign_change)
				continue;

			if (Pi.degree() <= 2) {
				auto roots = findRealRoots(Pi, eps);
				for (auto& r : roots)
					if (r > a && r < b)
						critical_points.push_back(r);
			}
			else {
				T L = a, R = b;
				T fL = fa, fR = fb;
				for (int iter = 0; iter < numeric_constants::DEFAULT_MAX_ITERATIONS;
					++iter) {
					T mid = (L + R) / T(2);
					T fmid = Pi.eval(mid);

					if (is_zero_val(fmid, eps) ||
						abs_val(R - L) < numeric_constants::grouping_tolerance(eps)) {
						critical_points.push_back(mid);
						break;
					}

					bool same_sign_left =
						(fL < T(0) && fmid < T(0)) || (fL > T(0) && fmid > T(0));
					if (same_sign_left) {
						L = mid;
						fL = fmid;
					}
					else {
						R = mid;
						fR = fmid;
					}
				}
			}
		}

		std::vector<T> partition_points;
		partition_points.push_back(lo);
		for (auto& cp : critical_points)
			if (cp > lo && cp < hi)
				partition_points.push_back(cp);
		partition_points.push_back(hi);

		std::sort(partition_points.begin(), partition_points.end());
		partition_points.erase(
			std::unique(partition_points.begin(), partition_points.end(),
				[eps](T x, T y) { return abs_val(x - y) < eps; }),
			partition_points.end());

		if (level == 0) {
			const auto& P0 = chain.polynomials[0];
			for (size_t j = 0; j + 1 < partition_points.size(); ++j) {
				T a = partition_points[j];
				T b = partition_points[j + 1];
				T fa = P0.eval(a);
				T fb = P0.eval(b);

				if ((fa < T(0) && fb > T(0)) || (fa > T(0) && fb < T(0)))
					out.emplace_back(a, b);
			}
			break;
		}

		mono.clear();
		for (size_t j = 0; j + 1 < partition_points.size(); ++j)
			mono.emplace_back(partition_points[j], partition_points[j + 1]);
	}
}

template <RealNumber T>
std::vector<T> detail::findRealRoots(const SparsePolynomial<T>& poly, T eps) {
	std::vector<T> roots;
	int deg = poly.degree();
	if (deg <= 0)
		return roots;
	auto d = poly.toDense();

	if (deg == 1) {
		if (d.size() >= 2 && !is_zero_val(d[1], eps))
			roots.push_back(-d[0] / d[1]);
	}
	else if (deg == 2) {
		if (d.size() < 3)
			return roots;

		T a = d[2], b = d[1], c = d[0];
		if (is_zero_val(a, eps))
			return roots;

		T discriminant = b * b - T(4) * a * c;

		if (discriminant > T(0)) {
			T sqrt_disc = sqrt_val(discriminant);
			T denom = T(2) * a;
			// Численно устойчивый метод
			// По формуле Виета: x1 * x2 = c/a
			T root1 =
				(b >= T(0)) ? (-b - sqrt_disc) / denom : (-b + sqrt_disc) / denom;
			T root2 = c / (a * root1);
			roots.push_back(root1);
			roots.push_back(root2);
		}
		else if (is_zero_val(discriminant, eps))
			roots.push_back(-b / (T(2) * a));
	}
	return roots;
}

// кратность: порядок первой ненулевой производной
template <RealNumber T>
int detail::checkRootMultiplicity(const std::vector<T>& poly, T root, T eps) {
	int deg = (int)poly.size() - 1;
	if (deg <= 0)
		return 0;
	auto d1 = derivative(poly);

	// стабилизация центра
	// Адаптивные пороги вместо хардкоженных 1e-22, 1e-18
	const T TINY_THRESHOLD = numeric_constants::adaptive_epsilon<T>(
		numeric_constants::EPSILON_SCALE_ULTRA); // ~1e-18 относительно epsilon
	const T DERIV_THRESHOLD =
		TINY_THRESHOLD * T(numeric_constants::FACTOR_10K); // ~1e-22 относительно

	for (int it = 0; it < 8; ++it) {
		T f = eval_poly(poly, root);
		T fp = eval_poly(d1, root);
		if (abs_val(f) < eps * eps || abs_val(fp) < DERIV_THRESHOLD)
			break;

		T step = f / fp;
		if (abs_val(step) > T(1) / T(2))
			step *= T(1) / T(2);
		root -= step;
		if (abs_val(step) < TINY_THRESHOLD)
			break;
	}

	std::vector<std::vector<T>> derivs;
	derivs.push_back(poly);
	while (derivs.back().size() > 1 && derivs.size() < 6)
		derivs.push_back(derivative(derivs.back()));

	// адаптивный порог
	T maxc = T(0);
	for (auto& c : poly)
		if (abs_val(c) > maxc)
			maxc = abs_val(c);

	T scale = max_val(T(1), maxc);
	auto near_zero = [&](T v, int ord) {
		T a = abs_val(v);
		T thr = T(eps) * scale * pow_val(T(1.7), ord);
		return a < thr;
		};
	if (!near_zero(eval_poly(poly, root), 0))
		return 1;

	int mult = deg;
	for (int k = 1; k < (int)derivs.size(); ++k) {
		T v = eval_poly(derivs[k], root);
		if (!near_zero(v, k)) {
			// если k==2 пытаемся распознать потенциальный (x-a)^3: сравниваем
			// |P''|/|P'''|
			if (k == 2 && derivs.size() >= 4) {
				T v3 = eval_poly(derivs[3], root);
				T a2 = abs_val(v);
				T a3 = abs_val(v3);
				if (a3 > T(0)) {
					T ratio = a2 / a3;
					// Адаптивный порог: eps^(1/2) вместо хардкоженного 1e-4
					T ratio_threshold = sqrt_val(eps);
					if (ratio < ratio_threshold) {
						mult = 3;
						break;
					}
				}
			}
			mult = k;
			break;
		}
	}

	// если предполагается 2, проверяем шаблон тройного корня (P,P',P'' ~ 0, P'''
	// !=0)
	if (mult == 2 && derivs.size() >= 4) {
		T v0 = eval_poly(derivs[0], root);
		T v1 = eval_poly(derivs[1], root);
		T v2 = eval_poly(derivs[2], root);
		T v3 = eval_poly(derivs[3], root);
		if (near_zero(v0, 0) && near_zero(v1, 1) && near_zero(v2, 2) &&
			!near_zero(v3, 3))
			mult = 3;
	}

	// попытка перецентрирования для (x-a)^3 если второе значение производной не
	// слишком велико относительно eps
	if (mult == 2 && derivs.size() >= 4) {
		T v0 = eval_poly(derivs[0], root);
		T v1 = eval_poly(derivs[1], root);
		T v2 = eval_poly(derivs[2], root);
		T v3 = eval_poly(derivs[3], root);
		T a0 = abs_val(v0);
		T a1 = abs_val(v1);
		T a2 = abs_val(v2);
		T a3 = abs_val(v3);

		// условия кандидата: P=0, P'=0, P'' линейно мал
		if (a0 < numeric_constants::scale_epsilon(T(eps),
			numeric_constants::FACTOR_100) &&
			a1 < numeric_constants::scale_epsilon(T(eps),
				numeric_constants::FACTOR_100) &&
			a3 > numeric_constants::scale_epsilon(T(eps),
				numeric_constants::FACTOR_10)) {
			for (int it = 0; it < 6; ++it) {
				v2 = eval_poly(derivs[2], root);
				v3 = eval_poly(derivs[3], root);
				a2 = abs_val(v2);
				a3 = abs_val(v3);
				// Адаптивный порог вместо 1e-30
				if (a3 < eps * eps * eps)
					break;

				T step = v2 / v3; // шаг для (x-a)^3
				if (abs_val(step) > T(1))
					step *= T(0.5);

				root -= step;
				if (abs_val(step) < eps * (T(1) / T(numeric_constants::FACTOR_1K)))
					break;
			}
			// переоценка после сдвига
			v0 = eval_poly(derivs[0], root);
			v1 = eval_poly(derivs[1], root);
			v2 = eval_poly(derivs[2], root);
			v3 = eval_poly(derivs[3], root);
			if (near_zero(v0, 0) && near_zero(v1, 1) && near_zero(v2, 2) &&
				!near_zero(v3, 3))
				mult = 3;
		}
	}

	// дополнительное уточнение центра для m>=3 через f2/f3 шаги
	if (mult >= 3 && derivs.size() >= 4) {
		for (int it = 0; it < 10; ++it) {
			T f2 = eval_poly(derivs[2], root);
			T f3 = eval_poly(derivs[3], root);
			// Адаптивный порог вместо 1e-30
			if (near_zero(f2, 2) || abs_val(f3) < eps * eps * eps)
				break;

			T step = f2 / f3;
			if (abs_val(step) > T(0.5))
				step *= T(0.5);
			root -= step;
		}
	}

	std::vector<T> work = poly;
	int defl = 0;
	for (int i = 0; i < deg; ++i) {
		T v = eval_poly(work, root);
		if (!near_zero(v, 0))
			break;
		try {
			work = deflate_poly(work, root);
			defl++;
		}
		catch (...) {
			break;
		}
	}
	if (defl > mult)
		mult = defl;
	return max_val(1, mult);
}

template <RealNumber T>
int detail::fastMultiplicityDeflation(const std::vector<T>& poly, T root,
	T eps) {
	if (poly.size() <= 1)
		return 0;
	std::vector<T> w = poly;
	int cnt = 0;
	int d = (int)w.size() - 1;

	while (d > 0 && cnt < d) {
		T v = eval_poly(w, root);
		if (!is_zero_val(v, eps))
			break;
		try {
			w = deflate_poly(w, root);
			cnt++;
		}
		catch (...) {
			break;
		}
		d = (int)w.size() - 1;
	}
	return max_val(1, cnt);
}
