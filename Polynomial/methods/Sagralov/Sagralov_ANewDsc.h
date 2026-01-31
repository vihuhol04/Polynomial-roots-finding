// Алгоритм ANewDsc, автор Мельхорн Сагралов
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#include "mathUtils.h"             
#include "polynomialUtils.h"       
#include "Sagralov_real_roots.h"   

#include <optional>               

// структура интервала для метода ANewDsc
// отличается от обычного интервала тем, что хранит число изменений знака, вычисленное по правилу Декарта для полиномиального преобразования
template <typename T>
struct Interval01 {
    T left;          // левая граница интервала
    T right;         // правая граница интервала
    int variability; // количество возможных корней внутри интервала

    // конструктор: по умолчанию количество возможных корней внутри интервала = -1 (ещё не вычислено)
    Interval01(T l, T r, int v = -1) : left(l), right(r), variability(v) {}

    // ширина интервала
    T width() const { return abs_val(right - left); }

    // середина интервала
    T midpoint() const { return (left + right) / T(2); }
};

namespace anewdsc_detail {
    // Декартово преобразование интервала [a, b]
    // преобразует исходный полином P(x) на интервале [a,b] в новый Q(y), где корни в [a,b] отображаются в [0, +бесконечность) 
    // используется для 01-теста
    template <typename T>
    std::vector<T> descartes_transform(const std::vector<T>& poly, T a, T b) {
        // преобразование Декарта:
        // x = (a*y + b) / (y + 1)
        // Q(y) = (y+1)^n * P((a*y + b) / (y + 1))

        int n = (int)poly.size() - 1; // степень полинома
        if (n < 0) return {}; 

        // итоговый полином степени n
        std::vector<T> result(n + 1, T(0));

        // Q(y) = сумма poly[i] * (a*y+b)^i * (y+1)^(n-i)
        for (int i = 0; i <= n; ++i) {
            // (a*y + b)^i
            std::vector<T> numerator(i + 1, T(0));
            numerator[0] = T(1);
            for (int k = 0; k < i; ++k) {
                for (int j = k + 1; j >= 1; --j)
                    numerator[j] = numerator[j] * b + numerator[j - 1] * a;
                numerator[0] *= b;
            }

            // (y + 1)^(n - i) 
            std::vector<T> denominator(n - i + 1, T(0));
            denominator[0] = T(1);
            for (int k = 0; k < n - i; ++k) {
                for (int j = k + 1; j >= 1; --j)
                    denominator[j] = denominator[j] + denominator[j - 1];
                // denominator[0] остаётся 1
            }

            // перемножаем обе части и добавляем в итоговый полином
            for (int p = 0; p <= i; ++p) 
                for (int q = 0; q <= n - i; ++q)
                    result[p + q] += poly[i] * numerator[p] * denominator[q] * T(p + q <= n);
        }
        return result;
    }

    // подсчёт изменений знака в последовательности коэффициентов
    // используется правило Декарта: количество положительных корней не больше числа изменений знака
    template <typename T>
    int count_sign_changes(const std::vector<T>& coeffs, T eps = 1e-12)
    {
        int changes = 0;
        bool have_sign = false; // были ли ненулевые коэффициенты ранее
        bool last_sign = false; // знак предыдущего коэффициента

        for (const auto& c : coeffs) {
            if (is_zero_val(c, eps)) continue;

            // знак текущего коэффициента
            bool current_sign = (c > T(0));
            // если знак изменился — увеличиваем счётчик
            changes = 1 + int(have_sign && current_sign != last_sign);
            // обновляем 
            last_sign = current_sign;
            have_sign = true;
        }
        return changes;
    }

    // 01-тест из метода ANewDsc
    // проверяет, сколько корней у полинома в интервале [a,b]:
    // 0 — нет корней
    // 1 — ровно один корень
    // -1 — неопределённо (несколько корней или возможны кратные)
    template <typename T>
    int test_01(const std::vector<T>& poly, T a, T b, T eps = 1e-10) {
        // если интервал слишком мал — считаем результат неопределённым
        if (abs_val(b - a) < T(eps)) return -1;

        // вычисляем значение полинома на концах интервала
        T fa = eval_poly(poly, a);
        T fb = eval_poly(poly, b);

        // если хотя бы на одной границе есть точный корень — считаем, что корень есть
        if (is_zero_val(fa, eps) || is_zero_val(fb, eps))
            return 1;

        // выполняем декартово преобразование к виду Q(y) на [0, +бесконечность)
        auto transformed = descartes_transform(poly, a, b);
        int v = count_sign_changes(transformed, eps);

        if (v == 0) return 0;
        if (v == 1) return 1;

        return -1;
    }

    // поиск допустимой точки внутри интервала, где значение полинома максимально далеко от нуля
    template <typename T>
    std::optional<T> find_admissible_point(const std::vector<T>& poly, T a, T b, int grid_size = 20, T eps = 1e-10) {
        T mid = (a + b) / T(2);
        T step = (b - a) / T(grid_size);

        // ищем точку внутри интервала, где |P(x)| максимально
        T max_val = T(0);
        T best_point = mid;

        for (int i = 0; i <= grid_size; ++i) {
            T x = a + step * T(i);
            T val = abs_val(eval_poly(poly, x));

            if (val > max_val) {
                max_val = val;
                best_point = x;
            }
        }

        // используем только точки, где значение достаточно далеко от 0
        T threshold = max_val / T(4);

        for (int i = grid_size / 2 - 2; i <= grid_size / 2 + 2; ++i) {
            if (i < 0 || i > grid_size) continue;

            T x = a + step * T(i);
            T val = abs_val(eval_poly(poly, x));

            if (val >= threshold) return x;
        }

        return mid; // если ничего лучшего нет — делим по середине
    }

    // тест Ньютона
    template <typename T>
    std::optional<Interval01<T>> newton_test(const std::vector<T>& poly, const Interval01<T>& interval, int N_param, T eps = 1e-10) {
        // применяем тест Ньютона только если внутри интервала более одного возможного корня
        if (interval.variability <= 1)
            return std::nullopt;

        T width = interval.width();
        T target_width = width / T(N_param); // новая предполагаемая ширина подынтервала

        auto deriv = derivative(poly);

        // берем три точки внутри интервала: четверть, середина, три четверти, чтобы сделать надежные итерации Ньютона
        std::vector<T> sample_points = { interval.left + width * T(0.25), interval.midpoint(), interval.left + width * T(0.75) };

        T avg_center = T(0);
        int valid_samples = 0;

        for (const auto& xi : sample_points) {
            T f_val = eval_poly(poly, xi);
            T fp_val = eval_poly(deriv, xi);

            // шаг Ньютона возможен только если f'(x) != 0
            if (!is_zero_val(fp_val, eps)) {
                T correction = f_val / fp_val;
                T candidate = xi - correction;

                // корректная точка Ньютона должна лежать в интервале
                avg_center += candidate * T(candidate >= interval.left && candidate <= interval.right);
                valid_samples = 1 * int(candidate >= interval.left && candidate <= interval.right);
            }
        }

        if (valid_samples == 0) return std::nullopt;

        // центрируем новый интервал вокруг среднего найденного значения
        T lambda = avg_center / T(valid_samples);
        T new_left = lambda - target_width / T(2);
        T new_right = lambda + target_width / T(2);

        // ограничиваем новыми границами внутри исходного интервала
        new_left = max_val(interval.left, new_left);
        new_right = min_val(interval.right, new_right);

        // снаружи корней быть не должно
        int test_left = test_01(poly, interval.left, new_left, eps);
        int test_right = test_01(poly, new_right, interval.right, eps);

        // подтверждение, что новый интервал — изолирующий
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

        return std::nullopt; // тест Ньютона не работает
    }
}

template <typename T>
std::vector<std::pair<T, T>> find_real_roots_ANewDsc(const std::vector<T>& input, T epsilon = 1e-8, int max_iterations = 1000) {
    if (input.empty() || input.size() == 1) return {};
    std::vector<T> poly = normalizePolynomial(input);

    // оценка границы корней по правилу Коши
    T maxc = T(0);
    for (size_t i = 0; i + 1 < poly.size(); ++i)
        maxc = max_val(maxc, abs_val(poly[i]));

    T lead = abs_val(poly.back());
    if (is_zero_val(lead))
        throw std::runtime_error("zero leading coefficient");

    T bound = T(1) + maxc / lead; // все действительные корни лежат в [-bound, bound]

    std::vector<Interval01<T>> active_intervals;
    active_intervals.emplace_back(-bound, bound, -1);

    // вектор уже изолированных интервалов
    std::vector<std::pair<T, T>> isolated_roots;
    int iteration = 0;

    while (!active_intervals.empty() && iteration < max_iterations) {
        ++iteration;

        // берем последний добавленный интервал 
        Interval01<T> current = active_intervals.back();
        active_intervals.pop_back();

        if (current.width() < T(epsilon)) {
            isolated_roots.emplace_back(current.left, current.right);
            continue;
        }

        // выполняем 01-тест
        int test_result = anewdsc_detail::test_01(poly, current.left, current.right, epsilon);

        if (test_result == 0) continue;

        if (test_result == 1) {
            // ровно один корень по варианту Декарта/границам, но проверяем ширину
            if (current.width() < T(epsilon)) {
                isolated_roots.emplace_back(current.left, current.right);
                continue;
            }
            test_result = -1;
        }

        // ускоряем поиск с помощью Ньютон-теста (если есть информация о количестве перемен знака)
        if (current.variability > 1) {
            int N_param = 4; // параметр, задающий, насколько узким мы хотим сделать новый интервал
            auto newton_result = anewdsc_detail::newton_test(poly, current, N_param, epsilon);

            if (newton_result) {
                active_intervals.push_back(*newton_result);
                continue;
            }
        }

        // если Ньютон-тест не помог — проводим бисекцию
        auto m_opt = anewdsc_detail::find_admissible_point(poly, current.left, current.right, 20, epsilon);
        T m = m_opt.value_or(current.midpoint());

        int left_test = anewdsc_detail::test_01(poly, current.left, m, epsilon);
        int right_test = anewdsc_detail::test_01(poly, m, current.right, epsilon);

        if (left_test != 0)
            active_intervals.emplace_back(current.left, m, -1);
        if (right_test != 0)
            active_intervals.emplace_back(m, current.right, -1);
    }

    if (isolated_roots.size() <= 1) return isolated_roots; 

    // сортируем для объединения перекрывающихся интервалов
    std::sort(isolated_roots.begin(), isolated_roots.end(), [](const std::pair<T, T>& a, const std::pair<T, T>& b) { return a.first < b.first; });

    std::vector<std::pair<T, T>> merged;
    merged.push_back(isolated_roots[0]);

    for (size_t i = 1; i < isolated_roots.size(); ++i) {
        auto& last = merged.back();
        auto& current = isolated_roots[i];

        // если текущий интервал начинается до конца последнего (с учётом epsilon) — считаем их перекрывающимися/соседними и объединяем
        if (current.first <= last.second + epsilon)
            last.second = max_val(last.second, current.second); 
        else
            // не перекрывается — просто добавляем
            merged.push_back(current);
    }
    return merged;
}
