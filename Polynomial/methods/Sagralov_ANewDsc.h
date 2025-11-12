#include "mathUtils.h"             
#include "polynomialUtils.h"       
#include "Sagralov_real_roots.h"   

#include <optional>               

// Структура интервала для метода ANewDsc
// Отличается от обычного интервала тем, что хранит "variability" — вариабельность (число изменений знака),
// вычисленное по правилу Декарта для полиномиального преобразования.
template <typename T>
struct Interval01
{
    T left;          // Левая граница интервала
    T right;         // Правая граница интервала
    int variability; // Вариабельность (количество возможных корней внутри интервала)

    // Конструктор: по умолчанию variability = -1 (ещё не вычислено)
    Interval01(T l, T r, int v = -1) : left(l), right(r), variability(v) {}

    // Ширина интервала
    T width() const { return abs_val(right - left); }

    // Середина интервала
    T midpoint() const { return (left + right) / T(2); }
};

namespace anewdsc_detail
{
    // -------------------------------
    // Декартово преобразование интервала [a, b]
    // Преобразует исходный полином P(x) на интервале [a,b] в новый Q(y),
    // где корни в [a,b] отображаются в [0, +∞). Используется для 01-теста.
    // -------------------------------
    template <typename T>
    std::vector<T> descartes_transform(const std::vector<T>& poly, T a, T b)
    {
        // Преобразование Декарта:
        // x = (a*y + b) / (y + 1)
        // Q(y) = (y+1)^n * P((a*y + b) / (y + 1))

        int n = (int)poly.size() - 1; // Степень полинома
        if (n < 0)
            return {}; // Пустой полином

        // Результирующий полином степени n
        std::vector<T> result(n + 1, T(0));

        // Разворачиваем формулу: Q(y) = Σ poly[i] * (a*y+b)^i * (y+1)^(n-i)
        for (int i = 0; i <= n; ++i)
        {
            // (a*y + b)^i — биномиальное разложение
            std::vector<T> numerator(i + 1, T(0));
            numerator[0] = T(1);
            for (int k = 0; k < i; ++k)
            {
                for (int j = k + 1; j >= 1; --j)
                    numerator[j] = numerator[j] * b + numerator[j - 1] * a;
                numerator[0] *= b;
            }

            // (y + 1)^(n - i) — биномиальное разложение
            std::vector<T> denominator(n - i + 1, T(0));
            denominator[0] = T(1);
            for (int k = 0; k < n - i; ++k)
            {
                for (int j = k + 1; j >= 1; --j)
                    denominator[j] = denominator[j] + denominator[j - 1];
                // denominator[0] остаётся 1
            }

            // Перемножаем обе части и добавляем в итоговый полином
            for (int p = 0; p <= i; ++p)
            {
                for (int q = 0; q <= n - i; ++q)
                {
                    if (p + q <= n)
                        result[p + q] += poly[i] * numerator[p] * denominator[q];
                }
            }
        }

        return result;
    }

    // Подсчёт изменений знака в последовательности коэффициентов
    // Используется правило Декарта: количество положительных корней не больше числа изменений знака.
    template <typename T>
    int count_sign_changes(const std::vector<T>& coeffs, T eps = 1e-12)
    {
        int changes = 0;
        bool have_sign = false; // Были ли ненулевые коэффициенты ранее
        bool last_sign = false; // Знак предыдущего коэффициента

        for (const auto& c : coeffs)
        {
            // Пропускаем нулевые коэффициенты
            if (is_zero_val(c, eps))
                continue;

            // Определяем знак текущего коэффициента
            bool current_sign = (c > T(0));

            // Если знак изменился — увеличиваем счётчик
            if (have_sign && current_sign != last_sign)
                ++changes;

            // Обновляем состояние
            last_sign = current_sign;
            have_sign = true;
        }

        return changes;
    }

    // Основной "01-тест" из метода ANewDsc.
    // Проверяет, сколько корней у полинома в интервале [a,b]:
    // 0 — нет корней
    // 1 — ровно один корень
    // -1 — неопределённо (несколько корней или возможны кратные)
    template <typename T>
    int test_01(const std::vector<T>& poly, T a, T b, T eps = 1e-10)
    {
        // Если интервал слишком мал — считаем результат неопределённым
        if (abs_val(b - a) < T(eps))
            return -1;

        // Вычисляем значение полинома на концах интервала
        T fa = eval_poly(poly, a);
        T fb = eval_poly(poly, b);

        // Если хотя бы на одной границе есть точный корень — считаем, что корень есть
        if (is_zero_val(fa, eps))
            return 1;
        if (is_zero_val(fb, eps))
            return 1;

        // Выполняем декартово преобразование к виду Q(y) на [0, +∞)
        auto transformed = descartes_transform(poly, a, b);

        // Считаем число изменений знака — это вариабельность
        int v = count_sign_changes(transformed, eps);

        // v == 0: корней нет на этом интервале
        if (v == 0)
            return 0;

        // v == 1: ровно один корень — интервал изолирующий
        if (v == 1)
            return 1;

        // Если знаки меняются >= 2 раз — может быть несколько корней
        // Надо сужать интервал (бисекция)
        return -1;
    }

    // Поиск "допустимой точки" внутри интервала, где значение полинома максимально далеко от нуля.
// Нужна для стабильной бисекции, чтобы не делить интервал точно по середине, если там f(x) ≈ 0.
    template <typename T>
    std::optional<T> find_admissible_point(const std::vector<T>& poly, T a, T b, int grid_size = 20, T eps = 1e-10)
    {
        T mid = (a + b) / T(2);
        T step = (b - a) / T(grid_size);

        // Ищем точку внутри интервала, где |P(x)| максимально — чтобы делить "умно"
        T max_val = T(0);
        T best_point = mid;

        for (int i = 0; i <= grid_size; ++i)
        {
            T x = a + step * T(i);
            T val = abs_val(eval_poly(poly, x));

            if (val > max_val)
            {
                max_val = val;
                best_point = x;
            }
        }

        // Используем только точки, где значение достаточно далеко от 0,
        // то есть избегаем «фальшивых» почти корней
        T threshold = max_val / T(4);

        for (int i = grid_size / 2 - 2; i <= grid_size / 2 + 2; ++i)
        {
            if (i < 0 || i > grid_size)
                continue;

            T x = a + step * T(i);
            T val = abs_val(eval_poly(poly, x));

            if (val >= threshold)
                return x; // Нашли осмысленную точку для деления
        }

        return mid; // Если ничего лучшего нет — делим по середине
    }

    // Ньютонов тест — ускорение сходимости:
    // Пробуем сделать один шаг Ньютона, чтобы быстро изолировать корень,
    // если интервал узкий и вариабельность > 1
    template <typename T>
    std::optional<Interval01<T>> newton_test(const std::vector<T>& poly, const Interval01<T>& interval, int N_param, T eps = 1e-10)
    {
        // Применяем тест Ньютона только если внутри интервала более одного возможного корня
        if (interval.variability <= 1)
            return std::nullopt;

        T width = interval.width();
        T target_width = width / T(N_param); // новая предполагаемая ширина подынтервала

        // Производная нужна для метода Ньютона
        auto deriv = derivative(poly);

        // Берём три точки внутри интервала: четверть, середина, три четверти
        // чтобы сделать надёжные итерации Ньютона
        std::vector<T> sample_points = {
            interval.left + width * T(0.25),
            interval.midpoint(),
            interval.left + width * T(0.75) };

        T avg_center = T(0);
        int valid_samples = 0;

        for (const auto& xi : sample_points)
        {
            T f_val = eval_poly(poly, xi);
            T fp_val = eval_poly(deriv, xi);

            // Шаг Ньютона возможен только если f'(x) != 0
            if (!is_zero_val(fp_val, eps))
            {
                T correction = f_val / fp_val;
                T candidate = xi - correction;

                // Корректная точка Ньютона должна лежать в интервале
                if (candidate >= interval.left && candidate <= interval.right)
                {
                    avg_center += candidate;
                    ++valid_samples;
                }
            }
        }

        if (valid_samples == 0)
            return std::nullopt; // Ничего не удалось улучшить

        // Центрируем новый интервал вокруг среднего найденного значения
        T lambda = avg_center / T(valid_samples);
        T new_left = lambda - target_width / T(2);
        T new_right = lambda + target_width / T(2);

        // Ограничиваем новыми границами внутри исходного интервала
        new_left = max_val(interval.left, new_left);
        new_right = min_val(interval.right, new_right);

        // Проверяем края: снаружи корней быть не должно
        int test_left = test_01(poly, interval.left, new_left, eps);
        int test_right = test_01(poly, new_right, interval.right, eps);

        // Подтверждение, что новый интервал — изолирующий
        if (test_left == 0 && test_right == 0)
        {
            int test_center = test_01(poly, new_left, new_right, eps);
            if (test_center >= 0)
                return Interval01<T>(new_left, new_right, test_center);
        }

        // Дополнительная проверка — делим интервал на крылья
        T delta = target_width / T(2);
        if (test_01(poly, interval.left + delta, interval.right, eps) == 0)
            return Interval01<T>(interval.left, interval.left + delta);

        if (test_01(poly, interval.left, interval.right - delta, eps) == 0)
            return Interval01<T>(interval.right - delta, interval.right);

        return std::nullopt; // Тест Ньютона не работает, нужен обычный спуск
    }
}

template <typename T>
std::vector<std::pair<T, T>> find_real_roots_ANewDsc(const std::vector<T>& input, T epsilon = 1e-8, int max_iterations = 1000)
{
    // Если полином пустой или константа — корней нет
    if (input.empty() || input.size() == 1)
        return {};

    // Нормализуем полином (удаляем ведущие нули, возможно деление на коэффициент и т.п.)
    // Ожидается, что normalizePolynomial вернёт вектора коэффициентов от свободного члена к старшим степеням.
    std::vector<T> poly = normalizePolynomial(input);

    // Оценка границы корней по неравенству типа Cauchy:
    // bound = 1 + max_{i < n} |a_i| / |a_n|
    T maxc = T(0);
    for (size_t i = 0; i + 1 < poly.size(); ++i)
        maxc = max_val(maxc, abs_val(poly[i]));

    T lead = abs_val(poly.back());
    if (is_zero_val(lead))
        throw std::runtime_error("zero leading coefficient"); // Неправильный вход: ведущий коэффициент не должен быть нулём

    T bound = T(1) + maxc / lead; // Все действительные корни лежат в [-bound, bound]

    // Начинаем из одного активного интервала [-bound, bound]
    std::vector<Interval01<T>> active_intervals;
    active_intervals.emplace_back(-bound, bound, -1);

    // Сюда будем собирать уже изолированные (узкие) интервалы, содержащие корни
    std::vector<std::pair<T, T>> isolated_roots;
    int iteration = 0;

    // Основной рабочий цикл: пока есть активные интервалы и не превышен лимит итераций
    while (!active_intervals.empty() && iteration < max_iterations)
    {
        ++iteration;

        // Берём последний добавленный интервал (стекоподобная стратегия)
        Interval01<T> current = active_intervals.back();
        active_intervals.pop_back();

        // Если интервал уже достаточно узкий — считаем его изолированным корнем
        if (current.width() < T(epsilon)) 
        {
            isolated_roots.emplace_back(current.left, current.right);
            continue;
        }

        // Выполняем 01-тест: возвращает 0 (нет корней), 1 (ровно один), -1 (неопределённо)
        int test_result = anewdsc_detail::test_01(poly, current.left, current.right, epsilon);

        if (test_result == 0)
            continue; // Корней нет в этом интервале -> пропускаем

        if (test_result == 1)
        {
            // Ровно один корень по варианту Декарта/границам, но проверяем ширину
            if (current.width() < T(epsilon))
            {
                // Достаточно узкий — добавляем в результаты
                isolated_roots.emplace_back(current.left, current.right);
                continue;
            }
            // Если интервал широк — нельзя принять его за окончательный, продолжаем дробление
            test_result = -1;
        }

        // Если test_result == -1 — неопределённость (возможно несколько корней или кратные)
        // Попробуем сначала ускорить поиск с помощью Ньютон-теста (если есть информация о вариабельности)
        if (current.variability > 1)
        {
            int N_param = 4; // Параметр, задающий, насколько узким мы хотим сделать новый интервал
            auto newton_result = anewdsc_detail::newton_test(poly, current, N_param, epsilon);

            if (newton_result)
            {
                // Если Ньютон-тест дал изолирующий подинтервал — используем его
                active_intervals.push_back(*newton_result);
                continue;
            }
        }

        // Если Ньютон-тест не помог или variability не > 1 — проводим бисекцию.
        // Для устойчивости сначала пытаемся найти "допустимую точку" для деления (avoid near-zero middle).
        auto m_opt = anewdsc_detail::find_admissible_point(poly, current.left, current.right, 20, epsilon);
        T m = m_opt.value_or(current.midpoint());

        // Проверяем, какие из подинтервалов могут содержать корни (по 01-тесту)
        int left_test = anewdsc_detail::test_01(poly, current.left, m, epsilon);
        int right_test = anewdsc_detail::test_01(poly, m, current.right, epsilon);

        // Добавляем в стек только те подинтервалы, которые потенциально содержат корни
        if (left_test != 0)
            active_intervals.emplace_back(current.left, m, -1);
        if (right_test != 0)
            active_intervals.emplace_back(m, current.right, -1);
    }

    // После обработки стек может содержать несколько найденных изолированных интервалов.
    // Теперь нужно объединить/дедуплицировать перекрывающиеся результаты.

    if (isolated_roots.size() <= 1)
        return isolated_roots; // Нечего объединять

    // Сортируем по левому концу для последовательного объединения перекрывающихся интервалов
    std::sort(isolated_roots.begin(), isolated_roots.end(),
        [](const std::pair<T, T>& a, const std::pair<T, T>& b)
        {
            return a.first < b.first;
        });

    std::vector<std::pair<T, T>> merged;
    merged.push_back(isolated_roots[0]);

    for (size_t i = 1; i < isolated_roots.size(); ++i)
    {
        auto& last = merged.back();
        auto& current = isolated_roots[i];

        // Если текущий интервал начинается до конца последнего (с учётом epsilon) —
        // считаем их перекрывающимися/соседними и объединяем
        if (current.first <= last.second + epsilon)
        {
            last.second = max_val(last.second, current.second); // расширяем правую границу
        }
        else
        {
            // Не перекрывается — просто добавляем
            merged.push_back(current);
        }
    }

    // Возвращаем объединённые интервалы-результаты
    return merged;
}