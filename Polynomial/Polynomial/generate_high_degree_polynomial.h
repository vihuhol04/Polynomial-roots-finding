#ifndef GEN_HIGH_DEGREE_POLY_FIXED_H
#define GEN_HIGH_DEGREE_POLY_FIXED_H

#include "mathUtils.h"
#include "polynomialUtils.h"

/*
Генератор многочлена высокой степени:
    Идеи и гарантии:
        1. Все вещественные корни лежат строго в [-1; 1]
        2. Пользователь указывает только количество кластеров и их радиусы
    Смещение от границ задавать не нужно, так как центры равномерно автоматически распределяются
    Радиусы могут автоматически уменьшиться, чтобы
        1. Кластеры остались в нужном промежутке
        2. Кластеры не пересеклись друг с другом
    Кластеры могут иметь разные радиусы
    Если радиусов меньше, чем кластеров, используется радиус по умолчанию
    Кратные корни и стандартные вещественные распределяются равномерно распределяются по свободным интервалам между кластерами
    Комплексные корни создаются сопряженными парами

На выходе имеем:
    Вектор коэффициентов
    Вектор вещественных корней с учетом кратности
    Вектор уникальных вещественных корней
    Вектор кратностей для каждого уникального корня
    Вектор сопряженных пар комплексных корней
*/
// Реализация: Павлова Анастасия, КМБО-01-22

// умножение на (x - r)
template<typename T>
void multiply_by_root(std::vector<T>& coeffs, const T& root) {
    std::vector<T> res(coeffs.size() + 1, T(0));
    for (size_t i = 0; i < coeffs.size(); ++i) {
        res[i] += coeffs[i] * (-root);
        res[i + 1] += coeffs[i];
    }
    coeffs.swap(res);
}

// умножение на квадратичный многочлен x^2 + b*x + c
template<typename T>
void multiply_by_quadratic(std::vector<T>& coeffs, const T& b, const T& c) {
    std::vector<T> res(coeffs.size() + 2, T(0));
    for (size_t i = 0; i < coeffs.size(); ++i) {
        res[i] += coeffs[i] * c;
        res[i + 1] += coeffs[i] * b;
        res[i + 2] += coeffs[i];
    }
    coeffs.swap(res);
}

// равномерно распределить N точек в [a,b]
// используется для генерации корней в кластере или свободных интервалах
template<typename T>
std::vector<T> evenly_spaced_in_interval(T a, T b, unsigned N) {
    std::vector<T> res;
    if (N == 0) return res;
    if (a > b) std::swap(a, b);
    if (N == 1) {
        res.push_back((a + b) / T(2));
        return res;
    }
    T step = (b - a) / T(N + 1);
    for (unsigned i = 1; i <= N; ++i) res.push_back(a + step * T(i));
    return res;
}

// вспомогательный аллокатор целых количеств по длинам интервалов 
// на вход: общее число элементов и массив длин интервалов
// на выход: сколько элементов отдать каждому интервалу (пропорционально длине)
static std::vector<unsigned> allocate_counts_proportional(unsigned total_items, const std::vector<double>& lengths) {
    std::vector<unsigned> counts(lengths.size(), 0);
    double total_len = 0.0;
    for (double L : lengths) total_len += L;
    if (lengths.empty() || total_len <= 0.0) {
        counts[0] = total_items * (!counts.empty());
        return counts;
    }
    std::vector<double> desired(lengths.size(), 0.0);
    std::vector<std::pair<double, unsigned>> fracs;
    unsigned sum = 0;
    for (size_t i = 0; i < lengths.size(); ++i) {
        double d = (lengths[i] / total_len) * static_cast<double>(total_items);
        unsigned f = static_cast<unsigned>(std::floor(d));
        counts[i] = f;
        sum += f;
        desired[i] = d - f; // fractional part
        fracs.emplace_back(desired[i], static_cast<unsigned>(i));
    }
    unsigned remain = (sum >= total_items) ? 0u : (total_items - sum);
    // распределяем остатки по наибольшим дробным частям
    std::sort(fracs.begin(), fracs.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    for (unsigned k = 0; k < remain; ++k) counts[fracs[k].second]++;
    return counts;
}

// Главная функция генерации
// Выходные параметры:
//  - coefficients: коэффициенты многочлена
//  - real_roots_repeated: все вещественные корни с учётом кратности (для умножения)
//  - unique_real_roots: уникальные вещественные корни (без повторений)
//  - real_root_multiplicities: кратности для unique_real_roots
//  - complex_roots_out: сгенерированные комплексные корни (включая сопряженные)

template<typename T>
void generate_high_degree_polynomial(
    unsigned P,
    unsigned num_complex_pairs,
    unsigned num_clusters,
    const std::vector<unsigned>& cluster_counts_in,   // простые корни в каждом кластере (может быть пустым)
    const std::vector<T>& cluster_radii_in,           // радиусы кластеров (может быть пустым)
    const std::vector<std::pair<unsigned, unsigned>>& multiplicity_groups, // (mult, count)
    T default_cluster_radius,
    bool normalize_coeffs,
    std::uint64_t seed,
    std::vector<T>& coefficients,
    std::vector<T>& real_roots_repeated,
    std::vector<T>& unique_real_roots,
    std::vector<unsigned>& real_root_multiplicities,
    std::vector<std::complex<T>>& complex_roots_out
) {
    // В начале generate_high_degree_polynomial_fixed<T>(...)
    if (P == 0) {
        throw std::invalid_argument("Степень многочлена P должна быть >= 1");
    }

    // Проверка радиусов
    for (size_t i = 0; i < cluster_radii_in.size(); ++i) {
        if (cluster_radii_in[i] <= T(0)) {
            throw std::invalid_argument("Радиус кластера должен быть положительным (cluster " + std::to_string(i) + ")");
        }
    }

    // Проверка суммарного радиуса кластеров
    T total_radius = std::accumulate(cluster_radii_in.begin(), cluster_radii_in.end(), T(0));
    if (total_radius >= T(1)) {
        throw std::invalid_argument("Сумма радиусов всех кластеров не должна превышать 1");
    }

    // инициализация случайного генератора
    std::random_device rd;
    std::uint64_t used_seed = (seed != 0) ? seed : static_cast<std::uint64_t>(rd());
    std::mt19937_64 gen(used_seed);

    // подготовка векторов кластеров и радиусов
    // если входные векторы не совпадают num_clusters по размеру, заполняем радиусами по умолчанию
    std::vector<unsigned> cluster_counts;
    if (num_clusters == 0) cluster_counts.clear();
    else if (cluster_counts_in.size() == num_clusters) cluster_counts = cluster_counts_in;
    else cluster_counts.assign(num_clusters, 1u);

    std::vector<T> cluster_radii;
    if (num_clusters == 0) cluster_radii.clear();
    else if (cluster_radii_in.size() == num_clusters) cluster_radii = cluster_radii_in;
    else cluster_radii.assign(num_clusters, default_cluster_radius);

    // копия групп кратностей
    std::vector<std::pair<unsigned, unsigned>> local_mult_groups = multiplicity_groups;

    // вычисления степеней (сколько уже занято)
    unsigned complex_degree = 2 * num_complex_pairs;
    if (complex_degree > P) throw std::invalid_argument("Недостаточная степень P для заданного числа пар комплексных корней");

    unsigned mult_degree = 0;
    unsigned mult_root_count = 0;
    for (auto& mc : local_mult_groups) {
        unsigned mult = mc.first; unsigned cnt = mc.second;
        if (mult < 1) continue;
        // защита от слишком большой кратности
        if (mult > P) throw std::invalid_argument("Слишком большая кратность в multiplicity_groups");
        mult_degree += mult * cnt;
        mult_root_count += cnt;
    }

    unsigned cluster_simple_roots = 0;
    for (auto c : cluster_counts) cluster_simple_roots += c;

    unsigned used_degree = complex_degree + mult_degree + cluster_simple_roots;

    // если не влезает — уменьшаем сначала кластерные простые корни, затем кратные группы
    if (used_degree > P) {
        unsigned overflow = used_degree - P;
        // уменьшаем cluster_counts по очереди, начиная с последних
        for (int i = static_cast<int>(cluster_counts.size()) - 1; i >= 0 && overflow > 0; --i) {
            unsigned rem = std::min<unsigned>(cluster_counts[i], overflow);
            cluster_counts[i] -= rem;
            overflow -= rem;
        }
        // пересчитываем
        cluster_simple_roots = 0; for (auto c : cluster_counts) cluster_simple_roots += c;
        used_degree = complex_degree + mult_degree + cluster_simple_roots;

        if (used_degree > P) {
            overflow = used_degree - P;
            // теперь уменьшаем multiplicity группы: удаляем целые корни (каждый удалённый корень уменьшает используемую степень на mult)
            for (auto& mc : local_mult_groups) {
                if (overflow == 0) break;
                unsigned mult = mc.first; unsigned& cnt = mc.second;
                if (cnt == 0) continue;
                // сколько корней этой группы нужно удалить (ceil(overflow / mult)) но не больше, чем cnt
                unsigned need = static_cast<unsigned>((overflow + mult - 1) / mult);
                unsigned rem = std::min<unsigned>(cnt, need);
                cnt -= rem;
                overflow = (rem * mult >= overflow) ? 0 : (overflow - rem * mult);
            }
            // финальная проверка
            cluster_simple_roots = 0; for (auto c : cluster_counts) cluster_simple_roots += c;
            mult_degree = 0; mult_root_count = 0; for (auto& mc : local_mult_groups) { mult_degree += mc.first * mc.second; mult_root_count += mc.second; }
            used_degree = complex_degree + mult_degree + cluster_simple_roots;
            if (used_degree > P) throw std::invalid_argument("Невозможно уместить все запрошенные корни/кратности в степень P (после попытки сокращения)");
        }
    }

    // оставшиеся простые корни, которые нужно разместить вне кластеров
    unsigned remaining_simple = P - (complex_degree + mult_degree + cluster_simple_roots);

    // размещение центров кластеров равномерно и корректировка радиусов
    std::vector<T> cluster_centers;
    cluster_centers.reserve(num_clusters);
    if (num_clusters > 0) {
        // идеальные центры
        double step = 2.0 / static_cast<double>(num_clusters + 1);
        for (unsigned i = 0; i < num_clusters; ++i) {
            double c = -1.0 + (i + 1) * step;
            cluster_centers.push_back(static_cast<T>(c));
        }
        // корректируем радиусы чтобы не выходили за границы и не перекрывались
        // сначала гарантируем, что радиусы не выходят за границы
        for (unsigned i = 0; i < num_clusters; ++i) {
            T max_by_border = std::min(cluster_centers[i] - T(-1), T(1) - cluster_centers[i]);
            if (cluster_radii[i] > max_by_border) cluster_radii[i] = max_by_border;
            if (cluster_radii[i] < T(0)) cluster_radii[i] = T(0);
        }
        // затем обеспечим отсутствие перекрытия: radius_i <= half distance to neighbors
        for (unsigned i = 0; i < num_clusters; ++i) {
            T left_half = (i == 0) ? (cluster_centers[i] - T(-1)) : (cluster_centers[i] - cluster_centers[i - 1]);
            if (i > 0) left_half /= T(2);
            T right_half = (i + 1 == num_clusters) ? (T(1) - cluster_centers[i]) : (cluster_centers[i + 1] - cluster_centers[i]);
            if (i + 1 < num_clusters) right_half /= T(2);
            T allowed = std::min(left_half, right_half);
            if (cluster_radii[i] > allowed) cluster_radii[i] = allowed;
        }
    }

    // свободные интервалы между кластерами 
    struct Interval { T a, b; };
    std::vector<Interval> free_intervals;
    if (num_clusters == 0) {
        free_intervals.push_back({ T(-1), T(1) });
    }
    else {
        // левый промежуток
        T left = T(-1);
        for (unsigned i = 0; i < num_clusters; ++i) {
            T right = cluster_centers[i] - cluster_radii[i];
            if (right > left) free_intervals.push_back({ left, right });
            left = cluster_centers[i] + cluster_radii[i];
        }
        if (left < T(1)) free_intervals.push_back({ left, T(1) });
    }
    // если по каким-то причинам нет свободных интервалов (очень редкий случай), делаем один общий
    if (free_intervals.empty()) free_intervals.push_back({ T(-1), T(1) });

    // формируем список уникальных корней, которые нужно разместить во free_intervals 
    // сначала создаём вектор пар (multiplicity, placeholder)
    std::vector<unsigned> mults_for_free; // кратности для каждой уникальной корни, которые пойдут в free_intervals
    for (auto& mc : local_mult_groups) {
        unsigned mult = mc.first; unsigned cnt = mc.second;
        for (unsigned k = 0; k < cnt; ++k) mults_for_free.push_back(mult);
    }
    // добавляем remaining_simple корней с кратностью 1
    for (unsigned k = 0; k < remaining_simple; ++k) mults_for_free.push_back(1u);

    unsigned M = static_cast<unsigned>(mults_for_free.size());

    // длины интервалов
    std::vector<double> free_lengths;
    for (auto& I : free_intervals) {
        double L = static_cast<double>(I.b - I.a);
        if (L < 0) L = 0;
        free_lengths.push_back(L);
    }

    // как распределить M уникальных корней по интервалам
    std::vector<unsigned> counts_per_interval = allocate_counts_proportional(M, free_lengths);

    // помещаем кратные/простые корни в интервалы 
    unique_real_roots.clear();
    real_root_multiplicities.clear();

    unsigned idx_mults = 0;
    for (size_t i = 0; i < free_intervals.size(); ++i) {
        unsigned k = counts_per_interval[i];
        if (k == 0) continue;
        auto pts = evenly_spaced_in_interval<T>(free_intervals[i].a, free_intervals[i].b, k);
        for (unsigned j = 0; j < k; ++j) {
            T r = pts[j];
            unique_real_roots.push_back(r);
            real_root_multiplicities.push_back(mults_for_free[idx_mults++]);
        }
    }

    // если ещё остались неприсвоенные (из-за округлений), добавим их в середину соответствующих интервалов
    while (idx_mults < M) {
        // добавляем в центр первой интервальной области
        unique_real_roots.push_back((free_intervals[0].a + free_intervals[0].b) / T(2));
        real_root_multiplicities.push_back(mults_for_free[idx_mults++]);
    }

    // добавляем корни из кластеров (они все простые, кратность = 1) 
    for (unsigned i = 0; i < num_clusters; ++i) {
        unsigned cnt = (i < cluster_counts.size()) ? cluster_counts[i] : 1u;
        if (cnt == 0) continue;
        T a = cluster_centers[i] - cluster_radii[i];
        T b = cluster_centers[i] + cluster_radii[i];
        if (a > b) std::swap(a, b);
        auto pts = evenly_spaced_in_interval<T>(a, b, cnt);
        for (auto& r : pts) {
            unique_real_roots.push_back(r);
            real_root_multiplicities.push_back(1u);
        }
    }

    // строим расширенный список real_roots_repeated по кратностям 
    real_roots_repeated.clear();
    for (size_t i = 0; i < unique_real_roots.size(); ++i) {
        T r = unique_real_roots[i];
        unsigned mult = real_root_multiplicities[i];
        for (unsigned k = 0; k < mult; ++k) real_roots_repeated.push_back(r);
    }

    // генерируем комплексные корни
    complex_roots_out.clear();
    complex_roots_out.reserve(2 * num_complex_pairs);
    std::uniform_real_distribution<double> real_dist(-0.95, 0.95);
    std::uniform_real_distribution<double> imag_dist(0.01, 0.5);
    for (unsigned i = 0; i < num_complex_pairs; ++i) {
        double a = real_dist(gen);
        double b = imag_dist(gen);
        complex_roots_out.emplace_back(static_cast<T>(a), static_cast<T>(b));
        complex_roots_out.emplace_back(static_cast<T>(a), static_cast<T>(-b));
    }

    // строим многочлен 
    coefficients.clear();
    coefficients.push_back(T(1));

    for (auto& r : real_roots_repeated) multiply_by_root(coefficients, r);
    for (size_t i = 0; i + 1 < complex_roots_out.size(); i += 2) {
        T a = complex_roots_out[i].real();
        T b = complex_roots_out[i].imag();
        multiply_by_quadratic(coefficients, T(-2) * a, a * a + b * b);
    }

    // нормировка
    if (normalize_coeffs && !coefficients.empty()) {
        T leading = coefficients.back();
        if (!is_zero_val(leading)) {
            for (auto& c : coefficients) c = c / leading;
        }
    }
}

// метод для печати полинома
template<typename T>
void print_polynomial(const std::vector<T>& coefficients, std::ostream& os = std::cout) {
    int degree = static_cast<int>(coefficients.size()) - 1;
    bool first_term = true;

    for (int i = degree; i >= 0; --i) {  // идем от старшего к младшему
        T coeff = coefficients[i];

        if (is_zero_val(coeff)) continue;

        // знак
        if (!first_term) {
            if (is_positive_val(coeff))
                os << " + ";
            else {
                os << " - ";
                coeff = -coeff;
            }
        }
        else {
            if (!is_positive_val(coeff)) {
                os << "-";
                coeff = -coeff;
            }
            first_term = false;
        }

        // печать коэффициента (пропускаем 1 для x^i, кроме свободного члена)
        bool skip_one = (i > 0 && is_one_val(coeff));
        if (!skip_one) os << coeff;

        // печать степени
        if (i >= 1) {
            os << "x";
            if (i > 1) os << "^" << i;
        }
    }

    if (first_term) os << "0"; // если все коэффициенты нулевые

    os << std::endl;
}

#endif // GEN_HIGH_DEGREE_POLY_FIXED_H