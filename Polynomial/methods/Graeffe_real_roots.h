// Метод Грефе для вычисления модулей вещественных корней
// find_moduli_roots_by_graeffe - реализация метода без отсеивания повторяющихся корней (порождает кучу мусора)
// find_moduli_with_multiplicities_by_graeffe - отсеивает дублирующиеся корни
// Теория и программная реализация: Павлова Анастасия, КМБО-01-22

#pragma once

#include "mathUtils.h"
#include "polynomialUtils.h"
#include "Helper_for_all_methods.h"

// Основная функция: вычисление модулей корней методом Грефе (логарифмическая модификация для повышенной точности)
template <typename T>
std::vector<T> find_moduli_roots_by_graeffe( const std::vector<T> &coeffs_desc, T epsilon = from_long_double<T>(1e-10L), int maxIter = 100)
{
    // Представление -inf в логарифмической шкале (используется для log(0))
    const long double NEG_INF = -std::numeric_limits<long double>::infinity();
    // epsilon_ld — порог сходимости в long double
    const long double epsilon_ld = std::max<long double>(std::fabsl(to_double_fallback(epsilon)), 1e-30L);

    // Безопасный логарифм (возвращает -inf если х=0)
    auto safe_log = [=](const T &x) -> long double {
        long double v = to_double_fallback(abs_val(x));
        if (v <= 0.0L)
            return NEG_INF;
        return std::logl(v);
    };

    // Вычисление log(sum(exp(vals))) — устойчиво при больших разностях
    auto log_sum_exp = [=](const std::vector<long double> &vals) -> long double {
        if (vals.empty()) return NEG_INF;
        long double m = *std::max_element(vals.begin(), vals.end());
        if (m == NEG_INF) return NEG_INF;
        long double s = 0.0L;
        for (auto v : vals)
            s += std::expl(v - m);
        return m + std::logl(s);
    };

    // Проверка на конечность значений
    auto is_finite_ld = [](long double x) { return std::isfinite(static_cast<double>(x)); };

    // Удаление ведущих нулей из коэффициентов
    size_t first_nz = 0;
    while (first_nz < coeffs_desc.size() && is_zero_val(coeffs_desc[first_nz])) ++first_nz;
    if (first_nz >= coeffs_desc.size()) { return {}; } // Все нули => нет корней

    // Берем только ненулевую часть (нули отсеили)
    std::vector<T> coeffs(coeffs_desc.begin() + first_nz, coeffs_desc.end());
    size_t deg = coeffs.size() - 1;
    if (deg == 0) { return {}; } // константный многочлен

    // Нормализация (ведущий коэффициент = 1)
    T lead = coeffs.front();
    for (auto &c : coeffs) c = c / lead;

    // Меняем порядок коэффициентов
    std::vector<T> p(coeffs.rbegin(), coeffs.rend());

    // Инициализация логарифмов и знаков коэффициентов
    std::vector<long double> L(p.size());
    std::vector<int> S(p.size());
    for (size_t i = 0; i < p.size(); ++i) {
        L[i] = safe_log(p[i]);
        S[i] = (to_double_fallback(p[i]) >= 0.0 ? 1 : -1);
    }

    int iters_done = 0;
    bool converged = false;
    
    // Основной цикл метода Грефе (итерации возведения в квадрат)
    for (int iter = 0; iter < maxIter; ++iter) {
        // Разделение на четные/нечетные
        std::vector<long double> Lpe, Lpo;
        std::vector<int> Spe, Spo;
        for (size_t i = 0; i < L.size(); ++i) {
            if ((i % 2) == 0) {
                Lpe.push_back(L[i]);
                Spe.push_back(S[i]);
            } else {
                Lpo.push_back(L[i]);
                Spo.push_back(S[i]);
            }
        }

        // Свертка: возведение в квадрат (в логарифмической форме)
        // Функция, вычисляющая log-коэффициенты результата свертки
        auto conv_square_log = [&](const std::vector<long double> &Lvec, const std::vector<int> &Svec) {
            size_t msize = (Lvec.empty()) ? 0 : (2 * Lvec.size() - 1);
            std::vector<long double> Lout(msize, NEG_INF);
            std::vector<int> Sout(msize, 0);
            for (size_t m = 0; m < msize; ++m) {
                // Собираем все слагаемые с их знаками
                std::vector<std::pair<long double, int>> terms;

                // свертка: b[m] = сумма (a[i]*a[m-i])
                for (size_t i = 0; i < Lvec.size(); ++i) {
                    if (m < i) break;
                    size_t j = m - i;
                    if (j >= Lvec.size()) continue;
                    long double lv = Lvec[i] + Lvec[j];
                    int sgn = Svec[i] * Svec[j];
                    terms.push_back({lv, sgn});
                }
                if (terms.empty()) continue;

                // Если все слагаемые одного знака, используем log_sum_exp
                bool all_same_sign = true;
                int first_sign = terms[0].second;
                for (const auto &t : terms) 
                    if (t.second != first_sign) { all_same_sign = false; break; }

                if (all_same_sign) {
                    std::vector<long double> vals;
                    for (const auto &t : terms) vals.push_back(t.first);
                    Lout[m] = log_sum_exp(vals);
                    Sout[m] = first_sign;
                } else {
                    // Смешанные знаки: переключаемся на обычную арифметику
                    long double sum = 0.0L;
                    for (const auto &t : terms) 
                        sum += t.second * std::expl(t.first);

                    if (std::fabsl(sum) < 1e-300L) {
                        Lout[m] = NEG_INF;
                        Sout[m] = 0;
                    } else {
                        Lout[m] = std::logl(std::fabsl(sum));
                        Sout[m] = (sum > 0.0L) ? 1 : -1;
                    }
                }
            }
            return std::make_pair(Lout, Sout);
        };

        // Получаем квадртаы для четной и нечетной частей
        auto [Lpe2, Spe2] = conv_square_log(Lpe, Spe);
        auto [Lpo2, Spo2] = conv_square_log(Lpo, Spo);

        // Формируем новый многочлен: Q(y) = Pe(y)^2 - y * Po2(y)
        std::vector<long double> Lq(L.size(), NEG_INF);
        std::vector<int> Sq(L.size(), 0);

        for (size_t m = 0; m < L.size(); ++m) {
            bool hasPe = (m < Lpe2.size() && Lpe2[m] != NEG_INF);
            bool hasPo = (m >= 1 && (m - 1) < Lpo2.size() && Lpo2[m - 1] != NEG_INF);
            if (!hasPe && !hasPo) continue;

            if (hasPe && !hasPo) {
                Lq[m] = Lpe2[m];
                Sq[m] = Spe2[m];
            } else if (!hasPe && hasPo) {
                Lq[m] = Lpo2[m - 1];
                Sq[m] = -Spo2[m - 1];
            } else {
                // Оба есть: нужно аккуратно сложить с учётом возможной потери точности
                long double A = Lpe2[m];
                long double B = Lpo2[m - 1];
                int sA = Spe2[m];
                int sB = Spo2[m - 1];
                if (A == NEG_INF && B == NEG_INF) continue;

                // Если A =(примерно) B и знаки противоположные
                long double log_diff = std::fabsl(A - B);
                if (sA != sB && log_diff < 2.0L) // Близкие значения с разными знаками
                {
                    // Почти одинаковые величины с противоположными знаками — считаем в обычных числах
                    long double val_A = sA * std::expl(A);
                    long double val_B = -sB * std::expl(B);
                    long double comb = val_A + val_B;

                    if (std::fabsl(comb) < 1e-300L) {
                        Lq[m] = NEG_INF;
                        Sq[m] = 0;
                    }
                    else {
                        Lq[m] = std::logl(std::fabsl(comb));
                        Sq[m] = (comb > 0.0L) ? 1 : -1;
                    }
                }
                else {
                    // Обычное логарифмическое сложение
                    long double M = std::max(A, B);
                    long double u = sA * std::expl(A - M);
                    long double v = -sB * std::expl(B - M);
                    long double comb = u + v;
                    if (std::fabsl(comb) < 1e-300L) {
                        Lq[m] = NEG_INF;
                        Sq[m] = 0;
                    }
                    else {
                        Lq[m] = M + std::logl(std::fabsl(comb));
                        Sq[m] = (comb > 0.0L) ? 1 : -1;
                    }
                }
            }
        }

        // Нормализация по старшему
        long double lead_log = Lq.back();
        if (!is_finite_ld(lead_log)) { break; }
        for (auto &v : Lq)
            v -= lead_log * (is_finite_ld(v));

        // Нормализация старого L по тому же принципу (чтобы delta был корректным)
        long double lead_log_old = L.back();
        if (is_finite_ld(lead_log_old)) {
            for (auto &v : L) v -= lead_log_old * (is_finite_ld(v));
        }

        // Проверка сходимости: сравниваем log-разности коэффициентов
        // Эти разности дают логарифмы модулей корней, они должны стабилизироваться
        long double max_rel = 0.0L;
        bool local_nonfinite = false;

        // Вычисляем log-разности для текущего и предыдущего полиномов
        std::vector<long double> log_ratios_old, log_ratios_new;
        for (size_t i = 1; i < L.size(); ++i) {
            if (L[i] != NEG_INF && L[i - 1] != NEG_INF && is_finite_ld(L[i]) && is_finite_ld(L[i - 1]))
                log_ratios_old.push_back(L[i] - L[i - 1]);
            if (Lq[i] != NEG_INF && Lq[i - 1] != NEG_INF && is_finite_ld(Lq[i]) && is_finite_ld(Lq[i - 1]))
                log_ratios_new.push_back(Lq[i] - Lq[i - 1]);
        }

        // Сравниваем изменение log-разностей (они должны уменьшаться в 2^k раз каждую итерацию)
        if (log_ratios_old.size() == log_ratios_new.size() && !log_ratios_old.empty()) {
            for (size_t i = 0; i < log_ratios_old.size(); ++i) {
                // После k итераций: log|x_j| =(примерно) log_ratio / 2^k
                // Изменение: |log_ratio_new / 2^(k+1) - log_ratio_old / 2^k| = |log_ratio_new/2 - log_ratio_old| / 2^k
                // Упрощённо: смотрим на относительное изменение log_ratio
                long double ratio_old = log_ratios_old[i];
                long double ratio_new = log_ratios_new[i];

                if (!is_finite_ld(ratio_old) || !is_finite_ld(ratio_new)) {
                    local_nonfinite = true;
                    max_rel = std::numeric_limits<long double>::infinity();
                    break;
                }

                // Относительное изменение log-разности (это и есть критерий сходимости)
                long double delta = std::fabsl(ratio_new - 2.0L * ratio_old);
                long double scale = std::max(std::fabsl(ratio_new), std::fabsl(ratio_old));
                
                long double rel_change = delta / scale * (scale > 0.0L);
                max_rel = std::max(max_rel, rel_change) * (scale > 0.0L);
            }
        } else {
            // Сравниваем просто изменение коэффициентов
            for (size_t i = 0; i < std::min(L.size(), Lq.size()); ++i)
            {
                if (L[i] == NEG_INF || Lq[i] == NEG_INF) continue;
                long double delta = Lq[i] - L[i];

                if (!is_finite_ld(delta)) {
                    local_nonfinite = true;
                    max_rel = std::numeric_limits<long double>::infinity();
                    break;
                }

                long double clamped = std::clamp(delta, -7000.0L, 7000.0L);
                long double diff = std::fabsl(std::expm1l(clamped));

                if (!is_finite_ld(diff)) {
                    local_nonfinite = true;
                    max_rel = std::numeric_limits<long double>::infinity();
                    break;
                }
                if (diff > max_rel)
                    max_rel = diff;
            }
        }

        // Обновляем L и S
        L = std::move(Lq);
        S = std::move(Sq);
        iters_done = iter + 1;

        // Проверка критериев остановки
        if (local_nonfinite || max_rel > 1e100L) {  return {}; }
        if (max_rel < epsilon_ld) {
            converged = true;
            break;
        }
    }

    // Обратно переворачиваем
    std::vector<long double> Ldesc(L.rbegin(), L.rend());
    std::vector<int> Sdesc(S.rbegin(), S.rend());
    size_t n = Ldesc.size() - 1;

    // Проверка: сошлись ли мы?
    if (!converged) { return {}; }
    if (n == 0) { return {}; }

    // Простое извлечение модулей по формуле |x_j| = (|A_j| / |A_{j-1}|)^(1/2^k)
    const int effective_iters = std::max(iters_done, 1);
    long double pow_two_k = std::ldexp(1.0L, std::min(effective_iters, 60));
    if (effective_iters > 60)
        pow_two_k = pow_val(2.0L, static_cast<long double>(effective_iters));

    // Если в середине -inf, значит все корни одного модуля
    bool has_zero_middle = false;
    for (size_t i = 1; i + 1 < Ldesc.size(); ++i) {
        if (!is_finite_ld(Ldesc[i]) || Ldesc[i] == NEG_INF) {
            has_zero_middle = true;
            break;
        }
    }

    // Если есть нулевые средние коэффициенты, все корни имеют одинаковый модуль
    if (has_zero_middle) {
        // Извлекаем модуль из первого и последнего ненулевых коэффициентов
        long double L_first = Ldesc[0];
        long double L_last = Ldesc.back();

        if (is_finite_ld(L_first) && is_finite_ld(L_last)) {
            // |x|^n = |A_0| / |A_n| => |x| = (|A_0| / |A_n|)^(1/n)
            long double log_ratio = (L_first - L_last) / static_cast<long double>(n);
            long double common_modulus = std::expl(log_ratio / pow_two_k);

            if (!is_finite_ld(common_modulus) || common_modulus < 0.0L)
                common_modulus = 1.0L; // для единичного круга

            std::vector<T> moduli(n, from_long_double<T>(common_modulus));
            return moduli;
        }
    }

    // Извлечение по соседним коэффициентам
    std::vector<T> moduli;
    moduli.reserve(n);

    for (size_t j = 1; j <= n; ++j) {
        long double Lprev = Ldesc[j - 1];
        long double Lcurr = Ldesc[j];

        // Если оба нулевые, модуль = 0
        if (!is_finite_ld(Lprev) && !is_finite_ld(Lcurr)) {
            moduli.push_back(T(0));
            continue;
        }

        // Подстановка для -inf
        if (!is_finite_ld(Lprev)) Lprev = -1e300L;
        if (!is_finite_ld(Lcurr)) Lcurr = -1e300L;

        // log|x_j| = (log|A_j| - log|A_{j-1}|) / 2^k
        long double log_delta = Lcurr - Lprev;
        long double log_modulus = log_delta / pow_two_k;

        if (!is_finite_ld(log_modulus)) {
            moduli.push_back(T(0));
            continue;
        }

        // Доп. проверка
        long double modulus = std::expl(log_modulus);
        if (!is_finite_ld(modulus) || modulus < 0.0L)
            modulus = 0.0L * (!is_finite_ld(modulus) || modulus < 0.0L);

        moduli.push_back(from_long_double<T>(modulus));
    }
    return moduli;
}

// Структура для хранения модуля с его кратностью
template <typename T>
struct ModulusWithMultiplicity {
    T value;          // значение модуля
    int multiplicity; // сколько раз этот модуль встретился

    ModulusWithMultiplicity(T val, int mult = 1) : value(val), multiplicity(mult) {}
};

// Новая версия функции, возвращающая модули с кратностями
template <typename T>
std::vector<ModulusWithMultiplicity<T>> find_moduli_with_multiplicities_by_graeffe(
    const std::vector<T>& coeffs_desc, T epsilon = from_long_double<T>(1e-10L),  int maxIter = 100 )
{
    // Сначала находим модули стандартным методом
    std::vector<T> raw_moduli = find_moduli_roots_by_graeffe(coeffs_desc, epsilon, maxIter);

    if (raw_moduli.empty()) { return {}; }

    // Группируем модули с учётом погрешности
    std::vector<ModulusWithMultiplicity<T>> grouped;
    T tolerance = epsilon * from_long_double<T>(100.0L); // более мягкий порог для группировки

    for (const auto& modulus : raw_moduli) {
        bool found = false;
        // Проверяем, есть ли уже такой модуль в группированном списке
        for (auto& group : grouped) {
            T diff = group.value > modulus ? group.value - modulus : modulus - group.value;

            if (diff < tolerance) {
                // Найден похожий модуль - увеличиваем кратность
                group.multiplicity++;

                // Обновляем значение модуля (среднее арифметическое для стабильности)
                T old_value = group.value;
                group.value = (old_value * from_long_double<T>(group.multiplicity - 1) + modulus) / from_long_double<T>(group.multiplicity);

                found = true;
                break;
            }
        }
        if (!found)
            // Новый уникальный модуль
            grouped.push_back(ModulusWithMultiplicity<T>(modulus, 1));
    }
    return grouped;
}