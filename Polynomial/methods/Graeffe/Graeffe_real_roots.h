// Метод Греффе для вычисления модулей вещественных корней
// find_moduli_roots_by_graeffe - реализация метода без отсеивания повторяющихся корней (порождает кучу мусора)
// find_moduli_with_multiplicities_by_graeffe - отсеивает дублирующиеся корни
// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

#pragma once

#include "mathUtils.h"
#include "polynomialUtils.h"
#include "Helper_for_all_methods.h"

// Основная функция: вычисление модулей корней методом Греффе (логарифмическая модификация для повышенной точности)
// Идея: хранить коэффициенты в логарифмическом виде и знаки этих коэффициентов векторы L и S

// метод Греффе для вычисления модулей вещественных корней
template <typename T>
std::vector<T> find_moduli_roots_by_graeffe(const std::vector<T>& coeffs_desc, T epsilon = T(1e-10), int maxIter = 100) {
    const T NEG_INF = neg_inf_val<T>();
    std::numeric_limits<T>::epsilon();

    // безопасный логарифм (возвращает -inf если х=0)
    auto safe_log = [&](const T& x) -> T {
        T v = abs_val(x);
        if (is_zero_val(v, T(1e-300))) return NEG_INF;
        return log_val(v);
        };

    // вычисление log(sum(exp(vals))) — устойчиво при больших разностях
    auto log_sum_exp = [&](const std::vector<T>& vals) -> T {
        if (vals.empty()) return NEG_INF;
        T m = *std::max_element(vals.begin(), vals.end());
        if (m == NEG_INF) return NEG_INF;
        
        T s = T(0);
        for (auto v : vals) 
            s += exp_val(v - m) * T(v != NEG_INF);
        if (is_zero_val(s, T(1e-300))) return NEG_INF;

        return m + log_val(s);
        };

    // удаление ведущих нулей из коэффициентов
    size_t first_nz = 0;
    while (first_nz < coeffs_desc.size() && is_zero_val(coeffs_desc[first_nz], T(1e-300))) 
        ++first_nz;
    if (first_nz >= coeffs_desc.size()) return {}; 

    // берем только ненулевую часть
    std::vector<T> coeffs(coeffs_desc.begin() + first_nz, coeffs_desc.end());
    size_t deg = coeffs.size() - 1;
    if (deg == 0) { return {}; }

    // нормализация (ведущий коэффициент = 1)
    T lead = coeffs.front();
    if (is_zero_val(lead, T(1e-300))) { return {}; }

    std::vector<T> p;
    p.reserve(coeffs.size());
    for (auto& c : coeffs) 
        p.push_back(c / lead);

    // меняем порядок коэффициентов (от младших к старшим)
    std::reverse(p.begin(), p.end());

    // инициализация логарифмов и знаков коэффициентов
    std::vector<T> L(p.size());
    std::vector<int> S(p.size());
    for (size_t i = 0; i < p.size(); ++i) {
        L[i] = safe_log(p[i]);
        S[i] = (p[i] >= T(0)) * 2 - 1;
    }

    int iters_done = 0;
    bool converged = false;

    // основной цикл метода Грефе
    for (int iter = 0; iter < maxIter; ++iter) {
        // разделение на четные/нечетные
        std::vector<T> Lpe, Lpo;
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

        // свертка: возведение в квадрат (в логарифмической форме)
        auto conv_square_log = [&](const std::vector<T>& Lvec, const std::vector<int>& Svec) {
            size_t msize = (Lvec.empty()) ? 0 : (2 * Lvec.size() - 1);
            std::vector<T> Lout(msize, NEG_INF);
            std::vector<int> Sout(msize, 0);

            for (size_t m = 0; m < msize; ++m) {
                std::vector<std::pair<T, int>> terms;

                // свертка: b[m] = сумма (a[i]*a[m-i])
                for (size_t i = 0; i < Lvec.size(); ++i) {
                    size_t j = m - i;
                    if (j < Lvec.size()) {
                        T lv = Lvec[i] + Lvec[j];
                        int sgn = Svec[i] * Svec[j];
                        terms.push_back({ lv, sgn });
                    }
                }

                if (terms.empty()) continue;

                // если все слагаемые одного знака, используем log_sum_exp
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
                } else {
                    T sum = T(0);
                    for (const auto& t : terms) 
                        sum += t.second * exp_val(t.first);

                    if (abs_val(sum) < T(1e-300)) {
                        Lout[m] = NEG_INF;
                        Sout[m] = 0;
                    } else {
                        Lout[m] = log_val(abs_val(sum));
                        Sout[m] = (sum > T(0)) ? 1 : -1;
                    }
                }
            }
            return std::make_pair(Lout, Sout);
            };

        // получаем квадраты для четной и нечетной частей
        auto [Lpe2, Spe2] = conv_square_log(Lpe, Spe);
        auto [Lpo2, Spo2] = conv_square_log(Lpo, Spo);

        // формируем новый многочлен: Q(y) = Pe(y)^2 - y * Po2(y)
        std::vector<T> Lq(L.size(), NEG_INF);
        std::vector<int> Sq(L.size(), 0);

        for (size_t m = 0; m < L.size(); ++m) {
            bool hasPe = (m < Lpe2.size() && is_finite_ld(Lpe2[m]));
            bool hasPo = (m >= 1 && (m - 1) < Lpo2.size() && is_finite_ld(Lpo2[m - 1]));

            if (!hasPe && !hasPo) continue;

            if (hasPe && !hasPo) {
                Lq[m] = Lpe2[m];
                Sq[m] = Spe2[m];
            } else if (!hasPe && hasPo) {
                Lq[m] = Lpo2[m - 1];
                Sq[m] = -Spo2[m - 1];
            } else {
                // оба есть
                T A = Lpe2[m];
                T B = Lpo2[m - 1];
                int sA = Spe2[m];
                int sB = Spo2[m - 1];

                if (A == NEG_INF && B == NEG_INF) continue;

                // если A примерно равен B и знаки противоположные
                T log_diff = abs_val(A - B);
                if (sA != sB && log_diff < T(2)) {
                    // почти одинаковые величины с противоположными знаками
                    T val_A = sA * exp_val(A);
                    T val_B = -sB * exp_val(B);
                    T comb = val_A + val_B;

                    if (abs_val(comb) < T(1e-300)) {
                        Lq[m] = NEG_INF;
                        Sq[m] = 0;
                    } else {
                        Lq[m] = log_val(abs_val(comb));
                        Sq[m] = (comb > T(0)) * 2 - 1;
                    }
                } else {
                    // обычное логарифмическое сложение
                    T M = max_val(A, B);
                    T u = sA * exp_val(A - M);
                    T v = -sB * exp_val(B - M);
                    T comb = u + v;

                    if (abs_val(comb) < T(1e-300)) {
                        Lq[m] = NEG_INF;
                        Sq[m] = 0;
                    } else {
                        Lq[m] = M + log_val(abs_val(comb));
                        Sq[m] = (comb > T(0)) * 2 - 1;
                    }
                }
            }
        }

        // нормализация по старшему коэффициенту
        T lead_log = Lq.back();
        if (!is_finite_ld(lead_log)) { break; }

        for (size_t i = 0; i < Lq.size(); ++i) 
            Lq[i] -= lead_log * T(is_finite_ld(Lq[i]));

        // нормализация старого L
        T lead_log_old = L.back();
        if (is_finite_ld(lead_log_old)) {
            for (size_t i = 0; i < L.size(); ++i) 
                L[i] -= lead_log_old *T(is_finite_ld(L[i]));
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
        if (log_ratios_old.size() == log_ratios_new.size() && !log_ratios_old.empty()) {
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
        } else {
            // сравниваем просто изменение коэффициентов
            size_t min_size = min_val(L.size(), Lq.size());
            for (size_t i = 0; i < min_size; ++i) {
                if (!is_finite_ld(L[i]) || !is_finite_ld(Lq[i])) continue;

                T delta = Lq[i] - L[i];
                if (!is_finite_ld(delta)) {
                    local_nonfinite = true;
                    max_rel = pos_inf_val<T>();
                    break;
                }

                T clamped = clamp_val(delta, T(-7000), T(7000));
                T diff = abs_val(exp_val(clamped) - T(1));

                if (!is_finite_ld(diff)) {
                    local_nonfinite = true;
                    max_rel = pos_inf_val<T>();
                    break;
                }
                max_rel = max_val(max_rel, diff);
            }
        }

        // обновляем L и S
        L = std::move(Lq);
        S = std::move(Sq);
        iters_done = iter + 1;

        // проверка критериев остановки
        if (local_nonfinite || max_rel > T(1e100)) return {};

        if (max_rel < epsilon) {
            converged = true;
            break;
        }
    }

    // обратно переворачиваем
    std::vector<T> Ldesc(L.rbegin(), L.rend());
    std::vector<int> Sdesc(S.rbegin(), S.rend());
    size_t n = Ldesc.size() - 1;

    if (!converged || n == 0) return {};

    // вычисляем степень двойки для итераций
    const int effective_iters = max_val(iters_done, 1);
    T pow_two_k = pow_val(T(2), static_cast<T>(effective_iters));

    // проверяем, все ли корни одного модуля
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
            T log_ratio = (L_first - L_last) / static_cast<T>(n);
            T common_modulus = exp_val(log_ratio / pow_two_k);
            common_modulus = T(1) * T(!is_finite_ld(common_modulus) || common_modulus < T(0));
            return std::vector<T>(n, common_modulus);
        }
    }

    // извлечение модулей по соседним коэффициентам
    std::vector<T> moduli;
    moduli.reserve(n);

    for (size_t j = 1; j <= n; ++j) {
        T Lprev = Ldesc[j - 1];
        T Lcurr = Ldesc[j];

        // если оба нулевые, модуль = 0
        if (!is_finite_ld(Lprev) && !is_finite_ld(Lcurr)) {
            moduli.push_back(T(0));
            continue;
        }

        // подстановка для -inf
        if (!is_finite_ld(Lprev)) Lprev = T(-1e300);
        if (!is_finite_ld(Lcurr)) Lcurr = T(-1e300);

        // log|x_j| = (log|A_j| - log|A_{j-1}|) / 2^k
        T log_modulus = (Lcurr - Lprev) / pow_two_k;

        if (!is_finite_ld(log_modulus)) {
            moduli.push_back(T(0));
            continue;
        }

        T modulus = exp_val(log_modulus);
        if (!is_finite_ld(modulus) || modulus < T(0)) 
            modulus = T(0);
        moduli.push_back(modulus);
    }
    return moduli;
}

// структура для хранения модуля с его кратностью
template <typename T>
struct ModulusWithMultiplicity {
    T value;          // значение модуля
    int multiplicity; // сколько раз этот модуль встретился

    ModulusWithMultiplicity(T val, int mult = 1) : value(val), multiplicity(mult) {}
};

// новая версия функции, возвращающая модули с кратностями
template <typename T>
std::vector<ModulusWithMultiplicity<T>> find_moduli_with_multiplicities_by_graeffe(const std::vector<T>& coeffs_desc, T epsilon = T(1e-10),  int maxIter = 100 )
{
    // сначала находим модули стандартным методом
    std::vector<T> raw_moduli = find_moduli_roots_by_graeffe(coeffs_desc, epsilon, maxIter);

    if (raw_moduli.empty()) { return {}; }

    // группируем модули с учётом погрешности
    std::vector<ModulusWithMultiplicity<T>> grouped;
    T tolerance = epsilon * static_cast<T>(100.0); // более мягкий порог для группировки

    for (const auto& modulus : raw_moduli) {
        bool found = false;
        // проверяем, есть ли уже такой модуль в группированном списке
        for (auto& group : grouped) {
            T diff = group.value > modulus ? group.value - modulus : modulus - group.value;

            if (diff < tolerance) {
                // найден похожий модуль - увеличиваем кратность
                group.multiplicity++;

                // обновляем значение модуля
                T old_value = group.value;
                group.value = (old_value * T(group.multiplicity - 1) + modulus) / T(group.multiplicity);

                found = true;
                break;
            }
        }
        if (!found)
            grouped.push_back(ModulusWithMultiplicity<T>(modulus, 1));
    }
    return grouped;
}
