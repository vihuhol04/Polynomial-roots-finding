#pragma once

#include <vector>
#include <algorithm>

#include "Graeffe_real_roots.h"
#include "NumericConcepts.h"
#include "mathUtils.h"

/**
 * @brief Модификация Юркша (ПРАВИЛЬНАЯ)
 *
 * Не масштабирует полином!
 * Использует повторный Graeffe + сглаживание отношений.
 */
template <RealNumber T>
std::vector<T> find_moduli_yurksch_wrapper(
    const std::vector<T>& coeffs_desc,
    int refinement_steps = 3,
    T epsilon = numeric_constants::adaptive_epsilon<T>(
        numeric_constants::EPSILON_SCALE_PRECISE),
    int maxIter = numeric_constants::DEFAULT_MAX_ITERATIONS,
    bool debug = false)
{
    if (coeffs_desc.size() < 2)
        return {};

    // 1. базовый Graeffe
    std::vector<T> moduli =
        find_moduli_roots_by_graeffe(coeffs_desc, epsilon, maxIter, debug);

    if (moduli.empty())
        return {};

    // 2. Юркш = итеративное сглаживание
    for (int step = 0; step < refinement_steps; ++step)
    {
        std::vector<T> refined = moduli;

        for (size_t i = 1; i < moduli.size(); ++i)
        {
            T prev = moduli[i - 1];
            T curr = moduli[i];

            if (prev > T(0) && curr > T(0))
            {
                // геометрическое сглаживание (аналог отношений)
                refined[i] = sqrt_val(prev * curr);
            }
        }

        // стабилизация (сортировка — важно для Греффе!)
        std::sort(refined.begin(), refined.end());

        moduli = refined;

        if (debug)
        {
            std::cerr << "[YURKSCH] step=" << step
                << " first=" << (double)moduli[0]
                << " last=" << (double)moduli.back()
                << std::endl;
        }
    }

    return moduli;
}