// Вспомогательные функции для метода Дженкинса-Трауба (CPOLY)
// Комплексный Горнер, комплексная дефляция, вычисление P(z) и P'(z)
// Поддерживает произвольную точность через mathUtils (float_precision и др.)

#pragma once

#include <complex>
#include <cstddef>
#include <vector>

#include "NumericConcepts.h"
#include "mathUtils.h"

namespace cpoly_detail {

template <typename T>
std::complex<T> horner_eval(const std::vector<std::complex<T>> &coeffs,
                            const std::complex<T> &z) {
	int n = int(coeffs.size()) - 1;
	std::complex<T> result = coeffs[n];
	for (int k = n - 1; k >= 0; --k)
		result = result * z + coeffs[k];
	return result;
}

template <typename T>
std::pair<std::complex<T>, std::complex<T>>
horner_eval_and_deriv(const std::vector<std::complex<T>> &coeffs,
                      const std::complex<T> &z) {
	int n = int(coeffs.size()) - 1;
	std::complex<T> p = coeffs[n];
	std::complex<T> dp(T(0), T(0));

	for (int k = n - 1; k >= 0; --k) {
		dp = dp * z + p;
		p = p * z + coeffs[k];
	}
	return {p, dp};
}

template <typename T>
std::vector<std::complex<T>>
complex_deflate(const std::vector<std::complex<T>> &coeffs,
                const std::complex<T> &alpha) {
	int n = int(coeffs.size()) - 1;
	if (n < 1)
		return {};

	std::vector<std::complex<T>> result(n);
	result[n - 1] = coeffs[n];
	for (int k = n - 2; k >= 0; --k)
		result[k] = coeffs[k + 1] + alpha * result[k + 1];
	return result;
}

template <typename T>
std::vector<std::complex<T>>
to_complex(const std::vector<T> &real_coeffs) {
	std::vector<std::complex<T>> result;
	result.reserve(real_coeffs.size());
	for (const auto &c : real_coeffs)
		result.emplace_back(c, T(0));
	return result;
}

} // namespace cpoly_detail
