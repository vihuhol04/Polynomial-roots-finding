// Аналог обычному polyprecision, но поддерживающий полиномы от нескольких переменных
// Реализация: Трудолюбов Никита, КМБО-01-22

#ifndef INC_MPOLYNOMIAL_H
#define INC_MPOLYNOMIAL_H

// Все необходимые включения стандартной библиотеки
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <sstream>
#include <numeric>
#include <cmath>
#include <stdexcept> // Для std::runtime_error

// Шаблонный класс для многочлена от нескольких переменных
template<class _TY>
class mpolynomial {
private:
    std::map<std::vector<unsigned int>, _TY> terms;

    void normalize() {
        for (auto it = terms.begin(); it != terms.end(); ) {
            if (it->second == _TY(0)) {
                it = terms.erase(it);
            }
            else {
                ++it;
            }
        }
    }
public:
    typedef _TY value_type;

    // Конструкторы
    mpolynomial() = default;
    mpolynomial(const _TY& constant) {
        if (constant != _TY(0)) {
            terms[{}] = constant;
        }
    }

    static mpolynomial variable(unsigned int var_index, unsigned int power = 1) {
        mpolynomial p;
        if (power > 0) {
            std::vector<unsigned int> exponents(var_index + 1, 0);
            exponents[var_index] = power;
            p.terms[exponents] = _TY(1);
        }
        return p;
    }

    // Методы доступа
    const std::map<std::vector<unsigned int>, _TY>& get_terms() const { return terms; }
    std::map<std::vector<unsigned int>, _TY>& get_terms() { return terms; }

    void add_term(const std::vector<unsigned int>& exponents, const _TY& coeff) {
        std::vector<unsigned int> exps = exponents;
        while (!exps.empty() && exps.back() == 0) {
            exps.pop_back();
        }

        terms[exps] += coeff;
        normalize();
    }

    bool empty() const { return terms.empty(); }

    size_t degree(size_t var_index) const {
        size_t max_deg = 0;
        if (empty()) return 0;
        for (const auto& term : terms) {
            if (var_index < term.first.size()) {
                if (term.first[var_index] > max_deg) {
                    max_deg = term.first[var_index];
                }
            }
        }
        return max_deg;
    }

    mpolynomial derivative(size_t var_index) const {
        mpolynomial result;
        for (const auto& term : terms) {
            if (var_index < term.first.size() && term.first[var_index] > 0) {
                _TY new_coeff = term.second * static_cast<_TY>(term.first[var_index]);
                std::vector<unsigned int> new_exponents = term.first;
                new_exponents[var_index]--;
                result.add_term(new_exponents, new_coeff);
            }
        }
        return result;
    }

    mpolynomial evaluate(const std::map<size_t, float_precision>& point) const {
        mpolynomial result;
        for (const auto& term : terms) {
            _TY new_coeff = term.second;
            std::vector<unsigned int> new_exponents = term.first;

            for (const auto& val_pair : point) {
                size_t var_idx = val_pair.first;
                const float_precision& value = val_pair.second;
                if (var_idx < new_exponents.size() && new_exponents[var_idx] > 0) {
                    for (unsigned int i = 0; i < new_exponents[var_idx]; ++i) {
                        new_coeff *= value;
                    }
                    new_exponents[var_idx] = 0;
                }
            }
            result.add_term(new_exponents, new_coeff);
        }
        return result;
    }

    bool is_zero() const {
        return terms.empty();
    }

    bool operator==(const mpolynomial& other) const { return this->terms == other.terms; }

    mpolynomial& operator+=(const mpolynomial& other) {
        for (const auto& term : other.terms) { this->add_term(term.first, term.second); }
        return *this;
    }

    mpolynomial& operator-=(const mpolynomial& other) {
        for (const auto& term : other.terms) { this->add_term(term.first, -term.second); }
        return *this;
    }

    mpolynomial& operator*=(const mpolynomial& other) {
        if (this->empty() || other.empty()) { terms.clear(); return *this; }
        mpolynomial result;
        for (const auto& term1 : terms) {
            for (const auto& term2 : other.terms) {
                _TY new_coeff = term1.second * term2.second;
                std::vector<unsigned int> new_exponents = term1.first;
                if (new_exponents.size() < term2.first.size()) { new_exponents.resize(term2.first.size(), 0); }
                for (size_t i = 0; i < term2.first.size(); ++i) { new_exponents[i] += term2.first[i]; }
                result.add_term(new_exponents, new_coeff);
            }
        }
        *this = result;
        return *this;
    }

    // Преобразование в строку
    std::string toString() const {
        if (empty()) return "0";
        std::stringstream ss;
        bool first_term = true;
        std::vector<std::pair<std::vector<unsigned int>, _TY>> sorted_terms(terms.begin(), terms.end());
        std::sort(sorted_terms.begin(), sorted_terms.end(), [](const auto& a, const auto& b) {
            return a.first > b.first;
            });

        for (const auto& term : sorted_terms) {
            const auto& exponents = term.first;
            const _TY& coeff = term.second;
            if (first_term) { if (coeff < _TY(0)) ss << "-"; }
            else { ss << (coeff < _TY(0) ? " - " : " + "); }
            _TY abs_coeff = (coeff < _TY(0)) ? -coeff : coeff;
            bool has_variables = false;
            for (unsigned int exp : exponents) { if (exp > 0) { has_variables = true; break; } }
            if (abs_coeff != _TY(1) || !has_variables) { ss << abs_coeff; }
            for (size_t i = 0; i < exponents.size(); ++i) {
                if (exponents[i] > 0) {
                    // Вот исправленная строка
                    if (abs_coeff == _TY(1) && has_variables && i > 0 && !ss.str().empty() && ss.str().back() != ' ') ss << "*";
                    ss << "x_" << i;
                    if (exponents[i] > 1) ss << "^" << exponents[i];
                }
            }
            first_term = false;
        }
        return ss.str();
    }
};

// Внешние операторы
template<class _TY>
mpolynomial<_TY> operator+(mpolynomial<_TY> lhs, const mpolynomial<_TY>& rhs) {
    lhs += rhs;
    return lhs;
}

template<class _TY>
mpolynomial<_TY> operator-(mpolynomial<_TY> lhs, const mpolynomial<_TY>& rhs) {
    lhs -= rhs;
    return lhs;
}

template<class _TY>
mpolynomial<_TY> operator*(mpolynomial<_TY> lhs, const mpolynomial<_TY>& rhs) {
    lhs *= rhs;
    return lhs;
}

template<class _TY>
std::ostream& operator<<(std::ostream& os, const mpolynomial<_TY>& p) {
    os << p.toString();
    return os;
}

template<class _TY>
mpolynomial<_TY> pow(const mpolynomial<_TY> base, int exp) {
    if (exp < 0) return mpolynomial<_TY>();
    if (exp == 0) return mpolynomial<_TY>(static_cast<_TY>(1));

    mpolynomial<_TY> result(static_cast<_TY>(1));
    mpolynomial<_TY> p = base;
    while (exp > 0) {
        if (exp % 2 == 1) result *= p;
        p *= p;
        exp /= 2;
    }
    return result;
}

#endif // INC_MPOLYNOMIAL_H