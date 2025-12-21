// Интервальная арифметика для 01-теста в методе ANewDsc
// Гарантированное содержание истинного значения

#pragma once

#include "mathUtils.h"

// Шаблон класса интервала [lower, upper]
// Хранит диапазон возможных значений числа
template <typename T>
class Interval
{
private:
    T lower; // нижняя граница интервала
    T upper; // верхняя граница интервала

public:
    // конструкторы
    Interval() : lower(0), upper(0) {} // по умолчанию интервал [0, 0]

    explicit Interval(T value) : lower(value), upper(value) {} // точечный интервал

    Interval(T l, T u) : lower(l), upper(u)
    {
        // гарантируем, что нижняя граница <= верхней
        if (l > u) {
            throw std::invalid_argument("Lower bound must be <= upper bound");
        }
    }

    // геттеры
    T getLower() const { return lower; } // начало
    T getUpper() const { return upper; } // конец
    T width() const { return upper - lower; }          // ширина интервала
    T midpoint() const { return (lower + upper) / T(2); } // центральная точка

    // проверки принадлежности
    bool contains(T value) const { return value >= lower && value <= upper; }

    // проверка, содержит ли интервал ноль
    bool containsZero() const { return lower <= T(0) && upper >= T(0); }

    // полностью положительный интервал
    bool isPositive() const { return lower > T(0); }

    // полностью отрицательный интервал
    bool isNegative() const { return upper < T(0); }

    // Арифметические операции
    // Все операции выполняются так, чтобы результат гарантированно содержал истинное значение результата операции над реальными числами.

    // сложение: [a,b] + [c,d] = [a+c, b+d]
    Interval operator+(const Interval& other) const { return Interval(lower + other.lower, upper + other.upper); }

    // вычитание: [a,b] - [c,d] = [a-d, b-c]
    Interval operator-(const Interval& other) const { return Interval(lower - other.upper, upper - other.lower); }

    // унарный минус: -[a,b] = [-b, -a]
    Interval operator-() const { return Interval(-upper, -lower); }

    // умножение: перебираем все комбинации концов интервалов
    Interval operator*(const Interval& other) const {
        T products[4] = {
            lower * other.lower,
            lower * other.upper,
            upper * other.lower,
            upper * other.upper
        };

        T min_prod = products[0];
        T max_prod = products[0];
        for (int i = 1; i < 4; ++i) {
            min_prod = min_val(min_prod, products[i]);
            max_prod = max_val(max_prod, products[i]);
        }

        return Interval(min_prod, max_prod);
    }

    // деление интервалов
    // деление запрещено, если делитель содержит 0
    Interval operator/(const Interval& other) const {
        if (other.containsZero()) {
            throw std::domain_error("Division by interval containing zero");
        }

        T quotients[4] = {
            lower / other.lower,
            lower / other.upper,
            upper / other.lower,
            upper / other.upper
        };

        T min_quot = quotients[0];
        T max_quot = quotients[0];
        for (int i = 1; i < 4; ++i) {
            min_quot = min_val(min_quot, quotients[i]);
            max_quot = max_val(max_quot, quotients[i]);
        }

        return Interval(min_quot, max_quot);
    }

    // операции с числовыми константами
    Interval operator+(T scalar) const { return Interval(lower + scalar, upper + scalar); }

    Interval operator-(T scalar) const { return Interval(lower - scalar, upper - scalar); }

    Interval operator*(T scalar) const {
        if (scalar >= T(0)) { return Interval(lower * scalar, upper * scalar); }
        else {
            // Если коэффициент отрицателен — границы меняются местами
            return Interval(upper * scalar, lower * scalar);
        }
    }

    Interval operator/(T scalar) const {
        if (scalar == T(0)) { throw std::domain_error("Division by zero"); }

        if (scalar > T(0)) { return Interval(lower / scalar, upper / scalar); }
        else { return Interval(upper / scalar, lower / scalar); }
    }

    // операции присваивания
    Interval& operator+=(const Interval& other) {
        *this = *this + other;
        return *this;
    }

    Interval& operator-=(const Interval& other) {
        *this = *this - other;
        return *this;
    }

    Interval& operator*=(const Interval& other) {
        *this = *this * other;
        return *this;
    }

    Interval& operator/=(const Interval& other) {
        *this = *this / other;
        return *this;
    }

    // Математические функции 
    // модуль интервала
    static Interval abs(const Interval& x) {
        if (x.lower >= T(0)) { return x; }
        else if (x.upper <= T(0)) { return -x; }
        else {
            // Если интервал пересекает 0
            return Interval(T(0), max_val(-x.lower, x.upper));
        }
    }

    // возведение в целую степень
    static Interval pow(const Interval& x, int n) {
        if (n == 0) { return Interval(T(1)); }
        if (n == 1) { return x; }
        if (n < 0) { return Interval(T(1)) / pow(x, -n); }

        // чётная степень: минимум может быть 0, если интервал пересекает 0
        if (n % 2 == 0) {
            if (x.lower >= T(0)) { return Interval(pow_val(x.lower, n), pow_val(x.upper, n)); }
            else if (x.upper <= T(0)) { return Interval(pow_val(x.upper, n), pow_val(x.lower, n)); }
            else {
                T max_pow = max_val(pow_val(x.lower, n), pow_val(x.upper, n));
                return Interval(T(0), max_pow);
            }
        } else {
            // Нечётная степень: монотонность сохранена
            return Interval(pow_val(x.lower, n), pow_val(x.upper, n));
        }
    }

    // Пересечение двух интервалов
    static Interval intersect(const Interval& a, const Interval& b) {
        T new_lower = max_val(a.lower, b.lower);
        T new_upper = min_val(a.upper, b.upper);

        if (new_lower > new_upper) {
            throw std::runtime_error("Intervals do not intersect");
        }

        return Interval(new_lower, new_upper);
    }

    // Объединение двух интервалов (наименьший интервал, содержащий оба)
    static Interval hull(const Interval& a, const Interval& b) {
        return Interval(min_val(a.lower, b.lower), max_val(a.upper, b.upper));
    }
};

// Внешние операции со скалярами (симметричные)
template <typename T>
Interval<T> operator+(T scalar, const Interval<T>& interval) {
    return interval + scalar;
}

template <typename T>
Interval<T> operator-(T scalar, const Interval<T>& interval) {
    return Interval<T>(scalar) - interval;
}

template <typename T>
Interval<T> operator*(T scalar, const Interval<T>& interval) {
    return interval * scalar;
}

template <typename T>
Interval<T> operator/(T scalar, const Interval<T>& interval) {
    return Interval<T>(scalar) / interval;
}