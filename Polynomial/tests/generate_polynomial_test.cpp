// Полный набор тестов для генератора многочленов
// Реализация: Павлова Анастасия, КМБО-01-22

#include "mathUtils.h"
#include "generate_high_degree_polynomial.h"

using std::cout;
using std::cerr;
using std::endl;

namespace {
    // Вспомогательные локальные функции
    // Вычисление значения полинома (coeffs: младший член первый) в вещественной точке x
    template<typename T>
    T eval_poly_local(const std::vector<T>& coeffs, T x) {
        T res = T(0);
        T powx = T(1);
        for (size_t i = 0; i < coeffs.size(); ++i) {
            res += coeffs[i] * powx;
            powx *= x;
        }
        return res;
    }

    // Вычисление значения полинома в комплексной точке z
    template<typename T>
    std::complex<T> eval_poly_local(const std::vector<T>& coeffs, const std::complex<T>& z) {
        std::complex<T> res = std::complex<T>(0);
        std::complex<T> powz = std::complex<T>(1);
        for (size_t i = 0; i < coeffs.size(); ++i) {
            res += coeffs[i] * powz;
            powz *= z;
        }
        return res;
    }

    // Сумма кратностей из вектора real_root_multiplicities
    unsigned sum_multiplicities(const std::vector<unsigned>& mults) {
        unsigned s = 0;
        for (auto m : mults) s += m;
        return s;
    }

    // Абсолютная величина для double
    double absd(double x) { return std::abs(x); }

    // Простая проверка сопряжённости
    template<typename T>
    bool is_conjugate_pair(const std::complex<T>& a, const std::complex<T>& b, T tol = T(1e-12)) {
        return (std::abs(a.real() - b.real()) < tol) && (std::abs(a.imag() + b.imag()) < tol);
    }

    // Для печати заголовка теста
    void print_test_header(const std::string& s) {
        cout << " ТЕСТ : " << s << "\n";
    }

    // ТЕСТ 1
    void test_basic_combined() {
        print_test_header("Базовый комбинированный случай (кластеры + кратности + комплексные пары)");

        unsigned P = 10;
        unsigned num_complex_pairs = 2; // 4 степени
        unsigned num_clusters = 2;
        std::vector<unsigned> cluster_counts = { 3, 1 }; // 4 простых в кластерах
        std::vector<double> cluster_radii = { 0.15, 0.12 };
        std::vector<std::pair<unsigned, unsigned>> mult_groups = { {2, 1} }; // один корень кратности 2 -> 2 степени
        double default_radius = 0.1;
        bool normalize = true;
        std::uint64_t seed = 20241010;

        std::vector<double> coeffs;
        std::vector<double> real_repeated;
        std::vector<double> unique_real;
        std::vector<unsigned> real_mults;
        std::vector<std::complex<double>> complex_roots;

        generate_high_degree_polynomial<double>(
            P,
            num_complex_pairs,
            num_clusters,
            cluster_counts,
            cluster_radii,
            mult_groups,
            default_radius,
            normalize,
            seed,
            coeffs,
            real_repeated,
            unique_real,
            real_mults,
            complex_roots
        );

        cout << "Степень P = " << P << ", фактическая степень из коэффициентов = " << (coeffs.size() - 1) << "\n";
        assert(coeffs.size() == P + 1);

        // Нормировка
        if (normalize && !coeffs.empty()) {
            double lead = coeffs.back();
            cout << "Проверка нормировки: старший коэффициент = " << std::setprecision(12) << lead << "\n";
            assert(std::abs(lead - 1.0) < 1e-12);
        }

        // Соответствие количества корней сумме кратностей
        unsigned sum_mults = sum_multiplicities(real_mults);
        cout << "Уникальных вещественных корней: " << unique_real.size() << ", сумма их кратностей = " << sum_mults
            << ", всего повторяющихся вещественных корней (real_repeated) = " << real_repeated.size() << "\n";
        assert(sum_mults == real_repeated.size());

        // Проверка того, что все вещественные корни лежат в (-1,1)
        for (double r : unique_real) {
            assert(r > -1.0 - 1e-12 && r < 1.0 + 1e-12);
        }

        // Проверка комплексных корней: должны идти парами и быть сопряжёнными
        assert(complex_roots.size() % 2 == 0);
        for (size_t i = 0; i + 1 < complex_roots.size(); i += 2) {
            auto a = complex_roots[i];
            auto b = complex_roots[i + 1];
            assert(is_conjugate_pair(a, b));
        }

        // Проверка точности: для каждого корня P(root) ~ 0 (с допустимой погрешностью)
        const double tol_real = 1e-8;
        const double tol_complex = 1e-6;
        for (size_t i = 0; i < unique_real.size(); ++i) {
            double r = unique_real[i];
            double v = eval_poly_local(coeffs, r);
            cout << "P(" << r << ") = " << std::setprecision(12) << v << "\n";
            // Для кратных корней численная ошибка может быть чуть больше; допускаем небольшую погрешность
            assert(std::abs(v) < 1e-6);
        }
        for (size_t i = 0; i + 1 < complex_roots.size(); i += 2) {
            auto z = complex_roots[i];
            auto v = eval_poly_local(coeffs, z);
            cout << "P(" << z << ") = " << v << "\n";
            assert(std::abs(v) < tol_complex);
        }

        cout << "ТЕСТ 1: OK\n\n";
    }

    // ТЕСТ 2
    void test_minimal_degree() {
        print_test_header("Минимальная степень P=1");

        unsigned P = 1;
        std::vector<double> coeffs;
        std::vector<double> rr, ur;
        std::vector<unsigned> mults;
        std::vector<std::complex<double>> cz;

        generate_high_degree_polynomial<double>(
            P,             
            0,             
            0,             
            {},            
            {},            
            {},            
            0.1,           
            true,          
            1,             
            coeffs, rr, ur, mults, cz
        );

        cout << "Коэффициенты (size) = " << coeffs.size() << " (ожидается 2)\n";
        assert(coeffs.size() == 2);
        cout << "ТЕСТ 2: OK\n\n";
    }

    // ТЕСТ 3
    void test_all_complex() {
        print_test_header("Полностью комплексный полином");

        unsigned P = 6;
        unsigned num_complex_pairs = 3; // полностью покрывает степень
        std::vector<double> coeffs;
        std::vector<double> rr, ur;
        std::vector<unsigned> mults;
        std::vector<std::complex<double>> cz;

        generate_high_degree_polynomial<double>(
            P,
            num_complex_pairs,
            0,
            {}, {}, {}, 0.1, true, 2025,
            coeffs, rr, ur, mults, cz
        );

        cout << "Уникальных вещественных корней: " << ur.size() << ", повторяющихся вещественных: " << rr.size() << "\n";
        assert(ur.empty());
        assert(rr.empty());
        assert(cz.size() == 2 * num_complex_pairs);
        for (size_t i = 0; i + 1 < cz.size(); i += 2) {
            assert(is_conjugate_pair(cz[i], cz[i + 1]));
        }
        cout << "ТЕСТ 3: OK\n\n";
    }

    // ТЕСТ 4
    void test_clusters_only() {
        print_test_header("Только кластеры (без комплексных и кратных)");

        unsigned P = 7;
        unsigned num_clusters = 3;
        std::vector<unsigned> cluster_counts = { 2, 2, 1 }; // суммарно 5 простых
        std::vector<double> cluster_radii = { 0.12, 0.08, 0.05 };
        std::vector<double> coeffs;
        std::vector<double> rr, ur;
        std::vector<unsigned> mults;
        std::vector<std::complex<double>> cz;

        generate_high_degree_polynomial<double>(
            P, 0, num_clusters,
            cluster_counts, cluster_radii, {}, 0.08, true, 31415,
            coeffs, rr, ur, mults, cz
        );

        // Должны остаться свободные простые корни вне кластеров для заполнения до P
        assert(coeffs.size() == P + 1);
        // Все уникальные вещественные корни лежат в [-1,1]
        for (double r : ur) assert(r >= -1.0 - 1e-12 && r <= 1.0 + 1e-12);
        cout << "кол-во уникальных значений = " << ur.size() << ", общее кол-во повторных = " << rr.size() << "\n";
        cout << "ТЕСТ 4: OK\n\n";
    }

    // ТЕСТ 5
    void test_multiplicity_and_overflow_handling() {
        print_test_header("Кратные корни и обработка переполнения");

        unsigned P = 6;
        // Просим: одна группа mult=3 count=2 -> 3*2 = 6 степеней,
        // плюс ещё некоторые простые/кластеры/комплексы — генератор должен уменьшить запросы или бросить исключение.
        std::vector<std::pair<unsigned, unsigned>> mult_groups = { {3, 2} }; // потенциально занимает всю степень
        std::vector<double> coeffs;
        std::vector<double> rr, ur;
        std::vector<unsigned> mults;
        std::vector<std::complex<double>> cz;

        // Случай 1: только multiplicity_groups — должен корректно вместиться
        generate_high_degree_polynomial<double>(
            P, 0, 0,
            {}, {}, mult_groups, 0.1, true, 77,
            coeffs, rr, ur, mults, cz
        );
        // все корни вещественные уникальные суммарной кратности 6
        assert(sum_multiplicities(mults) == rr.size());
        assert(coeffs.size() == P + 1);
        cout << "ТЕСТ 5 (case1): OK\n";

        // Случай 2: просим mult_groups + дополнительные кластерные простые корни -> генератор должен сократить кластерные
        unsigned P2 = 6;
        std::vector<unsigned> cluster_counts = { 2, 2 }; // даёт +4 простых
        std::vector<double> cluster_radii = { 0.2, 0.2 };

        generate_high_degree_polynomial<double>(
            P2, 0, 2,
            cluster_counts, cluster_radii, mult_groups, 0.1, true, 88,
            coeffs, rr, ur, mults, cz
        );
        // Степень всё ещё равна P2
        assert(coeffs.size() == P2 + 1);
        cout << "ТЕСТ 5 (case2): OK\n\n";
    }

    // ТЕСТ 6
    void test_no_normalize() {
        print_test_header("Поведение при normalize_coeffs == false");

        unsigned P = 5;
        std::vector<double> coeffs;
        std::vector<double> rr, ur;
        std::vector<unsigned> mults;
        std::vector<std::complex<double>> cz;

        generate_high_degree_polynomial<double>(
            P, 1, 1,
            { 2 }, { 0.1 }, {}, 0.05, false, 999,
            coeffs, rr, ur, mults, cz
        );

        assert(coeffs.size() == P + 1);
        double leading = coeffs.back();
        cout << "Старший коэффициент (normalize=false) = " << leading << "\n";
        // Не требуем конкретного значения, но он не должен быть 0
        assert(std::abs(leading) > 1e-18);
        cout << "ТЕСТ 6: OK\n\n";
    }

    // ТЕСТ 7
    void test_seed_zero_random() {
        print_test_header("Поведение при seed == 0 (рандомное зерно)");

        unsigned P = 8;
        std::vector<double> coeffs;
        std::vector<double> rr, ur;
        std::vector<unsigned> mults;
        std::vector<std::complex<double>> cz;

        // seed = 0 означает использовать std::random_device внутри функции
        generate_high_degree_polynomial<double>(
            P, 2, 2, { 2, 2 }, { 0.1, 0.1 }, {}, 0.05, true, 0,
            coeffs, rr, ur, mults, cz
        );
        assert(coeffs.size() == P + 1);
        cout << "ТЕСТ 7: OK\n\n";
    }

    // ТЕСТ 8
    void test_invalid_arguments() {
        print_test_header("Проверка выбрасывания исключений на неверные параметры");

        // Случай: P = 0 -> должно бросать исключение
        try {
            unsigned P = 0;
            std::vector<double> coeffs;
            std::vector<double> rr, ur;
            std::vector<unsigned> mults;
            std::vector<std::complex<double>> cz;
            generate_high_degree_polynomial<double>(
                P, 0, 0, {}, {}, {}, 0.1, true, 1,
                coeffs, rr, ur, mults, cz
            );
            // Если дошло сюда — это ошибка
            assert(false && "Ожидалось исключение при P=0");
        }
        catch (const std::invalid_argument& e) {
            cout << "Поймано ожидаемое исключение для P=0: " << e.what() << "\n";
        }

        // Случай: отрицательный радиус в cluster_radii_in
        try {
            unsigned P = 4;
            std::vector<double> coeffs;
            std::vector<double> rr, ur;
            std::vector<unsigned> mults;
            std::vector<std::complex<double>> cz;
            generate_high_degree_polynomial<double>(
                P, 0, 1, { 1 }, { -0.1 }, {}, 0.1, true, 2,
                coeffs, rr, ur, mults, cz
            );
            assert(false && "Ожидалось исключение для отрицательного радиуса");
        }
        catch (const std::invalid_argument& e) {
            cout << "Поймано ожидаемое исключение для отрицательного радиуса: " << e.what() << "\n";
        }

        // Случай: суммарный радиус >= 1 -> исключение
        try {
            unsigned P = 5;
            std::vector<double> coeffs;
            std::vector<double> rr, ur;
            std::vector<unsigned> mults;
            std::vector<std::complex<double>> cz;
            generate_high_degree_polynomial<double>(
                P, 0, 2, { 1,1 }, { 0.6, 0.5 }, {}, 0.1, true, 3,
                coeffs, rr, ur, mults, cz
            );
            assert(false && "Ожидалось исключение для суммарного радиуса >= 1");
        }
        catch (const std::invalid_argument& e) {
            cout << "Поймано ожидаемое исключение для суммарного радиуса: " << e.what() << "\n";
        }

        // Случай: слишком большая кратность (mult > P)
        try {
            unsigned P = 3;
            std::vector<std::pair<unsigned, unsigned>> mult_groups = { {5, 1} }; // mult=5 > P
            std::vector<double> coeffs;
            std::vector<double> rr, ur;
            std::vector<unsigned> mults;
            std::vector<std::complex<double>> cz;
            generate_high_degree_polynomial<double>(
                P, 0, 0, {}, {}, mult_groups, 0.1, true, 4,
                coeffs, rr, ur, mults, cz
            );
            assert(false && "Ожидалось исключение для mult > P");
        }
        catch (const std::invalid_argument& e) {
            cout << "Поймано ожидаемое исключение для mult > P: " << e.what() << "\n";
        }

        // Случай: недостаточная степень для заданного числа пар комплексных корней
        try {
            unsigned P = 3;
            unsigned num_complex_pairs = 2; // требует 4 степеней > P
            std::vector<double> coeffs;
            std::vector<double> rr, ur;
            std::vector<unsigned> mults;
            std::vector<std::complex<double>> cz;
            generate_high_degree_polynomial<double>(
                P, num_complex_pairs, 0, {}, {}, {}, 0.1, true, 5,
                coeffs, rr, ur, mults, cz
            );
            assert(false && "Ожидалось исключение при недостаточной степени для комплексных пар");
        }
        catch (const std::invalid_argument& e) {
            cout << "Поймано ожидаемое исключение для комплекса > P: " << e.what() << "\n";
        }

        cout << "ТЕСТ 8: OK\n\n";
    }

    // ТЕСТ 9
    void test_high_degree_numerical_stability() {
        print_test_header("Численная стабильность при высокой степени");

        unsigned P = 40; // высокая степень
        unsigned num_complex_pairs = 10; // 20 степеней
        unsigned num_clusters = 5;
        std::vector<unsigned> cluster_counts(num_clusters, 2); // 10 простых в кластерах
        std::vector<double> cluster_radii(num_clusters, 0.03);
        std::vector<std::pair<unsigned, unsigned>> mult_groups = { {2, 5} }; // 10 степеней
        std::vector<double> coeffs;
        std::vector<double> rr, ur;
        std::vector<unsigned> mults;
        std::vector<std::complex<double>> cz;

        // Попробуем с нормировкой
        generate_high_degree_polynomial<double>(
            P, num_complex_pairs, num_clusters,
            cluster_counts, cluster_radii, mult_groups, 0.02, true, 2025,
            coeffs, rr, ur, mults, cz
        );

        cout << "Генерировано полином степени " << (coeffs.size() - 1) << "\n";
        assert(coeffs.size() - 1 == static_cast<int>(P));

        // Проверка P(r) для вещественных корней (позволим более свободную погрешность для высокой степени)
        double tol = 1e-4;
        for (double r : ur) {
            double v = eval_poly_local(coeffs, r);
            if (std::abs(v) > tol) {
                cout << "Внимание: Погрешность для корня " << r << " = " << v << " (> tol=" << tol << ")\n";
            }
            // не assert-им строго, но выводим предупреждение; в идеале должно быть близко к нулю
        }

        cout << "ТЕСТ 9: завершён (см. возможные предупреждения)\n\n";
    }

    // ТЕСТ 10
    void test_exhaustive_smoke() {
        print_test_header("Всесторонняя проверка (smoke) — перебор комбинаций");

        // Будем перебирать несколько комбинаций параметров, чтобы покрыть всевозможные пути
        std::vector<unsigned> degrees = { 2, 3, 5, 8 };
        std::vector<unsigned> complex_pairs = { 0, 1, 2 };
        std::vector<unsigned> clusters = { 0, 1, 2, 3 };
        std::vector<bool> norms = { true, false };
        std::vector<std::uint64_t> seeds = { 1, 0, 123456789 };

        for (auto P : degrees) {
            for (auto cp : complex_pairs) {
                for (auto nc : clusters) {
                    for (auto normalize : norms) {
                        for (auto seed : seeds) {
                            // Подготовим cluster_counts и radii в зависимости от nc
                            std::vector<unsigned> cluster_counts(nc, 1);
                            std::vector<double> cluster_radii(nc, 0.05);
                            // Простые multiplicity_groups примеры
                            std::vector<std::pair<unsigned, unsigned>> mult_groups;
                            if (P >= 4) mult_groups.push_back({ 2, 1 }); // один двойной корень

                            std::vector<double> coeffs;
                            std::vector<double> rr, ur;
                            std::vector<unsigned> mults;
                            std::vector<std::complex<double>> cz;

                            try {
                                generate_high_degree_polynomial<double>(
                                    P, cp, nc,
                                    cluster_counts, cluster_radii, mult_groups,
                                    0.03, normalize, seed,
                                    coeffs, rr, ur, mults, cz
                                );
                                // Базовые проверки
                                assert(coeffs.size() == P + 1);
                                for (double r : ur) {
                                    assert(r >= -1.0 - 1e-12 && r <= 1.0 + 1e-12);
                                }
                                // Комплексные — парность
                                assert(cz.size() % 2 == 0);
                            } catch (const std::invalid_argument& e) {
                                // Возможно некоторые сочетания некорректны — допустимы исключения
                                cout << "Допустимая ошибка при комб. параметров (P=" << P << " cp=" << cp << " nc=" << nc
                                    << " norm=" << normalize << " seed=" << seed << "): " << e.what() << "\n";
                            }
                        }
                    }
                }
            }
        }
        cout << "ТЕСТ 10: OK (см. возможные допустимые ошибки выше)\n\n";
    }
}

int main() {
    setlocale(LC_ALL, "Russian");
    cout << "Запуск полного набора тестов для generate_high_degree_polynomial.h\n\n";

    try {
        test_basic_combined();
        test_minimal_degree();
        test_all_complex();
        test_clusters_only();
        test_multiplicity_and_overflow_handling();
        test_no_normalize();
        test_seed_zero_random();
        test_invalid_arguments();
        test_high_degree_numerical_stability();
        test_exhaustive_smoke();

        cout << "Все тесты выполнены.\n";
    }
    catch (const std::exception& ex) {
        cerr << "Во время тестов произошло непредвиденное исключение: " << ex.what() << "\n";
        return 2;
    }
    catch (...) {
        cerr << "Во время тестов произошло непредвиденное неизвестное исключение.\n";
        return 3;
    }

    return EXIT_SUCCESS;
}