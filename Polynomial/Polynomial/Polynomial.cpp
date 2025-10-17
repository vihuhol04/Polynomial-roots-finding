#include "Polynomial.h" 

// тестирование методов
// исправлено: генерация полиномов отныне по типам корней 
// рассеянные + кластеризованные + кратные + комплексные

/*int main() {
    setlocale(LC_ALL, "Russian");

    try {
        using fp_t = double;

        const int NUM_TESTS = 2;
        const double TOLERANCE = 1e-7;

        std::cout << std::fixed << std::setprecision(6);

        for (int degree = 6; degree <= 10; degree++) {

            for (int test = 0; test < NUM_TESTS; ++test) {
                std::vector<double> coefficients;
                std::vector<double> real_roots;
                std::vector<std::complex<double>> complex_pairs;

                unsigned N_complex_pairs = degree / 4;  // Одна комплексная пара = 2 корня
                unsigned N_scattered = degree - 2 * N_complex_pairs;

                generate_high_degree_polynomial<fp_t>(
                    degree,
                    N_scattered,
                    N_complex_pairs,
                    2,                // максимальная кратность = 2
                    -1, 1,            // корни в диапазоне [-1, 1]
                    0.5,              // кластеры с радиусом 0.5
                    coefficients,     // выход: коэффициенты
                    real_roots,       // выход: вещественные корни
                    complex_pairs     // выход: комплексные пары
                );

                std::vector<fp_t> found_graeffe = find_real_roots_by_graeffe(coefficients, TOLERANCE);
                //auto found_graeffe_complex = find_complex_roots_by_graeffe(coefficients, TOLERANCE);
                //auto found_Durand_Kerner = find_roots_by_Durand_Kerner(coefficients, TOLERANCE);
                std::vector<fp_t> found_cf = find_roots_by_Continued_Fractions(coefficients, TOLERANCE);
                auto found_sagralov = find_roots_by_Sagralov(coefficients, TOLERANCE);

                std::cout << "\n============================\n";
                std::cout << "Тест #" << test + 1 << " | Степень: " << degree << "\n";

                std::cout << "\nПолином:\n";
                print_polynomial(coefficients);

                std::cout << "Истинные корни:\n";

                std::cout << "  Вещественные:\n";
                for (fp_t root : real_roots) {
                    print_value(root);
                    std::cout << "\n";
                }

                std::cout << "  Комплексные пары:\n";
                for (const auto& z : complex_pairs) {
                    print_value(z);
                    std::cout << "\n";
                }

                std::cout << "Корни вещественные (Греф):\n";
                for (const auto& root : found_graeffe) {
                    print_value(root);
                    std::cout << " ";
                }
                std::cout << "\n";

                std::cout << "Корни комплексные (Греф):\n";
                for (const auto& root : found_graeffe_complex) {
                    print_value(root);
                    std::cout << " ";
                }
                std::cout << "\n";

                std::cout << "Корни комплексные (Дюран-Кернер):\n";
                for (const auto& root : found_Durand_Kerner) {
                    print_value(root);
                    std::cout << " ";
                }
                std::cout << "\n";

                std::cout << "Корни (Цепные дроби):\n";
                for (const auto& root : found_cf) {
                    print_value(root);
                    std::cout << " ";
                }
                std::cout << "\n";

                std::cout << "Интервалы (Сагралов):\n";
                for (const auto& interval : found_sagralov) {
                    print_value(interval);
                    std::cout << " ";
                }
                std::cout << "\n";

                

                fp_t abs_err = 0, rel_err = 0;
                int excess = 0, lost = 0;
            }

        }
        std::cout << "\n";

    }
    catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return 1;
    }*/

    // Тестирование QD алгоритма
    /*std::vector<int_precision> coeffs = {1, 4, -5};
    std::vector<float_precision> roots = find_roots_by_QD(coeffs);
    for (auto root : roots) {
        print_value(root);
        std::cout << std::endl;
    }
    
    // Тестирование CAD алгоритма
    using CoeffType = float_precision;
    using Poly = mpolynomial<CoeffType>;
    Poly x_0 = Poly::variable(0);
    Poly p = pow(x_0, 2) - Poly(CoeffType(9));
    find_roots_by_CAD(p, 1);
    
    // Тестирование бисекции
    std::vector<std::complex<double>> coeffs1 = { 1, -4, 5 };
    auto roots1 = find_roots_by_Durand_Kerner(coeffs1, 1e-5);
    std::cout << "Корни:\n";
    for (const auto& z : roots1) {
        print_value(z);
        std::cout << std::endl;
    }
    
    double TOLERANCE = 1e-5;
    std::vector<double> coefficients = { 1, -4, 5 };
    auto found_graeffe_complex = find_complex_roots_by_graeffe(coefficients, TOLERANCE);
    
    for (auto r : found_graeffe_complex) {
        print_value(r);
        std::cout << '\n';
    }

    std::vector<double> coeffs = { 1, -3, -7, 27, -18 }; 
    auto roots = find_moduli_roots_by_graeffe(coeffs);

    std::cout << "Корни (приближённые):\n";
    for (auto r : roots) {
        print_value(r);
        std::cout << "\n";
    }

    std::vector<double> coeffs = { 0.3, -0.6, 0.6 }; // 0.3 x^2 -0.6 x + 0.6
    std::complex<long double> center(0.0L, 0.0L);
    long double halfw = 2.0L;
    auto disks = CIsolate<double>(coeffs, center, halfw, 1e-6L);
    std::cout << "Found disks: " << disks.size() << '\n';
    for (auto& d : disks) {
        std::cout << "center: " << d.center << ", radius: " << d.radius << '\n';
    }

    

    return 0;
}*/

int main() 
{
    setlocale(LC_ALL, "Russian");
    using fp_t = double;

    std::vector<fp_t> coefficients;
    std::vector<fp_t> real_roots_repeated;
    std::vector<fp_t> unique_real_roots;
    std::vector<unsigned> real_root_multiplicities;
    std::vector<std::complex<fp_t>> complex_roots;

    generate_high_degree_polynomial<fp_t>(
        4,
        2,
        0,
        std::vector<unsigned>{  },
        std::vector<fp_t>{ },
        std::vector<std::pair<unsigned, unsigned>>{ {} },
        0.001,
        false,
        12345,
        coefficients,
        real_roots_repeated,
        unique_real_roots,
        real_root_multiplicities,
        complex_roots
    );
    
    print_polynomial(coefficients);
    print_roots_summary(unique_real_roots, real_root_multiplicities, complex_roots);

    fp_t TOLERANCE = 1e-4;

    std::reverse(coefficients.begin(), coefficients.end());
    auto found_graeffe = find_moduli_with_multiplicities_by_graeffe(coefficients, TOLERANCE);
    std::cout << "Корни вещественные (Греф):\n";
    for (const auto& root : found_graeffe) {
        std::cout << "root: " << root.value << "    mult: " << root.multiplicity << "\n";
    }
    std::cout << "\n";

    auto found_graeffe_complex = hosman_with_multiplicities(coefficients, found_graeffe, TOLERANCE);
    std::cout << "Корни комплексные (Хосман):\n";
    for (auto& r : found_graeffe_complex) {
        std::cout << r.real() << " " << r.imag() << "i\n";
    }

    /*std::vector<double> coefficients1 = {1, -4, -5};
    std::vector<double> found_graeffe1 = find_moduli_roots_by_graeffe(coefficients1, 1e-5);
    std::cout << "Корни вещественные (Греф):\n";
    for (const auto& root : found_graeffe1) {
        std::cout << root << "\n";
    }
    std::cout << "\n";
    
    std::vector<double> coefficients_1 = {50, -10, -9, 2, 1};
    print_polynomial(coefficients_1);
    auto found_graeffe_complex1 = hosman_modification_graeffe(coefficients_1, 1e-6);
    std::cout << "Корни комплексные (Хосман):\n";
    for (auto& r : found_graeffe_complex1) {
        std::cout << r.real() << " " << (r.imag() >= 0 ? "+" : "-")
            << std::abs(r.imag()) << "i\n";
    }*/

    // вызов Сагралова для вещественных корней
    // вызов Сагралова для комплексных корней

    return 0;
}
