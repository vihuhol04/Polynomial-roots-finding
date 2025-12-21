#include "Polynomial.h" 

// тестирование методов
// исправлено: генерация полиномов отныне по типам корней 
// рассеянные + кластеризованные + кратные + комплексные

using namespace std;

int main() 
{
    setlocale(LC_ALL, "Russian");
    /*using fp_t = float_precision;

    std::vector<fp_t> coefficients;
    std::vector<fp_t> real_roots_repeated;
    std::vector<fp_t> unique_real_roots;
    std::vector<unsigned> real_root_multiplicities;
    std::vector<std::complex<fp_t>> complex_roots;

    generate_high_degree_polynomial<fp_t>(
        10,
        0,
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
    fp_t TOLERANCE = 1e-8;

    cout << endl << "ИСХОДНЫЕ КОРНИ:" << endl;
    cout << "Вещественные корни: ";
    for (auto root : real_roots_repeated) {
        cout << root << " ";
    }
    cout << endl;

    cout << endl << "ВЫЧИСЛЕНИЕ МЕТОДОМ ГРЕФЕ" << endl;
    // меняем порядок коэффициентов
    vector<fp_t> coefficients_for_graeffe = coefficients;
    reverse(coefficients_for_graeffe.begin(), coefficients_for_graeffe.end());

    // ищем корни методом Грефе
    vector<fp_t> computed_moduli = find_moduli_roots_by_graeffe(coefficients_for_graeffe, TOLERANCE);
    cout << "Вычисленные модули корней: ";
    if (computed_moduli.empty()) {
        cout << "метод не нашел корней";
    }
    else {
        for (auto modulus : computed_moduli)
            cout << modulus << " ";
    }
    cout << endl << endl;*/

    std::vector<double> coeffs = { -6.0, 11.0, -6.0, 1.0 }; //{ -6.0, 11.0, -6.0, 1.0 };

    // Вызов метода Лагерра
    auto roots = find_roots_by_Continued_Fractions(coeffs, 1e-6);

    // Вывод результатов
    std::cout << "Корни полинома: " << std::endl;
    for (const auto& root : roots) {
        std::cout << "  " << root << std::endl;
    }

    return 0;
}
