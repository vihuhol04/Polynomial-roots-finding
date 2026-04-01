#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace wilf {
template <typename Real>
inline constexpr bool is_supported_real_v = std::is_floating_point_v<Real>;
enum class CoefficientOrder { ascending, descending };
namespace detail {
template <typename Real>
constexpr Real typed_constant(long double value) noexcept {
    static_assert(is_supported_real_v<Real>, "Шаблонный параметр типа должен быть вещественным");
    return static_cast<Real>(value);
}
template <typename Coefficient, typename Real>
inline constexpr bool is_supported_coefficient_v =
    std::is_same_v<Coefficient, Real> || std::is_same_v<Coefficient, std::complex<Real>>;
inline constexpr std::uint64_t default_random_seed = 123456789ULL;
inline constexpr std::size_t default_max_depth = 80;
inline constexpr std::size_t default_max_startup_doublings = 160;
inline constexpr std::size_t startup_attempts = 8;
inline constexpr std::size_t default_max_recenter_attempts = 3;
inline constexpr std::size_t sturm_sequence_extra_terms = 4;
template <typename Real>
struct numeric_constants {
    static constexpr Real default_shift_fraction = typed_constant<Real>(0.1L);
    static constexpr Real default_sign_safety = typed_constant<Real>(64.0L);
    static constexpr Real complex_trim_safety_factor = typed_constant<Real>(32.0L);
    static constexpr Real real_trim_safety_factor = typed_constant<Real>(64.0L);
    static constexpr Real corner_value_safety_factor = typed_constant<Real>(32.0L);
    static constexpr Real split_margin_fraction = typed_constant<Real>(0.1L);
    static constexpr Real lead_tolerance_safety_factor = typed_constant<Real>(64.0L);
    static constexpr Real startup_box_expansion_factor = typed_constant<Real>(2.0L);
};
}  
template <typename Real>
struct RootBox {
    static_assert(is_supported_real_v<Real>, "Шаблонный параметр типа должен быть вещественным");
    std::complex<Real> center{};
    Real width{}, height{}, error_bound{};
    std::size_t multiplicity{};
};
template <typename Real>
struct Options {
    static_assert(is_supported_real_v<Real>, "Шаблонный параметр типа должен быть вещественным");
    CoefficientOrder coefficient_order = CoefficientOrder::descending;
    std::uint64_t random_seed = detail::default_random_seed;
    std::size_t max_depth = detail::default_max_depth;
    std::size_t max_startup_doublings =
        detail::default_max_startup_doublings;
    std::size_t max_recenter_attempts =
        detail::default_max_recenter_attempts;
    Real shift_fraction =
        detail::numeric_constants<Real>::default_shift_fraction;
    Real sign_safety =
        detail::numeric_constants<Real>::default_sign_safety;
};
template <typename Real>
class Solver {
    static_assert(is_supported_real_v<Real>, "Шаблонный параметр типа должен быть вещественным");
public:
    using real_type = Real;
    using complex_type = std::complex<Real>;
    using root_box_type = RootBox<Real>;
    using options_type = Options<Real>;

    template <typename Coefficient>
    std::vector<root_box_type> localize(
        const std::vector<Coefficient>& coefficients,
        Real epsilon,
        const options_type& options = {}) const {
        // Внутри алгоритм всегда работает с комплексными
        // коэффициентами и рекурсивно делит стартовую область.
        static_assert(
            detail::is_supported_coefficient_v<Coefficient, Real>,
            "Коэффициенты должны иметь тип Real или std::complex<Real>");
        validate_options(options);
        auto polynomial = normalize_coefficients(
            to_complex_coefficients(coefficients),
            options.coefficient_order);
        const std::size_t degree = polynomial.size() - 1;

        if (degree == 0) {
            throw std::invalid_argument("Степень полинома должна быть положительной");
        }
        if (!std::isfinite(epsilon) || !(epsilon > zero)) {
            throw std::invalid_argument("Точность должна быть конечной и положительной");
        }

        std::mt19937_64 engine(options.random_seed);
        auto initial = initialize_box(polynomial, degree, options, engine);
        std::vector<BoxState> stack;
        stack.reserve(degree);
        stack.push_back(std::move(initial));
        std::vector<root_box_type> result;
        result.reserve(degree);
        const auto append_result = [&](const BoxState& state) {
            const Real width = state.box.width();
            const Real height = state.box.height();
            result.push_back(
                {
                    state.box.center(),
                    width,
                    height,
                    Real(0.5) * std::max(width, height),
                    state.multiplicity
                });
        };

        while (!stack.empty()) {
            BoxState entry = std::move(stack.back());
            stack.pop_back();

            if (entry.multiplicity == 0) continue;

            if (
                entry.depth >= options.max_depth ||
                std::max(entry.box.width(), entry.box.height()) <= epsilon) {
                append_result(entry);
                continue;
            }

            const auto refinement = refine_box(polynomial, entry, options, engine);
            if (!refinement) {
                append_result(entry);
                continue;
            }

            for (std::size_t index = refinement->size(); index > 0; --index) {
                const std::size_t child_index = index - 1;
                if ((*refinement)[child_index].multiplicity > 0) {
                    stack.push_back(std::move((*refinement)[child_index]));
                }
            }
        }

        std::size_t total = 0;
        for (const auto& box : result) total += box.multiplicity;
        if (total != degree) {
            throw std::runtime_error(
                "Внутренняя ошибка: сумма кратностей не совпала "
                "со степенью полинома");
        }

        std::sort(
            result.begin(),
            result.end(),
            [](const root_box_type& lhs, const root_box_type& rhs) {
                return std::make_tuple(
                           lhs.center.real(),
                           lhs.center.imag(),
                           lhs.multiplicity) <
                    std::make_tuple(
                        rhs.center.real(),
                        rhs.center.imag(),
                        rhs.multiplicity);
            });

        return result;
    }

    template <typename Coefficient>
    std::vector<complex_type> find_roots(
        const std::vector<Coefficient>& coefficients,
        Real epsilon,
        const options_type& options = {}) const {
        const auto boxes = localize(coefficients, epsilon, options);
        std::size_t total = 0;
        for (const auto& box : boxes) total += box.multiplicity;

        std::vector<complex_type> roots;
        roots.reserve(total);
        for (const auto& box : boxes) roots.insert(roots.end(), box.multiplicity, box.center);
        return roots;
    }

private:
    using ComplexPolynomial = std::vector<complex_type>;
    using RealPolynomial = std::vector<Real>;
    using constants = detail::numeric_constants<Real>;
    static constexpr std::size_t side_count = 4;
    static constexpr std::size_t child_count = 2;

    struct Box {
        Real x_min{}, x_max{}, y_min{}, y_max{};
        Real width() const noexcept { return x_max - x_min; }
        Real height() const noexcept { return y_max - y_min; }
        complex_type center() const noexcept {
            return {(x_min + x_max) / two, (y_min + y_max) / two};
        }
    };

    using SideProfile = std::vector<RealPolynomial>;
    using SideProfilePtr = std::shared_ptr<const SideProfile>;

    struct SideSlice {
        SideProfilePtr profile;
        Real t_begin{}, t_end{};
    };

    using Sides = std::array<SideSlice, side_count>;
    using CountResult = std::optional<std::size_t>;

    struct BoxState {
        Box box;
        std::size_t multiplicity{}, depth{};
        // Для каждого прямоугольника кэшируются подготовленные
        // профили сторон, чтобы не строить их заново.
        Sides sides{};
    };

    enum class SplitAxis { vertical, horizontal };

    struct Split {
        SplitAxis axis{};
        Real coordinate{}, lower_extent{}, upper_extent{};
    };

    using ChildStates = std::array<BoxState, child_count>;
    using ChildSides = std::array<Sides, child_count>;
    static constexpr Real zero = static_cast<Real>(0);
    static constexpr Real one = static_cast<Real>(1);
    static constexpr Real two = static_cast<Real>(2);
    static constexpr Real minus_one = static_cast<Real>(-1);

    static constexpr Real machine_epsilon() noexcept {
        return std::numeric_limits<Real>::epsilon();
    }

    static void validate_options(const options_type& options) {
        // Эти проверки отсекают параметры, при которых рекурсия и
        // численная эвристика теряют смысл.
        const auto require_positive = [](std::size_t value, const char* message) {
            if (value == 0) throw std::invalid_argument(message);
        };
        require_positive(options.max_depth, "Максимальная глубина должна быть положительной");
        require_positive(
            options.max_startup_doublings,
            "Максимальное число удвоений стартового квадрата должно быть положительным");
        if (!std::isfinite(options.shift_fraction) || options.shift_fraction < zero)
            throw std::invalid_argument("Доля смещения должна быть конечной и неотрицательной");
        if (!std::isfinite(options.sign_safety) || !(options.sign_safety > zero))
            throw std::invalid_argument("Запас устойчивости знака должен быть конечным и положительным");
    }

    template <typename Coefficient>
    static ComplexPolynomial to_complex_coefficients(
        const std::vector<Coefficient>& coefficients) {
        // Вещественный полином переводится в комплексную форму
        if constexpr (std::is_same_v<Coefficient, complex_type>) {
            return coefficients;
        } else {
            ComplexPolynomial result;
            result.reserve(coefficients.size());
            for (const Real value : coefficients) result.emplace_back(value, zero);
            return result;
        }
    }

    static ComplexPolynomial normalize_coefficients(
        const std::vector<complex_type>& coefficients,
        CoefficientOrder order) {
        if (coefficients.empty()) {
            throw std::invalid_argument("Вектор коэффициентов не должен быть пустым");
        }
        for (const complex_type& coefficient : coefficients) {
            if (!std::isfinite(coefficient.real()) || !std::isfinite(coefficient.imag())) {
                throw std::invalid_argument("Коэффициенты полинома должны быть конечными");
            }
        }

        ComplexPolynomial normalized = order == CoefficientOrder::ascending
            ? coefficients
            : ComplexPolynomial(coefficients.rbegin(), coefficients.rend());

        const std::size_t original_size = normalized.size();
        trim_polynomial(normalized, constants::complex_trim_safety_factor);
        if (normalized.size() != original_size) {
            throw std::invalid_argument(
                "Старший коэффициент равен нулю или слишком мал для указанной степени полинома");
        }
        if (normalized.size() < 2) {
            throw std::invalid_argument(
                "Степень полинома должна быть положительной, "
                "а старший коэффициент ненулевым");
        }
        normalize_polynomial(normalized);
        return normalized;
    }

    template <typename Polynomial>
    static void trim_polynomial(Polynomial& polynomial, Real safety_factor) {
        // Численный хвост обрезается относительно масштаба
        // коэффициентов, а не по абсолютному порогу.
        const Real scale = max_abs(polynomial);
        const Real tolerance = (scale == zero)
            ? zero
            : safety_factor * machine_epsilon() * scale;
        while (!polynomial.empty() && std::abs(polynomial.back()) <= tolerance) {
            polynomial.pop_back();
        }
    }

    template <typename Polynomial>
    static Real max_abs(const Polynomial& polynomial) {
        Real value = zero;
        for (const auto& coefficient : polynomial) {
            value = std::max(value, std::abs(coefficient));
        }
        return value;
    }

    template <typename Polynomial>
    static void normalize_polynomial(Polynomial& polynomial) {
        const Real scale = max_abs(polynomial);
        if (scale == zero) return;
        for (auto& coefficient : polynomial) coefficient /= scale;
    }

    static bool prepare_real_polynomial(RealPolynomial& polynomial) {
        trim_polynomial(polynomial, constants::real_trim_safety_factor);
        if (polynomial.empty()) return false;
        normalize_polynomial(polynomial);
        return true;
    }

    static void normalize_complex_polynomial(ComplexPolynomial& polynomial) {
        trim_polynomial(polynomial, constants::complex_trim_safety_factor);
        normalize_polynomial(polynomial);
    }

    static std::pair<Real, Real> evaluate_and_absolute_bound(
        const RealPolynomial& polynomial,
        Real x) {
        const Real ax = std::abs(x);
        Real value = zero;
        Real bound = zero;
        for (auto it = polynomial.rbegin(); it != polynomial.rend(); ++it) {
            value = value * x + *it;
            bound = bound * ax + std::abs(*it);
        }
        return {value, bound};
    }

    template <typename Polynomial, typename Value>
    static auto evaluate(const Polynomial& polynomial, const Value& x) {
        typename Polynomial::value_type value{};
        for (auto it = polynomial.rbegin(); it != polynomial.rend(); ++it) {
            value = value * x + *it;
        }
        return value;
    }

    static SideSlice segment(const SideSlice& slice, Real from, Real to) {
        const Real direction = (slice.t_end >= slice.t_begin) ? one : minus_one;
        return {slice.profile, slice.t_begin + direction * from, slice.t_begin + direction * to};
    }

    static Box box_from_center(const complex_type& center, Real side) {
        const Real half = side / two;
        return {
            center.real() - half,
            center.real() + half,
            center.imag() - half,
            center.imag() + half
        };
    }

    static Box scaled_about_center(const Box& box, Real factor) {
        const complex_type center = box.center();
        const Real half_width = box.width() * factor / two;
        const Real half_height = box.height() * factor / two;
        return {
            center.real() - half_width,
            center.real() + half_width,
            center.imag() - half_height,
            center.imag() + half_height
        };
    }

    static std::optional<Split> choose_split(
        const Box& box,
        Real width,
        Real height,
        const options_type& options,
        std::size_t attempt,
        std::mt19937_64& engine) {
        const SplitAxis axis =
            width >= height ? SplitAxis::vertical : SplitAxis::horizontal;
        const complex_type center = box.center();
        const Real length = axis == SplitAxis::vertical ? width : height;
        Real coordinate =
            axis == SplitAxis::vertical ? center.real() : center.imag();
        const Real radius = options.shift_fraction * length;
        if (attempt > 0 && radius > zero) {
            std::uniform_real_distribution<Real> distribution(-radius, radius);
            coordinate += distribution(engine);
        }

        const Real minimum =
            axis == SplitAxis::vertical ? box.x_min : box.y_min;
        const Real maximum =
            axis == SplitAxis::vertical ? box.x_max : box.y_max;
        const Real margin = length * constants::split_margin_fraction;
        coordinate = std::clamp(coordinate, minimum + margin, maximum - margin);
        if (!(coordinate > minimum && coordinate < maximum)) {
            return std::nullopt;
        }

        return Split{
            axis,
            coordinate,
            coordinate - minimum,
            maximum - coordinate
        };
    }

    static std::array<Box, child_count> split_children(const Box& box, const Split& split) {
        if (split.axis == SplitAxis::vertical) {
            return {{
                {box.x_min, split.coordinate, box.y_min, box.y_max},
                {split.coordinate, box.x_max, box.y_min, box.y_max}
            }};
        }

        return {{
            {box.x_min, box.x_max, box.y_min, split.coordinate},
            {box.x_min, box.x_max, split.coordinate, box.y_max}
        }};
    }

    static ChildSides split_child_sides(
        const Sides& parent,
        const SideSlice& cross_cut,
        const Split& split) {
        if (split.axis == SplitAxis::vertical) {
            const Real width = split.lower_extent + split.upper_extent;
            const Real height = std::abs(cross_cut.t_end - cross_cut.t_begin);
            return {{
                {{
                    segment(parent[0], zero, split.lower_extent),
                    segment(cross_cut, zero, height),
                    segment(parent[2], split.upper_extent, width),
                    parent[3]
                }},
                {{
                    segment(parent[0], split.lower_extent, width),
                    parent[1],
                    segment(parent[2], zero, split.upper_extent),
                    segment(cross_cut, height, zero)
                }}
            }};
        }

        const Real width = std::abs(cross_cut.t_end - cross_cut.t_begin);
        const Real height = split.lower_extent + split.upper_extent;
        return {{
            {{
                parent[0],
                segment(parent[1], zero, split.lower_extent),
                segment(cross_cut, width, zero),
                segment(parent[3], split.upper_extent, height)
            }},
            {{
                segment(cross_cut, zero, width),
                segment(parent[1], split.lower_extent, height),
                parent[2],
                segment(parent[3], zero, split.upper_extent)
            }}
        }};
    }

    static SideProfilePtr build_side_profile_from_transformed(ComplexPolynomial transformed) {
        // На параметризованной стороне p(z(t)) раскладывается на Re p(t) и
        // Im p(t),
        // а дальше по их последовательности Штурма можно считать
        // смены знака устойчиво.
        RealPolynomial real_part, imag_part;
        real_part.reserve(transformed.size()); imag_part.reserve(transformed.size());
        for (const auto& coefficient : transformed) {
            real_part.push_back(coefficient.real());
            imag_part.push_back(coefficient.imag());
        }
        const auto sequence = build_sturm_sequence(std::move(real_part), std::move(imag_part));
        if (!sequence) return {};
        return std::make_shared<SideProfile>(std::move(*sequence));
    }

    static std::optional<Sides> build_box_sides(
        const ComplexPolynomial& polynomial,
        const Box& box) {
        // На каждой стороне строится вещественный профиль, который
        // потом используется в подсчёте индекса.
        const auto build = [&](const complex_type& start, const complex_type& direction, Real length)
            -> std::optional<SideSlice> {
            const auto profile = build_side_profile_from_transformed(compose_linear(polynomial, start, direction));
            if (!profile) return std::nullopt;
            return SideSlice{profile, zero, length};
        };

        const auto bottom = build({box.x_min, box.y_min}, {one, zero}, box.width());
        const auto right = build({box.x_max, box.y_min}, {zero, one}, box.height()),
                   top = build({box.x_max, box.y_max}, {minus_one, zero}, box.width()),
                   left = build({box.x_min, box.y_max}, {zero, minus_one}, box.height());
        if (!bottom || !right || !top || !left) return std::nullopt;
        return {{{*bottom, *right, *top, *left}}};
    }

    static Real cauchy_radius(const ComplexPolynomial& polynomial) {
        // Оценка Коши даёт круг, в котором гарантированно лежат все
        // корни полинома.
        const Real lead = std::abs(polynomial.back());
        if (lead == zero) {
            throw std::invalid_argument(
                "Старший коэффициент должен быть ненулевым");
        }

        Real ratio = zero;
        for (std::size_t i = 0; i + 1 < polynomial.size(); ++i)
            ratio = std::max(ratio, std::abs(polynomial[i]) / lead);
        return one + ratio;
    }

    static complex_type random_center(std::mt19937_64& engine) {
        std::uniform_real_distribution<Real> distribution(zero, one);
        return {distribution(engine), distribution(engine)};
    }

    static Real random_positive_side(std::mt19937_64& engine) {
        std::uniform_real_distribution<Real> distribution(zero, one);
        const Real side = distribution(engine);
        return side <= std::sqrt(machine_epsilon()) ? one : side;
    }

    static Box cauchy_cover_box(Real radius, const complex_type& center = {}) {
        const Real reach = radius + std::max(std::abs(center.real()), std::abs(center.imag()));
        return box_from_center(center, constants::startup_box_expansion_factor * reach);
    }

    static std::array<complex_type, 9> startup_centers(
        Real radius,
        const options_type& options) {
        const Real shift = options.shift_fraction * std::max(one, radius);
        return {{
            {zero, zero},
            {shift, zero},
            {-shift, zero},
            {zero, shift},
            {zero, -shift},
            {shift, shift},
            {shift, -shift},
            {-shift, shift},
            {-shift, -shift}
        }};
    }

    static std::optional<BoxState> try_expand_box_to_cover_all_roots(
        const ComplexPolynomial& polynomial,
        std::size_t degree,
        const Box& initial_box,
        const options_type& options) {
        Box box = initial_box;
        for (std::size_t k = 0; k < options.max_startup_doublings; ++k) {
            if (const auto sides = build_box_sides(polynomial, box); sides) {
                // Если индекс на границе равен степени полинома,
                // стартовая область уже покрывает все корни.
                if (
                    const CountResult count =
                        count_zeros(polynomial, box, *sides, options);
                    count && *count == degree)
                    return BoxState{box, *count, 0, *sides};
            }
            box = scaled_about_center(box, constants::startup_box_expansion_factor);
        }
        return std::nullopt;
    }

    static BoxState initialize_box(
        const ComplexPolynomial& polynomial,
        std::size_t degree,
        const options_type& options,
        std::mt19937_64& engine) {
        // Стартовый квадрат ищется случайно, а затем при
        // необходимости расширяется до покрытия всех корней.
        for (std::size_t attempt = 0; attempt < detail::startup_attempts; ++attempt) {
            if (const auto initial = try_expand_box_to_cover_all_roots(
                    polynomial,
                    degree,
                    box_from_center(random_center(engine), random_positive_side(engine)),
                    options)) {
                return *initial;
            }
        }
        const Real radius = cauchy_radius(polynomial);
        for (const complex_type center : startup_centers(radius, options)) {
            if (const auto initial = try_expand_box_to_cover_all_roots(
                    polynomial,
                    degree,
                    cauchy_cover_box(radius, center),
                    options)) {
                return *initial;
            }
        }

        throw std::runtime_error("Не удалось построить начальный квадрат, содержащий все корни");
    }

    static bool corners_are_stable(const ComplexPolynomial& polynomial, const Box& box) {
        // Корень слишком близко к углу делает изменение аргумента на
        // соседних сторонах численно неоднозначным.
        // Исходный полином заранее нормализован по максимальному модулю коэффициента.
        const Real tolerance = constants::corner_value_safety_factor * machine_epsilon();

        for (const complex_type corner : std::array<complex_type, side_count>{{
                 {box.x_min, box.y_min},
                 {box.x_max, box.y_min},
                 {box.x_max, box.y_max},
                 {box.x_min, box.y_max}
             }}) {
            if (std::abs(evaluate(polynomial, corner)) <= tolerance) {
                return false;
            }
        }

        return true;
    }

    static CountResult count_zeros(
        const ComplexPolynomial& polynomial,
        const Box& box,
        const Sides& sides,
        const options_type& options) {
        // Здесь восстанавливается число корней внутри
        // прямоугольника по изменениям аргумента на его границе.
        if (!corners_are_stable(polynomial, box)) return std::nullopt;

        int sum = 0;
        for (const SideSlice& side : sides) {
            if (!side.profile) return std::nullopt;
            // Разность числа смен знака на концах стороны даёт вклад
            // этой стороны в индекс обхода.
            const auto v0 = sign_variations(*side.profile, side.t_begin, options);
            const auto v1 = sign_variations(*side.profile, side.t_end, options);
            if (!v0 || !v1) return std::nullopt;
            sum += *v1 - *v0;
        }

        return (sum >= 0 && (sum & 1) == 0)
            ? CountResult{static_cast<std::size_t>(sum / 2)}
            : std::nullopt;
    }

    static std::optional<ChildStates> refine_box(
        const ComplexPolynomial& polynomial,
        const BoxState& entry,
        const options_type& options,
        std::mt19937_64& engine) {
        // Один шаг рекурсии: строим один внутренний разрез и делим
        // прямоугольник на две дочерние области.
        const Box& box = entry.box;
        const Real width = box.width(), height = box.height();

        for (std::size_t attempt = 0; attempt <= options.max_recenter_attempts; ++attempt) {
            ChildStates refined{};
            const auto split = choose_split(box, width, height, options, attempt, engine);
            if (!split) continue;
            const auto children = split_children(box, *split);

            const auto cross_cut_profile =
                split->axis == SplitAxis::vertical
                ? build_side_profile_from_transformed(
                    compose_linear(
                        polynomial,
                        {split->coordinate, box.y_min},
                        {zero, one}))
                : build_side_profile_from_transformed(
                    compose_linear(
                        polynomial,
                        {box.x_min, split->coordinate},
                        {one, zero}));
            if (!cross_cut_profile) continue;

            const SideSlice cross_cut{
                cross_cut_profile,
                zero,
                split->axis == SplitAxis::vertical ? height : width
            };

            auto child_sides = split_child_sides(entry.sides, cross_cut, *split);

            std::size_t total = 0;
            bool good = true;

            for (std::size_t j = 0; j < children.size(); ++j) {
                const CountResult count =
                    count_zeros(polynomial, children[j], child_sides[j], options);
                if (!count || *count > entry.multiplicity) { good = false; break; }

                total += *count;
                refined[j] = {children[j], *count, entry.depth + 1, std::move(child_sides[j])};
            }

            if (good && total == entry.multiplicity) return refined;
        }

        return std::nullopt;
    }

    static ComplexPolynomial compose_linear(
        const ComplexPolynomial& polynomial,
        const complex_type& shift,
        const complex_type& direction) {
        // Здесь строится полином по t для композиции p(shift + direction * t).
        ComplexPolynomial result(1, complex_type(0, 0));
        ComplexPolynomial next;
        next.reserve(polynomial.size() + 1);

        for (auto it = polynomial.rbegin(); it != polynomial.rend(); ++it) {
            next.assign(result.size() + 1, complex_type(0, 0));
            for (std::size_t index = 0; index < result.size(); ++index) {
                next[index] += result[index] * shift;
                next[index + 1] += result[index] * direction;
            }
            next[0] += *it;
            result.swap(next);
        }

        normalize_complex_polynomial(result);
        return result;
    }

    static std::optional<SideProfile> build_sturm_sequence(
        RealPolynomial f1,
        RealPolynomial f2) {
        // Последовательность Штурма нужна для устойчивого подсчёта
        // смен знака на концах отрезка.
        if (!prepare_real_polynomial(f1) || !prepare_real_polynomial(f2)) return std::nullopt;

        const std::size_t max_length =
            f1.size() + f2.size() + detail::sturm_sequence_extra_terms;

        SideProfile sequence;
        sequence.reserve(max_length);
        sequence.push_back(std::move(f1));
        sequence.push_back(std::move(f2));

        while (sequence.back().size() > 1 && sequence.size() <= max_length) {
            auto remainder = negative_remainder(
                sequence[sequence.size() - 2],
                sequence.back());
            if (!remainder || !prepare_real_polynomial(*remainder)) return std::nullopt;
            sequence.push_back(std::move(*remainder));
        }

        return sequence.back().size() == 1
            ? std::optional<SideProfile>{std::move(sequence)}
            : std::nullopt;
    }

    static std::optional<RealPolynomial> negative_remainder(
        const RealPolynomial& dividend,
        RealPolynomial divisor) {
        // Для последовательности Штурма используется отрицательный остаток от деления.
        RealPolynomial remainder = dividend;
        trim_polynomial(remainder, constants::real_trim_safety_factor);
        trim_polynomial(divisor, constants::real_trim_safety_factor);

        if (divisor.empty()) return std::nullopt;

        const Real lead = divisor.back();
        const Real lead_tolerance = std::max(one, max_abs(divisor)) *
            constants::lead_tolerance_safety_factor * machine_epsilon();
        if (std::abs(lead) <= lead_tolerance) return std::nullopt;

        while (!remainder.empty() && remainder.size() >= divisor.size()) {
            const std::size_t shift = remainder.size() - divisor.size();
            const Real factor = remainder.back() / lead;
            for (std::size_t j = 0; j < divisor.size(); ++j)
                remainder[shift + j] -= factor * divisor[j];
            trim_polynomial(remainder, constants::real_trim_safety_factor);
        }

        for (Real& coefficient : remainder) coefficient = -coefficient;

        return remainder;
    }

    static std::optional<int> sign_variations(
        const std::vector<RealPolynomial>& sequence,
        Real x,
        const options_type& options) {
        int previous_sign = 0;
        int variations = 0;
        bool any_non_zero = false;

        for (const RealPolynomial& polynomial : sequence) {
            const auto [value, bound] = evaluate_and_absolute_bound(polynomial, x);
            const int sign = sign_of(value, bound, options);
            if (sign == 0) {
                // Почти нулевые значения пропускаются, чтобы шум не создавал ложные смены знака.
                continue;
            }

            any_non_zero = true;
            if (previous_sign != 0 && previous_sign != sign) ++variations;
            previous_sign = sign;
        }

        if (!any_non_zero) return std::nullopt;
        return variations;
    }

    static int sign_of(Real value, Real scale, const options_type& options) {
        const Real threshold = options.sign_safety * machine_epsilon() * std::max(one, scale);
        return std::abs(value) <= threshold ? 0 : (value > zero ? 1 : -1);
    }
};

template <typename Real, typename Coefficient>
std::vector<RootBox<Real>> localize_roots(
    const std::vector<Coefficient>& coefficients,
    Real epsilon,
    const Options<Real>& options = {}) {
    return Solver<Real>{}.localize(coefficients, epsilon, options);
}

template <typename Real, typename Coefficient>
std::vector<std::complex<Real>> find_roots(
    const std::vector<Coefficient>& coefficients,
    Real epsilon,
    const Options<Real>& options = {}) {
    return Solver<Real>{}.find_roots(coefficients, epsilon, options);
}

}  
