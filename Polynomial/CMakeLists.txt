cmake_minimum_required(VERSION 3.10)
project(PolynomialProject LANGUAGES CXX)

# Настройки компилятора
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Для Windows
if(WIN32)
    add_definitions(-DNOMINMAX -D_CRT_SECURE_NO_WARNINGS)
endif()

find_package(Threads REQUIRED)

# Создаем библиотеку PrecisionSource
add_library(PrecisionSource STATIC
    precisionSource/precisioncore.cpp
    precisionSource/mathUtils.h
    methods/Helper_for_all_methods.h
    precisionSource/polynomialUtils.h
    methods/Con_Frac.h
    methods/Con_Frac_Helper.h
    methods/Sagralov_real_roots.h
    methods/Durand-Kerner.h
     
    precisionSource/multipolyprecision.h
    methods/CAD.h
    methods/Graeffe_complex_roots.h
    methods/Sagralov_complex_roots.h)

target_include_directories(PrecisionSource PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/precisionSource
    ${CMAKE_CURRENT_SOURCE_DIR}/methods
)

target_link_libraries(PrecisionSource PUBLIC Threads::Threads)

enable_testing()

# Тестирование Грефе для вещественных корней
add_executable(graeffe_real_roots_test
    tests/graeffe_real_roots_test.cpp)

target_link_libraries(graeffe_real_roots_test PRIVATE PrecisionSource)

target_include_directories(graeffe_real_roots_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/precisionSource
    ${CMAKE_CURRENT_SOURCE_DIR}/methods
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial
)

add_test(NAME graeffe_real_roots_test COMMAND PrecisionSource)

# Тестирование Грефе для комплексных корней
add_executable(graeffe_complex_roots_test
    tests/graeffe_complex_roots_test.cpp)

target_link_libraries(graeffe_complex_roots_test PRIVATE PrecisionSource)

target_include_directories(graeffe_complex_roots_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/precisionSource
    ${CMAKE_CURRENT_SOURCE_DIR}/methods
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial
)

add_test(NAME graeffe_complex_roots_test COMMAND PrecisionSource)

# Тестирование генератора полиномов
add_executable(generate_polynomial_test
    tests/generate_polynomial_test.cpp
)

target_link_libraries(generate_polynomial_test PRIVATE PrecisionSource)

target_include_directories(generate_polynomial_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/precisionSource
    ${CMAKE_CURRENT_SOURCE_DIR}/methods
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial
)

add_test(NAME generate_polynomial_test COMMAND PrecisionSource)


# Подключаем подпроект Polynomial
add_subdirectory(Polynomial)