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
    precisionSource/multipolyprecision.h
    precisionSource/mathUtils.h
    precisionSource/polynomialUtils.h

    "methods/Continued_Fractions/Con_Frac.h"
    "methods/Continued_Fractions/Con_Frac_Helper.h"
    "methods/Sagralov/Sagralov_real_roots.h"
    methods/Durand-Kerner.h 
    methods/Helper_for_all_methods.h
    
    "methods/Graeffe/Graeffe_complex_roots.h"
    "methods/Sagralov/Sagralov_complex_roots.h" 
    "methods/Sagralov/Sagralov_detail.h"
    "methods/Sagralov/Sagralov_ANewDsc.h"
    "methods/Sagralov/IntervalArithmetic.h" )

target_include_directories(PrecisionSource PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/precisionSource
    ${CMAKE_CURRENT_SOURCE_DIR}/methods
)

target_link_libraries(PrecisionSource PUBLIC Threads::Threads)

enable_testing()

# Тестирование Грефе для вещественных корней
add_executable(graeffe_real_roots_test
    tests/graeffe_real_roots_test.cpp )

target_link_libraries(graeffe_real_roots_test PRIVATE PrecisionSource)

target_include_directories(graeffe_real_roots_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/precisionSource
    ${CMAKE_CURRENT_SOURCE_DIR}/methods
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial
)

add_test(NAME graeffe_real_roots_test COMMAND PrecisionSource)

# Тестирование Грефе для комплексных корней
add_executable(graeffe_complex_roots_test
    tests/graeffe_complex_roots_test.cpp )

target_link_libraries(graeffe_complex_roots_test PRIVATE PrecisionSource)

target_include_directories(graeffe_complex_roots_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/precisionSource
    ${CMAKE_CURRENT_SOURCE_DIR}/methods
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial
)

add_test(NAME graeffe_complex_roots_test COMMAND PrecisionSource)

# Тестирование генератора полиномов
add_executable(generate_polynomial_test 
    tests/generate_polynomial_test.cpp )

target_link_libraries(generate_polynomial_test PRIVATE PrecisionSource)

target_include_directories(generate_polynomial_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/precisionSource
    ${CMAKE_CURRENT_SOURCE_DIR}/methods
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial
)

# Тестирование ANewDsc для вещественных корней
add_executable(ANewDsc_real_roots_test
    tests/ANewDsc_real_roots_test.cpp )

target_link_libraries(ANewDsc_real_roots_test PRIVATE PrecisionSource)

target_include_directories(ANewDsc_real_roots_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/precisionSource
    ${CMAKE_CURRENT_SOURCE_DIR}/methods
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial
)

add_test(NAME ANewDsc_real_roots_test COMMAND PrecisionSource)

# Тестирование метода Сгралова для вещественных корней
add_executable(sagralov_real_roots_test
    tests/sagralov_real_roots_test.cpp )

target_link_libraries(sagralov_real_roots_test PRIVATE PrecisionSource)

target_include_directories(sagralov_real_roots_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/precisionSource
    ${CMAKE_CURRENT_SOURCE_DIR}/methods
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial
)

add_test(NAME sagralov_real_roots_test COMMAND PrecisionSource)

# Тестирование метода Сгралова для комплексных корней
add_executable(sagralov_complex_roots_test
    tests/sagralov_complex_roots_test.cpp )

target_link_libraries(sagralov_complex_roots_test PRIVATE PrecisionSource)

target_include_directories(sagralov_complex_roots_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/precisionSource
    ${CMAKE_CURRENT_SOURCE_DIR}/methods
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial
)

add_test(NAME sagralov_complex_roots_test COMMAND PrecisionSource)


# Подключаем подпроект Polynomial
add_subdirectory(Polynomial)