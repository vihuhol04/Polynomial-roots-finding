# Создаем исполняемый файл
add_executable(Polynomial
    "Polynomial.cpp"
    Polynomial.h
    
    "generate_high_degree_polynomial.h")

# Настройки для Windows
if(WIN32)
    target_link_libraries(Polynomial PRIVATE psapi)
    target_compile_definitions(Polynomial PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
endif()

# Линковка зависимостей
target_link_libraries(Polynomial PRIVATE PrecisionSource)

# Дополнительные настройки
target_include_directories(Polynomial PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../precisionSource
    ${CMAKE_CURRENT_SOURCE_DIR}/../methods
)