# Polynomial Root Finding Project

## Описание проекта

Данный проект посвящен поиску всех корней полиномов высокой степени с использованием библиотеки повышенной точности Arbitrary Precision.
Цель проекта - исследовать различные методы на численную устойчивость и сделать выводы об области их применения

### Ключевые компоненты

- **methods/** — Набор различных методов для тестирования
- **precisionSource/** — Библиотека arbitrary precision с перегруженными математическими функциями
- **tests/** — Набор тестов для проверки корректности и точности методов
- **documentation/** — Техническая документация проекта

### Авторы
Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru:
	1. methods/
		1.1 Graeffe_real_roots.h, Graeffe_complex_roots.h
		1.2 Con_Frac.h, Con_Frac_Helper.h
		1.3 Sagralov_real_roots.h, Sagralov_complex_roots.h
		1.4 Helper_for_all_methods.h
	2. /Polynomial
		2.1 generate_high_degree_polynomial.h
	3. /precisionSource
		3.1 polynomialUtils.h
		3.2 mathUtils.h
	4. /tests
		4.1 generate_polynomial_test.cpp
		4.2 graeffe_real_roots_test.cpp
		4.3 graeffe_complex_roots_test.cpp
	5. documentation
		5.1 Graeffe.md
		5.2 structure.md
