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
	- methods/
		- Graeffe_real_roots.h, Graeffe_complex_roots.h
		- Con_Frac.h, Con_Frac_Helper.h
		- Sagralov_real_roots.h, Sagralov_complex_roots.h
		- Helper_for_all_methods.h
	- /Polynomial
		- generate_high_degree_polynomial.h
	- /precisionSource
		- polynomialUtils.h
		- mathUtils.h
	- /tests
		- generate_polynomial_test.cpp
		- graeffe_real_roots_test.cpp
		- graeffe_complex_roots_test.cpp
	- documentation
		- Graeffe.md
		- structure.md
