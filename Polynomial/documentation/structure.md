Polynomial/
-├── CMakeLists.txt
-├── CMakePresets.json
-├── .gitignore
-├── Polynomial/
-│   └── generate_high_degree_polynomial.h (генератор)
-│   └── CMakeLists.txt 
-│   └── Polynomial.cpp	(общий вызов)
-│   └── Polynomial.h (подключение файлов и замер времени/памяти)
-└── precisionSource/ (библиотека повышенной точности)
-│   └── complexprecision.h
-│   └── fprecision.h
-│   └── fractionprecision.h
-│   └── intervalprecision.h
-│   └── iprecision.h
-│   └── mathprecision.h
-│   └── polyprecision.h
-│   └── precisioncore.cpp
-│   └── multipolyprecision.h (Аналог polyprecision для нескольких переменных)
-│   └── mathUtils.h (перегруженные математические функции)
-│   └── polynomialUtils.h (шаблонные методы для многочленов)
-└── methods/ (методы)
-│   └── CAD.h
-│   └── Con_Frac.h
-│   └── Con_Frac_Helper.h
-│   └── Durand-Kerner.h
-│   └── Graeffe_complex_roots.h
-│   └── Graeffe_real_roots.h
-│   └── Helper_for_all_methods.h
-│   └── Sagralov_complex_roots.h
-│   └── Sagralov_real_roots.h
-└── tests/ (тестирование методов)
│   └── generate_polynomial_test.cpp
│   └── graeffe_real_roots_test.cpp
└── documentation/ (документация)
│   └── Graeffe.md
│   └── README.md
│   └── structure.md
