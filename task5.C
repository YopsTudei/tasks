void task5() {
    using namespace TMath;
    
    // Создание гистограммы для метода Неймана
    TH1D *hist = new TH1D("hist", "histogram", 50, 0, 4);

    // Генератор случайных чисел ROOT
    TRandom3 *rnd = new TRandom3();

    // Параметры метода Неймана
    const int num_points = 10000;  // Количество случайных точек
    double A = 1.0 / exp(1.0);  // Площадь прямоугольника 1 / e
    int count_f = 0;  // Количество точек под кривой f(x)
    int count_fg = 0; // Количество точек для отношения f(r3) / g(r3)
    
    // Функции f(x) = exp(-x) * x^3 и g(x) = 1/3 * x^2
    auto f = [](double x) -> double {
        return exp(-x) * pow(x, 3);
    };
    auto g = [](double x) -> double {
        return (1.0 / 3.0) * pow(x, 2);
    };

    // Основной цикл для метода Неймана
    for (int i = 0; i < num_points; i++) {
        // Генерация случайных значений r1 (x) и r2 (y)
        double r1 = rnd->Uniform(0, 1);  // r1 в диапазоне [0, 1]
        double r2 = rnd->Uniform(0, 1 / exp(1.0));  // r2 в диапазоне [0, 1/e]

        // Проверка, если точка лежит под графиком f(r1)
        if (r2 <= f(r1)) {
            count_f++;
        }

        // Генерация r3 = r1^(1/3) и r4
        double r3 = pow(r1, 1.0 / 3.0);  // r3 = r1^(1/3)
        double r4 = rnd->Uniform(0, 3 / exp(1.0));  // r4 в диапазоне [0, 3/e]

        // Проверка отношения f(r3) / g(r3) с r4
        if (r4 <= f(r3) / g(r3)) {
            count_fg++;
        }

        // Заполнение гистограммы для r1
        hist->Fill(Sqrt(-Log(1 - r1)));  // Примерное заполнение, можете изменить
    }

    // Оценка интегралов
    double integral_f = ((double)count_f / num_points);
    double integral_fg = ((double)count_fg / num_points);  // Коррекция на площадь 3/e

    // Вывод результатов
    std::cout << "Приближенное значение интеграла для f(x) методом Неймана: " << integral_f << std::endl;
    std::cout << "Приближенное значение интеграла для отношения f(r3)/g(r3): " << integral_fg << std::endl;

    // Отрисовка гистограммы
    hist->Draw();
}
