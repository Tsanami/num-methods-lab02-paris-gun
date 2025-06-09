#ifndef ODESOLVER_HPP
#define ODESOLVER_HPP

#include <vector>
#include <stdexcept>




// Класс для решения системы ОДУ
class ODESolver {
private:
    double epsilon; // Точность для критерия сходимости
    double minStep; // Минимальный шаг
    double maxStep; // Максимальный шаг

public:
    // Конструктор с параметрами точности и шага
    ODESolver(double eps = 1e-6, double minSt = 1e-6, double maxSt = 1.0);

    // Правая часть системы ОДУ
    std::vector<double> derivative(const std::vector<double>& state, double t);

    // Метод Эйлера
    std::vector<double> euler(const std::vector<double>& state, double t, double dt);

    // Метод Рунге-Кутты 4-го порядка
    std::vector<double> rungeKutta4(const std::vector<double>& state, double t, double dt);

    // Адаптивный шаг интегрирования
    std::vector<std::vector<double>> solve(std::vector<double> initialState, double t0, double tEnd,
                                           double initialStep = 0.01, bool useRK4 = true);
};

#endif // ODESOLVER_HPP