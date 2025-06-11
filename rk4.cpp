#include "rk4.h"
#include "vector_operations.h"  // Подключаем перегрузки операторов
#include <iostream>
#include <vector>
#include <cmath>

using Vec = std::vector<double>;

void RK4::integrate(const Vec& y0, double t0, double tf, double dt, 
                    void (*f)(double, const Vec&, Vec&), 
                    std::vector<Vec>& results) {

    Vec y = y0;
    double t = t0;
    results.push_back(y);  // Добавляем начальное состояние

    while (t < tf) {
        double velocity = std::sqrt(y[2] * y[2] + y[3] * y[3]);
        if (velocity > 1000.0) {  // Ограничиваем скорость
            dt = std::min(dt, 0.05);  // Уменьшаем шаг
        }

        Vec k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size());

        f(t, y, k1);
        for (size_t i = 0; i < y.size(); i++) k1[i] *= dt;

        f(t + dt / 2.0, y + k1 / 2.0, k2);
        for (size_t i = 0; i < y.size(); i++) k2[i] *= dt;

        f(t + dt / 2.0, y + k2 / 2.0, k3);
        for (size_t i = 0; i < y.size(); i++) k3[i] *= dt;

        f(t + dt, y + k3, k4);
        for (size_t i = 0; i < y.size(); i++) k4[i] *= dt;

        for (size_t i = 0; i < y.size(); i++) {
            y[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;
        }

        // Проверка на физические ошибки (например, слишком большие скорости или NaN)
        if (std::isnan(y[0]) || std::isnan(y[1]) || std::isnan(y[2]) || std::isnan(y[3])) {
            std::cerr << "NaN detected, stopping simulation." << std::endl;
            break;  // Прекращаем вычисления при возникновении NaN
        }

        // Убедитесь, что высота не становится отрицательной
        if (y[1] <= 0.1) {
            std::cerr << "Projectile has hit the ground. Stopping simulation." << std::endl;
            break;
        }

        t += dt;
        if (y[0] > 0) {
            results.push_back(y);  // Сохраняем новый результат
        }
    }
}
