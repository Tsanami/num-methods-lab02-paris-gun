#include "rkf45.h"
#include "vector_operations.h"
#include <iostream>
#include <vector>
#include <cmath>

using Vec = std::vector<double>;

RKF45::RKF45(double tolerance) : tolerance_(tolerance) {}

void RKF45::integrate(const Vec& y0, double t0, double tf, double dt,
                      void (*f)(double, const Vec&, Vec&),
                      std::vector<Vec>& results) {
    Vec y = y0, y_new(y0.size()), yt(y0.size());
    double t = t0;
    results.clear();
    results.push_back(y);

    while (t < tf) {
        // 1) Вычисляем k1..k6
        Vec k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size()), k5(y.size()), k6(y.size());
        f(t, y, k1);
        for (size_t i = 0; i < y.size(); ++i) k1[i] *= dt;

        // k2
        for (size_t i = 0; i < y.size(); ++i) yt[i] = y[i] + k1[i] * 0.25;
        f(t + dt*0.25, yt, k2);
        for (auto &v:k2) v *= dt;

        // k3
        for (size_t i = 0; i < y.size(); ++i) yt[i] = y[i] + k2[i] * 0.25;
        f(t + dt*0.25, yt, k3);
        for (auto &v:k3) v *= dt;

        // k4
        for (size_t i = 0; i < y.size(); ++i) yt[i] = y[i] + k3[i] * 0.5;
        f(t + dt*0.5, yt, k4);
        for (auto &v:k4) v *= dt;

        // k5
        for (size_t i = 0; i < y.size(); ++i) yt[i] = y[i] + k4[i] * 0.75;
        f(t + dt*0.75, yt, k5);
        for (auto &v:k5) v *= dt;

        // k6
        for (size_t i = 0; i < y.size(); ++i) yt[i] = y[i] + k5[i];
        f(t + dt, yt, k6);
        for (auto &v:k6) v *= dt;

        // 2) Высчитываем шаг решения 4го порядка (для истории)
        for (size_t i = 0; i < y.size(); ++i) {
            y_new[i] = y[i] + ( k1[i]*2.0/9.0 + k3[i]*1.0/3.0 + k4[i]*4.0/9.0 + k6[i]*1.0/9.0 );
        }

        // 3) Оценка ошибки (разница 5го и 4го порядка)
        double error = 0.0;
        for (size_t i = 0; i < y.size(); ++i) {
            double y5 = y[i] + ( k1[i]*7.0/24.0 + k3[i]*1.0/4.0 + k4[i]*1.0/3.0 + k5[i]*1.0/8.0 );
            error += std::pow(y5 - y_new[i], 2);
        }
        error = std::sqrt(error);

        // 4) Адаптивный шаг: 

        // всегда продвигаем время, иначе застрянем
        t += dt;
        y = y_new;
        results.push_back(y);

        // проверяем приземление
        if (y[1] <= 0.01) {
            results.back()[1] = 0.0;
            break;
        }

        // корректируем dt
        if (error > tolerance_) {
            dt *= 0.5;
        } else if (error < tolerance_/4.0) {
            dt *= 1.5;
        }
        // ограничения
        dt = std::max(dt, 1e-4);
        dt = std::min(dt, 0.1);
    }
}
