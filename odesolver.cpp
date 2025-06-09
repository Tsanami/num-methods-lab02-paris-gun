#include "odesolver.hpp"
#include <iostream>
#include <cmath>
#include "atmosphere.h"




// Реализация методов класса ODESolver
ODESolver::ODESolver(double eps, double minSt, double maxSt)
        : epsilon(eps), minStep(minSt), maxStep(maxSt) {}

std::vector<double> ODESolver::derivative(const std::vector<double>& state, double t) {
    if (state.size() != 6) {
        throw std::invalid_argument("State vector must have 6 elements");
    }

    double x = state[0];
    double y = state[1];
    double z = state[2];
    double vx = state[3];
    double vy = state[4];
    double vz = state[5];

    // Высота (предполагаем, что z - высота)
    double h = z;

    // Атмосферные параметры
    double rho = density(h);
    double T = temperature(h);
    double a = soundSpeed(h);

    // Пример: сила сопротивления (пропорциональна квадрату скорости)
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    double dragCoeff = 0.1; // Коэффициент сопротивления
    double dragForce = 0.5 * rho * v * v * dragCoeff;

    // Ускорение (включая гравитацию и сопротивление)
    std::vector<double> deriv(6);
    deriv[0] = vx; // dx/dt = vx
    deriv[1] = vy; // dy/dt = vy
    deriv[2] = vz; // dz/dt = vz

    deriv[3] = -dragForce * vx / v; // dvx/dt = ax (сопротивление по x)
    deriv[4] = -dragForce * vy / v; // dvy/dt = ay (сопротивление по y)
    deriv[5] = -G0 - dragForce * vz / v; // dvz/dt = az (гравитация + сопротивление по z)
    return deriv;
}

std::vector<double> ODESolver::euler(const std::vector<double>& state, double t, double dt) {
    std::vector<double> deriv = derivative(state, t);
    std::vector<double> newState(6);

    for (int i = 0; i < 6; ++i) {
        newState[i] = state[i] + deriv[i] * dt;
    }

    return newState;
}

std::vector<double> ODESolver::rungeKutta4(const std::vector<double>& state, double t, double dt) {
    std::vector<double> k1 = derivative(state, t);

    std::vector<double> state2(6);
    for (int i = 0; i < 6; ++i) {
        state2[i] = state[i] + k1[i] * dt / 2.0;
    }
    std::vector<double> k2 = derivative(state2, t + dt/2.0);

    std::vector<double> state3(6);
    for (int i = 0; i < 6; ++i) {
        state3[i] = state[i] + k2[i] * dt / 2.0;
    }
    std::vector<double> k3 = derivative(state3, t + dt/2.0);

    std::vector<double> state4(6);
    for (int i = 0; i < 6; ++i) {
        state4[i] = state[i] + k3[i] * dt;
    }
    std::vector<double> k4 = derivative(state4, t + dt);

    std::vector<double> newState(6);
    for (int i = 0; i < 6; ++i) {
        newState[i] = state[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) * dt / 6.0;
    }

    return newState;
}

std::vector<std::vector<double>> ODESolver::solve(std::vector<double> initialState, double t0, double tEnd,
                                                  double initialStep, bool useRK4) {
    std::vector<std::vector<double>> solution;
    solution.push_back(initialState);

    double t = t0;
    double dt = initialStep;

    while (t < tEnd) {
        // Корректируем шаг, чтобы не выйти за tEnd
        if (t + dt > tEnd) {
            dt = tEnd - t;
        }

        std::vector<double> newState;
        try {
            if (useRK4) {
                newState = rungeKutta4(solution.back(), t, dt);
            } else {
                newState = euler(solution.back(), t, dt);
            }
        } catch (const std::exception& e) {
            std::cerr << "Error during integration at t = " << t << ": " << e.what() << std::endl;
            throw;
        }

        // Проверка сходимости (по изменению z-координаты)
        double deltaZ = std::abs(newState[2] - solution.back()[2]);

        if (deltaZ < epsilon || dt <= minStep) {
            // Шаг приемлем или достигнут минимальный шаг
            solution.push_back(newState);
            t += dt;

            // Увеличиваем шаг (но не больше максимального)
            if (deltaZ < epsilon / 2.0) {
                dt = std::min(dt * 1.5, maxStep);
            }
        } else {
            // Уменьшаем шаг (но не меньше минимального)
            dt = std::max(dt / 2.0, minStep);
            if (dt <= minStep) {
                std::cerr << "Warning: Minimum step size reached at t = " << t << std::endl;
            }
            continue; // Повторяем итерацию с новым шагом
        }
    }

    return solution;
}