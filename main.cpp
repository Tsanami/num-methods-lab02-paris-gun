#include <iostream>
#include "geodesy.h"
#include "atmosphere.h"
#include "odesolver.hpp"
#include "odevents.hpp"
#include <iostream>
#include <vector>
int main() {
    // Тест геодезии
    double lat = 50.0, lon = 14.0, h = 300.0;
    double x, y, z;
    geodeticToGeocentric(lat, lon, h, x, y, z);
    std::cout << "Geo: x=" << x << ", y=" << y << ", z=" << z << "\n";
    double lat2, lon2, h2;
    geocentricToGeodetic(x, y, z, lat2, lon2, h2);
    std::cout << "Back: lat=" << lat2 << ", lon=" << lon2 << ", h=" << h2 << "\n";

    // Тест атмосферы
    for (double alt : {0.0, 1000.0, 5000.0, 10000.0, 20000.0}) {
        std::cout << "h=" << alt << " m: T=" << temperature(alt)
                  << " K, p=" << pressure(alt) << " Pa, rho="
                  << density(alt) << " kg/m3, a=" << soundSpeed(alt) << " m/s\n";
    }
    // ODU Part
    try {
        // Начальные условия: [x, y, z, vx, vy, vz]
        std::vector<double> initialState = {0.0, 0.0, 10000.0, 200.0, 0.0, 0.0}; // 10 км высота, 200 м/с по x
        ODESolver solver(1e-3, 1e-3, 1.0); // Точность 1e-3, min шаг 1e-3, max шаг 1.0
        // Решение методом Рунге-Кутты 4-го порядка
        auto solutionRK4 = solver.solve(initialState, 0.0, 10.0, 0.01, true);
        // Решение методом Эйлера
        auto solutionEuler = solver.solve(initialState, 0.0, 10.0, 0.01, false);
        // Вывод результатов
        std::cout << "Time\tRK4 x\tRK4 z\tEuler x\tEuler z" << std::endl;
        double t = 0.0;
        double dt = 0.5; // Выводим каждые 0.5 секунды
        int step = static_cast<int>(dt / 0.01); // Шаг вывода
        for (size_t i = 0; i < solutionRK4.size(); i += step) {
            if (i >= solutionEuler.size()) break;

            std::cout << t << "\t"
                      << solutionRK4[i][0] << "\t" << solutionRK4[i][2] << "\t"
                      << solutionEuler[i][0] << "\t" << solutionEuler[i][2] << std::endl;
            t += dt;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    // Odevents part
    try {
        // 1. Создание интегратора с обработкой событий
        ODEIntegratorWithEvents integrator;

        // 2. Настройка обработчиков событий (новый синтаксис)
        integrator.setEventCallback([](const IntegrationEvent& event) {
            switch (event.type) {
                case IntegrationEvent::SURFACE_HIT:
                    std::cout << "\n=== Event got surface ==="
                              << "\ntime: " << event.time << " sec"
                              << "\ncoord: (" << event.state[0] << ", "
                              << event.state[1] << ", " << event.state[2] << ")"
                              << "\nspped: (" << event.state[3] << ", "
                              << event.state[4] << ", " << event.state[5] << ")\n";
                    break;

                case IntegrationEvent::MAX_ALTITUDE:
                    std::cout << "\n=== Event max Altitude ==="
                              << "\ntime: " << event.time << " sec"
                              << "\nheight: " << event.state[2] << " m\n";
                    break;
            }
        });

        // 3. Регистрация событий (явный вызов методов)
        integrator.addSurfaceEvent(0.0);    // Поверхность на z=0
        integrator.addMaxAltitudeEvent();   // Событие максимальной высоты

        // 4. Начальные условия [x, y, z, vx, vy, vz]
        std::vector<double> initialState = {0.0, 0.0, 10000.0, 200.0, 0.0, 0.0};

        // 5. Запуск интегрирования с обработкой событий
        std::cout << "Begin of integrating...\n";
        auto fullTrajectory = integrator.integrateWithEvents(
                initialState,   // Начальное состояние
                0.0,           // Начальное время
                100.0,         // Конечное время
                0.01,          // Шаг
                true           // Использовать RK4
        );

        // 6. Анализ результатов
        const auto& trajectorySegments = integrator.getTrajectoriesBetweenEvents();
        std::cout << "\nRes"
                  << "\nGot segment: " << trajectorySegments.size()
                  << "\nNum points: " << fullTrajectory.size()
                  << "\nLast point: t=" << fullTrajectory.back()[0]
                  << ", z=" << fullTrajectory.back()[2] << "\n";

    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
