#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include <vector>

// Атмосферные константы
static const double R_CONST = 8.31;           // Дж/(моль·К)
static const double MU = 0.0289644;            // кг/моль
static const double G0 = 9.80665;              // м/с^2
static const double T0 = 288.15;               // K
static const double P0 = 101325;               // Па
static const double LAPSE = -0.0065;           // K/м

// Функции атмосферы
double temperature(double h);
double pressure(double h);
double density(double h);
double soundSpeed(double h);

// Дополнительные функции для гравитации и сопротивления
double gravity(double h);                    // Функция для вычисления силы тяжести
double air_resistance(double vx, double vz, double rho, double Cd, double AREA); // Функция для вычисления сопротивления воздуха

#endif // ATMOSPHERE_H
