#include "atmosphere.h"
#include <cmath>
#include <iostream>

static const double GAMMA = 1.4;

double temperature(double h) {
    if (h < 11000.0) {
        return 288.15 - 0.0065 * h;
    } else if (h < 20000.0) {
        return 216.65;
    } else if (h < 32000.0) {
        return 216.65 + 0.001 * (h - 20000.0);
    } else if (h < 47000.0) {
        return 228.65 + 0.0028 * (h - 32000.0);
    } else if (h < 51000.0) {
        return 270.65;
    } else if (h < 71000.0) {
        return 270.65 - 0.0028 * (h - 51000.0);
    } else if (h < 84852.0) {
        return 214.65 - 0.002 * (h - 71000.0);
    } else {
        return 186.946; // выше 86 км
    }
}


double pressure(double h) {
    if (h < 11000.0) {
        double T = 288.15 - 0.0065 * h;
        return 101325.0 * std::pow(T / 288.15, 5.255877);
    } else if (h < 20000.0) {
        return 22632.06 * std::exp(-G0 * (h - 11000.0) / (R_CONST * 216.65 / MU));
    } else if (h < 32000.0) {
        double T = 216.65 + 0.001 * (h - 20000.0);
        return 5474.889 * std::pow(T / 216.65, -34.16319);
    } else if (h < 47000.0) {
        double T = 228.65 + 0.0028 * (h - 32000.0);
        return 868.0187 * std::pow(T / 228.65, -12.2011);
    } else if (h < 51000.0) {
        return 110.9063 * std::exp(-G0 * (h - 47000.0) / (R_CONST * 270.65 / MU));
    } else if (h < 71000.0) {
        double T = 270.65 - 0.0028 * (h - 51000.0);
        return 66.93887 * std::pow(T / 270.65, 12.2011);
    } else if (h < 84852.0) {
        double T = 214.65 - 0.002 * (h - 71000.0);
        return 3.956420 * std::pow(T / 214.65, 17.0816);
    } else {
        return 0.0;
    }
}


double density(double h) {
    if (h < 0.01) h = 0.0;  // Защита от отрицательных значений
    if (h > 100000.0) return 0.0;  // Ограничение на максимальную высоту

    double p = pressure(h);
    double T = temperature(h);

    if (T <= 0.0 || std::isnan(T)) {
        std::cerr << "Warning: invalid temperature at h = " << h << ", T = " << T << "\n";
        return 0.0;  // безопасный возврат
    }

    return p * MU / (R_CONST * T);
}




double soundSpeed(double h) {
    double T = temperature(h);
    return std::sqrt(GAMMA * R_CONST * T / MU);
}


// Modify gravity for high altitudes
double gravity(double h) {
    const double G = 6.67430e-11;
    const double M = 5.972e24;
    const double R = 6371000.0;
    double denom = R + h;

    if (denom <= 0.1) denom = R;  // Ensure no division by zero
    if (h > 100000.0) return 9.81;  // Gravitational acceleration becomes constant above 100 km

    return G * M / (denom * denom);  // Standard gravity model
}



// === Сила сопротивления воздуха ===
double air_resistance(double vx, double vz, double rho, double Cd, double AREA) {
    double velocity = std::sqrt(vx * vx + vz * vz);
    if (velocity < 1.0) return 0.01;  // минимальная скорость, чтобы избежать деления на 0
    return 0.5 * Cd * rho * velocity * velocity * AREA;
}
