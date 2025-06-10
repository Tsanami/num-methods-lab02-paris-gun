#include "atmosphere.h"
#include <cmath>

static const double GAMMA = 1.4;

double temperature(double h) {
    return T0 + LAPSE * h;
}

double pressure(double h) {
    double T = temperature(h);
    return P0 * std::pow(T / T0, -G0 * MU / (R_CONST * LAPSE));
}

double density(double h) {
    double p = pressure(h);
    double T = temperature(h);
    return p * MU / (R_CONST * T);
}

double soundSpeed(double h) {
    double T = temperature(h);
    return std::sqrt(GAMMA * R_CONST * T / MU);
}
