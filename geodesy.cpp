#include "geodesy.h"
#include <cmath>

void geodeticToGeocentric(double latDeg, double lonDeg, double h,
                           double &x, double &y, double &z) {
    double phi = latDeg * M_PI / 180.0;
    double lam = lonDeg * M_PI / 180.0;
    double sinPhi = std::sin(phi);
    double cosPhi = std::cos(phi);
    double cosLam = std::cos(lam);
    double sinLam = std::sin(lam);
    double N = A / std::sqrt(1 - E2 * sinPhi * sinPhi);
    x = (N + h) * cosPhi * cosLam;
    y = (N + h) * cosPhi * sinLam;
    z = (N * (1 - E2) + h) * sinPhi;
}

void geocentricToGeodetic(double x, double y, double z,
                           double &latDeg, double &lonDeg, double &h) {
    double lon = std::atan2(y, x);
    double p = std::sqrt(x * x + y * y);
    double phi = std::atan2(z, p * (1 - E2));
    double phiPrev;
    const double tol = 1e-12;
    do {
        phiPrev = phi;
        double sinPhi = std::sin(phi);
        double N = A / std::sqrt(1 - E2 * sinPhi * sinPhi);
        h = p / std::cos(phi) - N;
        phi = std::atan2(z, p * (1 - E2 * (N / (N + h))));
    } while (std::fabs(phi - phiPrev) > tol);
    latDeg = phi * 180.0 / M_PI;
    lonDeg = lon * 180.0 / M_PI;
}
