#include <iostream>
#include <cmath>
#include <vector>
#include <limits>

#include "geodesy.h"
#include "atmosphere.h"
#include "rk4.h"
#include "rkf45.h"
#include "vector_operations.h"

// === Настройки снаряда ===
const double Cd = 0.47;   // Коэффициент сопротивления
const double AREA = 0.01;  // Площадь поперечного сечения (м^2)
const double m = 1.0;     // Масса снаряда (кг)


// === Уравнения движения ===
// === Уравнения движения ===
void equations_of_motion(double /*t*/, const std::vector<double>& y, std::vector<double>& dydt) {
    if (dydt.size() != 4) dydt.resize(4);

    double x  = y[0];
    double z  = y[1];
    double vx = y[2];
    double vz = y[3];

    // Приземление
    if (z <= 0) {
        dydt = {0,0,0,0};
        return;
    }

    double rho = density(z);
    double g   = gravity(z);

    // правильный D0
    double D0 = 0.5 * Cd * rho * AREA;
    double v  = std::hypot(vx, vz);

    // Сопротивление: dv/dt = - (D0/m) * v * v_vec
    double ax_drag = - (D0/m) * v * vx;
    double az_drag = - (D0/m) * v * vz;

    dydt[0] = vx;
    dydt[1] = vz;
    dydt[2] = ax_drag;
    dydt[3] = -g + az_drag;
}





int main() {
    // Test Geodesy
    double lat = 50.0, lon = 14.0, h = 300.0;
    double x, y, z;
    geodeticToGeocentric(lat, lon, h, x, y, z);
    std::cout << "Geo: x=" << x << ", y=" << y << ", z=" << z << "\n";
    double lat2, lon2, h2;
    geocentricToGeodetic(x, y, z, lat2, lon2, h2);
    std::cout << "Back: lat=" << lat2 << ", lon=" << lon2 << ", h=" << h2 << "\n";

    // Atmosphere Test
    for (double alt : {0.0, 1000.0, 5000.0, 10000.0, 20000.0}) {
        std::cout << "h=" << alt << " m: T=" << temperature(alt)
                  << " K, p=" << pressure(alt) << " Pa, rho="
                  << density(alt) << " kg/m3, a=" << soundSpeed(alt) << " m/s\n";
    }

    // Initial Conditions
    double speed = 500.0;
    double angle = M_PI/4.0;  // 45 градусов
    std::vector<double> initial_conditions = {
        0.0,                       // x
        h,                         // z
        speed * std::cos(angle),   // vx
        speed * std::sin(angle)    // vz
    };


    double t0 = 0.0, tf = 300.0, dt = 0.1;
    std::vector<std::vector<double>> results;

    // Integrator Choice
    RK4 integrator; 
    //RKF45 integrator(1e-5); 
    integrator.integrate(initial_conditions, t0, tf, dt, equations_of_motion, results);

    // Output results
    std::cout << "\nTrajectory (x, z) over time:\n";
    for (const auto& result : results) {     
        std::cout << "x: " << result[0] << ", z: " << result[1]
                  << ", vx: " << result[2] << ", vz: " << result[3] << "\n";
    }

    return 0;
}

