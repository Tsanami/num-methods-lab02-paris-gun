#include <iostream>
#include "geodesy.h"
#include "atmosphere.h"

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
    return 0;
}
