#ifndef GEODESY_H
#define GEODESY_H

// Эллипсоидальные параметры
static const double A = 6378136.3;            // большая полуось, м
static const double F = 1.0 / 298.257223563;  // сжатие
static const double B = A * (1 - F);
static const double E2 = 1 - (B * B) / (A * A);  // первая эксцентриситет^2

// Преобразование геодезических ↔ геоцентрических координат
void geodeticToGeocentric(double latDeg, double lonDeg, double h,
                           double &x, double &y, double &z);
void geocentricToGeodetic(double x, double y, double z,
                           double &latDeg, double &lonDeg, double &h);

#endif // GEODESY_H
