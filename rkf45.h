#ifndef RKF45_H
#define RKF45_H

#include <vector>

class RKF45 {
public:
    RKF45(double tolerance = 1e-6);  // Конструктор с параметром для точности
    void integrate(const std::vector<double>& y0, double t0, double tf, double dt, 
                   void (*f)(double, const std::vector<double>&, std::vector<double>&),
                   std::vector<std::vector<double>>& results);
private:
    double tolerance_;
};

#endif  // RKF45_H
