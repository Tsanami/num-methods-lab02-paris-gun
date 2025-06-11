#ifndef RK4_H
#define RK4_H

#include <vector>

class RK4 {
public:
    void integrate(const std::vector<double>& y0, double t0, double tf, double dt, 
                   void (*f)(double, const std::vector<double>&, std::vector<double>&), 
                   std::vector<std::vector<double>>& results);
};

#endif // RK4_H
