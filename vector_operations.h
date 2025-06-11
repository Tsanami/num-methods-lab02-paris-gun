#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H

#include <vector>

inline std::vector<double> operator+(
    const std::vector<double>& a,
    const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] + b[i];
    return result;
}

inline std::vector<double> operator*(
    const std::vector<double>& a,
    double scalar) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] * scalar;
    return result;
}

inline std::vector<double> operator/(
    const std::vector<double>& a,
    double scalar) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] / scalar;
    return result;
}

#endif  // VECTOR_OPERATIONS_H
