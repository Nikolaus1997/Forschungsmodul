#include "limiter.h"
#include <algorithm>
#include <cmath>

void Limiter::setLimiterFunction(FunctionType type) {
    selectedFunction = type;
}

double Limiter::computeLimiter(double a, double b, double c,double h) {
    switch (selectedFunction) {
        case FunctionType::minmod:
            return minmod(a, b, c,h);
        case FunctionType::superbee:
            return superbee(a, b, c);
        case FunctionType::vanLeer:
            return vanLeer(a, b, c);
        case FunctionType::vanAlbada:
            return vanAlbada(a, b, c);
        default:
            throw std::runtime_error("Invalid limiter function selected");
    }
}

double Limiter::minmod(double a, double b, double c, double h) {
    // Implementation of the minmod limiter
    // Returns the minimum magnitude value with the same sign
    double sign_a = (a > 0) ? 1.0 : ((a < 0) ? -1.0 : 0.0);
    double sign_b = (b > 0) ? 1.0 : ((b < 0) ? -1.0 : 0.0);
    double sign_c = (c > 0) ? 1.0 : ((c < 0) ? -1.0 : 0.0);
    
    if(std::abs(a)>pow(h,2)){
        double sol = std::min({std::abs(a), std::abs(b), std::abs(c)});
        if(sol = std::abs(a))
            return sign_a*sol;
        if(sol = std::abs(b))
            return sign_b*sol;
        return sign_c*sol;
    }else if(std::abs(a)<=pow(h,2)){
        return a;
    }else{
        return 0.0;
    }
}

double Limiter::superbee(double a, double b, double c) {
    // Implementation of the superbee limiter
    // Returns the maximum of minmod(a, 2b) and minmod(2a, b)
    if (std::abs(b) < 1e-10) return 0.0;
    
    double r = a / b;
    return std::max(0.0, std::max(std::min(2.0 * r, 1.0), std::min(r, 2.0)));
}

double Limiter::vanLeer(double a, double b, double c) {
    // Implementation of the van Leer limiter
    // phi(r) = (r + |r|)/(1 + |r|)
    if (std::abs(b) < 1e-10) return 0.0;
    
    double r = a / b;
    return (r + std::abs(r)) / (1.0 + std::abs(r));
}

double Limiter::vanAlbada(double a, double b, double c) {
    // Implementation of the van Albada limiter
    // phi(r) = (r^2 + r)/(r^2 + 1)
    if (std::abs(b) < 1e-10) return 0.0;
    
    double r = a / b;
    double r2 = r * r;
    return (r2 + r) / (r2 + 1.0);
}