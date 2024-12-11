#include "flux.h"

void Flux::setFluxFunction(FunctionType type) {
    selectedFunction = type;
}

double Flux::compute(double u) {

    switch (selectedFunction) {
        case FunctionType::Burgers:
            return burgers(u);
        case FunctionType::Linear:
            return linear(u);
        case FunctionType::BuckleyLeverett:
            return buckleyLeverett(u);
        default:
            throw std::invalid_argument("Invalid flux function type");
    }
}

std::array<double, 2> Flux::compute(double u, double q, double m)
{

      switch (selectedFunction) {
            case FunctionType::Barenblatt:
                return barenBlattFlux(u,q, m);
            default:
                throw std::invalid_argument("Invalid flux function type");
        }
    
}

double Flux::burgers(double x) {
    // Example implementation for Burgers' flux
    return 0.5 * x * x;
}

double Flux::linear(double x) {
    // Example implementation for Linear flux
    double a = 1.0;
    return a*x;
}

double Flux::buckleyLeverett(double x) {
    // Example implementation for Buckley-Leverett flux
    return (x*x) / (x + (1-x)*(1-x));
}

std::array<double, 2> Flux::barenBlattFlux(double u, double q, double m)
{
    return {(m*pow(u,m-1)-sqrt(m*pow(u,m-1))*q),(-sqrt(m*pow(u,m-1)))};
}
