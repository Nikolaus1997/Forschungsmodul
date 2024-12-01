#include "flux.h"

void Flux::setFluxFunction(FunctionType type) {
    selectedFunction = type;
}

double Flux::compute(double x) {
    switch (selectedFunction) {
        case FunctionType::Burgers:
            return burgers(x);
        case FunctionType::Linear:
            return linear(x);
        case FunctionType::BuckleyLeverett:
            return buckleyLeverett(x);
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
