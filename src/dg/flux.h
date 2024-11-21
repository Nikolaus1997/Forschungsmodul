#pragma once

#include <stdexcept>
#include <functional>

class Flux {
public:
    enum class FunctionType { Burgers, Linear, BuckleyLeverett };

    // Constructor
    Flux() : selectedFunction(FunctionType::Burgers) {}

    // Set the function type
    void setFluxFunction(FunctionType type);

    // Compute the flux based on the selected function
    double compute(double x);

    // Individual flux functions
    double burgers(double x);
    double linear(double x);
    double buckleyLeverett(double x);

private:
    FunctionType selectedFunction; // Stores the currently selected function
};


