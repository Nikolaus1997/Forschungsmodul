#pragma once
#include <stdexcept>
#include <functional>

class NumericalFlux
{
    enum class FunctionType { upwind, downwind, lax, enquist };

    // Constructor
    NumericalFlux() : selectedFunction(FunctionType::upwind) {}

    // Set the function type
    void setNumFluxFunction(FunctionType type);

    // Compute the flux based on the selected function
    double computeNumFlux(double x_l, double x_r);

    // Individual flux functions
    double upwind(double x_l, double x_r);
    double downwind(double x_l, double x_r);
    double lax(double x_l, double x_r);
    double enquist(double x_l, double x_r);

private:
    FunctionType selectedFunction; // Stores the currently selected function
    friend class Computation;
};