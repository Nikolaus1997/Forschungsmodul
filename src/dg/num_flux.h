#pragma once
#include "flux.h"
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
    double computeNumFlux(double x_l, double x_r,Flux flux_);

    // Individual flux functions
    double upwind(double u_l, double u_r, Flux flux_);
    double downwind(double u_l, double u_r, Flux flux_);
    double lax(double u_l, double u_r, Flux flux_);
    double enquist(double u_l, double u_r, Flux flux_);
    double porousMedia(double u_l, double u_r,Flux flux_);

private:
    FunctionType selectedFunction; // Stores the currently selected function
    friend class Computation;
};