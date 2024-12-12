#pragma once
#include "flux.h"
#include "integration/quad.h"
#include <stdexcept>
#include <functional>
#include <array>
#include <memory>

class NumericalFlux
{
    enum class FunctionType { upwind, downwind, lax, enquist, porousMedia };

    // Constructor
    NumericalFlux() : selectedFunction(FunctionType::upwind) {}

    // Set the function type
    void setNumFluxFunction(FunctionType type);

    // Compute the flux based on the selected function
    double computeNumFlux(double x_l, double x_r,Flux flux_);
    std::array<double,2> computeNumFlux(double x_l, double x_r,double q_l, double q_r, double m,Flux flux_,const std::unique_ptr<Quadrature>& quad_);
    // Individual flux functions
    double upwind(double u_l, double u_r, Flux flux_);
    double downwind(double u_l, double u_r, Flux flux_);
    double lax(double u_l, double u_r, Flux flux_);
    double enquist(double u_l, double u_r, Flux flux_);
    std::array<double,2> porousMedia(double u_l, double u_r,double q_l, double q_r, double m,Flux flux_,const std::unique_ptr<Quadrature>& quad_);

private:
    FunctionType selectedFunction; // Stores the currently selected function
    friend class Computation;
};