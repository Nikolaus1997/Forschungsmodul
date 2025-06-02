#pragma once
#include "flux.h"
#include "integration/quad.h"
#include <stdexcept>
#include <functional>
#include <array>
#include <memory>
#include <cmath>
#include <iostream>

class NumericalFlux
{
    enum class FunctionType { upwind, downwind, lax, enquist, porousMedia, central, laxWen };

    // Constructor
    NumericalFlux() : selectedFunction(FunctionType::upwind) {}

    // Set the function type
    void setNumFluxFunction(FunctionType type);

    // Compute the flux based on the selected function
    double computeNumFlux(double x_l, double x_r,Flux flux_, double dt = 1.0, double meshWidth = 1.0);
    std::array<double,2> computeNumFlux(bool minus,double x_l, double x_r,double q_l, double q_r, double m,Flux flux_,const std::unique_ptr<Quadrature>& quad_, double u_mean, double u_mean_plus);
    // Individual flux functions
    double upwind(double u_l, double u_r, Flux flux_);
    double downwind(double u_l, double u_r, Flux flux_);
    double lax(double u_l, double u_r, Flux flux_, double dt, double meshWidth);
    double laxWen(double u_l, double u_r, Flux flux_, double dt, double meshWidth);
    double enquist(double u_l, double u_r, Flux flux_);
    bool almostEqual(double a, double b);
    std::array<double,2> porousMediaPlus(double u_l, double u_r,double q_l, double q_r, double m,Flux flux_,const std::unique_ptr<Quadrature>& quad_,double u_mean=.0, double u_mean_plus=0.0);
    std::array<double,2> porousMediaMinus(double u_l, double u_r,double q_l, double q_r, double m,Flux flux_,const std::unique_ptr<Quadrature>& quad_,double u_mean=.0, double u_mean_plus=0.0);      

private:
    FunctionType selectedFunction; // Stores the currently selected function
    friend class Computation;
};