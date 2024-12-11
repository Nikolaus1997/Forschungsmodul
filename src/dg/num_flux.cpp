#include "num_flux.h"



// Set the numerical flux function type
void NumericalFlux::setNumFluxFunction(FunctionType type)
{
    selectedFunction = type;
}

// Compute the numerical flux based on the selected function
double NumericalFlux::computeNumFlux(double u_l, double u_r,Flux flux_)
{
    switch (selectedFunction)
    {
    case FunctionType::upwind:
        return upwind(u_l, u_r,flux_);
    case FunctionType::downwind:
        return downwind(u_l, u_r,flux_);
    case FunctionType::lax:
        return lax(u_l, u_r,flux_);
    case FunctionType::enquist:
        return enquist(u_l, u_r,flux_);
    default:
        throw std::invalid_argument("Unknown numerical flux function type.");
    }
}

// Upwind flux function
double NumericalFlux::upwind(double u_l, double u_r,Flux flux_)
{
    // Example implementation: use left flux value
    return flux_.compute(u_l);
}

// Downwind flux function
double NumericalFlux::downwind(double u_l, double u_r, Flux flux_)
{
    // Example implementation: use right flux value
    return flux_.compute(u_r);
}

// Lax flux function
double NumericalFlux::lax(double u_l, double u_r, Flux flux_)
{
    // Example implementation: average of left and right
    return 0.5 * (flux_.compute(u_l) + flux_.compute(u_r))+0.5*(u_l-u_r);
}

// Enquist flux function
double NumericalFlux::enquist(double u_l, double u_r, Flux flux_)
{

    return (u_l > u_r) ? u_l : u_r;
}
