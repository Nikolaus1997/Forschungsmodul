#include "num_flux.h"

// Constructor implementation
NumericalFlux::NumericalFlux()
    : selectedFunction(FunctionType::upwind) // Default to upwind
{
}

// Set the numerical flux function type
void NumericalFlux::setNumFluxFunction(FunctionType type)
{
    selectedFunction = type;
}

// Compute the numerical flux based on the selected function
double NumericalFlux::computeNumFlux(double x_l, double x_r)
{
    switch (selectedFunction)
    {
    case FunctionType::upwind:
        return upwind(x_l, x_r);
    case FunctionType::downwind:
        return downwind(x_l, x_r);
    case FunctionType::lax:
        return lax(x_l, x_r);
    case FunctionType::enquist:
        return enquist(x_l, x_r);
    default:
        throw std::invalid_argument("Unknown numerical flux function type.");
    }
}

// Upwind flux function
double NumericalFlux::upwind(double x_l, double x_r)
{
    // Example implementation: use left flux value
    return x_l;
}

// Downwind flux function
double NumericalFlux::downwind(double x_l, double x_r)
{
    // Example implementation: use right flux value
    return x_r;
}

// Lax flux function
double NumericalFlux::lax(double x_l, double x_r)
{
    // Example implementation: average of left and right
    return 0.5 * (x_l + x_r);
}

// Enquist flux function
double NumericalFlux::enquist(double x_l, double x_r)
{
    // Example implementation: max of left and right
    return (x_l > x_r) ? x_l : x_r;
}
