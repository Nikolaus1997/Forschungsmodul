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

std::array<double, 2> NumericalFlux::computeNumFlux(double x_l, double x_r, double q_l, double q_r, double m, Flux flux_,const std::unique_ptr<Quadrature>& quad_)
{
    return porousMedia(x_l,x_r,q_l,q_r,m,flux_,quad_);
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

std::array<double,2> NumericalFlux::porousMedia(double u_l, double u_r, double q_l, double q_r, double m, Flux flux_,const std::unique_ptr<Quadrature>& quad_)
{   
    double g_plus =0.0,g_minus =0.0;
    double u_diff = u_r-u_l;
    double q_diff = q_r-q_l;
    double q_mean = 0.5*(q_r+q_l);
    double fraction_term =0.0;
    if(sqrt(pow(u_diff,2))<1e-18){
        fraction_term = -1.0*flux_.compute(u_l,0,m)[1];
    }else{
        g_plus= -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_r);
        g_minus= -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_l);
        fraction_term = (g_plus-g_minus)/u_diff;
    }
    double g_mean = 0.5*(g_plus+g_minus);
    double gamma = 0.5*fraction_term;
    return {-fraction_term*q_mean-gamma*q_diff,-g_mean+gamma*u_diff};
}
