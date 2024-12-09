#include <computation/initial_condition.h>
#include "initial_condition.h"

void InitialCondition::setInitialCondType(InitialCondType type){
    selectedFunction = type;
}


double InitialCondition::computeInitialCondition(double x, double a, double b)
{
    switch (selectedFunction){
        case InitialCondType::UnitStep:
            return unitStep( x, a, b);
        case InitialCondType::NegativeUnitStep:
            return negativeUnitStep(x,a,b);
        case InitialCondType::Sinus:
            return sinusFunc(x,a,b);    
        default:
            throw std::invalid_argument("Invalid initial condition type");
    }
}

double InitialCondition::unitStep(double x, double a, double b)
{
    if (x>=a and x<=b)
    {
        return 1.0;
    }else
    {
        return 0.0;
    }
    
}

double InitialCondition::negativeUnitStep(double x, double a, double b)
{
        if (x>=a and x<=b)
    {
        return -1.0;
    }else
    {
        return 0.0;
    }
}

double InitialCondition::sinusFunc(double x, double a, double b)
{
    return sin(2*x);
}
