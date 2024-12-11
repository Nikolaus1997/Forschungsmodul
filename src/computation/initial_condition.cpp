#include <computation/initial_condition.h>
#include "initial_condition.h"

void InitialCondition::setInitialCondType(InitialCondType type){
    selectedFunction = type;
}


double InitialCondition::computeInitialCondition(double x, double a, double b, double t, int m)
{
    if(m==0 and t ==0){
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
    }else{
        switch (selectedFunction)
        {
        case InitialCondType::Barenblatt:
        if(m==0)
            throw std::invalid_argument("The parameter m of the barenblatt function cant be 0");
            return barenBlatt(x,a,b,t,m);
        default:
            throw std::invalid_argument("Invalid initial condition type");
        }

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
//like in the papaer from ZhangWu
// never ever just use an int as multiplicator
double InitialCondition::barenBlatt(double x, double a, double b, double t, int m)
{
    double r = double(m);
    double k = 1/(r+1), factor = k*(r-1)/(2*r), timefactor = 1/(pow(t,2.0*k));
    double t_k = pow(t,-k), absx = pow(x,2);
    double brack = pow(std::max(1.0-absx*factor*timefactor,0.0),1/(r-1.0));
    //std::cout<<" k: "<<k<<" factor: "<<factor<<" absx: "<<absx<<" timefactor: "<<timefactor<<" brack: "<<brack<<"  "<<t_k <<std::endl;
    return t_k*brack;
}