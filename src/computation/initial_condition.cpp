#include <computation/initial_condition.h>
#include "initial_condition.h"

void InitialCondition::setInitialCondType(InitialCondType type){
    selectedFunction = type;
}


double InitialCondition::computeInitialCondition(double x, double a, double b, double t, double m)
{
    if(selectedFunction != InitialCondType::Barenblatt){
        switch (selectedFunction){
            case InitialCondType::UnitStep:
                return unitStep( x, a, b);
            case InitialCondType::NegativeUnitStep:
                return negativeUnitStep(x,a,b);
            case InitialCondType::Sinus:
                return sinusFunc(x,a,b); 
            case InitialCondType::exponential:
                return exp(-x*x);
            case InitialCondType::gaussian:
                return exp(-x*x*0.5)*1./(sqrt(2*M_1_PI)); 
            case InitialCondType::divorce:
                if(x>-M_PI and x<-M_PI/6. or x<M_PI and x>M_PI/6.)
                    return std::fabs(std::sin(x));
                else if (x>-M_PI/6. and x<M_PI/6.)
                    return 0.5;
                return 0.0;  
            default:
                throw std::invalid_argument("Invalid initial condition type");
        }
    }else{
        switch (selectedFunction)
        {
        case InitialCondType::Barenblatt:
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
    return sin(x)+2.;
}
//like in the papaer from ZhangWu
// never ever just use an int as multiplicator
double InitialCondition::barenBlatt(double x, double a, double b, double t, double m)
{
    double r = m;
    double k = 1./(r+1.), factor = k*(r-1.)/(2.*r), timefactor = 1./(pow(t,2.0*k));
    double t_k = pow(t,-k), absx = pow(fabs(x),2.);
    double brack = pow(std::max(1.0-absx*factor*timefactor,0.0),1./(r-1.0));
    //std::cout<<" k: "<<k<<" factor: "<<factor<<" absx: "<<absx<<" timefactor: "<<timefactor<<" brack: "<<brack<<"  "<<t_k <<std::endl;
    if(m>1)
        return t_k*brack;
    else
        return 1.0;
}