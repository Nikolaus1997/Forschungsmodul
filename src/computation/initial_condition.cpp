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
            default:
                throw std::invalid_argument("Invalid initial condition type");
        }
    }else{
        switch (selectedFunction)
        {
        case InitialCondType::Barenblatt:
            return barenBlatt(x,x,a,b,t,m);
        default:
            throw std::invalid_argument("Invalid initial condition type");
        }

    }

}
std::array<double,2> InitialCondition::computeInitialConditionGradient2D(double x, double y, double a, double b)
{
    switch (selectedFunction){
        case InitialCondType::Sinus:
            // Example: sinusFunc(x,y) = sin(pi*x/a)*sin(pi*y/b)
            {
                double dx = (M_PI / a) * cos(M_PI * x / a) * sin(M_PI * y / b);
                double dy = (M_PI / b) * sin(M_PI * x / a) * cos(M_PI * y / b);
                return {dx, dy};
            }

        case InitialCondType::exponential:
            // exp(-(x^2 + y^2))
            {
                double value = exp(-(x*x + y*y));
                double dx = -2.0 * x * value;
                double dy = -2.0 * y * value;
                return {dx, dy};
            }

        case InitialCondType::gaussian:
            // 2D normalized Gaussian: exp(-(x^2 + y^2)/2)/(2*pi)
            {
                double value = exp(-0.5*(x*x + y*y)) / (2.0*M_PI);
                double dx = -x * value;
                double dy = -y * value;
                return {dx, dy};
            }


        default:
            throw std::invalid_argument("Invalid initial condition type");
    }
}
double InitialCondition::computeInitialCondition2D(double x, double y, double a, double b, double t, double m)
{
    //if(selectedFunction != InitialCondType::Barenblatt){
        switch (selectedFunction){
            case InitialCondType::UnitStep:
                return unitStep( x,y, a, b);
            case InitialCondType::NegativeUnitStep:
                return negativeUnitStep(x,y,a,b);
            case InitialCondType::Sinus:
                return sinusFunc(x,y,a,b); 
            case InitialCondType::exponential:
                return exp(-(x*x + y*y));
            case InitialCondType::gaussian:
                return exp(-(x*x + y*y) * 0.5) * 1./(sqrt(2*M_1_PI)); 
            case InitialCondType::Barenblatt:
                return barenBlatt(x,y,a,b,t,m);
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
    return sin(x);
}
double InitialCondition::unitStep(double x,double y, double a, double b)
{
    if ((x>=a and x<=b) and (y>=a and y<=b))
    {
        return 1.0;
    }else
    {
        return 0.0;
    }
    
}

double InitialCondition::negativeUnitStep(double x,double y, double a, double b)
{
        if ((x>=a and x<=b) and (y>=a and y<=b))
    {
        return -1.0;
    }else
    {
        return 0.0;
    }
}

double InitialCondition::sinusFunc(double x,double y, double a, double b)
{
    return sin(x+y);
}
//like in the papaer from ZhangWu
// never ever just use an int as multiplicator
double InitialCondition::barenBlatt(double x, double y, double a, double b, double t, double m)
{
    double r = m;
    double k = 1./(r+1.), factor = k*(r-1.)/(2.*r), timefactor = 1./(pow(t,2.0*k));
    double t_k = pow(t,-k), absx = pow(sqrt(pow(x,2.)+pow(y,2.)),2.);
    double brack = pow(std::max(1.0-absx*factor*timefactor,0.0),1./(r-1.0));
    //std::cout<<" k: "<<k<<" factor: "<<factor<<" absx: "<<absx<<" timefactor: "<<timefactor<<" brack: "<<brack<<"  "<<t_k <<std::endl;
    if(m>1)
        return t_k*brack;
    else
        return 1.0;
}