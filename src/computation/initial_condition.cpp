#include <initial_condition.h>

InitialCondition::InitialCondition(std::array<double, 2> physicalSize, int nCells):
nCells_(nCells),u_(nCells_,(physicalSize[1]-physicalSize[0])/nCells)
{
}

Variable InitialCondition::unitStep(Variable x,double a,double b)
{
    for(int i = 0; i< x.size();i++)
    {
        if(x(i)>=a and x(i)<=b)
        {
            u_(i) = 1.0;
        }else 
        {
            u_(i) = 0.0;
        }
    }

    return u_;
}

Variable InitialCondition::negativeUnitStep(Variable x,double a,double b)
{
    for(int i = 0; i< x.size();i++)
    {
        if(x(i)>=a and x(i)<=b)
        {
            u_(i) = -1.0;
        }else 
        {
            u_(i) = 0.0;
        }
    }

    return u_;
}
