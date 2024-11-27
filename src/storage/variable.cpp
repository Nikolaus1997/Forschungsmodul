#include "variable.h"


Variable::Variable(const std::array<int,1> nCells,const std::array<double, 1>  meshWidth):
Array1D(nCells),meshWidth_(meshWidth)
{
}

const std::array<double, 1>  Variable::meshWidth()
{
    return meshWidth_;
}


