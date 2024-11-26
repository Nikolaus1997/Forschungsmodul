#include "variable.h"


Variable::Variable(int nCells, double meshWidt):
Array1D(nCells),meshWidth_(meshWidth)
{
}

const double Variable::meshWidth()
{
    return meshWidth_;
}


