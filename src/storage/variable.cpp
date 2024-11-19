#include "variable.h"

Variable::Variable(int size, double origin, double meshWidth):
Array1D(size),origin_(origin),meshWidth_(meshWidth)
{
}