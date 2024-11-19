#pragma once
#include "storage/array1d.h"
#include <array>

class Variable: public Array1D
{
public:
    //constructor
    Variable(int size, double origin, double meshWidth);
private:
    const double origin_;
    const double meshWidth_;
}