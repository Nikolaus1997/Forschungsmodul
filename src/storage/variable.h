#pragma once
#include "storage/array1d.h"
#include <array>

class Variable: public Array1D
{
public:
    //constructor
    Variable(int size, double meshWidth);
private:
    const double meshWidth_;
};