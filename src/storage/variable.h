#pragma once
#include "storage/array1d.h"
#include <array>

class Variable: public Array1D
{
public:
    //constructor
    Variable(int size,const double meshWidth);

    const double meshWidth();

private:
    const double meshWidth_;
};