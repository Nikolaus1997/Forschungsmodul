#pragma once
#include "storage/array1d.h"
#include <array>

class Variable: public Array1D
{
public:
    //constructor
    Variable(const std::array<int,1> size,const std::array<double, 1>  meshWidth);

    const std::array<double, 1>  meshWidth();

private:
    const std::array<double, 1> meshWidth_;
};