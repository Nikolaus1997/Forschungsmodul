#pragma once
#include "storage/array2d.h"
#include <array>
#include <iostream>
class Vandermonde: public Array2D
{
    public:
        Vandermonde(std::array<int,2> size, int nNodes);
        void printValues();
        void LprintValues();
        void LprimePrintValues();
    protected:
        Array2D VdM_;
        Array2D VdM_t_;
        Array2D L_;
        Array2D L_prime_;
    friend class Computation;
};