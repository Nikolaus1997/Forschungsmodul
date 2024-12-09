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
        
        Array2D &VdM();
        Array2D &VdMt();        

        double VdM(int i, int j) const;

        double &VdM(int i, int j);

        double VdMt(int i, int j) const;

        double &VdMt(int i, int j);

    protected:
        Array2D VdM_;
        Array2D VdM_t_;
        Array2D VdM1_;
        Array2D VdM2_;
        Array2D L_;
        Array2D L_prime_;
    friend class Computation;
};