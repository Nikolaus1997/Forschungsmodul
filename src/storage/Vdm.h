#pragma once
#include "storage/array2d.h"
#include "storage/array3d.h"
#include <array>
#include <iostream>
#include "array4d.h"
class Vandermonde: public Array3D
{
    public:
        Vandermonde(std::array<int,3> size, int nNodes);
        void printValues();
        void LprintValues();
        void LprimePrintValues();
        
        Array3D &VdM();
        Array3D &VdM1();
        Array3D &VdM2();
        Array3D &VdMt();      
        Array3D &VdMQ();    

        double VdM(int i, int j, int k) const;
        double &VdM(int i, int j, int k);
        

        double VdM1(int i, int j, int k) const;
        double &VdM1(int i, int j, int k);

        double VdM2(int i, int j, int k) const;
        double &VdM2(int i, int j, int k);

        double VdMQ(int i, int j, int k) const;
        double &VdMQ(int i, int j, int k);

        double VdMt(int i, int j, int k) const;
        double &VdMt(int i, int j, int k);

        double L(int i, int j) const;
        double &L(int i, int j);

        double L_prime(int i, int j) const;
        double &L_prime(int i, int j);

    protected:
        Array3D VdM_;
        Array3D VdM_t_;
        Array3D VdMQ_;
        Array3D VdM1_;
        Array4D VdMJ1_, VdMJ2_,VdMJ_,VdMJ_t_;
        Array3D VdM2_;
        Array2D L_;
        Array2D L_prime_;
    friend class Computation;
    friend class Grid;
};