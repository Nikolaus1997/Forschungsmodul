#pragma once
#include <storage/array2d.h>
#include <cmath>

class Projection
{
    public:
        Projection();
        void project(Array2D &u, const Array2D &x);
        std::pair<Array2D,Array2D> LR(Array2D &A);
        void makeProjection(Array2D &u, const Array2D x,int i, int order);
        Array2D MatMul(const Array2D &A, const Array2D &B);
        Array2D MakeTransPosed(const Array2D &A);
        Array2D MakeMonomBasis(const Array2D &u, int i, int order);
        Array2D backwardSubstitution(const Array2D &A, const Array2D &b);
        Array2D forwardSubstitution(const Array2D &A, const Array2D &b);

    private:

};
