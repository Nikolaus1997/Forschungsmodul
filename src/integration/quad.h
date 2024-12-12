#pragma once
#include <integration/basis.h>
#include <iostream>
#include <functional>
#include <storage/Vdm.h>


class Quadrature: public Basis
{
public:
    Quadrature(int N);
    double GaussLegendreQuad(std::function<double(double)> func, double a, double b);
    double IntGaussLegendreQuad(std::function<double(double)> func,int j ,double a, double b);
    double IntFluxGaussLegendreQuad(std::function<double(double)> func,int i,int j ,double a, double b,const Array2D& Vdm);
    double IntFluxQ(std::function<double(double)> func,int i,int j ,double a, double b,const Array2D& Vdm);
    double IntFluxU(std::function<double(double,double)> func,int i,int j ,double a, double b,const Array2D& Vdm, const Array2D& VdmQ);
};
