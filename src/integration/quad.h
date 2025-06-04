#pragma once
#include <integration/basis.h>
#include <iostream>
#include <functional>
#include <storage/Vdm.h>
#include <memory>


class Quadrature: public Basis
{
public:
    Quadrature(int N);
    double G(double u,double m);
    double IntJ_0(std::function<double(double)> func, int deg, double a, double b);
    double GaussLegendreQuad(std::function<double(double)> func, double a, double b);
    double IntGaussLegendreQuad(std::function<double(double)> func,int j ,double a, double b);
    double IntFluxGaussLegendreQuad(std::function<double(double)> func,int i,int j ,double a, double b,const Array2D& u);
    double IntFluxQ(std::function<double(double)> func,int i,int j ,double a, double b,const Array2D& u);
    double IntFluxU(std::function<double(double,double)> func,int i,int j ,double a, double b,const Array2D& u, const Array2D& q, const Array2D& source, bool isSource = false);
    double IntFluxU(std::function<double(double)> transportFunc, std::function<double(double, double)> func, int i, int j, double a, double b, const Array2D &u, const Array2D &q, const Array2D &source, bool isSource);
    double IntFluxU(int i, int j, double a, double b, const Array2D &u, const Array2D &source);
};
