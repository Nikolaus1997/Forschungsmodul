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
    double IntJ_0(const Array2D& u,double m, int i, int j);
    double surfaceInt2D(int j, double a, double b, const Array2D& faces);
    double volumeInt2D(std::function<double(double)> func, int deg, int iCell, int jCell, double a, double b, double ay, double by, const Array2D &elemId, const Array3D &u);
    double volumeInt2D(std::function<double(double,double)> func,int deg ,int iCell, int jCell, double a, double b, double ay, double by,const Array2D& elemId ,const Array3D &u,const Array3D &q);
    double GaussLegendreQuad(std::function<double(double)> func, double a, double b);
    double IntGaussLegendreQuad(std::function<double(double,double)> func,int j ,double a, double b, double ay, double by);
    double IntFluxGaussLegendreQuad(std::function<double(double)> func,int i,int j ,double a, double b,const Array2D& u);
    double IntFaceFluxQ(std::function<double(double)> func,int j ,double a, double b,const Array2D& u);
    double IntFluxQ(std::function<double(double)> func, int deg, int iCell, int jCell, double a, double b, const Array3D &u, Array2D &elemId);
    double IntFluxU(std::function<double(double, double)> func, int i, int j, double a, double b, const Array2D &u, const Array2D &q, const Array2D &source);
};
