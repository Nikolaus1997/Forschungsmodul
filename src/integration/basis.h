#pragma once
#include <tuple>
#include <storage/basis_storage.h>
#include "settings/settings.h"
#include <cmath>


class Basis 
{
public:
    Basis(int N);


    std::tuple<double, double> LegendrePolynomialAndDerivative(int N, double x);
    void LegendreGaussNodesAndWeights(int N);

    int N_;
    BasisStorage basis_;
};