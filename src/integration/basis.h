#pragma once
#include <array>
#include <storage/basis_storage.h>
#include "settings/settings.h"
#include <cmath>


class Basis 
{
public:
    Basis(int N);


    std::array<double,2> LegendrePolynomialAndDerivative(int N, double x);
    void LegendreGaussNodesAndWeights(int N);

    int N_;
    BasisStorage basis_;
};