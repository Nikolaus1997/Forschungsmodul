#pragma once
#include <tuple>
#include <storage/basis_storage.h>
#include "settings/settings.h"
#include <memory>


class Basis 
{
public:
    Basis(int N);
    std::tuple<double, double> LegendrePolynomialAndDerivative(int N, double x);
    std::tuple<std::vector<double>, std::vector<double>>  LegendreGaussNodesAndWeights(int N) const;

    Settings settings_;
    std::shared_ptr<BasisStorage> basis_;
};