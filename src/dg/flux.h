#pragma once

#include <stdexcept>
#include <functional>
#include <array>
#include <cmath>

class Flux {
public:
    enum class FunctionType { Burgers, Linear, BuckleyLeverett, Barenblatt};

    // Constructor
    Flux() : selectedFunction(FunctionType::Burgers) {}

    // Set the function type
    void setFluxFunction(FunctionType type);

    // Compute the flux based on the selected function
    double compute(double x);

    // Compute the flux based on the selected function
    std::array<double,2> compute(double x, double, double m );

    // Individual flux functions
    double burgers(double x);
    double linear(double x);
    double buckleyLeverett(double x);
    std::array<double,2> barenBlattFlux(double u, double q = 0, double m=0);

private:
    FunctionType selectedFunction; // Stores the currently selected function
    friend class OutputWriterParaview;
    friend class Computation;
};


