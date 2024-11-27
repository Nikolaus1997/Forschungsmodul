#pragma once

#include <stdexcept>
#include <functional>

class InitialCondition
{
public:
    enum class InitialCondType{ UnitStep, NegativeUnitStep};

    InitialCondition() : selectedFunction(InitialCondType::UnitStep) {}


    void setInitialCondType(InitialCondType type);

    double computeInitialCondition(double x, double a, double b);

    // Destructor (optional, for cleanup if needed)
    double unitStep(double x, double a, double b);

    double negativeUnitStep(double x, double a, double b);

    // Other member functions and data members can be added here
private:
    InitialCondType selectedFunction; 
};