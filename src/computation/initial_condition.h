#pragma once

#include <stdexcept>
#include <functional>
#include <cmath>
#include <iostream>

class InitialCondition
{
public:
    enum class InitialCondType{ UnitStep, NegativeUnitStep, Sinus,Barenblatt};

    InitialCondition() : selectedFunction(InitialCondType::UnitStep) {}


    void setInitialCondType(InitialCondType type);


    double computeInitialCondition(double x, double a, double b, double time = 0 , int m = 0);

    // Destructor (optional, for cleanup if needed)
    double unitStep(double x, double a, double b);

    double negativeUnitStep(double x, double a, double b);

    double sinusFunc(double x, double a, double b);

    double barenBlatt(double x, double a, double b, double time, int m);

    // Other member functions and data members can be added here
private:
    InitialCondType selectedFunction; 
    friend class Computation;
};