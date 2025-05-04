#pragma once

#include <stdexcept>
#include <functional>
#include <cmath>
#include <iostream>

class InitialCondition
{
public:
    enum class InitialCondType{ UnitStep, NegativeUnitStep, Sinus,Barenblatt, exponential, gaussian, linear, constant, none };

    InitialCondition() : selectedFunction(InitialCondType::UnitStep) {}


    void setInitialCondType(InitialCondType type);


    double computeInitialCondition(double x, double a, double b, double time = 0 , double m = 0.);

    double computeInitialCondition2D(double x, double y, double a, double b, double time = 0 , double m = 0.);

    // Destructor (optional, for cleanup if needed)
    double unitStep(double x, double a, double b);

    double negativeUnitStep(double x, double a, double b);

    double sinusFunc(double x, double a, double b);

    double barenBlatt(double x, double y, double a, double b, double time, double m);

    double unitStep(double x, double y, double a, double b);

    double negativeUnitStep(double x, double y, double a, double b);

    double sinusFunc(double x, double y, double a, double b);

    //double barenBlatt(double x, double y, double a, double b, double time, double m);

    // Other member functions and data members can be added here
private:
    InitialCondType selectedFunction; 
    friend class Computation;
};