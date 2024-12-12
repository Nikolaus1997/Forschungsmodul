#pragma once
#include <stdexcept>
#include <functional>
#include <array>
#include <memory>

class Limiter
{
    enum class FunctionType { minmod, superbee, vanLeer, vanAlbada};

    // Constructor
    Limiter() : selectedFunction(FunctionType::minmod) {}

    // Set the limiter type
    void setLimiterFunction(FunctionType type);

    // Compute the limiter based on the selected function
    double computeLimiter(double a, double b, double c, double h = 0);

    // Individual limiter functions
    double minmod(double a, double b, double c, double h);
    double superbee(double a, double b, double c);
    double vanLeer(double a, double b, double c);
    double vanAlbada(double a, double b, double c);

private:
    FunctionType selectedFunction; // Stores the currently selected function
    friend class Computation;
};