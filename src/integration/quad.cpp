#include "quad.h"


Quadrature::Quadrature(int N): Basis(N)
{
    // Initialization logic here
}

double Quadrature::GaussLegendreQuad(std::function<double(double)> func, double a, double b)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    
    for (int i = 0; i < length; i++)
    {
        double node = basis_.nodes(i);
        double weight = basis_.weights(i);
        // Transforming the node from [-1, 1] to [a, b]
        double transformedNode = 0.5 * (b - a) * node + 0.5 * (b + a);

        sol += weight * func(transformedNode);
    }

    // Scale by the length of the interval
    sol *= 0.5 * (b - a);
    return sol;
}