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

double Quadrature::IntGaussLegendreQuad(std::function<double(double)> func,int j ,double a, double b)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    
    for (int i = 0; i < length; i++)
    {
        double node = basis_.nodes(i);
        double weight = basis_.weights(i);
        // Transforming the node from [-1, 1] to [a, b]
        double transformedNode = 0.5 * (b - a) * node + 0.5 * (b + a);
        double L = LegendrePolynomialAndDerivative(j,node)[0];
        //std::cout<<"j "<<j<<" weight"<<weight<<" node: "<<node<<" transformedNode "<<transformedNode <<" "<<" func(trans) "<<func(transformedNode)<<" L: "<<L<<std::endl;
        sol += weight * func(transformedNode)*L;
    }

    // Scale by the length of the interval
    sol *= 0.5 * (b - a);
    return sol;
}

double Quadrature::IntFluxGaussLegendreQuad(std::function<double(double)> func,int i,int j ,double a, double b, Array2D Vdm)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    
    for (int k = 0; k < length; k++)
    {
        double node = basis_.nodes(k);
        double weight = basis_.weights(k);
        double evaluation=0.0;
        // Transforming the node from [-1, 1] to [a, b]
        double  L_prime = LegendrePolynomialAndDerivative(j,node)[1];
        for(int p = 0;p<Vdm.size()[1];p++)
        {
            evaluation += Vdm(i,p)*LegendrePolynomialAndDerivative(p,node)[0];
        }
        //std::cout<<"Eval: "<<evaluation<<"i: "<<i<<std::endl;
        sol += weight * func(evaluation)*L_prime;
    }

    // Scale by the length of the interval
    sol *= 0.5 * (b - a);
    return sol;
}
