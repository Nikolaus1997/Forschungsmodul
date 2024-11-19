#include <integration/basis.h>



Basis::Basis(int N)
{
}

std::tuple<double, double> Basis::LegendrePolynomialAndDerivative(int N, double x)
{
    double       L = 0.0;
    double L_prime = 0.0;
    if (N==0)
    {
        L       = 1.0;
        L_prime = 0.0;
        return {L, L_prime};
    }
    else if (N==1)
    {
        L       = x;
        L_prime = 1;
        return {L, L_prime};
    }
    else
    {   
        auto [L_1,L_1_prime] = LegendrePolynomialAndDerivative(N-1,x);
        auto [L_2,L_2_prime] = LegendrePolynomialAndDerivative(N-2,x);

        L = (2*N-1.0)/N*x*L_1-(N-1)/N*L_2;
        L_prime = L_2_prime+(2*N-1)*L_1;
    }
    return {L, L_prime};
}

std::tuple<std::vector<double>, std::vector<double>> Basis::LegendreGaussNodesAndWeights(int N) const
{
    int iterations = 4;
    if (N==0)
    {
        //TODO:: umschreiben sodass arrays übergeben werden somit kein shared pointer bevor er gebraucht wird.
    }
    
    return 0.0;
}