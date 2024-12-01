#include <integration/basis.h>
#include "basis.h"


Basis::Basis(int N):
N_(N),basis_(N_)
{
}

std::array<double,2> Basis::LegendrePolynomialAndDerivative(int N, double x)

{
        double L = 0.0, L_prime = 0.0;
        if (N == 0) {
            L = 1.0;
            L_prime = 0.0;
        } else if (N == 1) {
            L = x;
            L_prime = 1.0;
        } else {
            double L_2 = 1.0,L_2_prime = 0.0;  // P_0(x)
            double L_1 = x,L_1_prime = 1.0;    // P_1(x)
            for (int i = 2; i <= N; i++) {
                L = ((2 * i - 1) * x * L_1 - (i - 1) * L_2) / i;
                L_prime = L_2_prime+(2*i-1)*L_1;
                L_2 = L_1;
                L_1 = L;
                L_2_prime = L_1_prime;
                L_1_prime = L_prime;
            }
        }
        return {L, L_prime};
}
 void   Basis::LegendreGaussNodesAndWeights(int N)
{
        int iterations = 10;
        double dx = 0.0, tol = 1e-8;


        if (N==0)
        {
            basis_.weights(1) = 2.0;
            basis_.nodes(1) = 0;
        }else if (N==1)
        {
            basis_.weights(1) = 1.0;
            basis_.weights(2) = 1.0;
            basis_.nodes(1) = -std::sqrt(1/3);
            basis_.nodes(2) = - basis_.nodes(1);
        }else if(N!=0)
        {
            std::vector<double> nodes(N), weights(N);
            // Compute the roots (nodes) and weights
            for (int j = 0; j < floor((N + 1) / 2); j++) {
                // Initial guess using Chebyshev-Gauss nodes
                double x = -std::cos((2 * j + 1.0) / (2 * (N-1) + 2) * M_PI);

                // Refine root using Newton's method
                for (int k = 0; k < iterations; k++) {
                    auto [L, L_prime] = LegendrePolynomialAndDerivative(N, x);
                    dx = -L / L_prime;
                    x += dx;
                    if (std::abs(dx) <= tol * std::abs(x)) break;
                }

                // Store the root and its symmetric counterpart
                nodes[j] = x;
                nodes[N - j - 1] = -x;

                // Compute the weight for this node
                auto [L, L_prime] = LegendrePolynomialAndDerivative(N, x);
                double w = 2.0 / ((1 - x * x) * L_prime * L_prime);
                weights[j] = w;
                weights[N - j - 1] = w;
            }
        for (int i = 0; i <N; i++)
        {
            basis_.nodes(i+1) = nodes[i];
            basis_.weights(i+1) = weights[i];
        }
        }


if (N!=0)
{
basis_.weights(0) =  2.0/(N*(N+1));
basis_.weights(basis_.nodes_.size()[0]-1) = basis_.weights(0);
}
basis_.nodes(0) = -1.0;
basis_.nodes(basis_.nodes_.size()[0]-1)= 1.0;

}