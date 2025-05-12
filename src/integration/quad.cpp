#include "quad.h"


Quadrature::Quadrature(int N): Basis(N)
{
    // Initialization logic here
}

double Quadrature::GaussLegendreQuad(std::function<double(double)> func, double a, double b)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    
    for (int i = 1; i < length-1; i++)
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
    
    for (int i = 1; i < length-1; i++)
    {
        double node = basis_.nodes(i);
        double weight = basis_.weights(i);
        // Transforming the node from [-1, 1] to [a, b]
        double transformedNode = 0.5 * (b - a) * node + 0.5 * (b + a);
        double L = LegendrePolynomialAndDerivative(j,node)[0];
        //std::cout<<"j "<<j<<" weight: "<<weight<<" node: "<<node<<" transformedNode "<<transformedNode <<" "<<" func(trans) "<<func(transformedNode)<<" L: "<<L<<std::endl;
        sol += weight * func(transformedNode)*L;
    }

    // Scale by the length of the interval
    sol *= 0.5 * (b - a);
    //std::cout<<"sol: "<<sol<<" j "<<j<<std::endl;
    return sol;
}
// int i is the position index j is the polynomial index
double Quadrature::IntFluxGaussLegendreQuad(std::function<double(double)> func,int i,int j ,double a, double b,const Array2D& u)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    
    for (int k = 1; k < length-1; k++)
    {
        double node = basis_.nodes(k);
        double weight = basis_.weights(k);
        double  L_prime = LegendrePolynomialAndDerivative(j,node)[1];
        //std::cout<<"Eval: "<<func(evaluation)<<"i: "<<i<<" L_prime "<<L_prime<<std::endl;
        sol += weight * u(i,k)*L_prime;
    }

    // Scale by the length of the interval
    //sol *= 0.5 * (b - a);
    return sol;
}

double Quadrature::IntFluxQ(std::function<double(double)> func, int i, int j, double a, double b,const Array2D& u)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    
    for (int k = 1; k < length-1; k++)
    {
        double node = basis_.nodes(k);
        double weight = basis_.weights(k);

        double  L_prime = LegendrePolynomialAndDerivative(j,node)[1];
        double intermediate_sol= GaussLegendreQuad(func,0.0,u(i,k));
    
        //std::cout<<"FUNCEval: "<<func(evaluation)<<" i: "<<i<<" L_prime "<<L_prime<<std::endl;
        //std::cout<<"Evaluation "<< evaluation<<" -g(u): "<<intermediate_sol<<" i: "<<i<<" j "<<j<<" L_prime "<<L_prime<<std::endl;
        sol += weight * intermediate_sol*L_prime;
    }
    //std::cout<<"sol: "<<sol<<" i: "<<i<<" j: "<<j<<std::endl;
    return sol;
}
double Quadrature::IntFluxU(std::function<double(double, double)> func, int i, int j, double a, double b, 
                            const Array2D& u, const Array2D& q,const Array2D& source) {
    int length = basis_.weights_.size()[0];
    double sol1 = 0.0,sol2 = 0.0;

    for (int k = 1; k < length - 1; k++) {
        double node = basis_.nodes(k);
        double weight = basis_.weights(k);

        // Transform the node from [-1, 1] to [a, b]
        std::array<double,2> L = LegendrePolynomialAndDerivative(j, node);
        double intermediate_sol = func(u(i,k), q(i,k));
        sol1 +=weight * (-source(i,k)*L[0]);
        //std::cout<<" source "<<source(i,k)<<" u "<<u(i,k)<<" q "<<q(i,k)<<" L[0] "<<L[0]<<" L[1] "<<L[1]<<" i "<<i<<" k "<<k<<std::endl;
        sol2 += weight * (intermediate_sol * L[1]);
    }
    //std::cout<<"sol: "<<sol<<" i: "<<i<<" j: "<<j<<std::endl;
    // Scale by the length of the interval
    double sol = 0.5 * (b - a)*sol1 +sol2;//0.5 * (b - a)*sol1 + 
    return sol;
}

double Quadrature::G(double u, double m)
{
    return sqrt(m)*2./(m+1.)*sqrt(pow(u,m-1.))*u;
}

double Quadrature::IntJ_0(std::function<double(double)> func, int j, double a, double b)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    
    for (int k = 1; k < length-1; k++)
    {
        double node = basis_.nodes(k);
        double weight = basis_.weights(k);
        double evaluation=0.0;
        double transformedNode = 0.5 * (b - a) * node + 0.5 * (b + a);
        double  L_prime = LegendrePolynomialAndDerivative(j,node)[1];
        double intermediate_sol= func(transformedNode);
    
        //std::cout<<"FUNCEval: "<<func(evaluation)<<" i: "<<i<<" L_prime "<<L_prime<<std::endl;
        // if(i==8 or i ==9)
        //     std::cout<<"Evaluation "<< evaluation<<" -g(u): "<<intermediate_sol<<" i: "<<i<<" j "<<j<<" L_prime "<<L_prime<<std::endl;
        sol += weight * intermediate_sol*L_prime;
    }
    //std::cout<<"sol: "<<sol<<" i: "<<i<<" j: "<<j<<std::endl;
    return sol;
}
