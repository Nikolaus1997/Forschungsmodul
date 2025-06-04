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

double Quadrature::IntGaussLegendreQuad(std::function<double(double,double)> func,int j ,double a, double b, double ay, double by)
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
        for(int k = 1; k < length-1; k++)
        {
            double nodeY = basis_.nodes(k);
            double weightY = basis_.weights(k);
            // Transforming the node from [-1, 1] to [ay, by]
            double transformedNodeY = 0.5 * (by - ay) * nodeY + 0.5 * (by + ay); 
            double Ly = LegendrePolynomialAndDerivative(j,nodeY)[0];          
        //std::cout<<"j "<<j<<" weight: "<<weight<<" node: "<<node<<" transformedNode "<<transformedNode <<" "<<" L: "<<L<<std::endl;
        sol += weight * func(transformedNode,transformedNodeY)*L*Ly*weightY;
        //std::cout<<"i "<<i<<" k "<<k<<" weight: "<<weight<<" weighty: "<<weightY<<" node: "<<transformedNode<<" func "<<func(transformedNode,transformedNodeY) <<" "<<" transformednodey "<<transformedNodeY <<" L: "<<L<<std::endl; 
        //std::cout<<"j "<<j<<" weight: "<<weight<<" node: "<<node<<" transformedNode "<<transformedNode <<" "<<" L: "<<L<<std::endl;
        }   
    }

    // Scale by the length of the interval
    sol *= 0.5 * (b - a)*0.5 * (by - ay);
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
        double evaluation=0.0;
        double  L_prime = LegendrePolynomialAndDerivative(j,node)[1];
        //std::cout<<"Eval: "<<func(evaluation)<<"i: "<<i<<" L_prime "<<L_prime<<std::endl;
        sol += weight * func(u(i,k))*L_prime;
    }

    // Scale by the length of the interval
    //sol *= 0.5 * (b - a);
    return sol;
}

double Quadrature::IntFaceFluxQ(std::function<double(double)> func, int deg, double a, double b,const Array2D& u)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    
    for(int j = 0; j < u.size()[1]; j++)
    {
        for (int k = 1; k < length-1; k++)
        {
            double node = basis_.nodes(k);
            double weight = basis_.weights(k);
            double evaluation=0.0;

            double  L = LegendrePolynomialAndDerivative(deg,node)[0];
            double intermediate_sol= GaussLegendreQuad(func,0.0,u(k-1,j));
        
            //std::cout<<"FUNCEval: "<<func(evaluation)<<" i: "<<i<<" L_prime "<<L_prime<<std::endl;
            //std::cout<<"Evaluation "<< evaluation<<" -g(u): "<<intermediate_sol<<" i: "<<i<<" j "<<j<<" L_prime "<<L_prime<<std::endl;
            if(j == 0 or j == 1)  
                sol += -1.*weight *L*intermediate_sol*pow(-1.,j);
            else if(j == 2 or j==3)
                sol += weight *L*intermediate_sol;
            //sol += weight * intermediate_sol*L_prime;
        }
    }
    sol *= 0.5 * (b - a)*0.5 * (b - a);
    //std::cout<<"sol: "<<sol<<" i: "<<i<<" j: "<<j<<std::endl;
    return sol;
}

double Quadrature::IntFluxQ(std::function<double(double)> func, int deg, int iCell, int jCell, double a, double b,const Array3D& u, Array2D& elemId)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    int idx = elemId(iCell,jCell);
    for(int j = 1; j < length-1; j++)
    {
        double node = basis_.nodes(j);
        double weight = basis_.weights(j);

        std::array<double,2>  L = LegendrePolynomialAndDerivative(deg,node);
        for (int k = 1; k < length-1; k++)
        {
            double nodey = basis_.nodes(k);
            double weighty = basis_.weights(k);

            std::array<double,2>  Ly = LegendrePolynomialAndDerivative(deg,nodey);
            double intermediate_soly= GaussLegendreQuad(func,0.0,u(j-1,k-1,idx));
        
            //std::cout<<"FUNCEval: "<<func(evaluation)<<" i: "<<i<<" L_prime "<<L_prime<<std::endl;
            //std::cout<<"Evaluation "<< evaluation<<" -g(u): "<<intermediate_sol<<" i: "<<i<<" j "<<j<<" L_prime "<<L_prime<<std::endl;
                sol += (L[0]*Ly[1]+L[1]*Ly[0])*weighty *intermediate_soly*weight;
            //sol += weight * intermediate_sol*L_prime;
        }
    }
    sol *= 0.5 * (b - a);
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
        double evaluationU = 0.0;
        double evaluationQ = 0.0;

        // Transform the node from [-1, 1] to [a, b]
        std::array<double,2> L = LegendrePolynomialAndDerivative(j, node);
        double intermediate_sol = func(u(i,k), q(i,k));
        sol1 += weight * (-source(i,k)*L[0]);
        sol2 += weight * (-intermediate_sol * L[1]);
    }
    //std::cout<<"sol: "<<sol<<" i: "<<i<<" j: "<<j<<std::endl;
    // Scale by the length of the interval
    double sol = 0.5 * (b - a)*sol1 + sol2;
    return sol;
}

double Quadrature::G(double u, double m)
{
    return sqrt(m)*2./(m+1.)*sqrt(pow(u,m-1.))*u;
}

double Quadrature::IntJ_0(const Array2D& u, double m, int i, int j)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    
    for (int k = 1; k < length-1; k++)
    {
        double node = basis_.nodes(k);
        double weight = basis_.weights(k);
        double evaluation=0.0;

        double  L_prime = LegendrePolynomialAndDerivative(j,node)[1];
        double intermediate_sol= pow(u(i,k),m);
    
        //std::cout<<"FUNCEval: "<<func(evaluation)<<" i: "<<i<<" L_prime "<<L_prime<<std::endl;
        // if(i==8 or i ==9)
        //     std::cout<<"Evaluation "<< evaluation<<" -g(u): "<<intermediate_sol<<" i: "<<i<<" j "<<j<<" L_prime "<<L_prime<<std::endl;
        sol += weight * intermediate_sol*L_prime;
    }
    //std::cout<<"sol: "<<sol<<" i: "<<i<<" j: "<<j<<std::endl;
    return sol;
}


double Quadrature::surfaceInt2D(int deg ,double a, double b, const Array2D& faceFlux)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    
    for (int i = 0; i < faceFlux.size()[1]; i++)
    {

        for(int k = 1; k < length-1; k++)
        {
            double node = basis_.nodes(k);
            double weight = basis_.weights(k);
            // Transforming the node from [-1, 1] to [a, b]
            double transformedNode = 0.5 * (b - a) * node + 0.5 * (b + a);
            double L = LegendrePolynomialAndDerivative(deg,node)[0];                  
        if(i == 0 or i == 1)  
            sol += -1.*weight *L*faceFlux(k-1,i)*pow(-1.,deg);
        else if(i == 2 or i==3)
            sol += weight *L*faceFlux(k-1,i);
        }   
    }
    sol *= 0.5 * (b - a);
    return sol;
}

double Quadrature::volumeInt2D(std::function<double(double)> func,int deg ,int iCell, int jCell, double a, double b, double ay, double by,const Array2D& elemId ,const Array3D &u)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    int iElem = elemId(iCell,jCell);
    for (int i = 1; i < length-1; i++)
    {
        double node = basis_.nodes(i);
        double weight = basis_.weights(i);
        double evaluation=0.0;
        std::array<double,2> L = LegendrePolynomialAndDerivative(deg,node);
        for(int k = 1; k < length-1; k++)
        {
            double nodeY = basis_.nodes(k);
            double weightY = basis_.weights(k);
            // Transforming the node from [-1, 1] to [ay, by]
            double transformedNodeY = 0.5 * (by - ay) * nodeY + 0.5 * (by + ay); 
            std::array<double,2> Ly = LegendrePolynomialAndDerivative(deg,nodeY);          

            sol += weight * func(u(i-1,k-1,iElem))*(L[1]*Ly[0]+L[0]*Ly[1])*weightY;
        }   
        //sol += weight * func(u(i,k))*L_prime;
    }

    // Scale by the length of the interval
    sol *= 0.5 * (b - a) * 0.5 * (by - ay);
    return sol;
}
double Quadrature::volumeInt2D(std::function<double(double,double)> func,int deg ,int iCell, int jCell, double a, double b, double ay, double by,const Array2D& elemId ,const Array3D &u, const Array3D &q)
{
    int length = basis_.weights_.size()[0];
    double sol = 0.0;
    int iElem = elemId(iCell,jCell);
    for (int i = 1; i < length-1; i++)
    {
        double node = basis_.nodes(i);
        double weight = basis_.weights(i);
        double evaluation=0.0;
        std::array<double,2> L = LegendrePolynomialAndDerivative(deg,node);
        for(int k = 1; k < length-1; k++)
        {
            double nodeY = basis_.nodes(k);
            double weightY = basis_.weights(k);
            // Transforming the node from [-1, 1] to [ay, by]
            double transformedNodeY = 0.5 * (by - ay) * nodeY + 0.5 * (by + ay); 
            std::array<double,2> Ly = LegendrePolynomialAndDerivative(deg,nodeY);          

            sol += weight * func(u(i-1,k-1,iElem),q(i-1,k-1,iElem))*(L[1]*Ly[0]+L[0]*Ly[1])*weightY;
        }   
        //sol += weight * func(u(i,k))*L_prime;
    }

    // Scale by the length of the interval
    sol *= 0.5 * (b - a) * 0.5 * (by - ay);
    return sol;
}