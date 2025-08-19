#include "Vdm.h"

Vandermonde::Vandermonde(std::array<int, 3> size, int nNodes):Array3D(size),
                    VdM_(size),
                    VdMQ_(size),    
                    VdM1_(size),
                    VdM2_(size),
                    VdMJ_({size[0],size[1],size[2],2}),
                    VdMJ1_({size[0],size[1],size[2],2}),
                    VdMJ2_({size[0],size[1],size[2],2}),
                    VdMJ_t_({size[0],size[1],size[2],2}),
                    L_({nNodes,size[2]}),
                    L_prime_({nNodes,size[2]}), 
                    VdM_t_(size)
{
}

void Vandermonde::printValues()
{
    std::cout << "=== Array3D Contents ===\n";
    
    for (int k = 0; k < VdM_.size()[2]; ++k) { // over elements
        std::cout << "PPN:  " << k << ":\n";
        for (int i = 0; i < VdM_.size()[0]; ++i) { // over local nodes in i
            for (int j = 0; j < VdM_.size()[1]; ++j) { // over local nodes in j
                std::cout << "inCell  (" <<i<<", "<< j << "): ";
                std::cout << VdM_(i, j, k) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "----------------------\n";
    }
}

void Vandermonde::LprintValues()
{


        for (int i = 0; i < L_.size()[1]; i++)
            {
                std::cout<<"PP_N "<<i<<" :";
                for(int j= 0; j<L_.size()[0];j++){
                    std::cout<<L(j,i)<<" ";
                }
                std::cout<<";"<<std::endl;
            }
    
}

void Vandermonde::LprimePrintValues()
{

    for (int i = 0; i < L_prime_.size()[1]; i++)
        {
            std::cout<<"PP_N "<<i<<" :";
            for(int j= 0; j<L_prime_.size()[0];j++){
                std::cout<<L_prime(j,i)<<" ";
            }
            std::cout<<";"<<std::endl;
        }
    
}

Array3D &Vandermonde::VdM() 
{
    return VdM_;
}

Array3D &Vandermonde::VdM1() 
{
    return VdM1_;
}

Array3D &Vandermonde::VdM2() 
{
    return VdM2_;
}

Array3D &Vandermonde::VdMt()
{
    return VdM_t_;
}

Array3D &Vandermonde::VdMQ()
{
    return VdMQ_;
}

double Vandermonde::VdM(int i, int j, int k) const
{
    return VdM_(i,j,k);
}

double &Vandermonde::VdM(int i, int j, int k)
{
    return VdM_(i,j,k);
}

double Vandermonde::VdM1(int i, int j, int k) const
{
    return VdM1_(i,j,k);
}

double &Vandermonde::VdM1(int i, int j, int k)
{
    return VdM1_(i,j,k);
}

double Vandermonde::VdM2(int i, int j, int k) const
{
    return VdM2_(i,j,k);
}

double &Vandermonde::VdM2(int i, int j, int k)
{
    return VdM2_(i,j,k);
}

double Vandermonde::VdMQ(int i, int j, int k) const
{
    return VdMQ_(i,j,k);
}

double &Vandermonde::VdMQ(int i, int j, int k)
{
    return VdMQ_(i,j,k);
}

double Vandermonde::VdMt(int i, int j, int k) const
{
    return VdM_t_(i,j,k);
}

double &Vandermonde::VdMt(int i, int j, int k)
{
    return VdM_t_(i,j,k);
}

double Vandermonde::L(int i, int j) const
{
    return L_(i,j);
}

double &Vandermonde::L(int i, int j)
{
    // TODO: insert return statement here
    return L_(i,j);
}



double Vandermonde::L_prime(int i, int j) const
{
    return L_prime_(i,j);
}

double &Vandermonde::L_prime(int i, int j)
{
    // TODO: insert return statement here
    return L_prime_(i,j);
}