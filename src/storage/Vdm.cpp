#include "Vdm.h"

Vandermonde::Vandermonde(std::array<int, 2> size, int nNodes):Array2D(size),
VdM_(size),VdMQ_(size),VdM1_(size),VdM2_(size),L_({nNodes,size[1]}),L_prime_({nNodes,size[1]}),VdM_t_(size)
{
}

void Vandermonde::printValues()
{
    for (int i = 0; i < VdM_.size()[1]; i++)
    {
        std::cout<<"PP_N "<<i<<" :";
        for(int j= 0; j<VdM_.size()[0];j++){
            std::cout<<VdM_(j,i)<<" ";
        }
        std::cout<<";"<<std::endl;
    }
    
}

void Vandermonde::LprintValues()
{
    for (int i = 0; i < L_.size()[1]; i++)
    {
        std::cout<<"PP_N "<<i<<" :";
        for(int j= 0; j<L_.size()[0];j++){
            std::cout<<L_(j,i)<<" ";
        }
        std::cout<<";"<<std::endl;
    }
}

void Vandermonde::LprimePrintValues()
{
    for (int i = 0; i < VdM_.size()[1]; i++)
    {
        std::cout<<"PP_N "<<i<<" :";
        for(int j= 0; j<VdM_.size()[0];j++){
            std::cout<<L_prime_(j,i)<<" ";
        }
        std::cout<<";"<<std::endl;
    }
}

Array2D &Vandermonde::VdM() 
{
    return VdM_;
}

Array2D &Vandermonde::VdMt()
{
    return VdM_t_;
}

Array2D &Vandermonde::VdMQ()
{
    return VdMQ_;
}

double Vandermonde::VdM(int i, int j) const
{
    return VdM_(i,j);
}

double &Vandermonde::VdM(int i, int j)
{
    return VdM_(i,j);
}

double Vandermonde::VdMQ(int i, int j) const
{
    return VdMQ_(i,j);
}

double &Vandermonde::VdMQ(int i, int j)
{
    return VdMQ_(i,j);
}

double Vandermonde::VdMt(int i, int j) const
{
    return VdM_t_(i,j);
}

double &Vandermonde::VdMt(int i, int j)
{
    return VdM_t_(i,j);
}
