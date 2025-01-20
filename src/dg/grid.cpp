#include <dg/grid.h>
#include "grid.h"

Grid::Grid(std::array<int, 1>  nCells, std::array<double, 1>  meshWidth, int numberNodes):
nCells_(nCells), meshWidth_(meshWidth), u_      ({int(nCells_[0]*(numberNodes))},   meshWidth_),
                                        ut_      ({int(nCells_[0]*(numberNodes))},   meshWidth_),
                                        u1_     ({int(nCells_[0]*(numberNodes))},   meshWidth_),
                                        u2_     ({int(nCells_[0]*(numberNodes))},   meshWidth_),
                                        x_      ({int(nCells_[0]*(numberNodes))},   meshWidth_),
                                        faces_  ({nCells_[0]+1},   meshWidth_),
                                        rhs_    (nCells_,   meshWidth_),
                                        l2_error_({1},   meshWidth_),
                                        linf_error_({1},   meshWidth_)
{
}
const std::array<double,1> Grid::meshWidth() const
{
return meshWidth_;
}
const std::array<int,1> Grid::nCells() const
{
return nCells_;
}
Variable &Grid::u()
{
    return u_;
}

double Grid::u(int i) const
{
    return u_(i);
}

double &Grid::u(int i)
{
    return u_(i);
}

Variable &Grid::ut()
{
    return ut_;
}

double Grid::ut(int i) const
{
    return ut_(i);
}

double &Grid::ut(int i)
{
    return ut_(i);
}

const Variable &Grid::u1() const
{
    return u1_;
}

double Grid::u1(int i) const
{
    return u1_(i);
}

double &Grid::u1(int i)
{
    return u1_(i);
}

const Variable &Grid::u2() const
{
    return u2_;
}

double Grid::u2(int i) const
{
    return u2_(i);
}

double &Grid::u2(int i)
{
    return u2_(i);
}


const Variable &Grid::x() const
{
    return x_;
}

double Grid::x(int i) const
{
    return x_(i);
}

double &Grid::x(int i)
{
    return x_(i);
}

const Variable &Grid::faces() const
{

    return faces_;
}

double Grid::faces(int i) const
{
    return faces_(i);
}

double &Grid::faces(int i)
{

    return faces_(i);
}

const Variable &Grid::l2_error() const
{
    return l2_error_;
}

double Grid::l2_error(int i) const
{
    return l2_error_(i);
}

double &Grid::l2_error(int i)
{
    return l2_error_(i);
}

const Variable &Grid::linf_error() const
{
    return linf_error_;
}

double Grid::linf_error(int i) const
{
    return linf_error_(i);
}

double &Grid::linf_error(int i)
{
    return linf_error_(i);
}


const Variable &Grid::rhs() const
{

    return rhs_;
}

double Grid::rhs(int i) const
{
    return rhs_(i);
}

double &Grid::rhs(int i)
{
    return rhs_(i);
}

double Grid::dx() const
{
    return meshWidth_[0];
}

void Grid::fillSolution(Variable& x,const std::shared_ptr<Vandermonde> VdM)
{
    for(int i=0;i<VdM->size()[0];i++){
        for(int j=1;j<VdM->size()[1];j++){
                x(i*(VdM->size()[1]-1)+j-1) =0.0;
                for(int p=0;p<VdM->size()[1];p++){
                    x(i*(VdM->size()[1]-1)+j-1) +=VdM->VdM(i,p)*VdM->L(j,p);
                }
            }
    }
}


void Grid::fillDerivative(Variable& x,const std::shared_ptr<Vandermonde> VdM)
{
    for(int i=0;i<VdM->size()[0];i++){
        for(int j=1;j<VdM->size()[1];j++){
                x(i*(VdM->size()[1]-1)+j-1) =0.0;
                for(int p=0;p<VdM->size()[1];p++){
                    x(i*(VdM->size()[1]-1)+j-1) +=VdM->VdMt(i,p)*VdM->L(j,p);
                }
            }
    }
}
