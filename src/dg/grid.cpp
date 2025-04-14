#include <dg/grid.h>
#include "grid.h"

Grid::Grid(std::array<int, 1>  nCells, std::array<double, 1>  meshWidth, int numberNodes):
nCells_(nCells), meshWidth_(meshWidth), u_      ({int(nCells_[0]),(numberNodes+2)}),
                                        j_      ({int(nCells_[0]),(numberNodes+2)}),
                                        j_1_    ({int(nCells_[0]),(numberNodes+2)}),
                                        j_2_    ({int(nCells_[0]),(numberNodes+2)}),
                                        true_solution_({int(nCells_[0]),(numberNodes+2)}),
                                        q_      ({int(nCells_[0]),(numberNodes+2)}),
                                        ut_      ({int(nCells_[0]),(numberNodes+2)}),
                                        u1_     ({int(nCells_[0]),(numberNodes+2)}),
                                        u2_    ({int(nCells_[0]),(numberNodes+2)}),
                                        solution_({int(nCells_[0]*(numberNodes))},   meshWidth_),
                                        solutionJ_({int(nCells_[0]*(numberNodes))},   meshWidth_),
                                        derivative_({int(nCells_[0]*(numberNodes))},   meshWidth_),
                                        x_      ({int(nCells_[0]),(numberNodes+2)}),
                                        x_analyze_    ({300,1}),
                                        u_analyze_    ({300,1}),
                                        u_analyze_true_    ({300,1}),
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
Array2D &Grid::u()
{
    return u_;
}


Array2D &Grid::ut()
{
    return ut_;
}


const Array2D &Grid::u1() const
{
    return u1_;
}

double Grid::u1(int i, int j) const
{
    return u1_(i,j);
}

double &Grid::u1(int i, int j)
{
    return u1_(i,j);
}

const Array2D &Grid::u2() const
{
    return u2_;
}

double Grid::u2(int i, int j) const
{
    return u2_(i,j);
}

double &Grid::u2(int i, int j)
{
    return u2_(i,j);
}


const Array2D &Grid::x() const
{
    return x_;
}

double Grid::x(int i, int j) const
{
    return x_(i,j);
}

double &Grid::x(int i,int j)
{
    return x_(i,j);
}

Array2D &Grid::j()
{
    return j_;
}

double Grid::j(int i, int j) const
{
    return j_(i,j);
}

double &Grid::j(int i, int j)
{
    // TODO: insert return statement here
    return j_(i,j);
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

void Grid::fillArray(Array2D& x,const Array2D& VdM, const Array2D& L)
{
    for(int i=0;i<x.size()[0];i++){
        for(int j=0;j<x.size()[1];j++){
                x(i,j) =0.0;
                for(int p=0;p<VdM.size()[1];p++){
                    x(i,j) +=VdM(i,p)*L(j,p);
                    // if(abs(x(i,j))<1E-12)
                    //     x(i,j) = 0.0;
                }
            }
    }
}

void Grid::fillSolution(Variable& x,const Array2D& u)
{
    for(int i=0;i<u.size()[0];i++){
        for(int j=1;j<u.size()[1]-1;j++){
                    x(i*(u.size()[1]-2)+j-1) =u(i,j);
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
