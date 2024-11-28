#include <dg/grid.h>
#include "grid.h"

Grid::Grid(std::array<int, 1>  nCells, std::array<double, 1>  meshWidth, int numberNodes):
nCells_(nCells), meshWidth_(meshWidth), u_      ({nCells_[0]},   meshWidth_),
                                        u0_      ({nCells_[0]},   meshWidth_),
                                        x_      ({nCells_[0]*(numberNodes+1)+1},   meshWidth_),
                                        faces_  ({nCells_[0]+1},   meshWidth_),
                                        rhs_    (nCells_,   meshWidth_)
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
const Variable &Grid::u() const
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

const Variable &Grid::u0() const
{
    return u0_;
}

double Grid::u0(int i) const
{
    return u0_(i);
}

double &Grid::u0(int i)
{
    return u0_(i);
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
