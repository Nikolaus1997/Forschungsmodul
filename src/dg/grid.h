#pragma once

#include <storage/variable.h>


class Grid
{
public:
    Grid(std::array<int, 1>  nCells, std::array<double, 1>  meshWidth, int numberNodes);

    //get the mesh width, i.e. the length of a single cell in x and y direction 
    const std::array<double, 1> meshWidth() const;

    //get number of cells in each coordinate direction 
    const std::array<int, 1> nCells() const;


    const Variable &u() const;

    double u(int i) const;

    double &u(int i);

    const Variable &u2() const;

    double u2(int i) const;

    double &u2(int i);

    const Variable &u1() const;

    double u1(int i) const;

    double &u1(int i);

    const Variable &ut() const;

    double ut(int i) const;

    double &ut(int i);


    const Variable &x() const;

    double x(int i) const;

    double &x(int i);


    const Variable &faces() const;

    double faces(int i) const;

    double &faces(int i);


    const Variable &rhs() const;

    double rhs(int i) const;

    double &rhs(int i);
    

    double dx() const;

protected:
    const std::array<int, 1>        nCells_;
    const std::array<double, 1>     meshWidth_;
    Variable u_;
    Variable u2_;
    Variable u1_;
    Variable ut_;
    Variable x_;
    Variable faces_;
    Variable rhs_;
    friend class Computation;
    friend class OutputWriterParaview;

};