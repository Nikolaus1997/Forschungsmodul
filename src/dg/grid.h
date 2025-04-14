#pragma once

#include <storage/variable.h>
#include <memory>   
#include <vector>
#include <storage/Vdm.h>


class Grid
{
public:
    Grid(std::array<int, 1>  nCells, std::array<double, 1>  meshWidth, int numberNodes);

    //get the mesh width, i.e. the length of a single cell in x and y direction 
    const std::array<double, 1> meshWidth() const;

    //get number of cells in each coordinate direction 
    const std::array<int, 1> nCells() const;


    Array2D &u();


    const Array2D &u2() const;

    double u2(int i, int j) const;

    double &u2(int i, int j);

    const Array2D &u1() const;

    double u1(int i, int j) const;

    double &u1(int i, int j);

    Array2D &ut();

    double ut(int i) const;

    double &ut(int i);


    const Array2D &x() const;

    double x(int i, int j) const;

    double &x(int i, int j);

    Array2D &j() ;

    double j(int i, int j) const;

    double &j(int i, int j);

    const Variable &faces() const;

    double faces(int i) const;

    double &faces(int i);

    const Variable &l2_error() const;
    
    double l2_error(int i) const;
    
    double &l2_error(int i);

    const Variable &linf_error() const;

    double linf_error(int i) const;

    double &linf_error(int i);


    const Variable &rhs() const;

    double rhs(int i) const;

    double &rhs(int i);
    

    double dx() const;

    void fillArray(Array2D& x, const Array2D& VdM, const Array2D& L);
    void fillDerivative(Variable& x, std::shared_ptr<Vandermonde> VdM);
    void fillSolution(Variable& x,const Array2D& u);

protected:
    const std::array<int, 1>        nCells_;
    const std::array<double, 1>     meshWidth_;
    Variable solution_,solutionJ_;
    Variable derivative_;
    Array2D u_,u_analyze_, u_analyze_true_;
    Array2D j_,j_1_,j_2_;
    Array2D true_solution_;
    Array2D q_;
    Array2D u2_;
    Array2D u1_;
    Array2D ut_;
    Array2D x_, x_analyze_;
    Variable faces_;
    Variable l2_error_;
    Variable linf_error_;
    Variable rhs_;
    friend class Computation;
    friend class OutputWriterParaview;

};