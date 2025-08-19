#pragma once

#include <storage/variable.h>
#include <memory>   
#include <vector>
#include <cmath>
#include <storage/Vdm.h>
#include <storage/array3d.h>
#include <storage/array4d.h>
#include <storage/array2d.h>
#include <dg/num_flux.h>
#include <dg/flux.h>

class Grid
{
public:

    Grid(std::array<double, 2> physicalSizeX, std::array<double, 2> physicalSizeY, std::array<int, 2> nCells, std::array<double, 2> meshWidth, int numberNodes);

    // get the mesh width, i.e. the length of a single cell in x and y direction
    const std::array<double, 2> meshWidth() const;

    //get number of cells in each coordinate direction 
    const std::array<int, 2> nCells() const;


    Array3D &u();


    const Array3D &u2() const;

    double u2(int i, int j, int k) const;

    double &u2(int i, int j, int k);

    const Array3D &u1() const;

    double u1(int i, int j, int k) const;

    double &u1(int i, int j,int k);

    Array3D &ut();

    double ut(int i, int j, int k) const;

    double &ut(int i, int j, int k);


    const Array3D &x() const;

    double x(int i, int j, int k) const;

    double &x(int i, int j, int k);

    Array4D &j() ;

    double j(int i, int j, int k) const;

    double &j(int i, int j, int k);

    const Array2D&faces() const;

    double faces(int i, int j) const;

    double &faces(int i, int j);

    const Array1D &l2_error() const;
    
    double l2_error(int i) const;
    
    double &l2_error(int i);

    const Array1D &linf_error() const;

    double linf_error(int i) const;

    double &linf_error(int i);

    

    double dx() const;

    void fillArray(Array3D& x, const Array3D& VdM, const Array2D& L);
    void fillArray(Array4D &x, const Array4D &VdM, const Array2D &L);
    void fillFaces(Array3D &f, const Array3D &VdM, const Array2D &L);
    void fillFaces(Array3D &f, const Array4D &VdM, const Array2D &L);
    void prepareNodalDataForVisualization(const Array3D &VdM_, const std::unique_ptr<Quadrature> &quad);
    void fillDerivative(Array2D &x, std::shared_ptr<Vandermonde> VdM);
    void fillSolution(Array2D &x, const Array3D &u);
    void fillSolution(Array2D &x, const Array4D &u, int dim);
    
    Array3D getNodalSolutionComplete();

    double evaluatePolynomial(int cell_i, int cell_j, double xi, double eta, const Array3D VdM_, const std::unique_ptr<Quadrature> &quad) const;

    

protected:
    const std::array<int, 2>        nCells_;
    const std::array<double, 2>     meshWidth_;
    const std::array<double, 2>     physicalSizeX_,physicalSizeY_;
    Array2D solution_,solutionJ_;
    Array2D  derivative_;
    Array3D u_,u_analyze_, u_analyze_true_;
    Array4D j_,j_1_,j_2_,jt_;
    Array3D true_solution_, faceId, faceIdQ, faceId1, faceId2, face_dt;
    Array3D faceIdJ, faceIdJ1, faceIdJ2, face_dtJ;
    Array3D q_;
    Array3D u2_;
    Array3D u1_;
    Array3D ut_;
    Array3D nodal_solution_complete_;
    Array2D elemId;
    Array3D x_, x_analyze_,y_;
    Array2D faces_, faceFlux_, faceFluxJ_;
    Array1D l2_error_;
    Array1D linf_error_;
    friend class Computation;
    friend class OutputWriterParaview;
};