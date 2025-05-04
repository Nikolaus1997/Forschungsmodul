#include <dg/grid.h>
#include "grid.h"

Grid::Grid(std::array<int, 2>  nCells, std::array<double, 2>  meshWidth, int numberNodes):
nCells_(nCells), meshWidth_(meshWidth), u_      ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1])}),
                                        j_      ({int(nCells_[0]),int(nCells_[1]),(numberNodes+2)}),
                                        jt_     ({int(nCells_[0]),int(nCells_[1]),(numberNodes+2)}),
                                        j_1_    ({int(nCells_[0]),int(nCells_[1]),(numberNodes+2)}),
                                        j_2_    ({int(nCells_[0]),int(nCells_[1]),(numberNodes+2)}),
                                        true_solution_({int(nCells_[0]),int(nCells_[0]),(numberNodes+2)}),
                                        q_      ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1])}),
                                        ut_      ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1])}),
                                        u1_     ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1])}),
                                        u2_    ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1])}),
                                        elemId ({int(nCells_[0]),int(nCells_[1])}),
                                        faceFlux_({numberNodes,4}),
                                        faceFluxQ_({numberNodes,4}),
                                        faceIdQ({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        faceId({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        faceId1({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        faceId2({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        face_dt({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        solution_({int(nCells_[0]*(numberNodes)),int(nCells_[1]*(numberNodes))}),
                                        solutionJ_({int(nCells_[0]*(numberNodes)),int(nCells_[0]*(numberNodes))}),
                                        derivative_({int(nCells_[0]*(numberNodes)),int(nCells_[0]*(numberNodes))}),
                                        x_      ({(numberNodes+2),(numberNodes+2),nCells_[0]}),
                                        y_      ({(numberNodes+2),(numberNodes+2),nCells_[1]}),
                                        x_analyze_    ({300,1}),
                                        u_analyze_    ({300,1}),
                                        u_analyze_true_    ({300,1}),
                                        faces_  ({int(nCells_[0]+1),2}),
                                        l2_error_({1}),
                                        linf_error_({1})
{
}
const std::array<double,2> Grid::meshWidth() const
{
return meshWidth_;
}
const std::array<int,2> Grid::nCells() const
{
return nCells_;
}
Array3D &Grid::u()
{
    return u_;
}


Array3D &Grid::ut()
{
    return ut_;
}


const Array3D &Grid::u1() const
{
    return u1_;
}

double Grid::u1(int i, int j, int k) const
{
    return u1_(i,j,k);
}

double &Grid::u1(int i, int j, int k)
{
    return u1_(i,j,k);
}

const Array3D &Grid::u2() const
{
    return u2_;
}

double Grid::u2(int i, int j, int k) const
{
    return u2_(i,j,k);
}

double &Grid::u2(int i, int j, int k)
{
    return u2_(i,j,k);
}


const Array3D &Grid::x() const
{
    return x_;
}

double Grid::x(int i, int j, int k) const
{
    return x_(i,j,k);
}

double &Grid::x(int i,int j,int k)
{
    return x_(i,j,k);
}

Array3D &Grid::j()
{
    return j_;
}

double Grid::j(int i, int j, int k) const
{
    return j_(i,j,k);
}

double &Grid::j(int i, int j, int k)
{
    // TODO: insert return statement here
    return j_(i,j,k);
}

const Array2D &Grid::faces() const
{

    return faces_;
}

double Grid::faces(int i, int j) const
{
    return faces_(i,j);
}

double &Grid::faces(int i, int j)
{

    return faces_(i,j);
}

const Array1D &Grid::l2_error() const
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

const Array1D &Grid::linf_error() const
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



double Grid::dx() const
{
    return meshWidth_[0];
}

void Grid::fillArray(Array3D& x,const Array3D& VdM, const Array2D& L)
{

    int nNodes = L.size()[0]; // Number of interpolation points (per direction)
    // std::cout<<" FILL ARRAY "<<std::endl;
    // std::cout<<" SIZE "<<x.size()[0]<<" "<<x.size()[1]<<" "<<x.size()[2]<<std::endl;
    for (int iCell = 0; iCell < nCells_[0]; ++iCell) {
        for (int jCell = 0; jCell < nCells_[1]; ++jCell) {
            for (int node_i = 0; node_i < nNodes-2; ++node_i) {
                for (int node_j = 0; node_j < nNodes-2; ++node_j) {
                    int idx = elemId(iCell,jCell);
                    x(node_i, node_j, idx) = 0.0;
                    for (int l = 0; l < VdM.size()[2]; ++l) {
                        //std::cout<<" I "<<iCell<<" J "<<jCell<<" nodeI "<<node_i<<" nodeJ "<<node_j<<" l "<<l<<" idx "<<idx<<std::endl;
                        x(node_i, node_j, idx) += 
                            VdM(iCell, jCell, l) * L(node_i+1, l) * L(node_j+1, l);
                    }
                }
            }
        }
    }
}

void Grid::fillFaces(Array3D& f, const Array3D& VdM, const Array2D& L)
{
    //left Interface 0, right Interface 2, top Interface 3, bottom Interface 1
    int nNodes = L.size()[0]-2; // Number of interpolation points (per direction)
    for(int iCell = 0; iCell < nCells_[0]; ++iCell) {
        for (int jCell = 0; jCell < nCells_[1]; ++jCell) {
            for (int faceId = 0; faceId <4; ++faceId) {
                for (int node_j = 0; node_j < nNodes; ++node_j) {
                    int idx = elemId(iCell,jCell);
                    f(node_j, faceId, idx) = 0.0;
                    for (int l = 0; l < VdM.size()[2]; ++l) {
                        if(faceId==0 or faceId==1)
                            f(node_j, faceId, idx) += VdM(iCell, jCell, l) * L(0, l) * L(node_j+1, l);
                        else if(faceId==2 or faceId==3)
                            f(node_j, faceId, idx) += VdM(iCell, jCell, l) * L(node_j+1, l) * L(nNodes+1, l);
                    }
                }
            }
        }
    }

}

void Grid::fillSolution(Array2D& x,const Array3D& u)
{
    // std::cout<<" FILL SOLUTION "<<std::endl;
    // std::cout<<" SIZE "<<u.size()[0]<<" "<<u.size()[1]<<" "<<u.size()[2]<<std::endl;
    for (int iCell = 0; iCell < nCells_[0]; ++iCell) {
        for (int jCell = 0; jCell < nCells_[1]; ++jCell) {
            int idx = elemId(iCell,jCell);
            for (int node_i = 0; node_i < u.size()[0]; ++node_i) {
                for (int node_j = 0; node_j < u.size()[1]; ++node_j) {
                    x(node_i+iCell*(u.size()[0]), node_j+jCell*(u.size()[1])) = u(node_i, node_j, idx);
                }
            }
        }
    }   

}

void Grid::fillDerivative(Array2D& x,const std::shared_ptr<Vandermonde> VdM)
{
    for(int i=0;i<VdM->size()[0];i++){
        for(int k=0;k<VdM->size()[1];k++){
        for(int j=1;j<VdM->size()[2];j++){
                x(i*(VdM->size()[1]-1)+j-1,k*(VdM->size()[1]-1)+j-1) =0.0;
                for(int p=0;p<VdM->size()[1];p++){
                    x(i*(VdM->size()[1]-1)+j-1,k*(VdM->size()[1]-1)+j-1) +=VdM->VdMt(i,k,p)*VdM->L(j,p)*VdM->L(k,p);
                }
            }
        }
    }
}

