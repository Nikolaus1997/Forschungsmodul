#include <dg/grid.h>
#include "grid.h"

Grid::Grid(std::array<double, 2>  physicalSizeX,std::array<double, 2>  physicalSizeY,std::array<int, 2>  nCells, std::array<double, 2>  meshWidth, int numberNodes):
physicalSizeX_(physicalSizeX),physicalSizeY_(physicalSizeY),nCells_(nCells), meshWidth_(meshWidth), 
                                        u_      ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1])}),
                                        u1_     ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1])}),
                                        u2_     ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1])}),
                                        ut_     ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1])}),
                                        faceId({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        faceId1({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        faceId2({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        face_dt({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        j_      ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1]),2}),
                                        jt_     ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1]),2}),
                                        j_1_    ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1]),2}),
                                        j_2_    ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1]),2}),
                                        faceIdJ({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        faceIdJ1({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        faceIdJ2({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        face_dtJ({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),
                                        true_solution_({int(nCells_[0]),int(nCells_[0]),(numberNodes+2)}),
                                        q_      ({numberNodes,numberNodes,int(nCells_[0]*nCells_[1])}),
                                        elemId ({int(nCells_[0]),int(nCells_[1])}),
                                        faceFlux_({numberNodes,4}),
                                        faceFluxJ_({numberNodes,4}),
                                        faceIdQ({numberNodes,4,int(nCells_[0])*int(nCells_[1])}),

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
                                        linf_error_({1}),
                                        nodal_solution_complete_({numberNodes+2,numberNodes+2,int(nCells_[0]*nCells_[1])})
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

Array4D &Grid::j()
{
    return j_;
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

void Grid::fillArray(Array4D& x,const Array4D& VdM, const Array2D& L)
{
    int nNodes = L.size()[0]; // Number of interpolation points (per direction)
    // std::cout<<" FILL ARRAY "<<std::endl;
    // std::cout<<" SIZE "<<x.size()[0]<<" "<<x.size()[1]<<" "<<x.size()[2]<<std::endl;
    for(int i = 0; i <2; i++){
        for (int iCell = 0; iCell < nCells_[0]; ++iCell) {
            for (int jCell = 0; jCell < nCells_[1]; ++jCell) {
                for (int node_i = 0; node_i < nNodes-2; ++node_i) {
                    for (int node_j = 0; node_j < nNodes-2; ++node_j) {
                        int idx = elemId(iCell,jCell);
                        x(node_i, node_j, idx, i) = 0.0;
                        for (int l = 0; l < VdM.size()[2]; ++l) {
                            //std::cout<<" I "<<iCell<<" J "<<jCell<<" nodeI "<<node_i<<" nodeJ "<<node_j<<" l "<<l<<" idx "<<idx<<std::endl;
                            x(node_i, node_j, idx, i) += 
                                VdM(iCell, jCell, l, i) * L(node_i+1, l) * L(node_j+1, l);
                        }
                    }
                }
            }
        }
    }
}

void Grid::fillFaces(Array3D& f, const Array3D& VdM, const Array2D& L)
{   /**  _ 3 _
        |     |
       0|     |2
        |_ _ _|
           1
    */
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

void Grid::fillFaces(Array3D& f, const Array4D& VdM, const Array2D& L)
{   /**  _ 3 _
        |     |
       0|     |2
        |_ _ _|
           1
    */
    //left Interface 0, right Interface 2, top Interface 3, bottom Interface 1
    int nNodes = L.size()[0]-2; // Number of interpolation points (per direction)
    for(int iCell = 0; iCell < nCells_[0]; ++iCell) {
        for (int jCell = 0; jCell < nCells_[1]; ++jCell) {
            for (int faceId = 0; faceId <4; ++faceId) {
                for (int node_j = 0; node_j < nNodes; ++node_j) {
                    int idx = elemId(iCell,jCell);
                    f(node_j, faceId, idx) = 0.0;
                    for (int l = 0; l < VdM.size()[2]; ++l) {
                        if(faceId==0)
                            f(node_j, faceId, idx) += VdM(iCell, jCell, l,0) * L(0, l) * L(node_j+1, l);
                        else if(faceId==1)
                            f(node_j, faceId, idx) += VdM(iCell, jCell, l,1) * L(0, l) * L(node_j+1, l);
                        else if(faceId==2)
                            f(node_j, faceId, idx) += VdM(iCell, jCell, l,0) * L(node_j+1, l) * L(nNodes+1, l);
                        else if(faceId==3)
                            f(node_j, faceId, idx) += VdM(iCell, jCell, l,1) * L(node_j+1, l) * L(nNodes+1, l);
                    }
                }
            }
        }
    }

}

void Grid::prepareNodalDataForVisualization(const Array3D& VdM_, const std::unique_ptr<Quadrature>& quad)
{
    // --- LOGIC CORRECTION ---
    // The number of interior nodes is u_.size()[0], which is N-1.
    // The polynomial degree N is therefore u_.size()[0] + 1.
    const int polyDegree = u_.size()[0] ;
    const int nodesPerDim = polyDegree + 1;
    const int nCellsX = nCells_[0];
    const int nCellsY = nCells_[1];



    // Loop over each master cell to assemble its complete nodal data
    for (int j_cell = 0; j_cell < nCellsY; ++j_cell)
    {
        for (int i_cell = 0; i_cell < nCellsX; ++i_cell)
        {
            int cell_idx = this->elemId(i_cell, j_cell);

            // Loop over all (N+1)x(N+1) nodes of the target visualization grid
            for (int j_node = 0; j_node < nodesPerDim; ++j_node)
            {
                for (int i_node = 0; i_node < nodesPerDim; ++i_node)
                {
                    double value = 0.0;
                    
                    // The boundary checks must use `polyDegree`, not `polyDegree-1`.
                    bool is_interior_row = (j_node > 0 && j_node < polyDegree);
                    bool is_interior_col = (i_node > 0 && i_node < polyDegree);

                    if (is_interior_row && is_interior_col)
                    {
                        // Case 1: Interior Node. This logic is correct.
                        value = u_(i_node - 1, j_node - 1, cell_idx);
                    }
                    else if (is_interior_row && i_node == 0) // Left Face
                    {
                        // Case 2: Left Face. This logic is correct.
                        value = faceId(j_node - 1, 0, cell_idx);
                    }
                    else if (is_interior_row && i_node == polyDegree) // Right Face
                    {
                        // Case 3: Right Face. This logic is correct.
                        value = faceId(j_node - 1, 2, cell_idx);
                    }
                    else if (is_interior_col && j_node == 0) // Bottom Face
                    {
                        // Case 4: Bottom Face. This logic is correct.
                        value = faceId(i_node - 1, 1, cell_idx);
                    }
                    else if (is_interior_col && j_node == polyDegree) // Top Face
                    {
                        // Case 5: Top Face. This logic is correct.
                        value = faceId(i_node - 1, 3, cell_idx);
                    }
                    else
                    {
                        // --- LOGIC CORRECTION ---
                        // Case 6: Corner Node. The original loops were incorrect.
                        // We MUST call the (now corrected) evaluatePolynomial function
                        // to get the true value at the corner.
                        double xi  = quad->basis_.nodes(i_node);
                        double eta = quad->basis_.nodes(j_node);

                        value = this->evaluatePolynomial(i_cell, j_cell, xi, eta, VdM_, quad);
                    }
                    nodal_solution_complete_(i_node, j_node, cell_idx) = value;
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
void Grid::fillSolution(Array2D& x,const Array4D& u, int dim)
{
    // std::cout<<" FILL SOLUTION "<<std::endl;
    // std::cout<<" SIZE "<<u.size()[0]<<" "<<u.size()[1]<<" "<<u.size()[2]<<std::endl;
    for (int iCell = 0; iCell < nCells_[0]; ++iCell) {
        for (int jCell = 0; jCell < nCells_[1]; ++jCell) {
            int idx = elemId(iCell,jCell);
            for (int node_i = 0; node_i < u.size()[0]; ++node_i) {
                for (int node_j = 0; node_j < u.size()[1]; ++node_j) {
                    x(node_i+iCell*(u.size()[0]), node_j+jCell*(u.size()[1])) = u(node_i, node_j, idx,dim);
                }
            }
        }
    }   

}

Array3D Grid::getNodalSolutionComplete()
{
    return nodal_solution_complete_;
}

double Grid::evaluatePolynomial(int cell_i, int cell_j, double xi, double eta, const Array3D VdM_, const std::unique_ptr<Quadrature>& quad) const
{
    // This function requires access to the modal coefficients (VdM).
    // Let's assume they are stored in a member variable `VdM_`.
    // If they are not, this function will need access to them.
    // Placeholder for your modal coefficient array

    const int polyDegree = u_.size()[0]; // Assuming u_ is a square grid of size (N+1)x(N+1)
    const int nBasisFunctions = polyDegree + 1;
    double solutionValue = 0.0;

    // This loop must match your basis expansion
    for (int l_j = 0; l_j < nBasisFunctions; ++l_j) {

            double basis_val_x = quad->LegendrePolynomialAndDerivative(l_j, xi)[0];
            double basis_val_y = quad->LegendrePolynomialAndDerivative(l_j, eta)[0];
            
            // Note: VdM dimensions might be different, e.g., VdM(cell_idx, l)
            solutionValue += VdM_(cell_i, cell_j, l_j) * basis_val_x * basis_val_y;
        
    }
    return solutionValue;
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

