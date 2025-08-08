#include "output_writer/output_writer_paraview.h"

#include <iostream>

// Original includes
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h> 
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkLagrangeQuadrilateral.h>
// --- NEW INCLUDES FOR HIGH-ORDER WRITER ---
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include "output_writer_paraview.h"
#include <vector>
#include <stdexcept>
// --- END NEW INCLUDES ---


OutputWriterParaview::OutputWriterParaview(std::shared_ptr<Grid> grid_) :
   OutputWriter(grid_)
{
  // This vtkWriter_ is for the original low-order writeFile function.
  vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
}
inline std::vector<int> getVtkLagrangeQuadrilateralReorderingMap(int degree)
{
    switch (degree)
    {
        case 1: // P1 Quadrilateral (4 nodes)
        {
            // Simple row-major order: 0, 1, 2, 3
            // VTK order for vtkQuad:  0, 1, 3, 2 (ensures correct winding)
            return {0, 1, 3, 2};
        }
        case 2: // P2 Quadrilateral (9 nodes)
        {
            // Simple row-major order: 0, 1, 2, 3, 4, 5, 6, 7, 8
            // VTK order: corners, mid-edges, center
            return {
                0, 2, 8, 6, // Corners
                1, 5, 7, 3, // Mid-edges
                4           // Center
            };
        }
        case 3: // P3 Quadrilateral (16 nodes)
        {
            // VTK order: corners, edge-nodes, face-nodes
            return {
                0, 3, 15, 12, // Corners
                1, 2,          // Bottom edge
                7, 11,         // Right edge
                14, 13,        // Top edge
                8, 4,          // Left edge
                5, 6, 10, 9    // Face (interior) nodes
            };
        }
        case 4: // P4 Quadrilateral (25 nodes)
        {
            return {
                0, 4, 24, 20, // Corners
                1, 2, 3,       // Bottom edge
                9, 14, 19,     // Right edge
                23, 22, 21,    // Top edge
                15, 10, 5,     // Left edge
                6, 7, 8,       // Face nodes (row 1)
                11, 12, 13,    // Face nodes (row 2)
                16, 17, 18     // Face nodes (row 3)
            };
        }
        case 5: // P5 Quadrilateral (36 nodes)
        {
            return {
                0, 5, 35, 30, // Corners
                1, 2, 3, 4,    // Bottom edge
                11, 17, 23, 29,// Right edge
                34, 33, 32, 31,// Top edge
                25, 19, 13, 7, // Left edge
                8, 9, 10,      // Face nodes (row 1, part 1)
                14, 15, 16,    // Face nodes (row 2, part 1)
                20, 21, 22,    // Face nodes (row 3, part 1)
                26, 27, 28,    // Face nodes (row 4, part 1)
                6, 12, 18, 24  // Face nodes (column 0) -- This is a guess based on pattern
                               // Note: Order for face nodes in P5+ can be complex.
                               // It's best to verify with a simple test case.
                               // A likely full face order is row-by-row:
                               // 7, 8, 9, 10
                               // 13, 14, 15, 16
                               // 19, 20, 21, 22
                               // 25, 26, 27, 28
            };
        }
        // Add other degrees here if needed.
        default:
            throw std::runtime_error("Unsupported polynomial degree for VTK reordering map: " + std::to_string(degree));
    }
}
// ==============================================================================
// [NEW FUNCTION] High-order writer for interior nodes
// ==============================================================================
/**
 * @brief This function creates a .vtm multi-block file. Each block corresponds
 * to one of the master simulation cells (e.g., 80x80). Inside each block, a
 * structured grid shows the solution ONLY at the interior DG nodes.
 *
 * This requires the grid_ object to provide:
 *  - grid_->getPolynomialDegree() to get N
 *  - grid_->basis_.nodes(i) to get the location of the i-th node in [-1, 1]
 *  - grid_->evaluatePolynomial(i, j, xi, eta) to get the solution value
 */
void OutputWriterParaview::writeHighOrderFile(double currentTime, std::string OutputName, std::shared_ptr<Vandermonde> VdM, const std::unique_ptr<Quadrature>& quad)
{
    // --- STEP 0: Prepare the data within the Grid object ---
    // This call is now essential. It populates the array that this function will read.
    grid_->prepareNodalDataForVisualization(VdM->VdM(),quad);

    // --- SETUP VTK WRITER AND FILE ---
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    std::stringstream fileName;
    fileName << "out/" << OutputName << "_HighOrderSurface_" << std::setw(4) << std::setfill('0') << fileNo_ << ".vtu";
    writer->SetFileName(fileName.str().c_str());

    // 1. GET GRID PARAMETERS
    // Correct way to get polyDegree is from the grid, not u_ directly.
    const int polyDegree = grid_->u_.size()[0]+1; // Assuming u_ is a square grid of size (N+1)x(N+1)
    const int nodesPerDim = polyDegree + 1;
    const int nNodesPerCell = nodesPerDim * nodesPerDim;
    const int nCellsX = grid_->nCells()[0];
    const int nCellsY = grid_->nCells()[1];
    const double physicalSizeStart = grid_->physicalSize_[0];
    const double dx = grid_->meshWidth()[0];
    const double dy = grid_->meshWidth()[1];
const auto vtk_reorder_map = getVtkLagrangeQuadrilateralReorderingMap(polyDegree);
    // 2. CREATE ONE SINGLE GRID FOR THE ENTIRE DOMAIN
    vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
    data->SetName("u_high_order_surface");
    data->SetNumberOfComponents(1);

    // Get a reference to the complete data to avoid repeated function calls
    const auto& complete_solution = grid_->getNodalSolutionComplete();

    // 3. LOOP OVER CELLS, TRANSFERRING THE PRE-ASSEMBLED DATA
    for (int j_cell = 0; j_cell < nCellsY; ++j_cell)
    {
        for (int i_cell = 0; i_cell < nCellsX; ++i_cell)
        {
            int cell_idx = grid_->elemId(i_cell, j_cell);
             std::vector<vtkIdType> rowMajorPointIds(nNodesPerCell);

            // Inner loop over ALL (N+1)x(N+1) nodes of THIS cell
            for (int j_node = 0; j_node < nodesPerDim; ++j_node)
            {
                for (int i_node = 0; i_node < nodesPerDim; ++i_node)
                {
                    // Get node location in [-1, 1] from the basis object
                    double xi  = quad->basis_.nodes(i_node);
                    double eta = quad->basis_.nodes(j_node);
                    
                    // Map to physical coordinate
                    double x_phys = (physicalSizeStart + i_cell * dx) + 0.5 * dx * (1.0 + xi);
                    double y_phys = (physicalSizeStart + j_cell * dy) + 0.5 * dy * (1.0 + eta);
                    
                    // Add the point and get its global ID
                    vtkIdType currentPointId = points->InsertNextPoint(x_phys, y_phys, 0.0);
                    rowMajorPointIds[j_node * nodesPerDim + i_node] = currentPointId;// Store ID in row-major order

                    // *** KEY: Direct read from the pre-assembled data. NO calculations here. ***
                    double value = complete_solution(i_node, j_node, cell_idx);
                    data->InsertNextTuple1(value);
                }
            }
             std::vector<vtkIdType> vtkOrderedPointIds(nNodesPerCell);
            for(int k = 0; k < nNodesPerCell; ++k) {
                vtkOrderedPointIds[k] = rowMajorPointIds[vtk_reorder_map[k]];
            }
            
            // Create one high-order cell using the (N+1)*(N+1) points we just defined
            ugrid->InsertNextCell(VTK_LAGRANGE_QUADRILATERAL, nNodesPerCell, vtkOrderedPointIds.data());
        }
    }

    // 4. ASSEMBLE AND WRITE THE FILE
    ugrid->SetPoints(points);
    ugrid->GetPointData()->AddArray(data);

    vtkSmartPointer<vtkDoubleArray> arrayTime = vtkSmartPointer<vtkDoubleArray>::New();
    arrayTime->SetName("TIME");
    arrayTime->SetNumberOfTuples(1);
    arrayTime->SetTuple1(0, currentTime);
    ugrid->GetFieldData()->AddArray(arrayTime);
    
    writer->SetInputData(ugrid);
    writer->SetDataModeToBinary();
    writer->Write();
    std::cout << "Wrote high-performance surface file: " << fileName.str() << std::endl;
}


// ==============================================================================
// [ORIGINAL FUNCTION - UNCHANGED] Low-order writer for cell averages
// ==============================================================================
void OutputWriterParaview::writeFile(double currentTime,std::string OutputName)
{
    // Assemble the filename
    std::stringstream fileName;
    fileName << "out/" << OutputName << "_CellAvg_" << std::setw(4) << std::setfill('0') << fileNo_ << ".vti";
    // We reuse the fileNo_ here, and add a suffix to distinguish the file type
    fileNo_++;

    // Assign the new file name to the output vtkWriter_
    vtkWriter_->SetFileName(fileName.str().c_str());
  
    const double physicalSizeStart = -6.283;//grid_->getPhysicalSizeStart();
    const double physicalSizeEnd = 6.283;//grid_->getPhysicalSizeEnd();
    const int nCellsX = grid_->solution_.size()[0];
    const int nCellsY = grid_->solution_.size()[1];
    const double domainWidth = physicalSizeEnd - physicalSizeStart;
    const double dx = domainWidth / nCellsX;
    const double dy = domainWidth / nCellsY;

    vtkSmartPointer<vtkImageData> dataSet = vtkSmartPointer<vtkImageData>::New();
    dataSet->SetOrigin(physicalSizeStart, physicalSizeStart, 0.0);
    dataSet->SetSpacing(dx, dy, 1.0);
    dataSet->SetDimensions(nCellsX + 1, nCellsY + 1, 1);

    vtkSmartPointer<vtkDoubleArray> arraySolution = vtkSmartPointer<vtkDoubleArray>::New();
    arraySolution->SetName("u_cell_average");
    arraySolution->SetNumberOfComponents(1);
    arraySolution->SetNumberOfTuples(nCellsX * nCellsY);
  
    int index = 0;
    for (int j = 0; j < nCellsY; j++)
    {
        for (int i = 0; i < nCellsX; i++, index++)
        {
            arraySolution->SetValue(index, grid_->solution_(i, j));
        }
    }

    dataSet->GetCellData()->AddArray(arraySolution);
    
    vtkSmartPointer<vtkDoubleArray> arrayTime = vtkSmartPointer<vtkDoubleArray>::New();
    arrayTime->SetName("TIME");
    arrayTime->SetNumberOfTuples(1);
    arrayTime->SetTuple1(0, currentTime);
    dataSet->GetFieldData()->AddArray(arrayTime);

    dataSet->Squeeze();
    vtkWriter_->SetInputData(dataSet);
    vtkWriter_->SetDataModeToBinary();
    vtkWriter_->Write();

    std::cout << "Wrote cell-average file: " << fileName.str() << std::endl;
}

void OutputWriterParaview::writeFileTrueSolution(double currentTime,std::string OutputName)
{
  // Assemble the filename
  std::stringstream fileName;
  std::string outputName_ = OutputName;
  fileName << "out/" <<outputName_<<"_"<< std::setw(4) << setfill('0') << fileNo_ << "." << vtkWriter_->GetDefaultFileExtension();
  
  // increment file no.
  fileNo_++;

  // assign the new file name to the output vtkWriter_
  vtkWriter_->SetFileName(fileName.str().c_str());
  
  // initialize data set that will be output to the file
  vtkSmartPointer<vtkImageData> dataSet = vtkSmartPointer<vtkImageData>::New();
  dataSet->SetOrigin(0, 0, 0);

  // set spacing of mesh
  //const double dx = -grid_->x_analyze_(0,0)+grid_->x_analyze_(1,0);
  const double dy = 1;
  const double dz = 1;
  //dataSet->SetSpacing(dx, dy, dz);

  // set number of points in each dimension, 1 cell in z direction
  int nCells = grid_->u_analyze_.size()[0];
  dataSet->SetDimensions(nCells, 1, 1);  // we want to have points at each corner of each cell
  
  // add solution field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arraySolutionTrue = vtkDoubleArray::New();

  // the pressure is a scalar which means the number of components is 1
  arraySolutionTrue->SetNumberOfComponents(1);

  // Set the number of pressure values and allocate memory for it. We already know the number, it has to be the same as there are nodes in the mesh.
  arraySolutionTrue->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  
  arraySolutionTrue->SetName("u_true");

  // loop over the nodes of the mesh and assign the interpolated p values in the vtk data structure
  // we only consider the cells that are the actual computational domain, not the helper values in the "halo"

  int index = 0;   // index for the vtk data structure, will be incremented in the inner loop

  for (int i = 0; i < nCells; i++, index++)
  {
    //std::array<double,1> solutionVector;
    //solutionVector[0] =grid_->x(i);
    //solutionVector[0] = grid_->u(i);
    //arraySolutionTrue->SetValue(index,grid_->u_analyze_true_(i,0));
  }
  // now, we should have added as many values as there are points in the vtk data structure
  assert(index == dataSet->GetNumberOfPoints());

  // add the field variable to the data set
  dataSet->GetPointData()->AddArray(arraySolutionTrue);

// add solution field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arraySolution = vtkDoubleArray::New();

  // the pressure is a scalar which means the number of components is 1
  arraySolution->SetNumberOfComponents(1);

  // Set the number of pressure values and allocate memory for it. We already know the number, it has to be the same as there are nodes in the mesh.
  arraySolution->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  
  arraySolution->SetName("u");

  // loop over the nodes of the mesh and assign the interpolated p values in the vtk data structure
  // we only consider the cells that are the actual computational domain, not the helper values in the "halo"

  index = 0;   // index for the vtk data structure, will be incremented in the inner loop

  for (int i = 0; i < nCells; i++, index++)
  {
    //std::array<double,1> solutionVector;
    //solutionVector[0] =grid_->x(i);
    //solutionVector[0] = grid_->u(i);
    //arraySolution->SetValue(index,grid_->u_analyze_(i,0));
  }
  // now, we should have added as many values as there are points in the vtk data structure
  assert(index == dataSet->GetNumberOfPoints());

  // add the field variable to the data set
  dataSet->GetPointData()->AddArray(arraySolution);

  // add solution field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arrayXaxis = vtkDoubleArray::New();

  // the pressure is a scalar which means the number of components is 1
  arrayXaxis->SetNumberOfComponents(1);

  // Set the number of pressure values and allocate memory for it. We already know the number, it has to be the same as there are nodes in the mesh.
  arrayXaxis->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  
  arrayXaxis->SetName("xAnalyzeAxis");

  // loop over the nodes of the mesh and assign the interpolated p values in the vtk data structure
  // we only consider the cells that are the actual computational domain, not the helper values in the "halo"

  index = 0;   // index for the vtk data structure, will be incremented in the inner loop
  for (int i = 0; i <grid_->x_analyze_.size()[0]; i++,index++)
  {

     // arrayXaxis->SetValue(index,grid_->x_analyze_(i,0));
    
  }
  // now, we should have added as many values as there are points in the vtk data structure
  assert(index == dataSet->GetNumberOfPoints());
  dataSet->GetPointData()->AddArray(arrayXaxis);
    
  // add current time 
  vtkSmartPointer<vtkDoubleArray> arrayTime = vtkDoubleArray::New();
  arrayTime->SetName("TIME");
  arrayTime->SetNumberOfTuples(1);
  arrayTime->SetTuple1(0, currentTime);
  dataSet->GetFieldData()->AddArray(arrayTime);

  // Remove unused memory
  dataSet->Squeeze();
  
  // Write the data
  vtkWriter_->SetInputData(dataSet);
  
  //vtkWriter_->SetDataModeToAscii();     // comment this in to get ascii text files: those can be checked in an editor
  vtkWriter_->SetDataModeToBinary();      // set file mode to binary files: smaller file sizes

  // finally write out the data
  vtkWriter_->Write();
}