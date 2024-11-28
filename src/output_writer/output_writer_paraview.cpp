#include "output_writer/output_writer_paraview.h"

#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <iostream>

OutputWriterParaview::OutputWriterParaview(std::shared_ptr<Grid> grid_) :
   OutputWriter(grid_)
{
  // Create a vtkWriter_
  vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
}

void OutputWriterParaview::writeFile(double currentTime)
{
  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output_" << std::setw(4) << setfill('0') << fileNo_ << "." << vtkWriter_->GetDefaultFileExtension();
  
  // increment file no.
  fileNo_++;

  // assign the new file name to the output vtkWriter_
  vtkWriter_->SetFileName(fileName.str().c_str());
  
  // initialize data set that will be output to the file
  vtkSmartPointer<vtkImageData> dataSet = vtkSmartPointer<vtkImageData>::New();
  dataSet->SetOrigin(0, 0, 0);

  // set spacing of mesh
  const double dx = grid_->meshWidth()[0];
  const double dy = 1;
  const double dz = 1;
  dataSet->SetSpacing(dx, dy, dz);

  // set number of points in each dimension, 1 cell in z direction
  int nCells = grid_->u_.size()[0];
  dataSet->SetDimensions(nCells, 1, 1);  // we want to have points at each corner of each cell
  
  // add solution field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arraySolution = vtkDoubleArray::New();

  // the pressure is a scalar which means the number of components is 1
  arraySolution->SetNumberOfComponents(1);

  // Set the number of pressure values and allocate memory for it. We already know the number, it has to be the same as there are nodes in the mesh.
  arraySolution->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  
  arraySolution->SetName("solution");

  // loop over the nodes of the mesh and assign the interpolated p values in the vtk data structure
  // we only consider the cells that are the actual computational domain, not the helper values in the "halo"

  int index = 0;   // index for the vtk data structure, will be incremented in the inner loop

  for (int i = 0; i < nCells; i++, index++)
  {
    //std::array<double,1> solutionVector;
    //solutionVector[0] =grid_->x(i);
    //solutionVector[0] = grid_->u(i);
    arraySolution->SetValue(index,grid_->u(i));
  }

  // now, we should have added as many values as there are points in the vtk data structure
  assert(index == dataSet->GetNumberOfPoints());

  // add the field variable to the data set
  dataSet->GetPointData()->AddArray(arraySolution);

  // // add solution field variable
  // // ---------------------------
  // vtkSmartPointer<vtkDoubleArray> arrayXaxis = vtkDoubleArray::New();

  // // the pressure is a scalar which means the number of components is 1
  // arrayXaxis->SetNumberOfComponents(1);

  // // Set the number of pressure values and allocate memory for it. We already know the number, it has to be the same as there are nodes in the mesh.
  // arrayXaxis->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  
  // arrayXaxis->SetName("xAxis");

  // // loop over the nodes of the mesh and assign the interpolated p values in the vtk data structure
  // // we only consider the cells that are the actual computational domain, not the helper values in the "halo"

  // index = 0;   // index for the vtk data structure, will be incremented in the inner loop
  // nCells = grid_->x_.size()[0];
  // for (int i = 0; i < nCells; i++, index++)
  // {
  //   //std::array<double,1> solutionVector;
  //   //solutionVector[0] =grid_->x(i);
  //   //solutionVector[0] = grid_->u(i);
  //   arrayXaxis->SetValue(index,grid_->x(i));
  // }
  // // now, we should have added as many values as there are points in the vtk data structure
  // assert(index == dataSet->GetNumberOfPoints());
  // dataSet->GetPointData()->AddArray(arrayXaxis);
    
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