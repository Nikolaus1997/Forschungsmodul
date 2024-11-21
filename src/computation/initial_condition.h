#pragma once

#include <vector>
#include <storage/variable.h>

class InitialCondition
{
public:
    // Constructor
    InitialCondition(std::array<double,2> physicalSize, int nCells);

    // Destructor (optional, for cleanup if needed)
    Variable unitStep(Variable u);

    // Other member functions and data members can be added here

private:
    int nCells_; // Store the number of cells
    Variable u_;
};