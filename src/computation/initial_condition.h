#pragma once

#include <vector>

class InitialCondition
{
public:
    // Constructor
    InitialCondition(int nCells);

    // Destructor (optional, for cleanup if needed)
    auto unitStep();

    // Other member functions and data members can be added here

private:
    int nCells_; // Store the number of cells
};