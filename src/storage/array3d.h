#pragma once

#include <vector>
#include <array>

/**
 * This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 */
class Array3D
{
public:
    /**
    * constructor
    * @param size: number of cells
    */
    Array3D(std::array<int, 3> size);

    /**
    * get number of cells
    * @return number of cells
    */
    std::array<int, 3> size() const;

    /**
    * access the value at coordinate (i,j), declared not const, i.e. the value can be changed
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return reference to value at the grid cell (i,j)
    */
    double &operator()(int i, int j, int k);

    /**
    * get the value at coordinate (i,j), declared const, i.e. it is not possible to change the value
    * @param i: discretized position in x direction
    * @param j: discretized position in y direction
    * @return value at the grid cell (i,j)
    */
    double operator()(int i, int j, int k) const;

    void print() const;

protected:
    std::vector<double> data_;      //< storage array values, in row-major order
    const std::array<int, 3> size_; //< width, height of the domain
};