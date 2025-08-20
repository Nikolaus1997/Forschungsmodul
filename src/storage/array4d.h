#pragma once
#include <array>
#include <vector>
#include <iostream>
#include <cassert>

/**
 * @class Array4D
 * @brief A class representing a 4D array of doubles stored consecutively in memory.
 * 
 * The array can be accessed using four indices (i, j, k, l). Internally, data is
 * stored linearly in a std::vector<double>.
 */
class Array4D {
public:
    /**
     * @brief Construct a 4D array with given sizes in each dimension.
     * @param size Array of four integers representing the size in each dimension
     */
    Array4D(std::array<int, 4> size);

    /**
     * @brief Get the size of the array in each dimension.
     * @return std::array<int,4> with sizes
     */
    std::array<int,4> size() const;

    /**
     * @brief Access element (i,j,k,l) for modification.
     * @param i Index along first dimension
     * @param j Index along second dimension
     * @param k Index along third dimension
     * @param l Index along fourth dimension
     * @return Reference to the element at (i,j,k,l)
     */
    double &operator()(int i, int j, int k, int l);

    /**
     * @brief Access element (i,j,k,l) as const.
     * @param i Index along first dimension
     * @param j Index along second dimension
     * @param k Index along third dimension
     * @param l Index along fourth dimension
     * @return Value of the element at (i,j,k,l)
     */
    double operator()(int i, int j, int k, int l) const;

    /**
     * @brief Print the contents of the 4D array in a readable format.
     */
    void print() const;

private:
    std::array<int,4> size_;       ///< Sizes of the four dimensions
    std::vector<double> data_;     ///< Linear storage for the array
};

