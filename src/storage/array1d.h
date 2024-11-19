#pragma once

#include <vector>
#include <array>

/** This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 */

class Array1D
{
public:
  //! constructor
  Array1D(int size);

  //! get the size
  int size() const;

  //! access the value at coordinate (i,j), declared not const, i.e. the value can be changed
  double &operator()(int i);

  //! get the value at coordinate (i,j), declared const, i.e. it is not possible to change the value
  double operator()(int i) const;

protected:

  std::vector<double> data_;  //< storage array values, in row-major order
  const int size_;    //< width, height of the domain
};