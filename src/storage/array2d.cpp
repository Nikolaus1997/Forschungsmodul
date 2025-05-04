#include "storage/array2d.h"
#include <cassert>
#include "array2d.h"


Array2D::Array2D(std::array<int,2> size) :
  size_(size)
{
  // allocate data, initialize to 0
  data_.resize(size_[0]*size_[1], 0.0);
}

void Array2D::printValues()
{
  constexpr int width = 3; // width per number column

  for (int i = size_[1]-1; i >=0 ; --i)
  {
    std::cout << " [";
    for (int j = 0; j < size_[0]; ++j)
    {
      if(j!=0)
        std::cout << std::setw(width)<<"| " << operator()(j,i) << " ";
      else  
        std::cout << std::setw(width)<< operator()(j,i) << " ";
    }
    std::cout <<" ]"<< std::endl;
  }
}
//! get the size
std::array<int,2> Array2D::size() const
{
  return size_;
}


double &Array2D::operator()(int i, int j)
{
  const int index = j*size_[0] + i;

  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());

  return data_[index];
}


double Array2D::operator()(int i, int j) const

{
  const int index = j*size_[0] + i;

  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());

  return data_[index];

}