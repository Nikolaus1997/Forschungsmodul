#include "storage/array1d.h"
#include <cassert>
#include "array1d.h"

Array1D::Array1D(int size):size_(size)
{
    // allocate data, initialize to 0
    data_.resize(size_, 0.0);
}

int Array1D::size() const
{
    return size_;
}

double &Array1D::operator()(int i)
{
    const int index = i;

    // assert that indices are in range
    assert(0 <= i && i < size_);
    return data_[index];
}

double Array1D::operator()(int i) const
{
    const int index = i;

    // assert that indices are in range
    assert(0 <= i && i < size_);
    return data_[index];
}

void Array1D::printValues()
{
    for (double i : data_)
    {
        std::cout << i <<' ';
    }
    std::cout<<';'<<std::endl;
    
}
