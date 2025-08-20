#include "array4d.h"

Array4D::Array4D(std::array<int, 4> size)
    : size_(size)
{
    data_.resize(size_[0] * size_[1] * size_[2] * size_[3], 0.0);
}

std::array<int,4> Array4D::size() const {
    return size_;
}

double &Array4D::operator()(int i, int j, int k, int l) {
    const int index = i + size_[0]*j + size_[0]*size_[1]*k + size_[0]*size_[1]*size_[2]*l;

    #ifndef NDEBUG
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(0 <= k && k < size_[2]);
    assert(0 <= l && l < size_[3]);
    assert(index < (int)data_.size());
    #endif

    return data_[index];
}

double Array4D::operator()(int i, int j, int k, int l) const {
    const int index = i + size_[0]*j + size_[0]*size_[1]*k + size_[0]*size_[1]*size_[2]*l;

    #ifndef NDEBUG
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(0 <= k && k < size_[2]);
    assert(0 <= l && l < size_[3]);
    assert(index < (int)data_.size());
    #endif

    return data_[index];
}

void Array4D::print() const {
    std::cout << "=== Array4D Contents ===\n";
    for(int l = 0; l < size_[3]; ++l) {
        std::cout << "Slice l = " << l << ":\n";
        for(int k = 0; k < size_[2]; ++k) {
            std::cout << "  Element k = " << k << ":\n";
            for(int i = 0; i < size_[0]; ++i) {
                for(int j = 0; j < size_[1]; ++j) {
                    std::cout << (*this)(i,j,k,l) << " ";
                }
                std::cout << "\n";
            }
            std::cout << "----------------------\n";
        }
    }
}
