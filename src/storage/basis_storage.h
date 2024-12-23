#pragma once
#include <storage/variable.h>

class BasisStorage
{
public:
    //constructor
    BasisStorage(int N);
    
    const Array1D &weights() const;
    const Array1D &nodes() const;

    double weights(int i) const;
    double &weights(int i);

    double nodes(int i) const;
    double &nodes(int i);

    


    const int PP_N_;               // Polynomial degree of the basis 
    Array1D weights_;
    Array1D nodes_;
};