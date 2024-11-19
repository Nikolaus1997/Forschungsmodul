#pragma once
#include <storage/variable.h>

class BasisStorage
{
public:
    //constructor
    BasisStorage(int N);
    
    const Variable &weights() const;
    const Variable &nodes() const;

    double weights(int i) const;
    double &weights(int i);

    double nodes(int i) const;
    double &nodes(int i);

    

protected:
    const int PP_N;               // Polynomial degree of the basis 
    Variable weights_;
    Variable nodes_;
};