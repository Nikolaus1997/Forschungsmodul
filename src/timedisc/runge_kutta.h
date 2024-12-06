#pragma once

#include <vector>
#include <functional>
#include <storage/variable.h>

class LowStorageRungeKutta3 {
public:
    // Constructor
    LowStorageRungeKutta3(double dt);


    void advance(Variable u, Variable ut);

protected:
    double dt_; // Time step
};