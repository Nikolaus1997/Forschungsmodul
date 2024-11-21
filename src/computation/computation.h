#pragma once

#include "settings/settings.h"
#include "integration/basis.h"
#include "dg/flux.h"
#include "integration/quad.h"

#include <memory>
#include <cmath>
#include <algorithm>

/**
 * This class contains the main loop over all time steps of the simulation and all methods that are called in this loop.
*/

class Computation
{
    public:
        //initialize the computation object, parse the settings from file that is given as the only command line argument 
        void initialize(std::string filename);

        //run the whole simulation until tend 
        //void runSimulation();
    
    private:



        Settings settings_;
        std::unique_ptr<Quadrature> quad_;
        double dt_;    
        Flux flux_;

};