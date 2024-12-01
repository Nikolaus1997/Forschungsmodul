#pragma once

#include "settings/settings.h"
#include "integration/basis.h"
#include "storage/Vdm.h"
#include "dg/flux.h"
#include "integration/quad.h"
#include "dg/grid.h"
#include "computation/initial_condition.h"
#include "output_writer/output_writer_paraview.h"
#include "dg/num_flux.h"

#include <memory>
#include <cmath>
#include <algorithm>
#include <array>

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
        void fillX();
        void fillU();
        void fillUt();
        void initVdm();
        void eulerTimeStep();
        void fillFaces();
        void calcUdt();
        void fillNumFlux();
        double integralFlux( int i ,int j);
        double integralInit(double x, int j);
    
    private:
        std::array<double,1> meshWidth_;
        std::array<double,1> innerMeshWidth_;
        std::array<int, 1>  nCells_;
        Settings settings_;
        std::unique_ptr<Quadrature> quad_;
        std::shared_ptr<Grid> grid_;
        std::unique_ptr<OutputWriterParaview> outputWriterParaview_;        
        double dt_;    
        double b_;
        double a_;
        double nNodes;
        double initCondA_;
        double initCondB_;
        int PP_N_;
        Flux flux_;
        NumericalFlux gFlux_;
        InitialCondition initialCond_;
        std::shared_ptr<Vandermonde> VdM_;
};