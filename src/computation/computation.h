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
#include "analyze/timer.cpp"
#include "storage/array3d.h"
#include "storage/array2d.h"
#include "limiter/limiter.h"
#include "limiter/projection.h"

#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>
#include <array>
#include <unistd.h>
#include <thread> 

/**
 * This class contains the main loop over all time steps of the simulation and all methods that are called in this loop.
*/

class Computation
{
    public:
        //initialize the computation object, parse the settings from file that is given as the only command line argument 
        void initialize(std::string filename);

        //run the whole simulation until tend 
        void runSimulation();

        void fillX();
        void calcQ(Array2D& u);
        void initVdm();
        void initVdmJ();
        void eulerTimeStep();
        void rungeKutta();
        void fillFaces();
        void calcDt();
        void calcError(double currenTime);
        void calcUdt(Array2D& u,Array2D& q,Array2D& j, Array2D& VdM_t, double epsilon = 1.0);
        void calcUdt(Array2D& u, Array2D& VdM_t);
        void firstLimiter(Array2D& u);
        void secondLimiter(Array2D& u);
        void thirdLimiter(Array2D &u);
        void fillXanalyze(Array2D &x);
        void fillUanalyze(Array2D& u_analyze,const Array2D& x,const Array2D& Vdm);
    
    private:
        std::array<double,1> meshWidth_;
        std::array<double,1> innerMeshWidth_;
        std::array<int, 1>  nCells_;
        Settings settings_;
        std::unique_ptr<Quadrature> quad_;
        std::shared_ptr<Grid> grid_;
        std::unique_ptr<OutputWriterParaview> outputWriterParaview_;     
        bool useLimiter_;  
        int firstLimiterCalls_,secondLimiterCalls_; 
        double dt_;    
        double b_,m_;
        double a_;
        double nNodes;
        double initCondA_;
        double initCondB_;
        double epsilon_;
        int PP_N_;
        Projection proj_;
        Flux flux_;
        Limiter limiter_;
        NumericalFlux gFlux_;
        InitialCondition initialCond_;
        std::shared_ptr<Vandermonde> VdM_;
};