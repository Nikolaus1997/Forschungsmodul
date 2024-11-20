#include "computation.h"


//initialize the computation object, parse the settings from file that is given as the only command line argument 
void Computation::initialize(std::string filename)
{   
    //load settings
    settings_ = Settings();

    settings_.loadFromFile(filename);
    settings_.printSettings();
    bas_ = std::make_unique<Basis>(settings_.PP_N);
    bas_->LegendreGaussNodesAndWeights(settings_.PP_N);
    std::cout<< "Nodes: " <<std::endl;
    bas_->basis_.nodes_.printValues();
    std::cout<< "Weights: " <<std::endl;
    bas_->basis_.weights_.printValues();
    // meshWidth_[0] = settings_.physicalSize[0]/settings_.nCells[0];
    // meshWidth_[1] = settings_.physicalSize[1]/settings_.nCells[1];
    
    //initialize discretization
    std::cout<<"==========================================================================================================================="<<std::endl;
    // if (settings_.useDonorCell) {
    //     std::cout<<"Using DonorCell..."<<std::endl;
    //     discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    // }
    // else {
    //     std::cout<<"Using CentralDifferences..."<<std::endl;
    //     discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    // }
    
    // //initialize the pressure solver
    // if (settings_.pressureSolver == "SOR") {
    //     pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon,
    //                                             settings_.maximumNumberOfIterations, settings_.omega);
    // }
    // else if (settings_.pressureSolver == "GaussSeidel") {
    //     pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon,
    //                                             settings_.maximumNumberOfIterations);

    // } else {
    //     std::cout << "Solver not found!" << std::endl;
    // }
    //initialize output writer
    // outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
    // outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
}

//run the whole simulation until tend 
// void Computation::runSimulation()
// {
//     int t_i = 0;
//     double time = 0.0;
//     //main loop of the simulation over time
//     while(time<settings_.endTime){
//         t_i++;
//         applyBoundaryValues();
//         applyBoundaryValuesFandG();;
//         computeTimeStepWidth();

//         if (time + dt_>settings_.endTime){
//             dt_ = settings_.endTime-time;
//         }
//         time = time + dt_;
//         computePreliminaryVelocities();
//         computeRightHandSide();
//         computePressure();
//         computeVelocities();
//         //outputWriterText_->writeFile(time);
//         //outputWriterParaview_->writeFile(time);
//     }
// }



