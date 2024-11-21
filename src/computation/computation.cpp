#include "computation.h"


//initialize the computation object, parse the settings from file that is given as the only command line argument 
void Computation::initialize(std::string filename)
{   
    //load settings
    settings_ = Settings();

    flux_ = Flux();

    settings_.loadFromFile(filename);
    settings_.printSettings();
    if (settings_.PP_N==0)
    {
        
    }else
    {
        double nNodes = 2*settings_.PP_N-1.0;
        quad_ = std::make_unique<Quadrature>(nNodes);
        quad_->LegendreGaussNodesAndWeights(nNodes);
    }
    // quad_->basis_.weights_.printValues();
    // quad_->basis_.nodes_.printValues();

    // meshWidth_[0] = settings_.physicalSize[0]/settings_.nCells[0];
    // meshWidth_[1] = settings_.physicalSize[1]/settings_.nCells[1];
    
    //initialize discretization
    std::cout<<"==========================================================================================================================="<<std::endl;
    //initialize the flux function
    std::cout<<settings_.fluxFunction<<std::endl;
    if (settings_.fluxFunction == "linear") {
        flux_.setFluxFunction(Flux::FunctionType::Linear);
        std::cout<< "Choosing the linear flux function..."<<std::endl;
    }
    else if (settings_.fluxFunction == "burgers") {
        flux_.setFluxFunction(Flux::FunctionType::Burgers);
        std::cout<< "Choosing the burgers flux function..."<<std::endl;
    } 
    else if (settings_.fluxFunction == "buckley") {
        flux_.setFluxFunction(Flux::FunctionType::BuckleyLeverett);
        std::cout<< "Choosing the buckley flux function..."<<std::endl;
    } 
    else {
        std::cout << "flux function not found!" << std::endl;
    }
    double a = settings_.physicalSize[0];
    double b = settings_.physicalSize[1];
    std::cout<<"Integral Value: "<<quad_->GaussLegendreQuad([&](double x) { return flux_.compute(x); },a,b)<<std::endl;
    //initialize output writer
    // outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
    // outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
}




