#include "computation.h"


//initialize the computation object, parse the settings from file that is given as the only command line argument 
void Computation::initialize(std::string filename)
{   
    //load settings
    settings_ = Settings();
    settings_.loadFromFile(filename);
    settings_.printSettings();

    a_ = settings_.physicalSize[0];
    b_ = settings_.physicalSize[1];
    initCondA_ = settings_.initCondA;
    initCondB_ = settings_.initCondB;

    PP_N_= settings_.PP_N;
    nNodes = 2*settings_.PP_N-1.0;

    nCells_= settings_.nCells;
    meshWidth_[0] =  (b_-a_)/nCells_[0];
    innerMeshWidth_[0] = meshWidth_[0]/nNodes;

    flux_ = Flux();
    initialCond_ = InitialCondition();
    //basis_ = std::make_unique<Basis>(PP_N_);


    grid_ = std::make_shared<Grid>(nCells_,meshWidth_, nNodes);

    if (settings_.PP_N==0)
    {
        std::cout<< "Cant choose polynomial degree of 0 due to Legendre Polynomial calculation"<<std::endl;
    }else
    {
        quad_ = std::make_unique<Quadrature>(nNodes);
        quad_->LegendreGaussNodesAndWeights(nNodes);
    }


    // quad_->basis_.weights_.printValues();
    // quad_->basis_.nodes_.printValues();

    // meshWidth_[0] = settings_.physicalSize[0]/settings_.nCells[0];
    // meshWidth_[1] = settings_.physicalSize[1]/settings_.nCells[1];
    
    //initialize discretization
    std::cout <<R"(============================================================================================)"<<std::endl;
    //initialize the flux function
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

    //initialize the initialCondition
    if (settings_.initialCondition == "unitStep"){
        initialCond_.setInitialCondType(InitialCondition::InitialCondType::UnitStep);
        std::cout<< "Choosing the unit step as initial condition..."<<std::endl;
    }
    else if (settings_.initialCondition == "negativeUnitStep")
    {
        initialCond_.setInitialCondType(InitialCondition::InitialCondType::NegativeUnitStep);
        std::cout<< "Choosing the negative unit step as initial condition..."<<std::endl;
    }else {
        std::cout << "Initial Condition not set, choosing default" << std::endl;
    }
    
    fillFaces();
    fillX();
    fillU();

    //grid_->x_.printValues();
    //grid_->faces_.printValues();
    //grid_->u_.printValues();
    //std::cout<<grid_->x(4)<<std::endl;
    //std::cout<<"Integral Value: "<<quad_->GaussLegendreQuad([&](double x) { return initialCond_.computeInitialCondition(grid_->x(3),-0,0); },a,b)<<std::endl;
    //initialize output writer
    // outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(grid_);
    outputWriterParaview_->writeFile(1);
}

void Computation::fillX()
{
    double transformedNode = 0.0;
    int numberNodes = quad_->basis_.nodes_.size()[0];
    double mean = 0.0, diff  =0.0;
    std::cout<<numberNodes<<std::endl;
    for (int i = 0; i < grid_->faces_.size()[0]-1; i++)
    {
        mean = 0.5 * (grid_->faces_(i+1) + grid_->faces_(i));
        diff = 0.5 * (grid_->faces_(i+1) - grid_->faces_(i));
        for(int j = 1; j <numberNodes-1;j++)
        {
            transformedNode =  diff* quad_->basis_.nodes(j) + mean;
            if(i==0)
                grid_->x(i*(numberNodes-2)+j)     = transformedNode;
            if(i!=0)
                grid_->x(i*(numberNodes-1)+j)     = transformedNode;
        }
        grid_->x(i*(numberNodes-1)) = grid_->faces_(i);
    }   
    grid_->x(grid_->x_.size()[0]-1) = grid_->faces_(grid_->faces_.size()[0]-1);
}

void Computation::fillU() {
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <= PP_N_; j++) {
            if (j == 0) {
                // Compute mean value of the initial condition
                grid_->u0(i) += 1 / meshWidth_[0] *
                    quad_->GaussLegendreQuad([&](double x) {
                        return initialCond_.computeInitialCondition(x, initCondA_, initCondB_);
                    }, grid_->faces(i), grid_->faces(i+1));
            } else {
                // Compute higher-order Legendre coefficients
                double mean = 0.5 * (grid_->faces_(i+1) + grid_->faces_(i));
                double diff = 0.5 * (grid_->faces_(i+1) - grid_->faces_(i));

                // Iterate over Gauss-Legendre nodes for accuracy
                for (int k = 0; k < quad_->basis_.nodes_.size()[0]; k++) {
                    double transformedNode = diff * quad_->basis_.nodes(k) + mean;
                    auto [L, L_prime] = quad_->LegendrePolynomialAndDerivative(2*j+1, transformedNode);
                    //TODO: this is not right
                    grid_->u0(i) += L * (2*j+1)/meshWidth_[0] *
                        quad_->GaussLegendreQuad([&](double x) {
                            return integralInit(x,2*j+1);
                        }, grid_->faces(i), grid_->faces(i+1));
                }
            }
        }
    }
    quad_->basis_.nodes_.printValues();
    std::cout << nNodes << std::endl;
}

void Computation::fillFaces()
{
    for (int i = 0; i < nCells_[0]+1; i++)
    {
        grid_->faces(i) = a_+i*meshWidth_[0];
    }
    
}

double Computation::integralInit(double x, int j)
{
    auto [L, L_prime] = quad_->LegendrePolynomialAndDerivative(j,x);
    //quad_->basis_.weights_.printValues();
    double init =  initialCond_.computeInitialCondition(x,initCondA_,initCondB_);
    //std::cout<<init << std::endl;
    return init*L;
}
