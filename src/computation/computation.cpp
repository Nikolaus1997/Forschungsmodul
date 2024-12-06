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


    nNodes = 2*PP_N_-1;

    nCells_= settings_.nCells;
    dt_ = settings_.CFL*1/(nCells_[0]);
    meshWidth_[0] =  (b_-a_)/nCells_[0];
    innerMeshWidth_[0] = meshWidth_[0]/nNodes;
    std::array<int,2> t  ={nCells_[0],PP_N_+1};
    VdM_ = std::make_shared<Vandermonde>(t,nNodes+2);

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


    quad_->basis_.weights_.printValues();
    quad_->basis_.nodes_.printValues();

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

    // initialize the numerical flux function
    if (settings_.RiemannSolver == "upwind") {
        gFlux_.setNumFluxFunction(NumericalFlux::FunctionType::upwind);
        std::cout<< "Choosing the numerical upwind flux function..."<<std::endl;
    }
    else if (settings_.RiemannSolver == "downwind") {
        gFlux_.setNumFluxFunction(NumericalFlux::FunctionType::downwind);
        std::cout<< "Choosing the numerical downwind flux function..."<<std::endl;
    } 
    else if (settings_.RiemannSolver == "enquist") {
        gFlux_.setNumFluxFunction(NumericalFlux::FunctionType::lax);
        std::cout<< "Choosing the numerical Lax Friedrichs flux function..."<<std::endl;
    } 
    else {
        std::cout << "Numerical flux function not found!" << std::endl;
        gFlux_.setNumFluxFunction(NumericalFlux::FunctionType::upwind);
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
    }else if (settings_.initialCondition == "sinus")
    {
        initialCond_.setInitialCondType(InitialCondition::InitialCondType::Sinus);
        std::cout<< "Choosing the sinus as initial condition..."<<std::endl;
    }else {
        std::cout << "Initial Condition not set, choosing default" << std::endl;
    }
    Timer timer;
    timer.start();
    double time_ = 0.0;
    int iter = 0.0;
    fillFaces();
    int numberN = 1/dt_*0.1;
    fillX();
    initVdm();
    fillU();

while (time_<settings_.endTime and iter<settings_.maximumNumberOfIterations)
    {
        calcUdt();
        fillUt();
        eulerTimeStep();
        //rungeKutta();
        if(time_==0){
        outputWriterParaview_ = std::make_unique<OutputWriterParaview>(grid_);
        outputWriterParaview_->writeFile(time_,settings_.OutputName);
        }

        //fillU();

        time_=time_+dt_;
        //std::cout<<"CalcTimestep: "<<dt_<<std::endl;
        if(iter % numberN ==0)
            outputWriterParaview_->writeFile(time_,settings_.OutputName);
        iter++;
    }

    

    // std::cout<<"x: "<<std::endl;
    // grid_->x_.printValues();
    // std::cout<<"faces: "<<std::endl;
    // grid_->faces_.printValues();
    // std::cout<<"u: "<<std::endl;
    // grid_->u_.printValues();
    // std::cout<<"Vdm: "<<std::endl;
    // VdM_->printValues();
    // std::cout<<"L: "<<std::endl;
    // VdM_->LprintValues();
    //std::cout<<grid_->x(4)<<std::endl;
    //std::cout<<"Integral Value: "<<quad_->GaussLegendreQuad([&](double x) { return initialCond_.computeInitialCondition(grid_->x(3),-0,0); },a,b)<<std::endl;
    //initialize output writer
    // outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
    timer.stop();
    std::cout << "Elapsed time: " << timer.elapsedMilliseconds()/1000 << " s." << " nStates: "<<iter<< std::endl;

}

void Computation::fillX()
{
    double transformedNode = 0.0;
    double mean = 0.0, diff  =0.0;

    for (int i = 0; i < grid_->faces_.size()[0]-1; i++)
    {
        mean = 0.5 * (grid_->faces_(i+1) + grid_->faces_(i));
        diff = 0.5 * (grid_->faces_(i+1) - grid_->faces_(i));
            grid_->x(i) =mean;
    }   
}

void Computation::fillU() {
    for(int i=0;i<VdM_->VdM_.size()[0];i++){
        grid_->u(i) =0.0;
        for(int j=0;j<VdM_->VdM_.size()[1];j++){
            double L = quad_->LegendrePolynomialAndDerivative(j,quad_->basis_.nodes((nNodes+2)/2))[0];
                grid_->u(i) +=VdM_->L_((nNodes+2)/2,j)*VdM_-> VdM_(i,j);
            }
    }

}

void Computation::fillUt()
{
        for(int i=0;i<VdM_->VdM_.size()[0];i++){
            grid_->ut(i) = 0.0;
        for(int j=0;j<VdM_->VdM_.size()[1];j++){
                grid_->ut(i) +=VdM_->L_(int((nNodes+2)/2),j) *VdM_-> VdM_t_(i,j);
            }
    }
}

void Computation::initVdm()
{
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM_.size()[1]; j++) {
            // Iterate over Gauss-Legendre nodes for accuracy
            VdM_->VdM_(i,j) = quad_->IntGaussLegendreQuad([&](double x) {
                                return initialCond_.computeInitialCondition(x,initCondA_,initCondB_);
                            },j ,grid_->faces(i), grid_->faces(i+1))*(2*j+1)/meshWidth_[0];  

            for (int p = 0; p <quad_->basis_.nodes_.size()[0]; p++) {     
                        double L = quad_->LegendrePolynomialAndDerivative(j,quad_->basis_.nodes(p))[0] ;
                        VdM_->L_(p,j) = L;

            }
        }
    }
}

void Computation::eulerTimeStep()
{
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM_.size()[1]; j++) {
            VdM_->VdM_(i,j) += dt_*VdM_->VdM_t_(i,j);
            // if(sqrt(VdM_->VdM_(i,j)*VdM_->VdM_(i,j))<1E-16)
            //     VdM_->VdM_(i,j)= 0.0;
        }

        //std::cout<<grid_->ut(i)<<dt_<<std::endl;
        grid_->u(i) += dt_*grid_->ut(i);
        // if(sqrt(grid_->u(i)*grid_->u(i))<1E-16)
        //     grid_->u(i) = 0.0;
    }
}

void Computation::rungeKutta()
{
    calcUdt();
    fillUt();
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM_.size()[1]; j++) {
            VdM_->VdM1_(i,j) =VdM_->VdM_(i,j)+ dt_*VdM_->VdM_t_(i,j);
            if(sqrt(VdM_->VdM1_(i,j)*VdM_->VdM1_(i,j))<1E-16)
                VdM_->VdM1_(i,j)= 0.0;
        }

        //std::cout<<grid_->ut(i)<<dt_<<std::endl;
        grid_->u1(i) = grid_->u(i)+ dt_*grid_->ut(i);
        if(sqrt(grid_->u1(i)*grid_->u1(i))<1E-12)
            grid_->u1(i) = 0.0;
    }
    calcU1dt();
    fillUt();
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM_.size()[1]; j++) {
            VdM_->VdM2_(i,j) = VdM_->VdM1_(i,j)+dt_*-3/4*VdM_->VdM_t_(i,j);
            if(sqrt(VdM_->VdM2_(i,j)*VdM_->VdM2_(i,j))<1E-16)
                VdM_->VdM2_(i,j)= 0.0;
        }

        //std::cout<<grid_->ut(i)<<dt_<<std::endl;
        grid_->u2(i) = grid_->u1(i)+dt_*-3/4*grid_->ut(i);
        if(sqrt(grid_->u2(i)*grid_->u2(i))<1E-16)
            grid_->u2(i) = 0.0;
    }
    calcU2dt();
    fillUt();
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM_.size()[1]; j++) {
            VdM_->VdM_(i,j) =VdM_->VdM2_(i,j)+ dt_*1/3*VdM_->VdM_t_(i,j);
            if(sqrt(VdM_->VdM_(i,j)*VdM_->VdM_(i,j))<1E-12)
                VdM_->VdM_(i,j)= 0.0;
        }

        //std::cout<<grid_->ut(i)<<dt_<<std::endl;
        grid_->u(i) = grid_->u2(i)+ dt_*1/3*grid_->ut(i);
        if(sqrt(grid_->u(i)*grid_->u(i))<1E-16)
            grid_->u(i) = 0.0;
    }
}

void Computation::fillFaces()
{
    for (int i = 0; i < nCells_[0]+1; i++)
    {
        grid_->faces(i) = a_+i*meshWidth_[0];
    }
    
}

void Computation::calcUdt()
{
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM_.size()[1]; j++) {
            double flux_term =0.0;
            // Wrap around the grid for periodic boundary conditions
            if(i==0){
                double ul_i = VdM_->VdM_(i,j)*VdM_->L_(0,j);
                double ur_i = VdM_->VdM_(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM_->VdM_(nCells_[0]-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM_->VdM_(i+1,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i) + gFlux_.computeNumFlux(ur_i, ul_iplus);//* pow(-1, j) ;
            }else if (i==grid_->faces_.size()[0] - 2)
            {
                double ul_i = VdM_->VdM_(i,j)*VdM_->L_(0,j);
                double ur_i = VdM_->VdM_(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM_->VdM_(i-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM_->VdM_(0,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i)  + gFlux_.computeNumFlux(ur_i, ul_iplus);//* pow(-1, j) ;
            }else{
            // Compute the numerical flux
                double ul_i = VdM_->VdM_(i,j)*VdM_->L_(0,j);
                double ur_i = VdM_->VdM_(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM_->VdM_(i-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM_->VdM_(i+1,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i) + gFlux_.computeNumFlux(ur_i, ul_iplus);//* pow(-1, j) ;
            }
            //std::cout<<"flux  "<<flux_term<<std::endl;
            // Apply the formula for the update of VdM_t_* pow(-1, j) 
            double integ =integralFlux(i, j);
            // if(sqrt(flux_term*flux_term)<1E-16)
            //     flux_term=0.0;
            // if(sqrt(integ*integ)<1E-12)
            //     integ=0.0;
            //std::cout<<"integral "<<integ<<" fluxTerm "<<flux_term<<" j: "<<j<<std::endl;
            VdM_->VdM_t_(i, j) = integ*1/(2 * j + 1)*meshWidth_[0]- flux_term*1/(2 * j + 1)*meshWidth_[0];
        }
    }
}

void Computation::calcU1dt()
{
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM1_.size()[1]; j++) {
            double flux_term =0.0;
            // Wrap around the grid for periodic boundary conditions
            if(i==0){
                double ul_i = VdM_->VdM1_(i,j)*VdM_->L_(0,j);
                double ur_i = VdM_->VdM1_(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM_->VdM1_(nCells_[0]-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM_->VdM1_(i+1,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i) + gFlux_.computeNumFlux(ur_i, ul_iplus)* pow(-1, j) ;
            }else if (i==grid_->faces_.size()[0] - 2)
            {
                double ul_i = VdM_->VdM1_(i,j)*VdM_->L_(0,j);
                double ur_i = VdM_->VdM1_(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM_->VdM1_(i-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM_->VdM1_(0,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i)  + gFlux_.computeNumFlux(ur_i, ul_iplus)* pow(-1, j) ;
            }else{
            // Compute the numerical flux
                double ul_i = VdM_->VdM1_(i,j)*VdM_->L_(0,j);
                double ur_i = VdM_->VdM1_(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM_->VdM1_(i-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM_->VdM1_(i+1,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i) + gFlux_.computeNumFlux(ur_i, ul_iplus)* pow(-1, j) ;
            }
            //std::cout<<"flux  "<<flux_term<<std::endl;
            // Apply the formula for the update of VdM_t_* pow(-1, j) 
            double integ =integralFlux(i, j);
            if(sqrt(flux_term*flux_term)<1E-16)
                flux_term=0.0;
            if(sqrt(integ*integ)<1E-12)
                integ=0.0;
            //std::cout<<"integral "<<integ<<" fluxTerm "<<flux_term<<" j: "<<j<<std::endl;
            VdM_->VdM_t_(i, j) = integ*1/(2 * j + 1)*meshWidth_[0]- flux_term*1/(2 * j + 1)*meshWidth_[0];
        }
    }
}

void Computation::calcU2dt()
{
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM2_.size()[1]; j++) {
            double flux_term =0.0;
            // Wrap around the grid for periodic boundary conditions
            if(i==0){
                double ul_i = VdM_->VdM2_(i,j)*VdM_->L_(0,j);
                double ur_i = VdM_->VdM2_(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM_->VdM2_(nCells_[0]-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM_->VdM2_(i+1,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i) + gFlux_.computeNumFlux(ur_i, ul_iplus)* pow(-1, j) ;
            }else if (i==grid_->faces_.size()[0] - 2)
            {
                double ul_i = VdM_->VdM2_(i,j)*VdM_->L_(0,j);
                double ur_i = VdM_->VdM2_(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM_->VdM2_(i-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM_->VdM2_(0,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i)  + gFlux_.computeNumFlux(ur_i, ul_iplus)* pow(-1, j) ;
            }else{
            // Compute the numerical flux
                double ul_i = VdM_->VdM2_(i,j)*VdM_->L_(0,j);
                double ur_i = VdM_->VdM2_(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM_->VdM2_(i-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM_->VdM2_(i+1,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i) + gFlux_.computeNumFlux(ur_i, ul_iplus)* pow(-1, j) ;
            }
            //std::cout<<"flux  "<<flux_term<<std::endl;
            // Apply the formula for the update of VdM_t_* pow(-1, j) 
            double integ =integralFlux(i, j);
            if(sqrt(flux_term*flux_term)<1E-16)
                flux_term=0.0;
            if(sqrt(integ*integ)<1E-12)
                integ=0.0;
            //std::cout<<"integral "<<integ<<" fluxTerm "<<flux_term<<" j: "<<j<<std::endl;
            VdM_->VdM_t_(i, j) = integ*1/(2 * j + 1)*meshWidth_[0]- flux_term*1/(2 * j + 1)*meshWidth_[0];
        }
    }
}

//TODO this is not RIGHT???
double Computation::integralFlux(int i ,int j)
{
    double a = grid_->faces(i),b = grid_->faces(i+1);
    // double point = 0.0;
    // for (int k = 0; k < VdM_->VdM_.size()[1]; k++)
    // {
    //     point += VdM_->VdM_(i,k)*quad_->LegendrePolynomialAndDerivative(j,quad_->basis_.nodes(k))[0];
    // }
    
    return      quad_->IntFluxGaussLegendreQuad([&](double x) {
        //0.5 * (b - a) * x + 0.5 * (b + a)
            return flux_.compute(x);//grid_->u(i)*quad_->LegendrePolynomialAndDerivative(j,x)[0]);
            },i,j ,a, b,VdM_->VdM());
}

double Computation::integralInit(double x, int j)
{
    //auto [L, L_prime] = quad_->LegendrePolynomialAndDerivative(j,x);
    //quad_->basis_.weights_.printValues();
    double init =  initialCond_.computeInitialCondition(x,initCondA_,initCondB_);
    //std::cout<<init << std::endl;
    return init;
}

