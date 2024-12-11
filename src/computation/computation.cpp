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
    else if (settings_.RiemannSolver == "lax") {
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
    }else if (settings_.initialCondition == "barenblatt")
    {
        initialCond_.setInitialCondType(InitialCondition::InitialCondType::Barenblatt);
        std::cout<< "Choosing the barenblatt function as initial condition..."<<std::endl;        
    }else {
        std::cout << "Initial Condition not set, choosing default" << std::endl;
    }
    std::cout<<'-'<<std::flush;
    Timer timer;
    timer.start();
    double time_ = 0.0;
    int iter = 0.0;
    fillFaces();
    int numberN = 1/dt_*0.1;
    fillX();
    initVdm();
    fillU();
    double numberofIterations = settings_.endTime/dt_;
    if(time_<dt_){
        outputWriterParaview_ = std::make_unique<OutputWriterParaview>(grid_);
        outputWriterParaview_->writeFile(time_,settings_.OutputName);
    }

while (time_<settings_.endTime and iter<settings_.maximumNumberOfIterations)
    {
         calcUdt(VdM_->VdM());
         fillUt();
         eulerTimeStep();
        //rungeKutta5();

        //fillU();

        time_=time_+dt_;
        //std::cout<<"CalcTimestep: "<<dt_<<std::endl;
        if(iter % numberN ==0)
            outputWriterParaview_->writeFile(time_,settings_.OutputName);
        iter++;
        std::cout<<"\rCurrent Iteration: "<<iter<<" End Iter: "<<numberofIterations<< std::flush;
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
                //std::cout<<" IN FILL UT "<<VdM_->VdM_t_(i,j)<<std::endl;
            }
    }
}

void Computation::initVdm()
{
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM_.size()[1]; j++) {
            // Iterate over Gauss-Legendre nodes for accuracy
            if(InitialCondition::InitialCondType::Barenblatt==initialCond_.selectedFunction){
                VdM_->VdM_(i,j) = quad_->IntGaussLegendreQuad([&](double x) {
                            return initialCond_.computeInitialCondition(x,initCondA_,initCondB_,settings_.BarenblattTime,settings_.BarenblattM);
                            },j ,grid_->faces(i), grid_->faces(i+1))*(2*j+1)/meshWidth_[0]; 
            }else{
                VdM_->VdM_(i,j) = quad_->IntGaussLegendreQuad([&](double x) {
                            return initialCond_.computeInitialCondition(x,initCondA_,initCondB_);
                            },j ,grid_->faces(i), grid_->faces(i+1))*(2*j+1)/meshWidth_[0];  

            }


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
    calcUdt(VdM_->VdM_);
    fillUt();
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM_.size()[1]; j++) {
            VdM_->VdM1_(i,j) =VdM_->VdM_(i,j)+ dt_*VdM_->VdM_t_(i,j);
        }

        //std::cout<<grid_->ut(i)<<dt_<<std::endl;
        grid_->u1(i) = grid_->u(i)+ dt_*grid_->ut(i);
    }
    calcUdt(VdM_->VdM1_);
    fillUt();
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM_.size()[1]; j++) {
            VdM_->VdM2_(i,j) = VdM_->VdM1_(i,j)-dt_*3/4*VdM_->VdM_t_(i,j);
        }

        //std::cout<<grid_->ut(i)<<dt_<<std::endl;
        grid_->u2(i) = grid_->u1(i)+dt_*-3/4*grid_->ut(i);
    }
    calcUdt(VdM_->VdM2_);
    fillUt();
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM_.size()[1]; j++) {
            VdM_->VdM_(i,j) =VdM_->VdM2_(i,j)+ dt_*1/3*VdM_->VdM_t_(i,j);
        }

        //std::cout<<grid_->ut(i)<<dt_<<std::endl;
        grid_->u(i) = grid_->u2(i)+ dt_*1/3*grid_->ut(i);
    }
}

void Computation::fillFaces()
{
    for (int i = 0; i < nCells_[0]+1; i++)
    {
        grid_->faces(i) = a_+i*meshWidth_[0];
    }
    
}

void Computation::calcUdt(Array2D VdM)
{
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM.size()[1]; j++) {
            double flux_term =0.0;
            // Wrap around the grid for periodic boundary conditions
            if(i==0){
                double ul_i = VdM(i,j)*VdM_->L_(0,j);
                double ur_i = VdM(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM(nCells_[0]-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM(i+1,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,flux_)* pow(-1, j) + gFlux_.computeNumFlux(ur_i, ul_iplus,flux_) ;
            }else if (i==grid_->faces_.size()[0] - 2)
            {
                double ul_i = VdM(i,j)*VdM_->L_(0,j);
                double ur_i = VdM(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM(i-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM(0,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,flux_)* pow(-1, j)  + gFlux_.computeNumFlux(ur_i, ul_iplus,flux_);//* pow(-1, j) ;
            }else{
            // Compute the numerical flux
                double ul_i = VdM(i,j)*VdM_->L_(0,j);
                double ur_i = VdM(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM(i-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM(i+1,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,flux_)* pow(-1, j) + gFlux_.computeNumFlux(ur_i, ul_iplus,flux_);//* pow(-1, j) ;
            }
            // Apply the formula for the update of VdM_t_* pow(-1, j) 
            double integ =integralFlux(i, j,VdM);
            // if(sqrt(flux_term*flux_term)<1E-12)
            //     flux_term=0.0;
            // if(sqrt(integ*integ)<1E-12)
            //     integ=0.0;
            VdM_->VdM_t_(i,j) = integ*(2 * j + 1)/meshWidth_[0]- flux_term*(2 * j + 1)/meshWidth_[0];
        }
    }
}

//TODO this is not RIGHT???
double Computation::integralFlux(int i ,int j,Array2D VdM)
{
    double a = grid_->faces(i),b = grid_->faces(i+1);
    
    return      quad_->IntFluxGaussLegendreQuad([&](double x) {
                        return flux_.compute(x);
                             },i,j ,a, b,VdM);
}

double Computation::integralInit(double x, int j)
{
    //auto [L, L_prime] = quad_->LegendrePolynomialAndDerivative(j,x);
    //quad_->basis_.weights_.printValues();
    double init =  initialCond_.computeInitialCondition(x,initCondA_,initCondB_);
    //std::cout<<init << std::endl;
    return init;
}

void Computation::rungeKutta5()
{
    // Coefficients for 5-stage RK
    constexpr double a[5][5] = {
        {0, 0, 0, 0, 0},
        {1.0 / 5.0, 0, 0, 0, 0},
        {3.0 / 40.0, 9.0 / 40.0, 0, 0, 0},
        {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0, 0},
        {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0}
    };
    constexpr double b[5] = {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, 5103.0 / 18656.0};
    constexpr double c[5] = {0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0};

    // Temporary storage for each stage
    Array2D VdM_temp = Array2D(VdM_->VdM_.size());
    Array3D k = Array3D({5,VdM_->VdM_.size()[0], VdM_->VdM_.size()[1]}); // Stage derivatives
    Array1D u_temp = Array1D(grid_->u_.size());    

    // Loop through stages
    for (int stage = 0; stage < 5; ++stage) {
        // Compute intermediate state
        for (int i = 0; i < grid_->faces_.size()[0] - 1; ++i) {
            for (int j = 0; j < VdM_->VdM_.size()[1]; ++j) {
                double increment = 0.0;
                for (int s = 0; s < stage; ++s) {
                    increment += dt_ * a[stage][s] * k(s, i, j);  // Use k for stage derivatives
                }
                VdM_temp(i, j) = VdM_->VdM_(i, j) + increment;
            }
            double increment_u = 0.0;
            for (int s = 0; s < stage; ++s) {
                increment_u += dt_ * a[stage][s] * k(s, i, 0);  // Use k for u_t if applicable
            }
            u_temp(i) = grid_->u(i) + increment_u;
        }
        calcUdt(VdM_temp);
        fillUt();  // Ensure these are updating time derivatives
        for (int i = 0; i < grid_->faces_.size()[0] - 1; ++i) {
            for (int j = 0; j < VdM_->VdM_.size()[1]; ++j) {
                k(stage, i, j) = VdM_->VdM_t_(i, j);  // Store the time derivative for this stage
            }
            grid_->ut(i) = grid_->ut(i);  // Store u's time derivative appropriately
        }
    }

    for (int i = 0; i < grid_->faces_.size()[0] - 1; ++i) {
        for (int j = 0; j < VdM_->VdM_.size()[1]; ++j) {
            double update = 0.0;
            for (int stage = 0; stage < 5; ++stage) {
                update += b[stage] * k(stage, i, j);  // Use stage time derivatives
            }
            VdM_->VdM_(i, j) += dt_ * update;  // Final update
        }
        double update_u = 0.0;
        for (int stage = 0; stage < 5; ++stage) {
            update_u += b[stage] * grid_->ut(i);
        }
        grid_->u(i) += dt_ * update_u;
    }

}
