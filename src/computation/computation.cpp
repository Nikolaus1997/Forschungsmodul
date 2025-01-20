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


    nNodes = PP_N_;

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
    //TODO flesh out
    if(true){
        limiter_.setLimiterFunction(Limiter::FunctionType::minmod);
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
    else if (settings_.fluxFunction == "barenblatt") {
        flux_.setFluxFunction(Flux::FunctionType::Barenblatt);
        std::cout<< "Choosing the barenblatt flux function..."<<std::endl;
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
}

void Computation::runSimulation()
{
    std::cout<<'-'<<std::flush;
    Timer timer;
    timer.start();
    double time_ = 0.0;
    int iter = 0.0;
    fillFaces();
    fillX();
    initVdm();
    grid_->fillSolution(grid_->u(),VdM_);
    VdM_->LprintValues();
    VdM_->LprimePrintValues();
    if(settings_.BarenblattM!=0){
        dt_ = settings_.CFL*1/(nCells_[0]*nCells_[0])*1/double(settings_.BarenblattM);}
    else{
        dt_ = settings_.CFL*1/nCells_[0];
    }
    int numberN = 1/dt_*settings_.CFL;
    if(time_<dt_){
        outputWriterParaview_ = std::make_unique<OutputWriterParaview>(grid_);
        outputWriterParaview_->writeFile(time_,settings_.OutputName);
    }
    double errorTime =0.0;
    double numberofIterations = settings_.endTime/dt_;
    calcError(errorTime);
    while (time_<settings_.endTime)
        {
            if(settings_.BarenblattM==0){
                // calcUdt(VdM_->VdM_);
                // grid_->fillDerivative(grid_->ut(),VdM_);
                // eulerTimeStep();
                rungeKutta();
            }else{
                // calcQ(VdM_->VdM());
                // calcUdt(VdM_->VdM(),VdM_->VdMQ());
                // grid_->fillDerivative(grid_->ut_,VdM_);
                //applyLimiter(VdM_->VdMt());
                rungeKutta();
                //eulerTimeStep();
            }
            errorTime+=dt_;
            time_+=dt_;
            if(iter % numberN ==0){
                std::cout<<" TIME: "<<time_<<std::endl;
                outputWriterParaview_->writeFile(time_,settings_.OutputName);
                calcError(errorTime);
            }
            // calcDt();
            iter++;
            std::cout<<"\rCurrent Iteration: "<<iter<<" End Iter: "<<numberofIterations<< std::flush;
    }
    std::cout<<" TIME: "<<time_<<std::endl;
    outputWriterParaview_->writeFile(time_,settings_.OutputName);
    calcError(time_);
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
        for(int j = 1; j<=nNodes;j++){
            transformedNode = mean + diff*quad_->basis_.nodes(j);
            grid_->x(i*(nNodes)+j-1) = transformedNode;
        }
    }   
}

void Computation::calcDt(){
    double max = 0.0;
    double dx =0.0;
    double c = 0.0;
    for(int i=0;i<grid_->u_.size()[0];i++){
        double u = grid_->u(i);
        if(u==0)
            continue;
        dx = grid_->faces(i+1)-grid_->faces(i);
        c = flux_.compute(u,0.0,double(settings_.BarenblattM))[1]*flux_.compute(u,0.0,double(settings_.BarenblattM))[1];
        double dt = settings_.CFL*dx*dx/c*0.5;
        if(dt>max){
            max = dt;
        }
    }
    dt_ = max;
    //std::cout<<" DT "<<dt_<<std::endl;
}


void Computation::calcQ(const Array2D& VdM)
{

    double m = double(settings_.BarenblattM);
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM_.size()[1]; j++) {
            double flux_term =0.0;
            double ul_i = 0.0;
            double ur_i = 0.0;
            double ur_iminus  = 0.0;
            double ul_iplus =0.0;
            // Wrap around the grid for periodic boundary conditions
            if(i==0){
                for(int p=0; p<=PP_N_;p++){
                    ul_i += VdM(i,p)*VdM_->L_(0,p);
                    ur_i += VdM(i,p)*VdM_->L_(nNodes+1,p);
                    ur_iminus  += VdM(nCells_[0]-1,p)*VdM_->L_(nNodes+1,p);
                    ul_iplus += VdM(i+1,p)*VdM_->L_(0,p);
                }
            }else if (i==grid_->faces_.size()[0] - 2)
            {
                for(int p=0; p<=PP_N_;p++){
                    ul_i += VdM(i,p)*VdM_->L_(0,p);
                    ur_i += VdM(i,p)*VdM_->L_(nNodes+1,p);
                    ur_iminus  += VdM(i-1,p)*VdM_->L_(nNodes+1,p);
                    ul_iplus += VdM(0,p)*VdM_->L_(0,p);
                }
            }else{
            // Compute the numerical flux
                for(int p=0; p<=PP_N_;p++){
                    ul_i += VdM(i,p)*VdM_->L_(0,p);
                    ur_i += VdM(i,p)*VdM_->L_(nNodes+1,p);
                    ur_iminus  += VdM(i-1,p)*VdM_->L_(nNodes+1,p);
                    ul_iplus += VdM(i+1,p)*VdM_->L_(0,p);
                }
            }
            flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,0.,0.,m,flux_,quad_)[1]* pow(-1.0, double(j))
                    + gFlux_.computeNumFlux(ur_i, ul_iplus,0.,0.,m,flux_,quad_)[1] ;
            // Apply the formula for the update of VdM_t_* pow(-1, j) 
            double a = grid_->faces(i);
            double b = grid_->faces(i + 1);
            double integ =quad_->IntFluxQ([&](double x) {return flux_.compute(x, 0.0, m)[1];},i, j, a, b, VdM);
            // if(std::abs(flux_term)<1E-12)
            //     flux_term=0.0;
            // if(std::abs(integ)<1E-12)
            //     integ=0.0;
            VdM_->VdMQ_(i,j) = integ*(2.0 * double(j) + 1.0)/meshWidth_[0]- flux_term*(2.0 * double(j) + 1,0)/meshWidth_[0];
            //std::cout<<" IN CALCQ I "<<" Face i "<<grid_->faces(i)<<" Face i+1 "<<grid_->faces(i+1)<<" J "<<j<<" FLUX TERM "<<flux_term<<" INTEGRAL "<<integ<<" VDMQ "<<VdM_->VdMQ_(i,j)<<std::endl;
            //std::cout<<" IN CALCQ I " <<i<<" J "<<j<<" FLUX TERM "<<flux_term<<" INTEGRAL "<<integ<<" VDMQ "<<VdM_->VdMQ_(i,j)<<std::endl;
        }
    }
}

double Computation::integralQ(int i, int j, double m, const Array2D& VdM) {
    // Validate grid_ and indices
    if (!grid_) {
        throw std::runtime_error("grid_ is null");
    }
    if (i < 0 || i + 1 >= grid_->faces().size()[0]) {
        throw std::out_of_range("Index out of bounds for grid_->faces");
    }

    // Get integration bounds
    double a = grid_->faces(i);
    double b = grid_->faces(i + 1);

    // Validate quad_ pointer
    if (!quad_) {
        throw std::runtime_error("quad_ is null");
    }

    // Validate VdM dimensions
    if (VdM.size()[0] <= i || VdM.size()[1] <= j) {
        throw std::invalid_argument("VdM dimensions are invalid for indices i and j");
    }

    // Compute the integral
    return 0.0;
}

double Computation::integralU(int i, int j,double m, const Array2D &Vdm, const Array2D &VdmQ) {
    // Validate grid_ and indices
    if (!grid_) {
        throw std::runtime_error("grid_ is null");
    }
    if (i < 0 || i + 1 >= grid_->faces().size()[0]) {
        throw std::out_of_range("Index out of bounds for grid_->faces");
    }

    // Get integration bounds
    double a = grid_->faces(i);
    double b = grid_->faces(i + 1);

    // Validate quad_ pointer
    if (!quad_) {
        throw std::runtime_error("quad_ is null");
    }

    // Validate Vdm and VdmQ dimensions if applicable
    if (Vdm.size()[0] <= i || VdmQ.size()[0] <= i || Vdm.size()[1] <= j || VdmQ.size()[1] <= j) {
        throw std::invalid_argument("Vdm or VdmQ dimensions are invalid for indices i and j");
    }

    // Compute the integral
    return quad_->IntFluxU([&](double u, double q) { return flux_.compute(u, q, m)[0]; }, i, j, a, b, Vdm, VdmQ);
}

void Computation::initVdm() {
    // Validate critical pointers and sizes
    if (!grid_ || !VdM_ || !quad_) {
        throw std::runtime_error("grid_, VdM_, or quad_ is null");
    }
    if (grid_->faces_.size()[0] < 2) {
        throw std::invalid_argument("Insufficient grid faces for computation");
    }
    if (meshWidth_.empty() || meshWidth_[0] <= 0) {
        throw std::invalid_argument("Invalid meshWidth values");
    }

    // Iterate over grid cells and polynomial degrees
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j < VdM_->VdM_.size()[1]; j++) {
            // Determine integral lambda function
            double integral = 0.0;
            if (InitialCondition::InitialCondType::Barenblatt == initialCond_.selectedFunction) {
                integral = quad_->IntGaussLegendreQuad(
                    [&](double x) {
                        return initialCond_.computeInitialCondition(
                            x, initCondA_, initCondB_,
                            settings_.BarenblattTime, settings_.BarenblattM);
                    },
                    j, grid_->faces(i), grid_->faces(i + 1));
            } else {
                integral = quad_->IntGaussLegendreQuad(
                    [&](double x) {
                        return initialCond_.computeInitialCondition(x, initCondA_, initCondB_);
                    },
                    j, grid_->faces(i), grid_->faces(i + 1));
            }

            // Scale integral and store in VdM_
            VdM_->VdM_(i, j) = integral * (2.0 * double(j) + 1.0) / meshWidth_[0];

            // Compute and store Legendre polynomials
            for (int p = 0; p < quad_->basis_.nodes_.size()[0]; p++) {
                std::array<double,2> L = quad_->LegendrePolynomialAndDerivative(j, quad_->basis_.nodes(p));
                VdM_->L(p, j) = L[0];
                VdM_->L_prime(p, j) = L[1];
            }
        }
    }
}


void Computation::eulerTimeStep()
{
    for (int i = 0; i < VdM_->VdM_.size()[0]; i++) {
        for(int p = 0; p <=PP_N_; p++){
            VdM_->VdM_(i,p) += dt_*VdM_->VdM_t_(i,p);
        }
    }
    grid_->fillDerivative(grid_->ut(),VdM_);
    grid_->fillSolution(grid_->u(),VdM_);
}

void Computation::rungeKutta() {
    // Step 1: Compute the intermediate stage u^(1)
    if(settings_.BarenblattM==0){
        calcUdt(VdM_->VdM_); // Compute the time derivative for u^n
    }else
    {
        calcQ(VdM_->VdM_);
        calcUdt(VdM_->VdM(),VdM_->VdMQ());
    }
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <=PP_N_; j++) {
            // u^(1) = u^n + Δt * L_h(u^n)
            VdM_->VdM1_(i, j) = VdM_->VdM_(i, j) + dt_ * VdM_->VdM_t_(i, j);
        }
    }
    // Step 2: Compute the intermediate stage u^(2)
    if(settings_.BarenblattM==0){
        calcUdt(VdM_->VdM1_); // Compute the time derivative for u^n
    }else
    {
        calcQ(VdM_->VdM1_);
        calcUdt(VdM_->VdM1_,VdM_->VdMQ());
    }// Compute the time derivative for u^(1)
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <=PP_N_; j++) {
            // u^(2) = 3/4 * u^n + 1/4 * u^(1) + (1/4 * Δt) * L_h(u^(1))
            VdM_->VdM2_(i, j) = (3.0 / 4.0) * VdM_->VdM_(i, j) +
                                (1.0 / 4.0) * VdM_->VdM1_(i, j) +
                                (1.0 / 4.0) * dt_ * VdM_->VdM_t_(i, j);
        }
    }

    // Step 2: Compute the intermediate stage u^(2)
    if(settings_.BarenblattM==0){
        calcUdt(VdM_->VdM2_); // Compute the time derivative for u^n
    }else
    {
        calcQ(VdM_->VdM1_);
        calcUdt(VdM_->VdM2_,VdM_->VdMQ());
    } // Compute the time derivative for u^(2)
    // Ensure all updates to time derivatives are reflected
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <=PP_N_; j++) {
            // u^(n+1) = 1/3 * u^n + 2/3 * u^(2) + (2/3 * Δt) * L_h(u^(2))
            VdM_->VdM_(i, j) = (1.0 / 3.0) * VdM_->VdM_(i, j) +
                               (2.0 / 3.0) * VdM_->VdM2_(i, j) +
                               (2.0 / 3.0) * dt_ * VdM_->VdM_t_(i, j);
        }
    }
    grid_->fillDerivative(grid_->ut(),VdM_);
    grid_->fillSolution(grid_->u(),VdM_);
}

void Computation::fillFaces()
{
    for (int i = 0; i < nCells_[0]+1; i++)
    {
        grid_->faces(i) = a_+i*meshWidth_[0];
    }
    
}

void Computation::calcError(double currentTime)
{

    if(settings_.BarenblattM!=0){
        grid_->l2_error(0) = 0.0;
        for(int i =0;i<grid_->u_.size()[0];i++){
                grid_->l2_error(0) += pow((grid_->u(i)-initialCond_.computeInitialCondition(grid_->x(i),initCondA_,initCondB_,currentTime+1.0, settings_.BarenblattM)),2);        }
        
        grid_->linf_error(0) = 0.0;
        for(int i =0;i<grid_->u_.size()[0];i++){
            if(grid_->linf_error(0)<fabs(grid_->u(i)-initialCond_.computeInitialCondition(grid_->x(i),initCondA_,initCondB_,currentTime+1.0, settings_.BarenblattM))and fabs(grid_->x(i))>=1.5) 
                grid_->linf_error(0) = fabs(grid_->u(i)-initialCond_.computeInitialCondition(grid_->x(i),initCondA_,initCondB_,currentTime+1.0, settings_.BarenblattM));
        }
    }else{
        grid_->l2_error(0) = 0.0;
        for(int i =0;i<grid_->u_.size()[0];i++){
            grid_->l2_error(0) += pow(grid_->u(i)-sin(grid_->x(i)-currentTime),2.0);//initialCond_.computeInitialCondition(grid_->x(i),initCondA_,initCondB_),2.0);
                                
        }
        
        grid_->linf_error(0) = 0.0;
        for(int i =0;i<grid_->u_.size()[0];i++){
            if(grid_->linf_error(0)<fabs(grid_->u(i)-initialCond_.computeInitialCondition(grid_->x(i),initCondA_,initCondB_))) 
                grid_->linf_error(0) = fabs(grid_->u(i)-initialCond_.computeInitialCondition(grid_->x(i),initCondA_,initCondB_));
        }
    }

    std::cout<<"L2 Error: "<<sqrt(grid_->l2_error(0))<<" Linf Error: "<<grid_->linf_error(0)<<std::endl;
}

void Computation::calcUdt(const Array2D& VdM){
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <=PP_N_; j++) {
            double ul_i = 0.0;
            double ur_i = 0.0;
            double ur_iminus  =0.0;
            double ul_iplus = 0.0; 
            double flux_term =0.0;
            if(i==0){
            for(int p=0; p<=PP_N_;p++){
                ul_i += VdM(i,p)*VdM_->L_(0,p);
                ur_i += VdM(i,p)*VdM_->L_(nNodes+1,p);
                ur_iminus += VdM(nCells_[0]-1,p)*VdM_->L_(nNodes+1,p);
                ul_iplus += VdM(i+1,p)*VdM_->L_(0,p);
                }
            }else if (i==grid_->faces_.size()[0] - 2)
            {
                for(int p=0; p<=PP_N_;p++){
                    ul_i += VdM(i,p)*VdM_->L_(0,p);
                    ur_i += VdM(i,p)*VdM_->L_(nNodes+1,p);
                    ur_iminus += VdM(i-1,p)*VdM_->L_(nNodes+1,p);
                    ul_iplus += VdM(0,p)*VdM_->L_(0,p);
                }
            }else{
            //Compute the numerical flux
            for(int p=0; p<=PP_N_;p++){
                ul_i += VdM(i,p)*VdM_->L_(0,p);
                ur_i += VdM(i,p)*VdM_->L_(nNodes+1,p);
                ur_iminus += VdM(i-1,p)*VdM_->L_(nNodes+1,p);
                ul_iplus += VdM(i+1,p)*VdM_->L_(0,p);
                }
            }
            // Wrap around the grid for periodic boundary conditions
            flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,flux_)*VdM_->L_(0,j)  + gFlux_.computeNumFlux(ur_i, ul_iplus,flux_);
            double integ =quad_->IntFluxGaussLegendreQuad([&](double x) {return flux_.compute(x);}
                                                            ,i,j ,grid_->faces(i),grid_->faces(i+1),VdM);
            // if(std::abs(flux_term)<1E-12)
            //     flux_term=0.0;
            // if(std::abs(integ)<1E-12)
            //     integ=0.0;
            VdM_->VdM_t_(i,j) =integ*(2.0*double(j)+1.0)*1/meshWidth_[0]
                                - flux_term*1/meshWidth_[0]*(2.0*double(j)+1.0);
        }
    }
}

void Computation::applyLimiter(const Array2D &Vdm)
{
    for(int i=0;i<grid_->u_.size()[0];i++){
        double limit_l=0.0, limit_r=0.0;
        double u_r =0.0, u_l =0.0;
        if(i==0){
            for(int j = 0;j<VdM_->VdM_.size()[1];j++){
                u_l = Vdm(i,j)*VdM_->L_(0,j);
                u_r = Vdm(i,j)*VdM_->L_(nNodes,j);
            }
            limit_l = grid_->ut(i)-limiter_.computeLimiter(grid_->ut(i)-u_l,grid_->u(i)-grid_->ut(nCells_[0]-1),grid_->ut(i+1)-grid_->ut(i),meshWidth_[0]);
            limit_r = grid_->ut(i)+limiter_.computeLimiter(u_r-grid_->ut(i),grid_->u(i)-grid_->ut(nCells_[0]-1),grid_->ut(i+1)-grid_->ut(i),meshWidth_[0]);
        }else if(i==nCells_[0]-1){
            for(int j = 0;j<VdM_->VdM_.size()[1];j++){
                u_l = Vdm(i,j)*VdM_->L_(0,j);
                u_r = Vdm(i,j)*VdM_->L_(nNodes,j);
            }
            limit_l = grid_->ut(i)-limiter_.computeLimiter(grid_->ut(i)-u_l,grid_->ut(i)-grid_->ut(nCells_[0]-1),grid_->ut(0)-grid_->ut(i),meshWidth_[0]);
            limit_r = grid_->ut(i)+limiter_.computeLimiter(u_r-grid_->ut(i),grid_->ut(i)-grid_->ut(nCells_[0]-1),grid_->ut(0)-grid_->ut(i),meshWidth_[0]);
        }else{
            for(int j = 0;j<VdM_->VdM_.size()[1];j++){
                u_l = Vdm(i,j)*VdM_->L_(0,j);
                u_r = Vdm(i,j)*VdM_->L_(nNodes,j);
            }
            limit_l = grid_->ut(i)-limiter_.computeLimiter(grid_->ut(i)-u_l,grid_->ut(i)-grid_->ut(i-1),grid_->ut(i+1)-grid_->ut(i),meshWidth_[0]);
            limit_r = grid_->ut(i)+limiter_.computeLimiter(u_r-grid_->ut(i),grid_->ut(i)-grid_->ut(i-1),grid_->ut(i+1)-grid_->ut(i),meshWidth_[0]);
        }
        if(limit_l!=u_l or limit_r!=u_r){
            double u_temp=0.0;
            for(int k = 0;k<2;k++){
            for(int p=1; p<quad_->basis_.nodes_.size()[0]-1;p++){
                    u_temp+=Vdm(i,k)*VdM_->L_(p,k);
                }
            }
            grid_->ut(i) = u_temp;
        }
    }
}

void Computation::calcUdt(const Array2D& VdM,const Array2D& VdMQ)
{

    double m = double(settings_.BarenblattM);
    double integ=0.0;
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM.size()[1]; j++) {
            double flux_term =0.0;
            double ul_i =0.0;
            double ur_i =0.0;
            double ur_iminus  = 0.0;
            double ul_iplus = 0.0;
            double ql_i =0.0;
            double qr_i = 0.0;
            double qr_iminus =0.0;
            double ql_iplus = 0.0;
            // Wrap around the grid for periodic boundary conditions
            if(i==0){
                for(int p = 0;p<=PP_N_;p++){
                    ul_i += VdM(i,p)*VdM_->L_(0,p);
                    ur_i += VdM(i,p)*VdM_->L_(nNodes+1,p);
                    ur_iminus  += VdM(nCells_[0]-1,p)*VdM_->L_(nNodes+1,p);
                    ul_iplus += VdM(i+1,p)*VdM_->L_(0,p);
                    ql_i += VdMQ(i,p)*VdM_->L_(0,p);
                    qr_i += VdMQ(i,p)*VdM_->L_(nNodes+1,p);
                    qr_iminus += VdMQ(nCells_[0]-1,p)*VdM_->L_(nNodes+1,p);
                    ql_iplus += VdMQ(i+1,p)*VdM_->L_(0,p);
                }
            }else if (i==grid_->faces_.size()[0] - 2)
            {
            for(int p = 0;p<=PP_N_;p++){
                    ul_i += VdM(i,p)*VdM_->L_(0,p);
                    ur_i += VdM(i,p)*VdM_->L_(nNodes+1,p);
                    ur_iminus  += VdM(i-1,p)*VdM_->L_(nNodes+1,p);
                    ul_iplus += VdM(0,p)*VdM_->L_(0,p);
                    ql_i += VdMQ(i,p)*VdM_->L_(0,p);
                    qr_i += VdMQ(i,p)*VdM_->L_(nNodes+1,p);
                    qr_iminus += VdMQ(i-1,p)*VdM_->L_(nNodes+1,p);
                    ql_iplus += VdMQ(0,p)*VdM_->L_(0,p);
            }                                                      
            }else{
            // Compute the numerical flux
                for(int p = 0;p<=PP_N_;p++){            
                    ul_i += VdM(i,p)*VdM_->L_(0,p);
                    ur_i += VdM(i,p)*VdM_->L_(nNodes+1,p);
                    ur_iminus  += VdM(i-1,p)*VdM_->L_(nNodes+1,p);
                    ul_iplus += VdM(i+1,p)*VdM_->L_(0,p);
                    ql_i += VdMQ(i,p)*VdM_->L_(0,p);
                    qr_i += VdMQ(i,p)*VdM_->L_(nNodes+1,p);
                    qr_iminus += VdMQ(i-1,p)*VdM_->L_(nNodes+1,p);
                    ql_iplus += VdMQ(i+1,p)*VdM_->L_(0,p);
                }
            }
            flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,qr_iminus,ql_i,m,flux_,quad_)[0]* pow(-1.0, double(j))
                                            +gFlux_.computeNumFlux(ur_i, ul_iplus,qr_i,ql_iplus,m,flux_,quad_)[0]   ;
            integ =quad_->IntFluxU([&](double u, double q) { return flux_.compute(u, q, m)[0]; }, i, j,
                                            grid_->faces(i), grid_->faces(i+1), VdM, VdMQ);
            // Apply the formula for the update of VdM_t_* pow(-1, j) 
            // if(std::abs(flux_term)<1E-12)
            //     flux_term=0.0;
            // if(std::abs(integ)<1E-12)
            //     integ=0.0;
            VdM_->VdM_t_(i,j) = integ*(2.0 * double(j) + 1.0)/meshWidth_[0]- flux_term*(2.0 * double(j) + 1.0)/meshWidth_[0];
        }
    }   
}




double Computation::integralInit(double x, int j)
{

    double init =  initialCond_.computeInitialCondition(x,initCondA_,initCondB_);
    
    return init;
}

