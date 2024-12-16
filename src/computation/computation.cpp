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
        if(settings_.BarenblattM==0){
            //calcUdt(VdM_->VdM());
            //fillUt();
            //eulerTimeStep();
            rungeKutta5();
        }else{
            // calcQ(VdM_->VdM());
            // calcUdt(VdM_->VdM(),VdM_->VdMQ());
            // applyLimiter(VdM_->VdMt());
            rungeKutta();
            //eulerTimeStep();
            //
        }


        //fillU();

        time_=time_+dt_;
        //std::cout<<"CalcTimestep: "<<dt_<<std::endl;
        if(iter % numberN ==0){
            std::cout<<" TIME: "<<time_<<std::endl;
            outputWriterParaview_->writeFile(time_,settings_.OutputName);
        }

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
            for(int p=1; p<quad_->basis_.nodes_.size()[0]-1;p++)
            {
                grid_->u(i) +=VdM_->L_(p,j)*VdM_-> VdM_(i,j);
            }
        }
        grid_->u(i)*=1/(nNodes);
    }

}

void Computation::fillUt()
{
    for(int i=0;i<VdM_->VdM_.size()[0];i++){
            grid_->ut(i) = 0.0;
        for(int j=0;j<VdM_->VdM_.size()[1];j++){
            //double counter = 0.0;
            for(int p=1; p<quad_->basis_.nodes_.size()[0]-1;p++)
            {
                grid_->ut(i) +=VdM_->L_(p,j) *VdM_-> VdM_t_(i,j);
                //std::cout<<" IN FILL UT "<<VdM_->VdM_t_(i,j)<<std::endl;
                //counter += 1.0; 
            }
            //
        }
        grid_->ut(i) *= 1/nNodes;
    }
}

void Computation::calcQ(const Array2D& VdM)
{

    double m = double(settings_.BarenblattM);
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM_->VdM_.size()[1]; j++) {
            double flux_term =0.0;
            // Wrap around the grid for periodic boundary conditions
            if(i==0){
                double ul_i = VdM(i,j)*VdM_->L_(0,j);
                double ur_i = VdM(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM(nCells_[0]-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM(i+1,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,0.,0.,m,flux_,quad_)[1]* pow(-1, j) 
                                + gFlux_.computeNumFlux(ur_i, ul_iplus,0.,0.,m,flux_,quad_)[1] ;
            }else if (i==grid_->faces_.size()[0] - 2)
            {
                double ul_i = VdM(i,j)*VdM_->L_(0,j);
                double ur_i = VdM(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM(i-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM(0,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,0.,0.,m,flux_,quad_)[1]* pow(-1, j)
                                 + gFlux_.computeNumFlux(ur_i, ul_iplus,0.,0.,m,flux_,quad_)[1] ;
            }else{
            // Compute the numerical flux
                double ul_i = VdM(i,j)*VdM_->L_(0,j);
                double ur_i = VdM(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM(i-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM(i+1,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,0.,0.,m,flux_,quad_)[1]* pow(-1, j)
                                 + gFlux_.computeNumFlux(ur_i, ul_iplus,0.,0.,m,flux_,quad_)[1] ;
            }
            // Apply the formula for the update of VdM_t_* pow(-1, j) 
            double integ =integralQ(i,j,double(settings_.BarenblattM),VdM);
            // if(sqrt(flux_term*flux_term)<1E-12)
            //     flux_term=0.0;
            // if(sqrt(integ*integ)<1E-12)
            //     integ=0.0;
            VdM_->VdMQ_(i,j) = integ*(2 * j + 1)/meshWidth_[0]- flux_term*(2 * j + 1)/meshWidth_[0];
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
    return quad_->IntFluxQ(
        [&](double x) {
            auto result = flux_.compute(x, 0.0, m);
            if (result.size() <= 1) {
                throw std::runtime_error("flux_.compute did not return enough components");
            }
            return result[1];
        },
        i, j, a, b, VdM);
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
            VdM_->VdM_(i, j) = integral * (2 * j + 1) / meshWidth_[0];

            // Compute and store Legendre polynomials
            for (int p = 0; p < quad_->basis_.nodes_.size()[0]; p++) {
                double L = quad_->LegendrePolynomialAndDerivative(j, quad_->basis_.nodes(p))[0];
                VdM_->L_(p, j) = L;
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
        for (int j = 0; j < VdM_->VdM_.size()[1]; j++) {
            // u^(1) = u^n + Δt * L_h(u^n)
            VdM_->VdM1_(i, j) = VdM_->VdM_(i, j) + dt_ * VdM_->VdM_t_(i, j);
        }
    }

    // Update u1 for grid
    for (int i = 0; i < VdM_->VdM_.size()[0]; i++) {
        grid_->u1(i) = 0.0;
        for (int j = 0; j < VdM_->VdM_.size()[1]; j++) {
            for (int p = 1; p < quad_->basis_.nodes_.size()[0] - 1; p++) {
                grid_->u1(i) += VdM_->L_(p, j) * VdM_->VdM1_(i, j);
            }
        }
        grid_->u1(i) *= 1.0 / nNodes;
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
        for (int j = 0; j < VdM_->VdM_.size()[1]; j++) {
            // u^(2) = 3/4 * u^n + 1/4 * u^(1) + (1/4 * Δt) * L_h(u^(1))
            VdM_->VdM2_(i, j) = (3.0 / 4.0) * VdM_->VdM_(i, j) +
                                (1.0 / 4.0) * VdM_->VdM1_(i, j) +
                                (1.0 / 4.0) * dt_ * VdM_->VdM_t_(i, j);
        }
    }

    // Update u2 for grid
    for (int i = 0; i < VdM_->VdM_.size()[0]; i++) {
        grid_->u2(i) = 0.0;
        for (int j = 0; j < VdM_->VdM_.size()[1]; j++) {
            for (int p = 1; p < quad_->basis_.nodes_.size()[0] - 1; p++) {
                grid_->u2(i) += VdM_->L_(p, j) * VdM_->VdM2_(i, j);
            }
        }
        grid_->u2(i) *= 1.0 / nNodes;
    }

    // Step 2: Compute the intermediate stage u^(2)
    if(settings_.BarenblattM==0){
        calcUdt(VdM_->VdM2_); // Compute the time derivative for u^n
    }else
    {
        calcQ(VdM_->VdM1_);
        calcUdt(VdM_->VdM2_,VdM_->VdMQ());
    } // Compute the time derivative for u^(2)
    fillUt(); // Ensure all updates to time derivatives are reflected
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j < VdM_->VdM_.size()[1]; j++) {
            // u^(n+1) = 1/3 * u^n + 2/3 * u^(2) + (2/3 * Δt) * L_h(u^(2))
            VdM_->VdM_(i, j) = (1.0 / 3.0) * VdM_->VdM_(i, j) +
                               (2.0 / 3.0) * VdM_->VdM2_(i, j) +
                               (2.0 / 3.0) * dt_ * VdM_->VdM_t_(i, j);
        }
    }

    // Update grid's u for the final stage
    for (int i = 0; i < VdM_->VdM_.size()[0]; i++) {
        double u_temp = grid_->u(i);
        double u_next = 0.0;
        for (int j = 0; j < VdM_->VdM_.size()[1]; j++) {
            for (int p = 1; p < quad_->basis_.nodes_.size()[0] - 1; p++) {
                u_next += VdM_->L_(p, j) * VdM_->VdM_(i, j);
            }
        }
        u_next *= 1.0 / nNodes;

        // Final update: u^(n+1) = 1/3 * u^n + 2/3 * u^(2) + 2/3 * u_next
        grid_->u(i) =  u_next;
    }
}

void Computation::fillFaces()
{
    for (int i = 0; i < nCells_[0]+1; i++)
    {
        grid_->faces(i) = a_+i*meshWidth_[0];
    }
    
}

void Computation::calcUdt(const Array2D& VdM){
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
            for(int p=1; p<quad_->basis_.nodes_.size()[0]-1;p++){
                for(int k = 0;k<2;k++){
                    u_temp+=Vdm(i,k)*VdM_->L_(p,k);
                }
                u_temp*=1/nNodes;
            }
            grid_->ut(i) = u_temp;
        }
    }
}

void Computation::calcUdt(const Array2D& VdM,const Array2D& VdMQ)
{

        double m = double(settings_.BarenblattM);
        for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        for (int j = 0; j <VdM.size()[1]; j++) {
            double flux_term =0.0;
            // Wrap around the grid for periodic boundary conditions
            if(i==0){
                double ul_i = VdM(i,j)*VdM_->L_(0,j);
                double ur_i = VdM(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM(nCells_[0]-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM(i+1,j)*VdM_->L_(0,j);
                double ql_i = VdMQ(i,j)*VdM_->L_(0,j);
                double qr_i = VdMQ(i,j)*VdM_->L_(nNodes,j);
                double qr_iminus = VdMQ(nCells_[0]-1,j)*VdM_->L_(nNodes,j);
                double ql_iplus = VdMQ(i+1,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,qr_iminus,ql_i,m,flux_,quad_)[0]* pow(-1, j)  
                                    +gFlux_.computeNumFlux(ur_i, ul_iplus,qr_i,ql_iplus,m,flux_,quad_)[0] ;
            }else if (i==grid_->faces_.size()[0] - 2)
            {
                double ul_i = VdM(i,j)*VdM_->L_(0,j);
                double ur_i = VdM(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM(i-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM(0,j)*VdM_->L_(0,j);
                double ql_i = VdMQ(i,j)*VdM_->L_(0,j);
                double qr_i = VdMQ(i,j)*VdM_->L_(nNodes,j);
                double qr_iminus = VdMQ(i-1,j)*VdM_->L_(nNodes,j);
                double ql_iplus = VdMQ(0,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,qr_iminus,ql_i,m,flux_,quad_)[0]* pow(-1, j)  
                                    +gFlux_.computeNumFlux(ur_i, ul_iplus,qr_i,ql_iplus,m,flux_,quad_)[0] ;
            }else{
            // Compute the numerical flux
                double ul_i = VdM(i,j)*VdM_->L_(0,j);
                double ur_i = VdM(i,j)*VdM_->L_(nNodes,j);
                double ur_iminus  = VdM(i-1,j)*VdM_->L_(nNodes,j);
                double ul_iplus = VdM(i+1,j)*VdM_->L_(0,j);
                double ql_i = VdMQ(i,j)*VdM_->L_(0,j);
                double qr_i = VdMQ(i,j)*VdM_->L_(nNodes,j);
                double qr_iminus = VdMQ(i-1,j)*VdM_->L_(nNodes,j);
                double ql_iplus = VdMQ(i+1,j)*VdM_->L_(0,j);
                flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,qr_iminus,ql_i,m,flux_,quad_)[0]* pow(-1, j)  
                                    +gFlux_.computeNumFlux(ur_i, ul_iplus,qr_i,ql_iplus,m,flux_,quad_)[0] ;
            }
            // Apply the formula for the update of VdM_t_* pow(-1, j) 
            double integ =integralU(i,j,m,VdM,VdMQ);
            // if(sqrt(flux_term*flux_term)<1E-12)
            //     flux_term=0.0;
            // if(sqrt(integ*integ)<1E-12)
            //     integ=0.0;
            VdM_->VdM_t_(i,j) = integ*(2 * j + 1)/meshWidth_[0]- flux_term*(2 * j + 1)/meshWidth_[0];
        }
    }   
}



double Computation::integralFlux(int i ,int j,const Array2D& Vdm)
{
    double a = grid_->faces(i),b = grid_->faces(i+1);
    return      quad_->IntFluxGaussLegendreQuad([&](double x) {
                        return flux_.compute(x);
                             },i,j ,a, b,Vdm);
}

double Computation::integralInit(double x, int j)
{

    double init =  initialCond_.computeInitialCondition(x,initCondA_,initCondB_);
    
    return init;
}

void Computation::rungeKutta5()
{
    // Coefficients for 5-stage RK (Cash-Karp coefficients)
    constexpr double a[5][5] = {
        {0, 0, 0, 0, 0},
        {1.0/5.0, 0, 0, 0, 0},
        {3.0/40.0, 9.0/40.0, 0, 0, 0},
        {44.0/45.0, -56.0/15.0, 32.0/9.0, 0, 0},
        {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0, 0}
    };
    constexpr double b[5] = {9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, 5103.0/18656.0};
    constexpr double c[5] = {0, 1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0};

    // Temporary storage for each stage
    Array2D VdM_temp = Array2D(VdM_->VdM_.size());
    Array3D k = Array3D({5, VdM_->VdM_.size()[0], VdM_->VdM_.size()[1]}); // Stage derivatives
    Array1D u_temp = Array1D(grid_->u_.size());

    // Loop through stages
    for (int stage = 0; stage < 5; ++stage) {
        // Compute intermediate state
        for (int i = 0; i < grid_->faces_.size()[0] - 1; ++i) {
            for (int j = 0; j < VdM_->VdM_.size()[1]; ++j) {
                double increment = 0.0;
                for (int s = 0; s < stage; ++s) {
                    increment += a[stage][s] * k(s, i, j);
                }
                VdM_temp(i, j) = VdM_->VdM_(i, j) + dt_ * increment;
            }
            
            // Update u_temp correctly using the basis functions
            u_temp(i) = 0.0;
            for (int j = 0; j < VdM_->VdM_.size()[1]; ++j) {
                for (int p = 1; p < quad_->basis_.nodes_.size()[0] - 1; p++) {
                    u_temp(i) += VdM_->L_(p, j) * VdM_temp(i, j);
                }
            }
            u_temp(i) *= 1.0 / nNodes;
        }

        // Calculate time derivatives for this stage
        if(settings_.BarenblattM == 0) {
            calcUdt(VdM_temp);
        } else {
            calcQ(VdM_temp);
            calcUdt(VdM_temp, VdM_->VdMQ());
        }
        fillUt();

        // Store the derivatives for this stage
        for (int i = 0; i < grid_->faces_.size()[0] - 1; ++i) {
            for (int j = 0; j < VdM_->VdM_.size()[1]; ++j) {
                k(stage, i, j) = VdM_->VdM_t_(i, j);
            }
        }
    }

    // Final update
    for (int i = 0; i < grid_->faces_.size()[0] - 1; ++i) {
        for (int j = 0; j < VdM_->VdM_.size()[1]; ++j) {
            double update = 0.0;
            for (int stage = 0; stage < 5; ++stage) {
                update += b[stage] * k(stage, i, j);
            }
            VdM_->VdM_(i, j) += dt_ * update;
        }

        // Update u using the final VdM_ values and basis functions
        double u_next = 0.0;
        for (int j = 0; j < VdM_->VdM_.size()[1]; ++j) {
            for (int p = 1; p < quad_->basis_.nodes_.size()[0] - 1; p++) {
                u_next += VdM_->L_(p, j) * VdM_->VdM_(i, j);
            }
        }
        grid_->u(i) = u_next * (1.0 / nNodes);
    }
}