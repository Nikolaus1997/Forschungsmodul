#include "computation.h"


//initialize the computation object, parse the settings from file that is given as the only command line argument 
void Computation::initialize(std::string filename)
{   
    //load settings
    settings_ = Settings();
    settings_.loadFromFile(filename);
    settings_.printSettings();

    epsilon_ = settings_.Epsilon;
    m_ = double(settings_.BarenblattM);
    a_ = settings_.physicalSize[0];
    b_ = settings_.physicalSize[1];
    initCondA_ = settings_.initCondA;
    initCondB_ = settings_.initCondB;

    firstLimiterCalls_ = 0;
    secondLimiterCalls_=0;

    PP_N_= settings_.PP_N;

    nNodes = PP_N_;

    nCells_= settings_.nCells;
    dt_ = settings_.CFL*1/(nCells_[0]);
    meshWidth_[0] =  (b_-a_)/double(nCells_[0]);

    innerMeshWidth_[0] = meshWidth_[0]/(nNodes+2);

    //initialize the Vandermonde Matrix
    std::array<int,2> t  ={nCells_[0],PP_N_+1};
    VdM_ = std::make_shared<Vandermonde>(t,nNodes+2);

    //initialize Projection Operator
    Projection proj_ = Projection();

    //initialize the flux function;
    flux_ = Flux();
    //initialize the initial Condition
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
    if(settings_.useLimiter=="true"){
        useLimiter_ = true;
    }else{
        useLimiter_ = false;
    }
    // fillXanalyze(grid_->x_analyze_);
    // grid_->x_analyze_.printValues();
    // Array2D x_ = Array2D({1,3});
    // Array2D u_ = Array2D({1,3});
    // x_(0,0) = 0.0; x_(0,1) = 1.0; x_(0,2) = 2.0;// x_(0,3) = 3.0;
    // u_(0,0) = 6.0; u_(0,1) = 0.0; u_(0,2) = 0.0;// u_(0,3) = 20.0;
    // proj_.makeProjection(u_,x_,0,2);
    // double sol = limiter_.computeLimiter(-.011,-.000025,-3.,.1);
    // std::cout<<" LIM "<<sol<<std::endl;
    
}

void Computation::runSimulation()
{
    std::cout<<'-'<<std::flush;
    Timer timer;
    timer.start();
    double time_ = 0.0;
    int iter = 0;
    fillFaces();
    fillX();
    //grid_->x_.printValues();
    fillXanalyze(grid_->x_analyze_);
    initVdm();
    grid_->fillArray(grid_->u(),VdM_->VdM_,VdM_->L_);
    grid_->fillSolution(grid_->solution_,grid_->u());
    initVdmJ();
    grid_->fillArray(grid_->j(),VdM_->VdMJ_,VdM_->L_);
    grid_->fillSolution(grid_->solutionJ_,grid_->j());
    VdM_->LprintValues();
    VdM_->LprimePrintValues();
    if(flux_.getFluxFunction()==Flux::FunctionType::Barenblatt){
        dt_ = settings_.CFL*1/(nCells_[0]*nCells_[0])*1/double(settings_.BarenblattM);}
    else{
        dt_ = settings_.CFL*1/nCells_[0];
    }
    double numberN = settings_.nWriteState/dt_*settings_.CFL;
    if(time_<dt_){
        grid_->fillSolution(grid_->solution_,grid_->u_);
        grid_->fillSolution(grid_->derivative_,grid_->ut_);
        outputWriterParaview_ = std::make_unique<OutputWriterParaview>(grid_);
        outputWriterParaview_->writeFile(time_,settings_.OutputName);

    }
    double errorTime =0.0;
    double numberofIterations = settings_.endTime/dt_;
    calcError(errorTime);
    // // grid_->x_.printValues();
    // grid_->u_.printValues();
    if(useLimiter_){
        std::cout<<"Using Limiter"<<std::endl;
    firstLimiter(grid_->u());
    secondLimiter(grid_->u());
    }
    calcError(errorTime);
    // grid_->u_.printValues();
    while (time_<settings_.endTime)
        {
            if(flux_.getFluxFunction()!=Flux::FunctionType::Barenblatt){
                if(settings_.timeStepping=="euler" or settings_.timeStepping=="Euler"){
                    eulerTimeStep();

                }else if(settings_.timeStepping=="RK" or settings_.timeStepping=="rungeKutta" or settings_.timeStepping=="RungeKutta"){
                    rungeKutta();
                }
            }else{
                if(settings_.timeStepping=="euler" or settings_.timeStepping=="Euler"){
                    calcDt();
                    eulerTimeStep();

                    secondLimiter(grid_->u());
                    firstLimiter(grid_->u());

                }else if(settings_.timeStepping=="RK" or settings_.timeStepping=="rungeKutta" or settings_.timeStepping=="RungeKutta"){
                    calcDt();
                    if(time_+dt_>settings_.endTime)
                        dt_ = settings_.endTime-time_;
                    rungeKutta();                
                }
            }
            errorTime+=dt_;
            time_+=dt_;
            iter++;
            if(iter%int(settings_.nWriteState)==0){
                grid_->fillSolution(grid_->solution_,grid_->u_);
                grid_->fillSolution(grid_->derivative_,grid_->ut_);
                outputWriterParaview_->writeFile(time_,settings_.OutputName);
                outputWriterParaview_->writeFileTrueSolution(time_,settings_.OutputName+"TrueSolution");
                std::cout<<"Write State TIME: "<<time_<<std::endl;
                calcError(time_); 
            }

            // calcDt();

            std::cout<<"\rCurrent Time: "<<time_<<" End Time: "<<settings_.endTime<< std::flush;
    }
    std::cout<<" TIME: "<<time_<<std::endl;
    grid_->fillSolution(grid_->solution_,grid_->u_);
    grid_->fillSolution(grid_->derivative_,grid_->ut_);
    calcError(time_);
    outputWriterParaview_->writeFile(time_,settings_.OutputName);
    outputWriterParaview_->writeFileTrueSolution(time_,settings_.OutputName+"TrueSolution");
    //calcError(time_);
    timer.stop();
    std::cout << "Elapsed time: " << timer.elapsedMilliseconds()/1000 << " s." << " nStates: "<<iter<<" firstLimiterCalls "<<firstLimiterCalls_<<" secondlimiterCalls "<<secondLimiterCalls_ <<std::endl;
}


void Computation::calcDt(){
    double max = 0.0;
    double dx =0.0;
    double c = 0.0;
    for(int i=0;i<grid_->u_.size()[0];i++){
        for(int j = 0; j<grid_->u_.size()[1];j++){
            double u_temp = grid_->u_(i,j);
            // if(abs(u_temp)<1E-11)
            //     continue;
            double m_ = double(settings_.BarenblattM);
            if(initialCond_.selectedFunction==InitialCondition::InitialCondType::Barenblatt){
                c = m_*pow(u_temp,m_-1.0);
            }else{
            c =     1.;
            }
            if(c>max){
                max = c;
                }
        }
    }
    dx = meshWidth_[0];    
    dt_ = settings_.CFL*dx*dx/max*0.5;
    if(dt_>settings_.dt)
        dt_ = settings_.dt;
    //std::cout<<" DT "<<dt_<<std::endl;
}


void Computation::calcQ(const Array2D& u)
{
    double m = double(settings_.BarenblattM);
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        double flux_term =0.0,ul_i = 0.0,ur_i = 0.0,ur_iminus  = 0.0,ul_iplus =0.0;
        double u_mean = 0.0,u_mean_plus = 0.0, u_mean_minus = 0.0;
        // Wrap around the grid for periodic boundary conditions
        if(i==0){
                ul_i = u(i,0);
                ur_i = u(i,nNodes+1);
                ur_iminus  = u(nCells_[0]-1,nNodes+1);
                ul_iplus = u(i+1,0);
        }else if (i==grid_->faces_.size()[0] - 2)
        {
                ul_i = u(i,0);
                ur_i = u(i,nNodes+1);
                ur_iminus  = u(i-1,nNodes+1);
                ul_iplus = u(0,0);
        }else{
        // Compute the numerical flux
                ul_i = u(i,0);
                ur_i = u(i,nNodes+1);
                ur_iminus  = u(i-1,nNodes+1);
                ul_iplus = u(i+1,0);
                for(int k =0; k<u.size()[1];k++){
                    u_mean += u(i,k);
                    u_mean_plus += u(i+1,k);
                    u_mean_minus += u(i-1,k);
                }
                u_mean = u_mean/(u.size()[1]);
                u_mean_plus = u_mean_plus/(u.size()[1]);  
                u_mean_minus = u_mean_minus/(u.size()[1]);
        }

        double gFlux_minus = -gFlux_.computeNumFlux(true,ur_iminus,ul_i,0.,0.,m,flux_,quad_,u_mean_minus,u_mean)[1];
        //std::cout<<" FROM CALCQ "<<" g_jminus "<<" i "<<i<<" gFlux_minus "<<gFlux_minus<<std::endl;

        double gFlux_plus  =  gFlux_.computeNumFlux(false,ur_i, ul_iplus,0.,0.,m,flux_,quad_,u_mean,u_mean_plus)[1] ;
        //std::cout<<" FROM CALCQ "<<" g_jplus "<<" gFlux_plus "<<gFlux_plus<<std::endl;
        for (int j = 0; j <=PP_N_; j++) {
            double l = j;
            // Compute the numerical flux
            //std::cout<<" FROM CALCQ "<<std::endl;
            flux_term = gFlux_plus+gFlux_minus*VdM_->L(0,j);
            // Apply the formula for the update of VdM_t_* pow(-1, j) 
            double a = grid_->faces(i);
            double b = grid_->faces(i + 1);
            double integ =quad_->IntFluxQ([&](double x) {return flux_.compute(x, 0.0, m)[1];},i, j, a, b, u);
            if(std::abs(flux_term)<1E-12)
                flux_term=0.0;
            if(std::abs(integ)<1E-12)
                integ=0.0;
            VdM_->VdMQ_(i,j) = (integ-flux_term)*(2. * l + 1.)/meshWidth_[0];
            // std::cout<<" IN CALCQ I "<<" i "<<i<<" Face i "<<grid_->faces(i)<<" Face i+1 "<<grid_->faces(i+1)<<" J "<<j<<" FLUX TERM "<<flux_term<<" INTEGRAL "<<integ<<" VDMQ "<<VdM_->VdMQ_(i,j)<<std::endl;
            // std::cout<<" IN CALCQ I " <<i<<" J "<<j<<" gFlux_minus "<<gFlux_minus<<" gFlux_plus "<<gFlux_plus<<" L "<<VdM_->L(0,j)<<std::endl;  
        }
    }
    grid_->fillArray(grid_->q_,VdM_->VdMQ_,VdM_->L_);
    // std::cout<<" Q "<<std::endl;
    // grid_->q_.printValues();
}

void Computation::eulerTimeStep()
{   
    if(flux_.getFluxFunction()!=Flux::FunctionType::Barenblatt){
        calcUdt(grid_->u_,VdM_->VdM_t_);
    }else{  
        calcQ(grid_->u_);
        calcUdt(grid_->u_,grid_->q_, VdM_->VdMJ_t_, epsilon_);
    }
    grid_->fillArray(grid_->jt_,VdM_->VdMJ_t_,VdM_->L_);
    calcUdt(grid_->jt_,VdM_->VdM_t_);
    grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
    for(int i = 0; i<grid_->u_.size()[0];i++){
        for(int j = 0; j<grid_->u_.size()[1];j++){
            grid_->u_(i,j) += dt_*grid_->ut_(i,j);
        }
    }
}

void Computation::rungeKutta() {
    // Step 1: Compute the intermediate stage u^(1)
    if(flux_.getFluxFunction()!=Flux::FunctionType::Barenblatt){
        calcUdt(grid_->u_,VdM_->VdM_t_); // Compute the time derivative for u^n
    }else
    {
        calcQ(grid_->u_);
        calcUdt(grid_->u_,grid_->q_,VdM_->VdMJ_t_,epsilon_);
    }
    grid_->fillArray(grid_->jt_,VdM_->VdMJ_t_,VdM_->L_);
    for (int i = 0; i < grid_->j_.size()[0]; i++) {
        for (int j = 0; j <grid_->j_.size()[1]; j++) {
            // u^(1) = u^n + Δt * L_h(u^n)
            grid_->j_1_(i, j)  = grid_->j_(i,j)+ dt_ * grid_->jt_(i, j);
        }
    }
    for(int i = 0; i<VdM_->VdM_.size()[0];i++){
        for(int j = 0; j<VdM_->VdM_.size()[1];j++){
            VdM_->VdMJ1_(i,j) = VdM_->VdMJ_(i,j)+ dt_ * VdM_->VdMJ_t_(i,j);
        }
    }

    if(flux_.getFluxFunction()==Flux::FunctionType::Barenblatt)
        calcUdt(grid_->j_1_,VdM_->VdM_t_); // Compute the time derivative for u^n
    grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
    for (int i = 0; i < grid_->j_.size()[0]; i++) {
        for (int j = 0; j <grid_->j_.size()[1]; j++) {
            // u^(1) = u^n + Δt * L_h(u^n)
            grid_->u1_(i, j)  = grid_->u_(i,j)+ dt_ * grid_->ut_(i, j);
        }
    }
    for(int i = 0; i<VdM_->VdM_.size()[0];i++){
        for(int j = 0; j<VdM_->VdM_.size()[1];j++){
            VdM_->VdM1_(i,j) = VdM_->VdM_(i,j)+ dt_ * VdM_->VdM_t_(i,j);
        }
    }

    if(useLimiter_){
    //    std::cout<<" U1 "<<std::endl;
    //    grid_->u1_.printValues();
    //secondLimiter(grid_->j_1_);
    firstLimiter(grid_->u1_);
    //    std::cout<<" AFTER first LIMITER U1 "<<std::endl;
    //    grid_->u1_.printValues();
    secondLimiter(grid_->u1_);
    }
    //    std::cout<<" AFTER LIMITER U1 "<<std::endl;
    //    grid_->u1_.printValues();
    // Step 2: Compute the intermediate stage u^(2)
    if(flux_.getFluxFunction()!=Flux::FunctionType::Barenblatt){
        calcUdt(grid_->u1_,VdM_->VdM_t_); // Compute the time derivative for u^n
    }else
    {
        calcQ(grid_->u1_);
        calcUdt(grid_->u1_,grid_->q_,VdM_->VdMJ_t_,epsilon_);
    }// Compute the time derivative for u^(1)

    grid_->fillArray(grid_->jt_,VdM_->VdMJ_t_,VdM_->L_);
    for (int i = 0; i < grid_->j_.size()[0]; i++) {
        for (int j = 0; j <grid_->j_.size()[1]; j++) {
            // u^(2) = 3/4 * u^n + 1/4 * u^(1) + (1/4 * Δt) * L_h(u^(1))
            grid_->j_2_(i, j) = (3.0 / 4.0) * grid_->j_(i, j) +
                                (1.0 / 4.0) * grid_->j_1_(i, j) +
                                (1.0 / 4.0) * dt_ * grid_->jt_(i, j);
        }
    }
    for(int i = 0;i<VdM_->VdM_.size()[0];i++){
        for(int j = 0;j<VdM_->VdM_.size()[1];j++){
            VdM_->VdMJ2_(i,j) = (3.0 / 4.0) * VdM_->VdMJ_(i,j) +
                                (1.0 / 4.0) * VdM_->VdMJ1_(i,j) +
                                (1.0 / 4.0) * dt_ * VdM_->VdMJ_t_(i,j);
        }
    }
    if(flux_.getFluxFunction()==Flux::FunctionType::Barenblatt)
        calcUdt(grid_->j_2_,VdM_->VdM_t_); // Compute the time derivative for u^n
    grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
    for (int i = 0; i < grid_->j_.size()[0]; i++) {
        for (int j = 0; j <grid_->j_.size()[1]; j++) {
            // u^(2) = 3/4 * u^n + 1/4 * u^(1) + (1/4 * Δt) * L_h(u^(1))
            grid_->u2_(i, j) = (3.0 / 4.0) * grid_->u_(i, j) +
                                (1.0 / 4.0) * grid_->u1_(i, j) +
                                (1.0 / 4.0) * dt_ * grid_->ut_(i, j);
        }
    }
    for(int i = 0;i<VdM_->VdM_.size()[0];i++){
        for(int j = 0;j<VdM_->VdM_.size()[1];j++){
            VdM_->VdM2_(i,j) = (3.0 / 4.0) * VdM_->VdM_(i,j) +
                                (1.0 / 4.0) * VdM_->VdM1_(i,j) +
                                (1.0 / 4.0) * dt_ * VdM_->VdM_t_(i,j);
        }
    } 
    if(useLimiter_){
    //  std::cout<<" U2 "<<std::endl;
    //  grid_->u2_.printValues();
    //secondLimiter(grid_->u2_); 
    firstLimiter(grid_->u2_);
    secondLimiter(grid_->u2_); 
    }
    //  std::cout<<" AFTER LIMITER U2 "<<std::endl;
    //  grid_->u2_.printValues();
    // Step 2: Compute the intermediate stage u^(2)
    if(flux_.getFluxFunction()!=Flux::FunctionType::Barenblatt){
        calcUdt(grid_->u2_,VdM_->VdM_t_); // Compute the time derivative for u^n
    }else
    {
        calcQ(grid_->u2_);
        calcUdt(grid_->u2_,grid_->q_,VdM_->VdMJ_t_);
    } // Compute the time derivative for u^(2)
    // Ensure all updates to time derivatives are reflected
    grid_->fillArray(grid_->jt_,VdM_->VdMJ_t_,VdM_->L_);
    for (int i = 0; i < grid_->u_.size()[0]; i++) {
        for (int j = 0; j <grid_->u_.size()[1]; j++) {
            // u^(n+1) = 1/3 * u^n + 2/3 * u^(2) + (2/3 * Δt) * L_h(u^(2))
            grid_->j_(i, j) = (1.0 / 3.0) *  grid_->j_(i, j) +
                               (2.0 / 3.0) *  grid_->j_2_(i, j) +
                               (2.0 / 3.0) * dt_ *  grid_->jt_(i, j);
        }
    }
    for(int i = 0;i<VdM_->VdM_.size()[0];i++){
        for(int j = 0;j<VdM_->VdM_.size()[1];j++){
            VdM_->VdMJ_(i,j) = (1.0 / 3.0) *  VdM_->VdMJ_(i,j) +
                               (2.0 / 3.0) *  VdM_->VdMJ2_(i,j) +
                               (2.0 / 3.0) * dt_ *  VdM_->VdMJ_t_(i,j);
        }
    }

    if(flux_.getFluxFunction()!=Flux::FunctionType::Barenblatt){
        calcUdt(grid_->j_,VdM_->VdM_t_); // Compute the time derivative for u^n
    }
    grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
    for (int i = 0; i < grid_->u_.size()[0]; i++) {
        for (int j = 0; j <grid_->u_.size()[1]; j++) {
            // u^(n+1) = 1/3 * u^n + 2/3 * u^(2) + (2/3 * Δt) * L_h(u^(2))
            grid_->u_(i, j) = (1.0 / 3.0) *  grid_->u_(i, j) +
                               (2.0 / 3.0) *  grid_->u2_(i, j) +
                               (2.0 / 3.0) * dt_ *  grid_->ut_(i, j);
        }
    }
    for(int i = 0;i<VdM_->VdM_.size()[0];i++){
        for(int j = 0;j<VdM_->VdM_.size()[1];j++){
            VdM_->VdM_(i,j) = (1.0 / 3.0) *  VdM_->VdM_(i,j) +
                               (2.0 / 3.0) *  VdM_->VdM2_(i,j) +
                               (2.0 / 3.0) * dt_ *  VdM_->VdM_t_(i,j);
        }
    }

    if(useLimiter_){
    //   std::cout<<" U "<<std::endl;
    //   grid_->u_.printValues();
    //secondLimiter(grid_->j_);
    firstLimiter(grid_->u_);
    secondLimiter(grid_->u_);
    //   std::cout<<" AFTER LIMITER U "<<std::endl;
    //   grid_->u_.printValues();   
    } 
}



void Computation::firstLimiter(Array2D &u)
{
    for(int i=0;i<grid_->u_.size()[0];i++){
        double limit_l=0.0, limit_r=0.0;
        double u_r =0.0, u_l =0.0,u_mean = 0.0,u_mean_plus = 0.0,u_mean_minus = 0.0;
        bool check = false;
        u_r = u(i,nNodes+1);
        u_l = u(i,0);
        if(i==0){
            for(int k = 0; k<u.size()[1];k++){
                u_mean += u(i,k);
                u_mean_plus += u(i+1,k);
                u_mean_minus += u(nCells_[0]-1,k);
            }
        }else if(i==nCells_[0]-1){
            for(int k = 0; k<u.size()[1];k++){
                u_mean += u(i,k);
                u_mean_plus += u(0,k);
                u_mean_minus += u(i-1,k);
            }
            }else{
            for(int k = 0; k<u.size()[1];k++){
                u_mean += u(i,k);
                u_mean_plus += u(i+1,k);
                u_mean_minus += u(i-1,k);
            }
                    }
        u_mean = u_mean/(u.size()[1]);
        u_mean_plus = u_mean_plus/(u.size()[1]); 
        u_mean_minus = u_mean_minus/(u.size()[1]);
              
        double a = u_r -u_mean;
        double b = u_mean - u_mean_minus;
        double c = u_mean_plus - u_mean;    
        limit_r = u_mean+limiter_.computeLimiter(a,b,c,meshWidth_[0]);
        double a_ = u_mean - u_l;
        double b_ = u_mean - u_mean_minus;
        double c_ = u_mean_plus - u_mean;
        limit_l = u_mean-limiter_.computeLimiter(a_,b_,c_,meshWidth_[0]);
        
        if(std::abs(limit_r-u_r)>1E-12){
            // for(int k = 1; k<u.size()[1]-1;k++){
            //     u(i,k) = u_mean;
            // }
            proj_.makeProjection(u,grid_->x_,i,2);
            check = true;
        }else
        if(std::abs(limit_l-u_l)>1E-12){
            // for(int k = 1; k<u.size()[1]-1;k++){
            //     u(i,k) = u_mean;
            // }
            proj_.makeProjection(u,grid_->x_,i,2);
            check = true;
        }    
        if(check){
            double mean_i= 0.0,mean_iminus = 0.0,mean_iplus = 0.0;   
            for(int j = 0;j<u.size()[1];j++){
                mean_i += u(i,j);
                if((i>0))   
                    mean_iminus += u(i-1,j);
                if((i<nCells_[0]-1))
                    mean_iplus += u(i+1,j);
            }
            mean_i = mean_i/(u.size()[1]);
            mean_iminus = mean_iminus/(u.size()[1]);
            mean_iplus = mean_iplus/(u.size()[1]);
            double slope = (u(i,nNodes+1)-u(i,0))/meshWidth_[0];
            for(int j = 0;j<u.size()[1];j++){
                u(i,j)= mean_i+(grid_->x_(i,j)-(grid_->x(i,0)+grid_->x(i,nNodes+1))/2.0)
                            *limiter_.computeLimiter(slope,(mean_i-mean_iminus)*2./double(meshWidth_[0]),(mean_iplus-mean_i)*2./double(meshWidth_[0]),double(meshWidth_[0]));
                    }
        firstLimiterCalls_++;
        }
    }

}

void Computation::secondLimiter(Array2D &u)
{
    for(int i=0;i<u.size()[0];i++){
        for(int j=0;j<u.size()[1];j++){
            if(u(i,j)<0){
                proj_.makeProjection(u,grid_->x_,i,2);

              double x_j = (grid_->x_(i,0)+grid_->x_(i,nNodes+1))/2.0;
              double mean = 0.0;

            for(int k = 0;k<u.size()[1];k++){
                mean += u(i,k);
            }
            mean = mean/(u.size()[1]);
            //std::cout<<" PROJECTION i "<<i<<" "<<u(i,0)<<" "<<u(i,1)<<" "<<u(i,2)<<" "<<u(i,3)<<" "<<u(i,4)<<" mean "<<mean<<" x_j "<<x_j<<std::endl;
                    if(u(i,nNodes+1)<0.){
                        for(int k = 0;k<u.size()[1];k++){
                            if(abs(grid_->x(i,k)-double(x_j)-meshWidth_[0]/2.)<1E-14){
                                u(i,k)=0.0;
                            }else{

                            u(i,k)=(1.0-2.0/double(meshWidth_[0])*(grid_->x(i,k)-x_j))*mean;
                            }
                            if(mean<0)
                                u(i,k) =0.;
                            //std::cout<<"RIGHT "<<" ij "<<i<<" "<<j<<" u "<<u(i,k)<<" mean "<<mean<<" x_j "<<x_j<<" x "<<grid_->x(i,k)<<" "<<grid_->x(i,k)-x_j-meshWidth_[0]/2.<<std::endl; 
                        }
                        secondLimiterCalls_++;
                        break;
                    }else
                    if(u(i,0)<0.){
                        for(int k = 0;k<u.size()[1];k++){
                            if(abs(grid_->x(i,k)-x_j+meshWidth_[0]/2)<1E-14){
                                u(i,k)=0.0;
                            }else{

                            u(i,k)=(1.0+2.0/meshWidth_[0]*(grid_->x(i,k)-x_j))*mean;
                            }
                            if(mean<0)
                                u(i,k) =0.;
                            //std::cout<<"LEFT "<<" ij "<<i<<" "<<j<<" u "<<u(i,k)<<" mean "<<mean<<" x_j "<<x_j<<" x "<<grid_->x(i,k)<<" "<<grid_->x(i,k)-x_j-meshWidth_[0]/2.<<std::endl; 

                        }
                        secondLimiterCalls_++;
                        break;
                }
            }
        }
    }



}



void Computation::calcUdt(const Array2D& u,const Array2D& q, Array2D& VdM_t, double epsilon)
{
    double flux_term =0.0,integ=0.0;
    double ul_i =0.0,ur_i =0.0,ur_iminus  = 0.0,ul_iplus = 0.0;
    double ql_i =0.0, qr_i = 0.0, qr_iminus =0.0, ql_iplus = 0.0;
    double m = double(settings_.BarenblattM);
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        double u_mean = 0.0, u_mean_plus = 0.0, u_mean_minus = 0.0;
            ul_i =0.0;
            ur_i =0.0;
            ur_iminus  = 0.0;
            ul_iplus = 0.0;
            ql_i =0.0;
            qr_i = 0.0;
            qr_iminus =0.0;
            ql_iplus = 0.0;
            // Wrap around the grid for periodic boundary conditions
            if(i==0){
                ul_i        = u(i,0);
                ur_i        = u(i,nNodes+1);
                ur_iminus   = u(nCells_[0]-1,nNodes+1);
                ul_iplus    = u(i+1,0);
                ql_i        = q(i,0);
                qr_i        =  q(i,nNodes+1);
                qr_iminus   = q(q.size()[0]-1,nNodes+1);
                ql_iplus    = q(i+1,0);
            }else if (i==grid_->faces_.size()[0] - 2)
            {
                ul_i        = u(i,0);
                ur_i        = u(i,nNodes+1);
                ur_iminus   = u(i-1,nNodes+1);
                ul_iplus    = u(0,0);
                ql_i        = q(i,0);
                qr_i        =  q(i,nNodes+1);
                qr_iminus   = q(i-1,nNodes+1);
                ql_iplus    = q(0,0);                                              
            }else{
            // Compute the numerical flux
                ul_i        = u(i,0);
                ur_i        = u(i,nNodes+1);
                ur_iminus   = u(i-1,nNodes+1);
                ul_iplus    = u(i+1,0);
                ql_i        = q(i,0);
                qr_i        =  q(i,nNodes+1);
                qr_iminus   = q(i-1,nNodes+1);
                ql_iplus    = q(i+1,0); 
                for(int k = 0; k<u.size()[1];k++){
                    u_mean += u(i,k);
                    u_mean_plus += u(i+1,k);
                    u_mean_minus += u(i-1,k);
                }
                u_mean = u_mean/(u.size()[1]);
                u_mean_plus = u_mean_plus/(u.size()[1]);  
                u_mean_minus = u_mean_minus/(u.size()[1]);
            }
           double g_minus =  -gFlux_.computeNumFlux(true,ur_iminus,ul_i,qr_iminus,ql_i,m,flux_,quad_,u_mean_minus,u_mean)[0];
           double g_plus  =   gFlux_.computeNumFlux(false,ur_i, ul_iplus,qr_i,ql_iplus,m,flux_,quad_,u_mean,u_mean_plus)[0]; 
        for (int j = 0; j <=PP_N_; j++) {
            double l = double(j);
            flux_term = g_minus*VdM_->L_(0,j)+g_plus;
            //std::cout<<" FROM CALC UUUU "<<" i "<<i<<std::endl;
            if(abs(flux_term)<1E-12)
                flux_term = 0.0;
            integ =quad_->IntFluxU([&](double u, double q) { return flux_.compute(u, q, m)[0]; }, i, j,
                                            grid_->faces(i), grid_->faces(i+1),u, q,grid_->j_);
            if(abs(integ)<1E-12)
                integ = 0.0;
            // Apply the formula for the update of VdM_t_* pow(-1, j) 
            VdM_t(i,j) = 1.0/epsilon*(integ- flux_term)*double((2.0 * l + 1.0)/meshWidth_[0]);///(meshWidth_[0]);
            //std::cout<<" IN CALCUDT I " <<i<<" J "<<j<<" FLUX TERM "<<flux_term<<" INTEGRAL "<<integ<<" VDMT "<<VdM_->VdM_t_(i,j)<<" meshWidth "<<double(meshWidth_[0])<<std::endl;
        }
    }   
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
    fillUanalyze(grid_->u_analyze_,grid_->x_analyze_,VdM_->VdM_);

    if(flux_.getFluxFunction()==Flux::FunctionType::Barenblatt){
        for(int i = 0; i<grid_->u_.size()[0];i++){
            for(int j=0; j<grid_->u_.size()[1];j++){
                grid_->true_solution_(i,j) =initialCond_.computeInitialCondition(grid_->x_(i,j),initCondA_,initCondB_,currentTime+settings_.BarenblattTime, settings_.BarenblattM);     
            }
        }
        grid_->l2_error(0) = 0.0;
        for(int i =0;i<grid_->u_analyze_.size()[0];i++){
            grid_->u_analyze_true_(i,0) = initialCond_.computeInitialCondition(grid_->x_analyze_(i,0),initCondA_,initCondB_,currentTime+settings_.BarenblattTime, settings_.BarenblattM);
            grid_->l2_error(0) += pow(grid_->u_analyze_(i,0)- grid_->u_analyze_true_(i,0),2.);     
        }
        
        grid_->linf_error(0) = 0.0;

        for(int i =0;i<grid_->u_analyze_.size()[0];i++){
            if(grid_->linf_error(0)<abs(grid_->u_analyze_(i,0)-initialCond_.computeInitialCondition(grid_->x_analyze_(i,0),initCondA_,initCondB_,currentTime+settings_.BarenblattTime, settings_.BarenblattM))) 
                grid_->linf_error(0) = abs(grid_->u_analyze_(i,0)-initialCond_.computeInitialCondition(grid_->x_analyze_(i,0),initCondA_,initCondB_,currentTime+settings_.BarenblattTime, settings_.BarenblattM));
        }
    }
    else if(settings_.initialCondition =="sinus"){
        for(int i = 0; i<grid_->u_.size()[0];i++){
            for(int j=0; j<grid_->u_.size()[1];j++){
                grid_->true_solution_(i,j) =sin(grid_->x_(i,j));     
            }
        }
        grid_->l2_error(0) = 0.0;
        grid_->linf_error(0) = 0.0;
        for(int i =0;i<grid_->u_analyze_.size()[0];i++){
            double sol = initialCond_.computeInitialCondition(grid_->x_analyze_(i,0),initCondA_,initCondB_);
            //std::cout<<" SOL "<<sol<<" "<<grid_->u_analyze_(i,0)<<std::endl;
            grid_->l2_error(0) += pow(grid_->u_analyze_(i,0)-sol ,2);     
        if(grid_->linf_error(0)<abs(grid_->u_analyze_(i,0)-sol))
                grid_->linf_error(0) = abs(grid_->u_analyze_(i,0)-sin(grid_->x_analyze_(i,0)));
        }
    }
    std::cout<<"L2 Error: "<<sqrt(grid_->l2_error(0)*1./double(grid_->u_analyze_.size()[0]))<<" Linf Error: "<<grid_->linf_error(0)<<std::endl;
}

void Computation::fillX()
{
    double transformedNode = 0.0;
    double mean = 0.0, diff  =0.0;
    for (int i = 0; i < grid_->faces_.size()[0]-1; i++)
    {
        for(int j=0 ; j<grid_->x_.size()[1];j++){
            mean = 0.5 * (grid_->faces_(i+1) + grid_->faces_(i));
            diff = 0.5 * (grid_->faces_(i+1) - grid_->faces_(i));
            transformedNode = mean + diff*quad_->basis_.nodes(j);
            grid_->x(i,j) = transformedNode;
        }
    }   
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

void Computation::initVdmJ(){
    for(int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        double flux_term_left = pow(grid_->u_(i,0),m_);
        double flux_term_right = pow(grid_->u_(i,nNodes+1),m_);
        for (int j = 0; j < VdM_->VdM_.size()[1]; j++) {
            // Compute the integral for the j-th polynomial degree
            double integral = quad_->IntJ_0(grid_->u_,m_,i,j);
            double flux_term = -flux_term_left*VdM_->L(0,j) + flux_term_right;

            VdM_->VdMJ_(i,j) = (-integral + flux_term) * (2.0 * double(j) + 1.0) / meshWidth_[0];
            }
        }
}

void Computation::calcUdt(const Array2D& u_, Array2D& VdM_t){
    Flux uflux_;
    uflux_.setFluxFunction(Flux::FunctionType::Linear);
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
            double ul_i = 0.0;
            double ur_i = 0.0;
            double ur_iminus  =0.0;
            double ul_iplus = 0.0; 
            double flux_term =0.0;
            if(i==0){
                ul_i        = u_(i,0);
                ur_i        = u_(i,nNodes+1);
                ur_iminus   = u_(nCells_[0]-1,nNodes+1);
                ul_iplus    = u_(i+1,0);
            }else if (i==grid_->faces_.size()[0] - 2)
            {
                ul_i        = u_(i,0);
                ur_i        = u_(i,nNodes+1);
                ur_iminus   = u_(i-1,nNodes+1);
                ul_iplus    = u_(0,0);
            }else{
            //Compute the numerical flux
                ul_i        = u_(i,0);
                ur_i        = u_(i,nNodes+1);
                ur_iminus   = u_(i-1,nNodes+1);
                ul_iplus    = u_(i+1,0);
            }
            for (int j = 0; j <=PP_N_; j++) {
            // Wrap around the grid for periodic boundary conditions
            flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,uflux_,dt_,meshWidth_[0])*VdM_->L_(0,j)  + gFlux_.computeNumFlux(ur_i, ul_iplus,uflux_,dt_,meshWidth_[0]);
            double integ =quad_->IntFluxGaussLegendreQuad([&](double x) {return uflux_.compute(x);}
                                                             ,i,j ,grid_->faces(i),grid_->faces(i+1),u_);

            VdM_t(i,j) =integ*(2.0*double(j)+1.0)*1/meshWidth_[0] - flux_term*1/meshWidth_[0]*(2.0*double(j)+1.0);
        }
    }
}

void Computation::fillXanalyze(Array2D &x)
{
    double stepsize_ = (settings_.initCondB-settings_.initCondA)/double(grid_->x_analyze_.size()[0]-1);
    for(int i=0;i<grid_->x_analyze_.size()[0];i++){
        x(i,0) = settings_.initCondA+i*stepsize_;
    }
}

void Computation::fillUanalyze(Array2D &u_analyze, const Array2D &x, const Array2D &Vdm)
{
    for(int i =0;i<x.size()[0];i++){
        u_analyze(i,0) = 0.0;
        for(int k=0;k<grid_->faces().size()[0]-1;k++){
            if(x(i,0)>=grid_->faces(k) and x(i,0)<grid_->faces(k+1)){
                double transformedNode =  (x(i,0) - 0.5 * (grid_->faces(k+1) + grid_->faces(k)))/(meshWidth_[0]*0.5);
                for(int j = 0; j<Vdm.size()[1];j++){
                    u_analyze(i,0) += Vdm(k,j)*quad_->LegendrePolynomialAndDerivative(j,transformedNode)[0];
                }
                break;
            }
        }

    }

}
