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
    innerMeshWidth_[0] = meshWidth_[0]/(nNodes+2);

    //initialize the Vandermonde Matrix
    std::array<int,2> t  ={nCells_[0],PP_N_+1};
    VdM_ = std::make_shared<Vandermonde>(t,nNodes+2);


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
    grid_->fillArray(grid_->u(),VdM_->VdM_,VdM_->L_);
    grid_->fillSolution(grid_->solution_,grid_->u());
    VdM_->LprintValues();
    VdM_->LprimePrintValues();
    if(settings_.BarenblattM!=0){
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
    //calcError(errorTime);
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
                    // if(iter<80){
                    //     std::cout<<"U"<<std::endl;
                    //     grid_->u_.printValues();
                    //     std::cout<<"Q"<<std::endl;
                    //     grid_->q_.printValues();
                    // }
                    eulerTimeStep();

                    applyLimiter(grid_->u());
                    // if(iter<80){
                    //     std::cout<<" AFTER LIMIT U"<<std::endl;
                    //     grid_->u_.printValues();
                    //     std::cout<<" AFTER LIMIT Q"<<std::endl;
                    //     grid_->q_.printValues();
                    // }

                }else if(settings_.timeStepping=="RK" or settings_.timeStepping=="rungeKutta" or settings_.timeStepping=="RungeKutta"){
                    rungeKutta();                 
                }
            }
            errorTime+=dt_;
            time_+=dt_;
            if(iter % int(numberN) ==0){
                grid_->fillSolution(grid_->solution_,grid_->u_);
                grid_->fillSolution(grid_->derivative_,grid_->ut_);
                outputWriterParaview_->writeFile(time_,settings_.OutputName);
                //calcError(errorTime);
            }

            // calcDt();
            iter++;
            std::cout<<"\rCurrent Iteration: "<<iter<<" End Iter: "<<numberofIterations<< std::flush;
    }
    std::cout<<" TIME: "<<time_<<std::endl;
    grid_->fillSolution(grid_->solution_,grid_->u_);
    grid_->fillSolution(grid_->derivative_,grid_->ut_);
    outputWriterParaview_->writeFile(time_,settings_.OutputName);
    //calcError(time_);
    timer.stop();
    std::cout << "Elapsed time: " << timer.elapsedMilliseconds()/1000 << " s." << " nStates: "<<iter<< std::endl;
}


void Computation::calcDt(){
    double max = 0.0;
    double dx =0.0;
    double c = 0.0;
    // for(int i=0;i<grid_->u_.size()[0];i++){
    //     double u = grid_->u(i);
    //     if(u==0)
    //         continue;
    //     dx = grid_->faces(i+1)-grid_->faces(i);
    //     c = flux_.compute(u,0.0,double(settings_.BarenblattM))[1]*flux_.compute(u,0.0,double(settings_.BarenblattM))[1];
    //     double dt = settings_.CFL*dx*dx/c*0.5;
    //     if(dt>max){
    //         max = dt;
    //     }
    // }
    // dt_ = max;
    //std::cout<<" DT "<<dt_<<std::endl;
}


void Computation::calcQ(const Array2D& u)
{
    double m = double(settings_.BarenblattM);
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
        double flux_term =0.0,ul_i = 0.0,ur_i = 0.0,ur_iminus  = 0.0,ul_iplus =0.0;
        double u_mean = 0.0,u_mean_plus = 0.0;
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
                ul_iplus = u(i,0);
                for(int k = 1; k<u.size()[1]-1;k++){
                    u_mean += u(i,k);
                    u_mean_plus += u(i+1,k);
                }
                u_mean = u_mean/(u.size()[1]-2);
                u_mean_plus = u_mean_plus/(u.size()[1]-2);  
        }

        for (int j = 0; j <=PP_N_; j++) {
            double l = j;
            // Compute the numerical flux
            //std::cout<<" FROM CALCQ "<<std::endl;
            flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,0.,0.,m,flux_,quad_,u_mean,u_mean_plus)[1]*VdM_->L_(0,j)
                    + gFlux_.computeNumFlux(ur_i, ul_iplus,0.,0.,m,flux_,quad_,u_mean,u_mean_plus)[1] ;
            // Apply the formula for the update of VdM_t_* pow(-1, j) 
            double a = grid_->faces(i);
            double b = grid_->faces(i + 1);
            double integ =quad_->IntFluxQ([&](double x) {return flux_.compute(x, 0.0, m)[1];},i, j, a, b, u);
            // if(std::abs(flux_term)<1E-6)
            //     flux_term=0.0;
            // if(std::abs(integ)<1E-6)
            //     integ=0.0;
            VdM_->VdMQ_(i,j) = (integ-flux_term)*(2. * l + 1.)/meshWidth_[0];
            //std::cout<<" IN CALCQ I "<<" Face i "<<grid_->faces(i)<<" Face i+1 "<<grid_->faces(i+1)<<" J "<<j<<" FLUX TERM "<<flux_term<<" INTEGRAL "<<integ<<" VDMQ "<<VdM_->VdMQ_(i,j)<<std::endl;
            //std::cout<<" IN CALCQ I " <<i<<" J "<<j<<" FLUX TERM "<<flux_term<<" INTEGRAL "<<integ<<" VDMQ "<<VdM_->VdMQ_(i,j)<<std::endl;
        }
    }
    grid_->fillArray(grid_->q_,VdM_->VdMQ_,VdM_->L_);
}


void Computation::eulerTimeStep()
{   
    if(flux_.getFluxFunction()!=Flux::FunctionType::Barenblatt){
        calcUdt(grid_->u_);
    }else{  
        calcQ(grid_->u_);
        calcUdt(grid_->u_,grid_->q_);
    }
    grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
    for(int i = 0; i<grid_->u_.size()[0];i++){
        for(int j = 0; j<grid_->u_.size()[1];j++){
            grid_->u_(i,j) += dt_*grid_->ut_(i,j);
        }
    }
}

void Computation::rungeKutta() {
    // Step 1: Compute the intermediate stage u^(1)
    if(settings_.BarenblattM==0){
        calcUdt(grid_->u_); // Compute the time derivative for u^n
    }else
    {
        calcQ(grid_->u_);
        calcUdt(grid_->u_,grid_->q_);
    }
    grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
    for (int i = 0; i < grid_->u_.size()[0]; i++) {
        for (int j = 0; j <grid_->u_.size()[1]; j++) {
            // u^(1) = u^n + Δt * L_h(u^n)
            grid_->u1_(i, j)  = grid_->u_(i,j)+ dt_ * grid_->ut_(i, j);
        }
    }
    applyLimiter(grid_->u1_);
    // Step 2: Compute the intermediate stage u^(2)
    if(settings_.BarenblattM==0){
        calcUdt(grid_->u1_); // Compute the time derivative for u^n
    }else
    {
        calcQ(grid_->u1_);
        calcUdt(grid_->u1_,grid_->q_);
    }// Compute the time derivative for u^(1)

    grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
    for (int i = 0; i < grid_->u_.size()[0]; i++) {
        for (int j = 0; j <grid_->u_.size()[1]; j++) {
            // u^(2) = 3/4 * u^n + 1/4 * u^(1) + (1/4 * Δt) * L_h(u^(1))
            grid_->u2_(i, j) = (3.0 / 4.0) * grid_->u_(i, j) +
                                (1.0 / 4.0) * grid_->u1_(i, j) +
                                (1.0 / 4.0) * dt_ * grid_->ut_(i, j);
        }
    }
    applyLimiter(grid_->u2_);
    // Step 2: Compute the intermediate stage u^(2)
    if(settings_.BarenblattM==0){
        calcUdt(grid_->u2_); // Compute the time derivative for u^n
    }else
    {
        calcQ(grid_->u2_);
        calcUdt(grid_->u2_,grid_->q_);
    } // Compute the time derivative for u^(2)
    // Ensure all updates to time derivatives are reflected
    grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
    for (int i = 0; i < grid_->u_.size()[0]; i++) {
        for (int j = 0; j <grid_->u_.size()[1]; j++) {
            // u^(n+1) = 1/3 * u^n + 2/3 * u^(2) + (2/3 * Δt) * L_h(u^(2))
            grid_->u_(i, j) = (1.0 / 3.0) *  grid_->u_(i, j) +
                               (2.0 / 3.0) *  grid_->u2_(i, j) +
                               (2.0 / 3.0) * dt_ *  grid_->ut_(i, j);
        }
    }
    applyLimiter(grid_->u_);
}



void Computation::applyLimiter(Array2D &u)
{
    for(int i=0;i<grid_->u_.size()[0];i++){
        double limit_l=0.0, limit_r=0.0;
        double u_r =0.0, u_l =0.0,u_mean = 0.0,u_mean_plus = 0.0,u_mean_minus = 0.0;
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
        
        int midpoint = 0;
        if(PP_N_==1){
            midpoint = 1;
        }else if(PP_N_%2==0){
            midpoint = int(PP_N_/2);
        }else{
            midpoint = int(PP_N_/2)+1;
        }

        for(int k = 0; k<u.size()[1];k++){
            if(u(i,k)<0 ){
                if(k>=0){
                    u(i,k) = (1.0-2.0/meshWidth_[0]*innerMeshWidth_[0])*u(i,k);
                }else if(k<midpoint){
                    u(i,k) = (1.0+2.0/meshWidth_[0]*innerMeshWidth_[0])*u(i,k);
                }else{
                    u(i,k) = 0;    
                }
            }
            // }else if(u(i,k)<0 and u_mean>0){
            //     if(k<midpoint){
            //         u(i,k) = (1-2/meshWidth_[0]*innerMeshWidth_[0]*(k-midpoint))*u_mean;
            //     }else{
            //         u(i,k) = (1+2/meshWidth_[0]*innerMeshWidth_[0]*(k-midpoint))*u_mean;
            //     }
            // }  
        }           
        double a = u_r -u_mean;
        double b = u_mean - u_mean_minus;
        double c = u_mean_plus - u_mean;    
        limit_r = u_mean+limiter_.computeLimiter(a,b,c,innerMeshWidth_[0]);
        a = u_mean - u_l;
        b = u_mean - u_mean_minus;
        c = u_mean_plus - u_mean;
        limit_l = u_mean-limiter_.computeLimiter(a,b,c,innerMeshWidth_[0]);
        
        if(std::abs(limit_r-u_r)>1E-8){
            for(int k = 1; k<u.size()[1]-1;k++){
                u(i,k) = u_mean;
            }
        }
        if(std::abs(limit_l-u_l)>1E-8){
            for(int k = 1; k<u.size()[1]-1;k++){
                u(i,k) = u_mean;
            }
        }    
        if(u_r<0) {
                u(i,nNodes+1) = (1-2)*u_mean;
                }
        if(u_l<0){
                u(i,0) =(1-2)*u_mean;
                }
    }

}

void Computation::calcUdt(const Array2D& u,const Array2D& q)
{
    double flux_term =0.0,integ=0.0;
    double ul_i =0.0,ur_i =0.0,ur_iminus  = 0.0,ul_iplus = 0.0;
    double ql_i =0.0, qr_i = 0.0, qr_iminus =0.0, ql_iplus = 0.0;
    double m = double(settings_.BarenblattM);
    double u_mean = 0.0, u_mean_plus = 0.0;
    for (int i = 0; i < grid_->faces_.size()[0] - 1; i++) {
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
                for(int k = 1; k<u.size()[1]-1;k++){
                    u_mean += u(i,k);
                    u_mean_plus += u(i+1,k);
                }
                u_mean = u_mean/(u.size()[1]-2);
                u_mean_plus = u_mean_plus/(u.size()[1]-2);  
            }
        for (int j = 0; j <=PP_N_; j++) {
            double l = double(j);
            //std::cout<<" FROM CALC UUUU "<<std::endl;
            flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,qr_iminus,ql_i,m,flux_,quad_,u_mean,u_mean_plus)[0]*VdM_->L_(0,j)
                                            +gFlux_.computeNumFlux(ur_i, ul_iplus,qr_i,ql_iplus,m,flux_,quad_,u_mean,u_mean_plus)[0]   ;
            integ =quad_->IntFluxU([&](double u, double q) { return flux_.compute(u, q, m)[0]; }, i, j,
                                            grid_->faces(i), grid_->faces(i+1),u, q);
            // Apply the formula for the update of VdM_t_* pow(-1, j) 
            VdM_->VdM_t_(i,j) = (integ- flux_term)*double((2.0 * l + 1.0)/meshWidth_[0]);///(meshWidth_[0]);
            //std::cout<<" IN CALCUDT I " <<i<<" J "<<j<<" FLUX TERM "<<flux_term<<" INTEGRAL "<<integ<<" VDMT "<<VdM_->VdM_t_(i,j)<<" direct "<<double(meshWidth_[0])<<std::endl;
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

    // if(settings_.BarenblattM!=0){
    //     grid_->l2_error(0) = 0.0;
    //     for(int i =0;i<grid_->u_.size()[0];i++){
    //             grid_->l2_error(0) += pow((grid_->u(i)-initialCond_.computeInitialCondition(grid_->x(i),initCondA_,initCondB_,currentTime+1.0, settings_.BarenblattM)),2);        }
        
    //     grid_->linf_error(0) = 0.0;
    //     for(int i =0;i<grid_->u_.size()[0];i++){
    //         if(grid_->linf_error(0)<fabs(grid_->u(i)-initialCond_.computeInitialCondition(grid_->x(i),initCondA_,initCondB_,currentTime+1.0, settings_.BarenblattM))and fabs(grid_->x(i))>=1.5) 
    //             grid_->linf_error(0) = fabs(grid_->u(i)-initialCond_.computeInitialCondition(grid_->x(i),initCondA_,initCondB_,currentTime+1.0, settings_.BarenblattM));
    //     }
    // }else{
    //     grid_->l2_error(0) = 0.0;
    //     for(int i =0;i<grid_->u_.size()[0];i++){
    //         grid_->l2_error(0) += pow(grid_->u(i)-sin(grid_->x(i)-currentTime),2.0);//initialCond_.computeInitialCondition(grid_->x(i),initCondA_,initCondB_),2.0);
                                
    //     }
        
    //     grid_->linf_error(0) = 0.0;
    //     for(int i =0;i<grid_->u_.size()[0];i++){
    //         if(grid_->linf_error(0)<fabs(grid_->u(i)-initialCond_.computeInitialCondition(grid_->x(i),initCondA_,initCondB_))) 
    //             grid_->linf_error(0) = fabs(grid_->u(i)-initialCond_.computeInitialCondition(grid_->x(i),initCondA_,initCondB_));
    //     }
    // }

    std::cout<<"L2 Error: "<<sqrt(grid_->l2_error(0))<<" Linf Error: "<<grid_->linf_error(0)<<std::endl;
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

void Computation::calcUdt(const Array2D& u_){
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
            flux_term = -gFlux_.computeNumFlux(ur_iminus,ul_i,flux_)*VdM_->L_(0,j)  + gFlux_.computeNumFlux(ur_i, ul_iplus,flux_);
            double integ =quad_->IntFluxGaussLegendreQuad([&](double x) {return flux_.compute(x);}
                                                             ,i,j ,grid_->faces(i),grid_->faces(i+1),u_);
            VdM_->VdM_t_(i,j) =integ*(2.0*double(j)+1.0)*1/meshWidth_[0] - flux_term*1/meshWidth_[0]*(2.0*double(j)+1.0);
        }
    }
}