#include "computation.h"


//initialize the computation object, parse the settings from file that is given as the only command line argument 
void Computation::initialize(std::string filename)
{   
    //load settings
    settings_ = Settings();
    settings_.loadFromFile(filename);
    settings_.printSettings();

    epsilon_ = settings_.Epsilon;
    std::cout<<"Epsilon: "<<epsilon_<<std::endl;
    beta_ = settings_.Beta;
    delta_ = settings_.Delta;

    m_ = double(settings_.BarenblattM);
    aX_ = settings_.physicalSizeX[0];
    bX_ = settings_.physicalSizeX[1];
    aY_ = settings_.physicalSizeY[0];
    bY_ = settings_.physicalSizeY[1];
    initCondA_ = settings_.initCondA;
    initCondB_ = settings_.initCondB;

    firstLimiterCalls_ = 0;
    secondLimiterCalls_=0;

    PP_N_= settings_.PP_N;

    nNodes = PP_N_;

    nCells_= settings_.nCells;
    dt_ = settings_.CFL*1/(nCells_[0]);
    meshWidth_[0] =  (bX_-aX_)/double(nCells_[0]);
    meshWidth_[1] =  (bY_-aY_)/double(nCells_[1]);
    innerMeshWidth_[0] = meshWidth_[0]/(nNodes+2);
    innerMeshWidth_[1] = meshWidth_[1]/(nNodes+2);
    //initialize the Vandermonde Matrix
    std::array<int,3> t  ={nCells_[0],nCells_[1],PP_N_+1};
    VdM_ = std::make_shared<Vandermonde>(t,nNodes+2);

    //initialize Projection Operator
    Projection proj_ = Projection();

    //initialize the flux function;
    flux_ = Flux();
    //initialize the initial Condition
    initialCond_ = InitialCondition();
    //basis_ = std::make_unique<Basis>(PP_N_);


    grid_ = std::make_shared<Grid>(settings_.physicalSizeX,settings_.physicalSizeY,nCells_,meshWidth_, nNodes);

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
    else if (settings_.RiemannSolver == "enquist") {
        gFlux_.setNumFluxFunction(NumericalFlux::FunctionType::enquist);
        std::cout<< "Choosing the numerical EnquistOsher flux function..."<<std::endl;
    } 
    else if (settings_.RiemannSolver == "barenblatt") {
        gFlux_.setNumFluxFunction(NumericalFlux::FunctionType::porousMedia);
        std::cout<< "Choosing the numerical BarenBlattFlux flux function..."<<std::endl;
    } 
    else if (settings_.RiemannSolver == "central") {
        gFlux_.setNumFluxFunction(NumericalFlux::FunctionType::central);
        std::cout<< "Choosing the numerical Central flux function..."<<std::endl;
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
    }else if (settings_.initialCondition == "gaussian")
        {
            initialCond_.setInitialCondType(InitialCondition::InitialCondType::gaussian);
            std::cout<< "Choosing the gaussian as initial condition..."<<std::endl;
        }   
    else if (settings_.initialCondition == "exponential")
        {
            initialCond_.setInitialCondType(InitialCondition::InitialCondType::exponential);
            std::cout<< "Choosing the exponential as initial condition..."<<std::endl;
        }
    else {
        std::cout << "Initial Condition not set, choosing default" << std::endl;
    }
    if(settings_.useLimiter=="true"){
        useLimiter_ = true;
    }else{
        useLimiter_ = false;
    }
    std::cout <<R"(============================================================================================)"<<std::endl;

    // //grid_->solution_.printValues();
    // outputWriterParaview_ = std::make_unique<OutputWriterParaview>(grid_);
    // outputWriterParaview_->writeFile(0.0,settings_.OutputName);
}

void Computation::runSimulation()
{
    std::cout<<'-'<<std::flush;
    Timer timer;
    timer.start();
    double time_ = 0.0;
    int iter = 0;
    fillFaces();
    //std::cout<< "Faces: "<<std::endl;
    //grid_->faces_.printValues();
    //grid_->elemId.printValues();

   fillX();
//    std::cout<< "X: "<<std::endl;
    // grid_->x_.print();
//    std::cout<< "Y : "<<std::endl;  
//    grid_->y_.print();
   initVdm();
   initVdmJ();
   //VdM_->printValues();
   VdM_->LprintValues();
   // VdM_->LprimePrintValues();
   grid_->fillArray(grid_->u(),VdM_->VdM_,VdM_->L_);
   grid_->fillFaces(grid_->faceId,VdM_->VdM_,VdM_->L_);
   //std::cout<< "U: "<<std::endl;
   grid_->fillArray(grid_->j_,VdM_->VdMJ_,VdM_->L_);
   grid_->fillFaces(grid_->faceIdJ,VdM_->VdMJ_,VdM_->L_);

   //grid_->faceId.print();
   
   //grid_->u().print();
   grid_->fillSolution(grid_->solution_,grid_->u());
//    std::cout<< "Solution: "<<std::endl;
//    for(int iCell=0; iCell<grid_->nCells()[0]; iCell++){

//         for(int iNode=0; iNode<nNodes; iNode++){
//             for(int jNode=0; jNode<nNodes; jNode++){
//             std::cout<<"U: "<<grid_->u()(iNode,jNode,iCell)<<" ";
//             std::cout<<"Solution: ";
//             std::cout<<grid_->solution_(iNode,jNode*iCell)<<" ";
//             }
        
//         }
//     }
    // fillFaces();
    // fillX();
    // //grid_->x_.printValues();
    // fillXanalyze(grid_->x_analyze_);
    // initVdm();
    // grid_->fillArray(grid_->u(),VdM_->VdM_,VdM_->L_);
    // grid_->fillSolution(grid_->solution_,grid_->u());
    // initVdmJ();
    // grid_->fillArray(grid_->j(),VdM_->VdMJ_,VdM_->L_);
    // grid_->fillSolution(grid_->solutionJ_,grid_->j());
    // VdM_->LprintValues();
    // VdM_->LprimePrintValues();
    calcDt();
    double numberN = settings_.nWriteState/dt_*settings_.CFL;
    if(time_<dt_){
        grid_->fillSolution(grid_->solution_,grid_->j_,0);
        grid_->fillSolution(grid_->derivative_,grid_->ut_);
        outputWriterParaview_ = std::make_unique<OutputWriterParaview>(grid_);
        outputWriterParaview_->writeFile(time_,settings_.OutputName+"dX",grid_->solution_);
        grid_->fillSolution(grid_->solution_,grid_->j_,1);
        outputWriterParaview_->writeFile(time_,settings_.OutputName+"dY",grid_->solution_);

    }
    double errorTime =0.0;
    double numberofIterations = settings_.endTime/dt_;
    // calcError(errorTime);
    // // grid_->x_.printValues();
    // grid_->u_.printValues();
    // if(useLimiter_){
    //     std::cout<<"Using Limiter"<<std::endl;
    // firstLimiter(grid_->u());
    // //firstLimiter(grid_->j());
    // secondLimiter(grid_->u());
    // }
    // calcError(errorTime);
    // grid_->u_.printValues();
    while (time_<settings_.endTime)
        {
            if(flux_.getFluxFunction()!=Flux::FunctionType::Barenblatt){
                if(settings_.timeStepping=="euler" or settings_.timeStepping=="Euler"){
                    //std::cout<<"EULER TIME STEP"<<std::endl;
                    eulerTimeStep();

                }else if(settings_.timeStepping=="RK" or settings_.timeStepping=="rungeKutta" or settings_.timeStepping=="RungeKutta"){
                    rungeKutta();
                }
            }else{
                if(settings_.timeStepping=="euler" or settings_.timeStepping=="Euler"){
                    calcDt();
                    if(time_+dt_>settings_.endTime)
                        dt_ = settings_.endTime-time_;
                    eulerTimeStep();
                    // firstLimiter(grid_->u());
                    // secondLimiter(grid_->u());


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
                //grid_->fillSolution(grid_->solution_,grid_->j_);

                //grid_->prepareNodalDataForVisualization(VdM_->VdM_,quad_);

                //outputWriterParaview_->writeHighOrderFile(time_,settings_.OutputName, VdM_, quad_);
                outputWriterParaview_->writeFile(time_,settings_.OutputName,grid_->solution_);
                //std::cout<<"Write State TIME: "<<time_<<std::endl;
                //calcError(time_); 
            }

            // calcDt();

            std::cout<<"\rCurrent Time: "<<time_<<" End Time: "<<settings_.endTime<< std::flush;
    }
    std::cout<<" TIME: "<<time_<<std::endl;
    grid_->fillSolution(grid_->solution_,grid_->u_);
    grid_->fillSolution(grid_->derivative_,grid_->ut_);

    outputWriterParaview_->writeFile(time_,settings_.OutputName, grid_->solution_);
    outputWriterParaview_->writeFileTrueSolution(time_,settings_.OutputName+"TrueSolution");

    timer.stop();
    std::cout << "Elapsed time: " << timer.elapsedMilliseconds()/1000 << " s." << " nStates: "<<iter<<" firstLimiterCalls "<<firstLimiterCalls_<<" secondlimiterCalls "<<secondLimiterCalls_ <<std::endl;
}


void Computation::calcDt(){
    double max = 0.0;
    double maxj = 0.0;
    double dx =0.0;
    double c = 0.0;
    
    for(int iCell = 0; iCell<grid_->elemId.size()[0]; iCell++){
        for(int jCell = 0; jCell<grid_->elemId.size()[1]; jCell++){
            int idx = grid_->elemId(iCell,jCell);
            for(int node_i=0; node_i<grid_->u_.size()[0];node_i++){
                for(int node_j=0; node_j<grid_->u_.size()[1];node_j++){
                    double currentValue = grid_->u_(node_i,node_j,idx);
                    if(abs(currentValue)>max){
                        max = abs(currentValue);
                    }
                }
            }
            for(int faceIdx=0; faceIdx<grid_->faceId.size()[1];faceIdx++){
                for(int node_i=0; node_i<grid_->faceId.size()[0];node_i++){
                    double currentValue = grid_->faceId(node_i,faceIdx,idx);
                    if(abs(currentValue)>max){
                        max = abs(currentValue);
                    }
                }
            }
        }
    }
    for (int iCell = 0; iCell < grid_->elemId.size()[0]; iCell++) {
            for (int jCell = 0; jCell <grid_->elemId.size()[1]; jCell++) {
                // u^(1) = u^n + Δt * L_h(u^n)
                int idx = grid_->elemId(iCell,jCell);
                for(int node_i=0; node_i<grid_->u_.size()[0];node_i++){
                    for(int node_j=0; node_j<grid_->u_.size()[1];node_j++){
                        double currentValue1 = grid_->j_(node_i,node_j,idx,0);
                        double currentValue2 = grid_->j_(node_i,node_j,idx,1);
                        if(abs(currentValue1)>maxj){
                            maxj = abs(currentValue1);
                        }
                        if(abs(currentValue2)>maxj){
                            maxj = abs(currentValue2);
                        }
                    }
                }
                for(int faceIdx=0; faceIdx<grid_->faceId.size()[1];faceIdx++){
                    for(int node_i=0; node_i<grid_->faceId.size()[0];node_i++){
                        double currentValue = grid_->faceIdJ(node_i,faceIdx,idx); 
                        if(abs(currentValue)>maxj){
                            maxj = abs(currentValue);
                        }
                    }
                }
            }
        }
    dx = meshWidth_[0];    
    dt_ = settings_.CFL/(PP_N_+1)*(dx*dx/maxj*0.5*epsilon_+dx/max*sqrt(epsilon_));
    if(dt_>settings_.dt)
        dt_ = settings_.dt;
    //std::cout<<" DT "<<dt_<<std::endl;
}



void Computation::eulerTimeStep()
{   
    if(epsilon_==0.0){
        calcUdt(grid_->u_,grid_->faceId,VdM_->VdM_t_);
        grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
        grid_->fillFaces(grid_->face_dt,VdM_->VdM_t_,VdM_->L_);
        //std::cout<<" U "<<std::endl;
    }else{  
        calcJdt(grid_->u_,grid_->j_,grid_->faceId,VdM_->VdMJ_t_, epsilon_); // Compute the time derivative for j
        grid_->fillArray(grid_->jt_,VdM_->VdMJ_t_,VdM_->L_);
        grid_->fillFaces(grid_->face_dtJ,VdM_->VdMJ_t_,VdM_->L_);
        
        //firstLimiter(grid_->jt_);
        for (int iCell = 0; iCell < grid_->elemId.size()[0]; iCell++) {
            for (int jCell = 0; jCell <grid_->elemId.size()[1]; jCell++) {
                // u^(1) = u^n + Δt * L_h(u^n)
                int idx = grid_->elemId(iCell,jCell);
                for(int node_i=0; node_i<grid_->u_.size()[0];node_i++){
                    for(int node_j=0; node_j<grid_->u_.size()[1];node_j++){
                        grid_->j_(node_i,node_j,idx,0)  = grid_->j_(node_i,node_j,idx,0)+ dt_ * grid_->jt_(node_i, node_j,idx,0);
                        grid_->j_(node_i,node_j,idx,1)  = grid_->j_(node_i,node_j,idx,1)+ dt_ * grid_->jt_(node_i, node_j,idx,1);
                    }
                }
                for(int faceIdx=0; faceIdx<grid_->faceId.size()[1];faceIdx++){
                    for(int node_i=0; node_i<grid_->faceId.size()[0];node_i++){
                        grid_->faceIdJ(node_i,faceIdx,idx) =grid_->faceIdJ(node_i,faceIdx,idx)+ dt_*grid_->face_dtJ(node_i,faceIdx,idx); 
                    }
                }
            }
        }
        calcUdt(grid_->j_,grid_->faceIdJ,VdM_->VdM_t_); // Compute the time derivative for u^n
        grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
        grid_->fillFaces(grid_->face_dt,VdM_->VdM_t_,VdM_->L_);
    }
    for(int iCell = 0; iCell<grid_->elemId.size()[0]; iCell++){
        for(int jCell = 0; jCell<grid_->elemId.size()[1]; jCell++){
            int idx = grid_->elemId(iCell,jCell);
            for(int node_i=0; node_i<grid_->u_.size()[0];node_i++){
                for(int node_j=0; node_j<grid_->u_.size()[1];node_j++){
                    grid_->u_(node_i,node_j,idx) += dt_*grid_->ut_(node_i,node_j,idx);
                }
            }
            for(int faceIdx=0; faceIdx<grid_->faceId.size()[1];faceIdx++){
                for(int node_i=0; node_i<grid_->faceId.size()[0];node_i++){
                    grid_->faceId(node_i,faceIdx,idx) += dt_*grid_->face_dt(node_i,faceIdx,idx); 
                }
            }
        }
    }

}
void Computation::rungeKutta(){
    // Step 1: Compute the intermediate stage u^(1)
    if(epsilon_==0.0){
        calcUdt(grid_->u_,grid_->faceId,VdM_->VdM_t_); // Compute the time derivative for u^n
        grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
        grid_->fillFaces(grid_->face_dt,VdM_->VdM_t_,VdM_->L_);
    }else
    {
        calcJdt(grid_->u_,grid_->j_,grid_->faceId,VdM_->VdMJ_t_, epsilon_); // Compute the time derivative for j
        grid_->fillArray(grid_->jt_,VdM_->VdMJ_t_,VdM_->L_);
        grid_->fillFaces(grid_->face_dtJ,VdM_->VdMJ_t_,VdM_->L_);
        
        //firstLimiter(grid_->jt_);
        for (int iCell = 0; iCell < grid_->elemId.size()[0]; iCell++) {
            for (int jCell = 0; jCell <grid_->elemId.size()[1]; jCell++) {
                // u^(1) = u^n + Δt * L_h(u^n)
                int idx = grid_->elemId(iCell,jCell);
                for(int node_i=0; node_i<grid_->u_.size()[0];node_i++){
                    for(int node_j=0; node_j<grid_->u_.size()[1];node_j++){
                        grid_->j_1_(node_i,node_j,idx,0)  = grid_->j_(node_i,node_j,idx,0)+ dt_ * grid_->jt_(node_i, node_j,idx,0);
                        grid_->j_1_(node_i,node_j,idx,1)  = grid_->j_(node_i,node_j,idx,1)+ dt_ * grid_->jt_(node_i, node_j,idx,1);
                    }
                }
                for(int faceIdx=0; faceIdx<grid_->faceId.size()[1];faceIdx++){
                    for(int node_i=0; node_i<grid_->faceId.size()[0];node_i++){
                        grid_->faceIdJ1(node_i,faceIdx,idx) =grid_->faceIdJ(node_i,faceIdx,idx)+ dt_*grid_->face_dtJ(node_i,faceIdx,idx); 
                    }
                }
            }
        }
        calcUdt(grid_->j_1_,grid_->faceIdJ1,VdM_->VdM_t_); // Compute the time derivative for u^n
        grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
        grid_->fillFaces(grid_->face_dt,VdM_->VdM_t_,VdM_->L_);
    }

    //firstLimiter(grid_->jt_);

    //firstLimiter(grid_->jt_);
    for (int iCell = 0; iCell < grid_->elemId.size()[0]; iCell++) {
        for (int jCell = 0; jCell <grid_->elemId.size()[1]; jCell++) {
            // u^(1) = u^n + Δt * L_h(u^n)
            int idx = grid_->elemId(iCell,jCell);
            for(int node_i=0; node_i<grid_->u_.size()[0];node_i++){
                for(int node_j=0; node_j<grid_->u_.size()[1];node_j++){
                    grid_->u1_(node_i,node_j,idx)  = grid_->u_(node_i,node_j,idx)+ dt_ * grid_->ut_(node_i, node_j,idx);
                }
            }
            for(int faceIdx=0; faceIdx<grid_->faceId.size()[1];faceIdx++){
                for(int node_i=0; node_i<grid_->faceId.size()[0];node_i++){
                    grid_->faceId1(node_i,faceIdx,idx) =grid_->faceId(node_i,faceIdx,idx)+ dt_*grid_->face_dt(node_i,faceIdx,idx); 
                }
            }
        }
    }

    if(epsilon_==0.0){
        calcUdt(grid_->u1_,grid_->faceId1,VdM_->VdM_t_); // Compute the time derivative for u^n
        grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
        grid_->fillFaces(grid_->face_dt,VdM_->VdM_t_,VdM_->L_);
    }else
    {
        calcJdt(grid_->u1_,grid_->j_1_,grid_->faceId1,VdM_->VdMJ_t_, epsilon_); // Compute the time derivative for j
        grid_->fillArray(grid_->jt_,VdM_->VdMJ_t_,VdM_->L_);
        grid_->fillFaces(grid_->face_dtJ,VdM_->VdMJ_t_,VdM_->L_);
        
        //firstLimiter(grid_->jt_);
        for (int iCell = 0; iCell < grid_->elemId.size()[0]; iCell++) {
            for (int jCell = 0; jCell <grid_->elemId.size()[1]; jCell++) {
                // u^(1) = u^n + Δt * L_h(u^n)
                int idx = grid_->elemId(iCell,jCell);
                for(int node_i=0; node_i<grid_->u_.size()[0];node_i++){
                    for(int node_j=0; node_j<grid_->u_.size()[1];node_j++){
                        grid_->j_2_(node_i,node_j,idx,0)  = (3.0 / 4.0) *grid_->j_(node_i,node_j,idx,0)+(1.0 / 4.0) *grid_->j_1_(node_i,node_j,idx,0) 
                                                                                                    + (1.0 / 4.0) *dt_ * grid_->jt_(node_i, node_j,idx,0);
                        grid_->j_2_(node_i,node_j,idx,1)  = (3.0 / 4.0) *grid_->j_(node_i,node_j,idx,1)+(1.0 / 4.0) *grid_->j_1_(node_i,node_j,idx,1) 
                                                                                                    + (1.0 / 4.0) *dt_ * grid_->jt_(node_i, node_j,idx,1);
                    }
                }
                for(int faceIdx=0; faceIdx<grid_->faceId.size()[1];faceIdx++){
                    for(int node_i=0; node_i<grid_->faceId.size()[0];node_i++){
                        grid_->faceIdJ2(node_i,faceIdx,idx) =(3.0 / 4.0) *grid_->faceIdJ(node_i,faceIdx,idx)+(1.0 / 4.0) *grid_->faceIdJ1(node_i,faceIdx,idx) 
                                                                + (1.0 / 4.0)*dt_ * grid_->face_dtJ(node_i,faceIdx,idx); 
                    }
                }
            }
        }
        calcUdt(grid_->j_2_,grid_->faceIdJ2,VdM_->VdM_t_); // Compute the time derivative for u^n
        grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
        grid_->fillFaces(grid_->face_dt,VdM_->VdM_t_,VdM_->L_);
    }

    //firstLimiter(grid_->jt_);
    for (int iCell = 0; iCell < grid_->elemId.size()[0]; iCell++) {
        for (int jCell = 0; jCell <grid_->elemId.size()[1]; jCell++) {
            // u^(1) = u^n + Δt * L_h(u^n)
            int idx = grid_->elemId(iCell,jCell);
            for(int node_i=0; node_i<grid_->u_.size()[0];node_i++){
                for(int node_j=0; node_j<grid_->u_.size()[1];node_j++){
                    grid_->u2_(node_i,node_j,idx)  = (3.0 / 4.0) * grid_->u_(node_i,node_j,idx)+(1.0 / 4.0) * grid_->u1_(node_i,node_j,idx) 
                                                                + (1.0 / 4.0)*dt_ * grid_->ut_(node_i, node_j,idx);
                }
            }
            for(int faceIdx=0; faceIdx<grid_->faceId.size()[1];faceIdx++){
                for(int node_i=0; node_i<grid_->faceId.size()[0];node_i++){
                    grid_->faceId2(node_i,faceIdx,idx) = (3.0 / 4.0) *grid_->faceId(node_i,faceIdx,idx)+(1.0 / 4.0) *grid_->faceId1(node_i,faceIdx,idx) 
                                                                + (1.0 / 4.0)*dt_ * grid_->face_dt(node_i,faceIdx,idx); 
                }
            }
        }
    }
    if(epsilon_==0.0){
        calcUdt(grid_->u2_,grid_->faceId2,VdM_->VdM_t_); // Compute the time derivative for u^n
        grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
        grid_->fillFaces(grid_->face_dt,VdM_->VdM_t_,VdM_->L_);
    }else
    {
        calcJdt(grid_->u2_,grid_->j_2_,grid_->faceId2,VdM_->VdMJ_t_, epsilon_); // Compute the time derivative for j
        grid_->fillArray(grid_->jt_,VdM_->VdMJ_t_,VdM_->L_);
        grid_->fillFaces(grid_->face_dtJ,VdM_->VdMJ_t_,VdM_->L_);
        
        //firstLimiter(grid_->jt_);
        for (int iCell = 0; iCell < grid_->elemId.size()[0]; iCell++) {
            for (int jCell = 0; jCell <grid_->elemId.size()[1]; jCell++) {
                // u^(1) = u^n + Δt * L_h(u^n)
                int idx = grid_->elemId(iCell,jCell);
                for(int node_i=0; node_i<grid_->u_.size()[0];node_i++){
                    for(int node_j=0; node_j<grid_->u_.size()[1];node_j++){
                        grid_->j_(node_i,node_j,idx,0)  = (1.0 / 3.0) *grid_->j_(node_i,node_j,idx,0)+(2.0 / 3.0) *grid_->j_2_(node_i,node_j,idx,0) 
                                                                                                    + (2.0 / 3.0)*dt_ * grid_->jt_(node_i, node_j,idx,0);
                        grid_->j_(node_i,node_j,idx,1)  = (1.0 / 3.0) *grid_->j_(node_i,node_j,idx,1)+(2.0 / 3.0) *grid_->j_2_(node_i,node_j,idx,1) 
                                                                                                    + (2.0 / 3.0)*dt_ * grid_->jt_(node_i, node_j,idx,1);
                    }
                }
                for(int faceIdx=0; faceIdx<grid_->faceId.size()[1];faceIdx++){
                    for(int node_i=0; node_i<grid_->faceId.size()[0];node_i++){
                        grid_->faceIdJ(node_i,faceIdx,idx) =(1.0 / 3.0) *grid_->faceIdJ(node_i,faceIdx,idx)+(2.0 / 3.0) *grid_->faceIdJ2(node_i,faceIdx,idx) 
                                                                + (2.0 / 3.0) *dt_ * grid_->face_dtJ(node_i,faceIdx,idx); 
                    }
                }
            }
        }
        calcUdt(grid_->j_,grid_->faceIdJ,VdM_->VdM_t_); // Compute the time derivative for u^n
        grid_->fillArray(grid_->ut_,VdM_->VdM_t_,VdM_->L_);
        grid_->fillFaces(grid_->face_dt,VdM_->VdM_t_,VdM_->L_);
    }

    //firstLimiter(grid_->jt_);
    for (int iCell = 0; iCell < grid_->elemId.size()[0]; iCell++) {
        for (int jCell = 0; jCell <grid_->elemId.size()[1]; jCell++) {
            // u^(1) = u^n + Δt * L_h(u^n)
            int idx = grid_->elemId(iCell,jCell);
            for(int node_i=0; node_i<grid_->u_.size()[0];node_i++){
                for(int node_j=0; node_j<grid_->u_.size()[1];node_j++){
                    grid_->u_(node_i,node_j,idx)  = (1.0 / 3.0) *grid_->u_(node_i,node_j,idx)+(2.0 / 3.0) *grid_->u2_(node_i,node_j,idx) 
                                                                + (2.0 / 3.0)*dt_ * grid_->ut_(node_i, node_j,idx);
                }
            }
            for(int faceIdx=0; faceIdx<grid_->faceId.size()[1];faceIdx++){
                for(int node_i=0; node_i<grid_->faceId.size()[0];node_i++){
                    grid_->faceId(node_i,faceIdx,idx) =(1.0 / 3.0) *grid_->faceId(node_i,faceIdx,idx)+(2.0 / 3.0) *grid_->faceId2(node_i,faceIdx,idx) 
                                                                + (2.0 / 3.0)*dt_ * grid_->face_dt(node_i,faceIdx,idx); 
                }
            }
        }
    }
}


void Computation::fillFaces()
{
    for (int i = 0; i < grid_->faces_.size()[0]; i++)
    {
        for(int k = 0; k<grid_->faces_.size()[1];k++){
            if(k==0){
                grid_->faces_(i,k) = aX_+i*meshWidth_[0];
            }else if(k==1){
                grid_->faces_(i,k) = aY_+i*meshWidth_[1];
            }
        }
    }

    for(int i = 0; i < grid_->elemId.size()[0]; i++)
    {
        for(int k = 0; k<grid_->elemId.size()[1];k++){
            grid_->elemId(i,k) = i+k*grid_->elemId.size()[0];
        }
    }
}


void Computation::fillX()
{
    double transformedNode = 0.0;
    double mean = 0.0, diff  =0.0;
 
    for(int l=0; l<grid_->x_.size()[2]; l++){
        for (int i = 0; i < grid_->x_.size()[0]; i++)
        {
            for(int j = 0; j < grid_->x_.size()[1]; j++)
            {
                for(int k = 0; k < grid_->faces_.size()[2]+1; k++)
                {
                    if(k==0){
                        mean = 0.5 * (grid_->faces_(l+1,k) + grid_->faces_(l,k));
                        diff = 0.5 * (grid_->faces_(l+1,k) - grid_->faces_(l,k));
                        transformedNode = mean + diff*quad_->basis_.nodes(i);
                        grid_->x_(i,j,l) = transformedNode;
                    }else if(k==1){
                        mean = 0.5 * (grid_->faces_(l+1,k) + grid_->faces_(l,k));
                        diff = 0.5 * (grid_->faces_(l+1,k) - grid_->faces_(l,k));
                        transformedNode = mean + diff*quad_->basis_.nodes(j);
                        grid_->y_(i,j,l) = transformedNode;
                    }
                }
            }
        }
    }   
}


void Computation::initVdm() {
    // Iterate over grid cells and polynomial degrees
    for(int deg = 0; deg < VdM_->VdM_.size()[2];deg++){
        for (int i = 0; i < VdM_->VdM_.size()[0]; i++) {
            for (int k = 0; k < VdM_->VdM_.size()[1]; k++) {
                // Determine integral lambda function
                double integral = 0.0;
                    integral = quad_->IntGaussLegendreQuad(
                        [&](double x, double y) {
                            return initialCond_.computeInitialCondition2D(x,y, initCondA_, initCondB_,settings_.BarenblattTime, settings_.BarenblattM);
                        },
                        deg, grid_->faces(i,0), grid_->faces(i+1,0), grid_->faces(k,1), grid_->faces(k+1,1));
                
                // Scale integral and store in VdM_
                VdM_->VdM_(i, k,deg) = integral * (2.0 * double(deg) + 1.0) / meshWidth_[0]* (2.0 * double(deg) + 1.0) / meshWidth_[0];
                // Compute and store Legendre polynomials
                for (int p = 0; p < quad_->basis_.nodes_.size()[0]; p++) {
                    std::array<double,2> L = quad_->LegendrePolynomialAndDerivative(deg, quad_->basis_.nodes(p));
                    VdM_->L(p, deg) = L[0];
                    VdM_->L_prime(p, deg) = L[1];
                }
            }
        }
    }
}

void Computation::initVdmJ(){
    for (int iCell = 0; iCell < grid_->elemId.size()[0] ; iCell++) {
        for(int jCell = 0; jCell < grid_->elemId.size()[1]; jCell++) {
            for(int deg = 0; deg <=PP_N_; deg++) {
                double surface_intX = quad_->IntSurfaceLegendreGaussX(
                    [&](double x,double y) {
                        return initialCond_.computeInitialCondition2D(x,y, initCondA_, initCondB_,settings_.BarenblattTime, settings_.BarenblattM);
                    },
                    deg, grid_->faces(iCell,0), grid_->faces(iCell+1,0), grid_->faces(jCell,1), grid_->faces(jCell+1,1),m_);
                double surface_intY = quad_->IntSurfaceLegendreGaussY(
                    [&](double x,double y) {
                        return initialCond_.computeInitialCondition2D(x,y, initCondA_, initCondB_,settings_.BarenblattTime, settings_.BarenblattM);
                    },
                    deg, grid_->faces(iCell,0), grid_->faces(iCell+1,0), grid_->faces(jCell,1), grid_->faces(jCell+1,1),m_);
                double integ = quad_->IntGaussLegendreQuadDeriv(
                    [&](double x, double y) {
                        return initialCond_.computeInitialCondition2D(x,y, initCondA_, initCondB_,settings_.BarenblattTime, settings_.BarenblattM);
                    },
                    deg, grid_->faces(iCell,0), grid_->faces(iCell+1,0), grid_->faces(jCell,1), grid_->faces(jCell+1,1));                
                    VdM_->VdMJ_(iCell,jCell,deg,0) = (integ-surface_intX)*1/meshWidth_[0]*(2.0*double(deg)+1.0)*1/meshWidth_[0]*(2.0*double(deg)+1.0);
                    VdM_->VdMJ_(iCell,jCell,deg,1) = (integ-surface_intY)*1/meshWidth_[0]*(2.0*double(deg)+1.0)*1/meshWidth_[0]*(2.0*double(deg)+1.0);
            }
        }
    }
}

void Computation::calcJdt(const Array3D& u, const Array4D& j, Array3D& faceId ,Array4D& VdMJ_t, double epsilon){
    for (int iCell = 0; iCell < grid_->elemId.size()[0] ; iCell++) {
        for(int jCell = 0; jCell < grid_->elemId.size()[1]; jCell++) {

            double left = grid_->faces(iCell,0);
            double right = grid_->faces(iCell+1,0);
            double bottom = grid_->faces(jCell,1);
            double top = grid_->faces(jCell+1,1);

            gFlux_.fillFluxArray(flux_,faceId,grid_->faceFluxJ_,grid_->elemId,iCell,jCell,dt_,meshWidth_[0]);

            for (int deg = 0; deg <=PP_N_; deg++) {
                double flux_termX = quad_->surfaceInt2DX(deg,grid_->faces_(iCell,1),grid_->faces_(iCell+1,1),grid_->faceFluxJ_,m_);
                double flux_termY = quad_->surfaceInt2DY(deg,grid_->faces_(jCell,0),grid_->faces_(jCell+1,0),grid_->faceFluxJ_,m_);

                double integX =quad_->volumeInt2DX([&](double x) {return flux_.compute(x);},deg,iCell,jCell,
                                            left,right,bottom,top,grid_->elemId,u,m_) - quad_->volumeInt2DJ(0,deg,iCell,jCell,left,right,bottom,top,grid_->elemId,j);
                double integY =quad_->volumeInt2DY([&](double x) {return flux_.compute(x);},deg,iCell,jCell,
                                            left,right,bottom,top,grid_->elemId,u,m_) - quad_->volumeInt2DJ(1,deg,iCell,jCell,left,right,bottom,top,grid_->elemId,j);                                       
                
                VdMJ_t(iCell,jCell,deg,0) =1.0/epsilon*(integX-flux_termX)*1/meshWidth_[0]*(2.0*double(deg)+1.0)*1/meshWidth_[1]*(2.0*double(deg)+1.0);
                VdMJ_t(iCell,jCell,deg,1) =1.0/epsilon*(integY-flux_termY)*1/meshWidth_[0]*(2.0*double(deg)+1.0)*1/meshWidth_[1]*(2.0*double(deg)+1.0);
            }
        }
    }
}

void Computation::calcUdt(const Array4D& j, Array3D& faceId ,Array3D& VdM_t){
    for (int iCell = 0; iCell < grid_->elemId.size()[0] ; iCell++) {
        for(int jCell = 0; jCell < grid_->elemId.size()[1]; jCell++) {
            gFlux_.fillFluxArray(flux_,faceId,grid_->faceFlux_,grid_->elemId,iCell,jCell,dt_,meshWidth_[0]);

            for (int deg = 0; deg <=PP_N_; deg++) {
                double flux_term = quad_->surfaceInt2D(deg,grid_->faces_(iCell,0),grid_->faces_(iCell+1,0),grid_->faceFlux_);

                double integ =quad_->volumeInt2D([&](double x) {return flux_.compute(x);},deg,iCell,jCell,
                                             grid_->faces(iCell,0),grid_->faces(iCell+1,0),grid_->faces(jCell+1,1),grid_->faces(jCell,1),grid_->elemId,j);
                       
                VdM_t(iCell,jCell,deg) =(integ-flux_term)*1/meshWidth_[0]*(2.0*double(deg)+1.0)*1/meshWidth_[1]*(2.0*double(deg)+1.0);

            }
        }
    }
}

void Computation::calcUdt(const Array3D& u_, Array3D& faceId ,Array3D& VdM_t){
    for (int iCell = 0; iCell < grid_->elemId.size()[0] ; iCell++) {
        for(int jCell = 0; jCell < grid_->elemId.size()[1]; jCell++) {
            gFlux_.fillFluxArray(flux_,faceId,grid_->faceFlux_,grid_->elemId,iCell,jCell,dt_,meshWidth_[0]);

            for (int deg = 0; deg <=PP_N_; deg++) {
                double flux_term = quad_->surfaceInt2D(deg,grid_->faces_(iCell,0),grid_->faces_(iCell+1,0),grid_->faceFlux_);

                double integ =quad_->volumeInt2D([&](double x) {return flux_.compute(x);},deg,iCell,jCell,
                                            grid_->faces(iCell,0),grid_->faces(iCell+1,0),grid_->faces(jCell+1,1),grid_->faces(jCell,1),grid_->elemId,u_);
                       
                VdM_->VdM_t_(iCell,jCell,deg) =(integ-flux_term)*1/meshWidth_[0]*(2.0*double(deg)+1.0)*1/meshWidth_[0]*(2.0*double(deg)+1.0);

            }
        }
    }
}

// void Computation::fillXanalyze(Array2D &x)
// {
//     double stepsize_ = (settings_.initCondB-settings_.initCondA)/double(grid_->x_analyze_.size()[0]-1);
//     for(int i=0;i<grid_->x_analyze_.size()[0];i++){
//         x(i,0) = settings_.initCondA+i*stepsize_;
//     }
// }

// void Computation::fillUanalyze(Array2D &u_analyze, const Array2D &x, const Array2D &Vdm)
// {
//     for(int i =0;i<x.size()[0];i++){
//         u_analyze(i,0) = 0.0;
//         for(int k=0;k<grid_->faces().size()[0]-1;k++){
//             if(x(i,0)>=grid_->faces(k) and x(i,0)<grid_->faces(k+1)){
//                 double transformedNode =  (x(i,0) - 0.5 * (grid_->faces(k+1) + grid_->faces(k)))/(meshWidth_[0]*0.5);
//                 for(int j = 0; j<Vdm.size()[1];j++){
//                     u_analyze(i,0) += Vdm(k,j)*quad_->LegendrePolynomialAndDerivative(j,transformedNode)[0];
//                 }
//                 break;
//             }
//         }
//     }


// void Computation::calcError(double currentTime)
// {
//     fillUanalyze(grid_->u_analyze_,grid_->x_analyze_,VdM_->VdM_);

//     if(flux_.getFluxFunction()==Flux::FunctionType::Barenblatt){
//         for(int i = 0; i<grid_->u_.size()[0];i++){
//             for(int j=0; j<grid_->u_.size()[1];j++){
//                 grid_->true_solution_(i,j) =initialCond_.computeInitialCondition(grid_->x_(i,j),initCondA_,initCondB_,currentTime+settings_.BarenblattTime, settings_.BarenblattM);     
//             }
//         }
//         grid_->l2_error(0) = 0.0;
//         for(int i =0;i<grid_->u_analyze_.size()[0];i++){
//             grid_->u_analyze_true_(i,0) = initialCond_.computeInitialCondition(grid_->x_analyze_(i,0),initCondA_,initCondB_,currentTime+settings_.BarenblattTime, settings_.BarenblattM);
//             grid_->l2_error(0) += pow(grid_->u_analyze_(i,0)- grid_->u_analyze_true_(i,0),2.);     
//         }
        
//         grid_->linf_error(0) = 0.0;

//         for(int i =0;i<grid_->u_analyze_.size()[0];i++){
//             if(grid_->linf_error(0)<abs(grid_->u_analyze_(i,0)-initialCond_.computeInitialCondition(grid_->x_analyze_(i,0),initCondA_,initCondB_,currentTime+settings_.BarenblattTime, settings_.BarenblattM))) 
//                 grid_->linf_error(0) = abs(grid_->u_analyze_(i,0)-initialCond_.computeInitialCondition(grid_->x_analyze_(i,0),initCondA_,initCondB_,currentTime+settings_.BarenblattTime, settings_.BarenblattM));
//         }
//     }
//     else if(settings_.initialCondition =="sinus"){
//         for(int i = 0; i<grid_->u_.size()[0];i++){
//             for(int j=0; j<grid_->u_.size()[1];j++){
//                 grid_->true_solution_(i,j) =sin(grid_->x_(i,j));     
//             }
//         }
//         grid_->l2_error(0) = 0.0;
//         grid_->linf_error(0) = 0.0;
//         for(int i =0;i<grid_->u_analyze_.size()[0];i++){
//             double sol = initialCond_.computeInitialCondition(grid_->x_analyze_(i,0),initCondA_,initCondB_);
//             //std::cout<<" SOL "<<sol<<" "<<grid_->u_analyze_(i,0)<<std::endl;
//             grid_->l2_error(0) += pow(grid_->u_analyze_(i,0)-sol ,2);     
//         if(grid_->linf_error(0)<abs(grid_->u_analyze_(i,0)-sol))
//                 grid_->linf_error(0) = abs(grid_->u_analyze_(i,0)-sin(grid_->x_analyze_(i,0)));
//         }
//     }
//     std::cout<<"L2 Error: "<<sqrt(grid_->l2_error(0)*1./double(grid_->u_analyze_.size()[0]))<<" Linf Error: "<<grid_->linf_error(0)<<std::endl;
// }

// void Computation::calcQ(const Array3D& u, Array2D& faceIdQ)
// {
//     double m = double(settings_.BarenblattM);
//     for (int iCell = 0; iCell < grid_->elemId.size()[0]; iCell++) {
//         for (int jCell = 0; jCell < grid_->elemId.size()[0]; jCell++) {
//         double flux_term =0.0,ul_i = 0.0,ur_i = 0.0,ur_iminus  = 0.0,ul_iplus =0.0;
//         double u_mean = 0.0,u_mean_plus = 0.0, u_mean_minus = 0.0;
//         // Wrap around the grid for periodic boundary conditions
//         // void fillFluxArray(bool QTrue,const std::unique_ptr<Quadrature> &quad_, Array3D& u, Flux flux, Array3D &faceId, Array3D &faceIdQ, Array2D &faceFluxQ
//         //     , const Array2D &elemId, int iCell, int jCell, double dt, double meshWidth, double m);
//         gFlux_.fillFluxArray(true,quad_,u,flux_,grid_->faceId,grid_->faceIdQ,grid_->faceFlux_,grid_->faceFluxQ_,grid_->elemId,iCell,jCell,m_);
//         int idx = grid_->elemId(iCell,jCell);
//         //std::cout<<" FROM CALCQ "<<" g_jplus "<<" gFlux_plus "<<gFlux_plus<<std::endl;
//         for (int l = 0; l <=PP_N_; l++) {
//             double deg = double(l); 
//             // Compute the numerical flux
//             //std::cout<<" FROM CALCQ "<<std::endl;
//             flux_term = quad_->IntFaceFluxQ([&](double x) {return flux_.compute(x,0.0,m)[1];},deg,grid_->faces(iCell,0),grid_->faces(iCell+1,0),faceIdQ);
//             // Apply the formula for the update of VdM_t_* pow(-1, j) 
//             double a = grid_->faces(iCell,0);
//             double b = grid_->faces(iCell + 1,0);
//             double integ =quad_->IntFluxQ([&](double x) {return flux_.compute(x, 0.0, m)[1];},deg,iCell, jCell, a, b, u,grid_->elemId);
//             // if(std::abs(flux_term)<1E-12)
//             //     flux_term=0.0;
//             // if(std::abs(integ)<1E-12)
//             //     integ=0.0;
//             VdM_->VdMQ_(iCell,jCell,deg) = (integ-flux_term)*(2. * deg + 1.)/meshWidth_[0]*(2. * deg + 1.)/meshWidth_[0];
//             // std::cout<<" IN CALCQ I "<<" i "<<i<<" Face i "<<grid_->faces(i)<<" Face i+1 "<<grid_->faces(i+1)<<" J "<<j<<" FLUX TERM "<<flux_term<<" INTEGRAL "<<integ<<" VDMQ "<<VdM_->VdMQ_(i,j)<<std::endl;
//             // std::cout<<" IN CALCQ I " <<i<<" J "<<j<<" gFlux_minus "<<gFlux_minus<<" gFlux_plus "<<gFlux_plus<<" L "<<VdM_->L(0,j)<<std::endl;  
//         }
//     }
//     grid_->fillArray(grid_->q_,VdM_->VdMQ_,VdM_->L_);
//     grid_->fillFaces(grid_->faceIdQ,VdM_->VdMQ_,VdM_->L_);
//     // std::cout<<" Q "<<std::endl;
//     // grid_->q_.printValues();
//     }
// }
