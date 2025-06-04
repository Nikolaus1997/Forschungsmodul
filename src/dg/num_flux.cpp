#include "num_flux.h"



// Set the numerical flux function type
void NumericalFlux::setNumFluxFunction(FunctionType type)
{
    selectedFunction = type;
}

void NumericalFlux::fillFluxArray(Flux flux, Array3D &faceId, Array2D &faceFlux, const Array2D &elemId, int iCell, int jCell, double dt, double meshWidth)
{
    int iElem_0, iElem_2, iElem_3, iElem_1;
    int iElem = elemId(iCell,jCell);
    if(iCell > 0 and iCell < elemId.size()[0]-1){
        iElem_0 = elemId(iCell-1,jCell);
        iElem_2 = elemId(iCell+1,jCell);
    }
    else if(iCell == 0){
        iElem_0 = elemId(elemId.size()[0]-1,jCell);
        iElem_2 = elemId(iCell+1,jCell);
    }
    else if(iCell == elemId.size()[0]-1){
        iElem_0 = elemId(iCell-1,jCell);
        iElem_2 = elemId(0,jCell);
    }
    if(jCell> 0 and jCell < elemId.size()[1]-1){
        iElem_1 = elemId(iCell,jCell-1);
        iElem_3 = elemId(iCell,jCell+1);
    }else if(jCell == 0){
        iElem_1 = elemId(iCell,elemId.size()[1]-1);
        iElem_3 = elemId(iCell,jCell+1);
    }else if(jCell == elemId.size()[1]-1){
        iElem_1 = elemId(iCell,jCell-1);
        iElem_3 = elemId(iCell,0);
    }
    for(int i = 0; i < faceId.size()[0]; ++i) {
        for (int j = 0; j < faceId.size()[1]; ++j) {
            if(j == 0){
                faceFlux(i,j) = computeNumFlux(faceId(i,j,iElem_0),faceId(i,j,iElem),flux,dt,meshWidth);
            }
            else if(j == 1){
                faceFlux(i,j) = computeNumFlux(faceId(i,j,iElem_1),faceId(i,j,iElem),flux,dt,meshWidth);
            }
            else if(j == 2){
                faceFlux(i,j) = computeNumFlux(faceId(i,j,iElem),faceId(i,j,iElem_2),flux,dt,meshWidth);
            }
            else if(j == 3){
                faceFlux(i,j) = computeNumFlux(faceId(i,j,iElem),faceId(i,j,iElem_3),flux,dt,meshWidth);
            }
        }

    }
}
// Compute the numerical flux based on the selected function
double NumericalFlux::computeNumFlux(double u_l, double u_r,Flux flux_, double dt, double meshWidth)
{
    switch (selectedFunction)
    {
    case FunctionType::upwind:
        return upwind(u_l, u_r,flux_);
    case FunctionType::downwind:
        return downwind(u_l, u_r,flux_);
    case FunctionType::lax:
        return lax(u_l, u_r,flux_, dt, meshWidth);
    case FunctionType::enquist:
        return enquist(u_l, u_r,flux_);
    case FunctionType::central:
        return 0.5*(flux_.compute(u_l) + flux_.compute(u_r));
    default:
        throw std::invalid_argument("Unknown numerical flux function type.");
    }
}

std::array<double, 2> NumericalFlux::computeNumFlux(bool minus,double x_l, double x_r, double q_l, double q_r, double m, Flux flux_,const std::unique_ptr<Quadrature>& quad_, double u_mean, double u_mean_plus)
{
    if(minus)
        return porousMediaMinus(x_l,x_r,q_l,q_r,m,flux_,quad_, u_mean,u_mean_plus);
    else
        return porousMediaPlus(x_l,x_r,q_l,q_r,m,flux_,quad_, u_mean,u_mean_plus);
}

// Upwind flux function
double NumericalFlux::upwind(double u_l, double u_r,Flux flux_)
{
    // Example implementation: use left flux value
    return flux_.compute(u_l);
}

// Downwind flux function
double NumericalFlux::downwind(double u_l, double u_r, Flux flux_)
{
    // Example implementation: use right flux value
    return flux_.compute(u_r);
}

// Lax flux function
double NumericalFlux::lax(double u_l, double u_r, Flux flux_, double dt, double meshWidth)
{
    // Example implementation: average of left and right
    return 0.5*(flux_.compute(u_l) + flux_.compute(u_r)) + 0.5 * (u_r - u_l);
    //return 0.0;
}

// Enquist flux function
double NumericalFlux::enquist(double u_l, double u_r, Flux flux_)
{

    if(abs(u_l) > abs(u_r)) 
        return u_l;
    else
        return u_r;
    return 0.0;
}

void NumericalFlux::fillFluxArray(bool QTrue,const std::unique_ptr<Quadrature>& quad_,const Array3D& u,Flux flux, Array3D& faceId, Array3D& faceIdQ, Array2D& faceFlux, Array2D& faceFluxQ, const Array2D& elemId,int iCell, int jCell, double m)
{   
    int iElem_0, iElem_2, iElem_3, iElem_1;
    int iElem = elemId(iCell,jCell);
    if(iCell > 0 and iCell < elemId.size()[0]-1){
        iElem_0 = elemId(iCell-1,jCell);
        iElem_2 = elemId(iCell+1,jCell);
    }
    else if(iCell == 0){
        iElem_0 = elemId(elemId.size()[0]-1,jCell);
        iElem_2 = elemId(iCell+1,jCell);
    }
    else if(iCell == elemId.size()[0]-1){
        iElem_0 = elemId(iCell-1,jCell);
        iElem_2 = elemId(0,jCell);
    }
    if(jCell> 0 and jCell < elemId.size()[1]-1){
        iElem_1 = elemId(iCell,jCell-1);
        iElem_3 = elemId(iCell,jCell+1);
    }else if(jCell == 0){
        iElem_1 = elemId(iCell,elemId.size()[1]-1);
        iElem_3 = elemId(iCell,jCell+1);
    }else if(jCell == elemId.size()[1]-1){
        iElem_1 = elemId(iCell,jCell-1);
        iElem_3 = elemId(iCell,0);
    }

    for (int j = 0; j < faceId.size()[1]; ++j) {
        for(int i = 0; i < faceId.size()[0]; ++i) {        
            if(j == 0){
                double u_mean = 0.0, u_mean_plus = 0.0;
                for(int k = 0; k<u.size()[0];k++){
                    u_mean += u(k,i,iElem_0);
                    u_mean_plus += u(k,i,iElem);
                }
                u_mean += u(i,j,iElem_0)+u(i,j,iElem_0);
                u_mean_plus += u(i,j,iElem)+u(i,j,iElem);
                u_mean = u_mean/(u.size()[0]+2);
                u_mean_plus = u_mean_plus/(u.size()[0]+2);
                if(QTrue){
                        faceFluxQ(i,j) = porousMediaMinus( faceId(i,j,iElem_0),faceId(i,j,iElem), 0.0,0.0,
                                                            m, flux,quad_, u_mean, u_mean_plus)[1];
                }else{
                    faceFlux(i,j) = porousMediaMinus(faceId(i,j,iElem_0),faceId(i,j,iElem), faceIdQ(i,j,iElem_0),
                                                    faceIdQ(i,j,iElem), m, flux,quad_, u_mean, u_mean_plus)[0];
                }
            }
            else if(j == 1){
                double u_mean = 0.0, u_mean_plus = 0.0;
                for(int k = 0; k<u.size()[1];k++){
                    u_mean += u(i,k,iElem_1);
                    u_mean_plus += u(i,k,iElem);
                }
                u_mean += u(i,j,iElem_1)+u(i,j,iElem_1);
                u_mean_plus += u(i,j,iElem)+u(i,j,iElem);
                u_mean = u_mean/(u.size()[1]+2);
                u_mean_plus = u_mean_plus/(u.size()[1]+2);
                if(QTrue){
                    faceFluxQ(i,j) = porousMediaMinus(faceId(i,j,iElem_1),faceId(i,j,iElem),0.0,0.0,
                                                         m, flux,quad_, u_mean, u_mean_plus)[1];
                }
                else{
                    faceFlux(i,j) = porousMediaMinus(faceId(i,j,iElem_1),faceId(i,j,iElem), faceIdQ(i,j,iElem_1),
                                                    faceIdQ(i,j,iElem), m, flux,quad_, u_mean, u_mean_plus)[0];
                }
            }
            else if(j == 2){
                double u_mean = 0.0, u_mean_plus = 0.0;
                for(int k = 0; k<u.size()[0];k++){
                    u_mean += u(k,i,iElem_2);
                    u_mean_plus += u(k,i,iElem);
                }
                u_mean_plus += u(i,j,iElem_2)+u(i,j,iElem_2);
                u_mean += u(i,j,iElem)+u(i,j,iElem);
                u_mean = u_mean/(u.size()[0]+2);
                u_mean_plus = u_mean_plus/(u.size()[0]+2);
                if(QTrue){
                    faceFluxQ(i,j) = porousMediaPlus(faceId(i,j,iElem),faceId(i,j,iElem_2), 0.0,0.0,
                                                        m, flux,quad_, u_mean, u_mean_plus)[1];
                }
                else{
                    faceFlux(i,j) = porousMediaPlus(faceId(i,j,iElem),faceId(i,j,iElem_2), faceIdQ(i,j,iElem),
                                                    faceIdQ(i,j,iElem_2), m, flux,quad_, u_mean, u_mean_plus)[0];
                }
            }
            else if(j == 3){
                double u_mean = 0.0, u_mean_plus = 0.0;
                for(int k = 0; k<u.size()[1];k++){
                    u_mean_plus += u(i,k,iElem_3);
                    u_mean += u(i,k,iElem);
                }
                u_mean_plus += u(i,j,iElem_3)+u(i,j,iElem_3);
                u_mean += u(i,j,iElem)+u(i,j,iElem);
                u_mean = u_mean/(u.size()[1]+2);
                u_mean_plus = u_mean_plus/(u.size()[1]+2);
                if(QTrue){
                    faceFluxQ(i,j) = porousMediaPlus(faceId(i,j,iElem),faceId(i,j,iElem_3), 0.0,0.0,
                                                    m, flux,quad_, u_mean, u_mean_plus)[1];
                }
                else{
                    faceFlux(i,j) = porousMediaPlus(faceId(i,j,iElem),faceId(i,j,iElem_3), faceIdQ(i,j,iElem),
                                                    faceIdQ(i,j,iElem_3), m, flux,quad_, u_mean, u_mean_plus)[0];
                }
            }
        }
    }
}

std::array<double,2> NumericalFlux::porousMediaPlus(double u_l, double u_r, double q_l, double q_r, double m, Flux flux_,const std::unique_ptr<Quadrature>& quad_, double u_mean, double u_mean_plus)
{   
    double g_plus =0.0,g_minus =0.0;
    double u_diff = u_r-u_l;
    double q_diff = q_r-q_l;
    double q_mean = 0.5*(q_r+q_l);
    double fraction_term =0.0;
    double gamma = 0.0;
    double u_mean_diff = u_mean_plus-u_mean;
    double g_mean_plus = 0.0,g_mean_minus = 0.0;
    double g_diff=0.0;
    if((std::abs(u_diff))<1E-14){
        //std::cout<<"IM IN THE IF"<<std::endl;
        fraction_term = -1.0*flux_.compute(u_l,0,m)[1];
        // g_plus= quad_->G(u_r,m);//-1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_r);
        // g_minus= quad_->G(u_l,m);//-1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_l);
        g_plus= -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_r);
        g_minus= -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_l);

        //std::cout<<"IN IF fraction term "<<fraction_term<<" UDIFF "<< u_diff<<" G_plus "<<g_plus<<" u_R "<<u_r<<" G_minus "<<g_minus<<" U_l"<<u_l <<std::endl;
    }else{
        // g_plus= quad_->G(u_r,m);//-1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_r);
        // g_minus= quad_->G(u_l,m);//-1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_l);

        g_plus= -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_r);
        g_minus= -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_l);

        //std::cout<<"IN ELSE fraction term "<<fraction_term<<" UDIFF "<< u_diff<<" G_plus "<<g_plus<<" u_R "<<u_r<<" G_minus "<<g_minus<<" U_l"<<u_l<<" u_mean_plus "<< u_mean_plus<< " u_mean "<<u_mean <<std::endl;

        if(abs(u_r-u_mean_plus)<1E-14 ) {
            u_mean_plus = u_r;       
        }
        if(abs(u_l-u_mean)<1E-14 ) {
            u_mean = u_l;
        }
        g_diff = g_plus-g_minus;
        fraction_term = (g_diff)/u_diff;
    }
    if((std::abs(u_mean_diff))<1E-14){
        gamma = -0.5*flux_.compute(u_mean,0,m)[1];
    }else{
        // g_mean_plus = quad_->G(u_mean_plus,m);//-1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_mean_plus);
        // g_mean_minus = quad_->G(u_mean,m);//-1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_mean);
 
        g_mean_plus = -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_mean_plus);
        g_mean_minus = -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_mean);

        gamma = 0.5*(g_mean_plus-g_mean_minus)/(u_mean_plus-u_mean);
    }
    double g_mean = 0.5*(g_plus+g_minus);
    // std::cout<< " G MEAN "<<g_mean<<" GAMMA "<<gamma<<" sol0 "<<-fraction_term*q_mean-gamma*q_diff<<" sol1 "<<gamma*u_diff<<" u_meandiff "<<u_mean_diff<<std::endl;
    // std::cout<<" IN NUMERICAL FLUXPLUS"<<" u_l "<<u_l<<" u_r "<<u_r <<" udiff " << u_diff<<" g_mean "<<g_mean<<std::endl;
    // std::cout<<" IN NUMERICAL FLUXPLUS"<<" q_l "<<q_l<<" q_r "<<q_r <<" q_diff "<<q_diff<<" q_mean "<<q_mean<<" g_plus "<<g_plus<<" g_minus "<<g_minus<<" gamma "<<gamma<<std::endl;
    return {-fraction_term*q_mean-gamma*q_diff,-g_mean+gamma*u_diff};
}

std::array<double, 2> NumericalFlux::porousMediaMinus(double u_l, double u_r, double q_l, double q_r, double m, Flux flux_, const std::unique_ptr<Quadrature> &quad_, double u_mean, double u_mean_plus)
{
    double g_plus =0.0,g_minus =0.0;
    double u_diff = u_r-u_l;
    double q_diff = q_r-q_l;
    double q_mean = 0.5*(q_r+q_l);
    double fraction_term =0.0;
    double gamma = 0.0;
    double u_mean_diff = u_mean_plus-u_mean;
    double g_mean_plus = 0.0,g_mean_minus = 0.0;
    double g_diff=0.0;
    if((std::abs(u_diff))<1E-14 ){
        //std::cout<<"IM IN THE IF"<<std::endl;
        fraction_term = -1.0*flux_.compute(u_r,0,m)[1];
        // g_plus= quad_->G(u_r,m);//-1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_r);
        // g_minus= quad_->G(u_l,m);//-1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_l);
        g_plus= -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_r);
        g_minus= -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_l);


        //std::cout<<"IN IF fraction term "<<fraction_term<<" UDIFF "<< u_diff<<" G_plus "<<g_plus<<" u_R "<<u_r<<" G_minus "<<g_minus<<" U_l"<<u_l <<std::endl;
    }else{
        // g_plus= quad_->G(u_r,m);//-1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_r);
        // g_minus= quad_->G(u_l,m);//-1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_l);

        g_plus= -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_r);
        g_minus= -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_l);


        //std::cout<<"IN ELSE fraction term "<<fraction_term<<" UDIFF "<< u_diff<<" G_plus "<<g_plus<<" u_R "<<u_r<<" G_minus "<<g_minus<<" U_l"<<u_l<<" u_mean_plus "<< u_mean_plus<< " u_mean "<<u_mean <<std::endl;

        if(abs(u_r-u_mean_plus)<1E-14 ) {
            u_mean_plus = u_r;       
        }
        if(abs(u_l-u_mean)<1E-14 ) {
            u_mean = u_l;
        }
        g_diff = g_plus-g_minus;
        fraction_term = (g_diff)/u_diff;
        //gamma = 0.5*fraction_term;
    }
    if((std::abs(u_mean_diff))<1E-14){
        gamma = -0.5*flux_.compute(u_mean_plus,0,m)[1];
    }else{
        // g_mean_plus = quad_->G(u_mean_plus,m);//-1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_mean_plus);
        // g_mean_minus = quad_->G(u_mean,m);//-1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_mean);
 
        g_mean_plus = -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_mean_plus);
        g_mean_minus = -1.0*quad_->GaussLegendreQuad([&](double x) {return flux_.compute(x,0,m)[1];},0.0,u_mean);

        gamma = 0.5*(g_mean_plus-g_mean_minus)/(u_mean_plus-u_mean);
    }

    double g_mean = 0.5*(g_plus+g_minus);
    // std::cout<< " G MEAN "<<g_mean<<" GAMMA "<<gamma<<" sol0 "<<-fraction_term*q_mean-gamma*q_diff<<" sol1 "<<gamma*u_diff<<" u_meandiff "<<u_mean_diff<<std::endl;
    // std::cout<<" IN NUMERICAL FLUXMINUS"<<" u_l "<<u_l<<" u_r "<<u_r <<" udiff " << u_diff<<" g_mean "<<g_mean<<std::endl;
    // std::cout<<" IN NUMERICAL FLUXMINUS"<<" q_l "<<q_l<<" q_r "<<q_r <<" q_diff "<<q_diff<<" q_mean "<<q_mean<<" g_plus "<<g_plus<<" g_minus "<<g_minus<<" gamma "<<gamma<<std::endl;
    return {-fraction_term*q_mean-gamma*q_diff,-g_mean+gamma*u_diff};
}
