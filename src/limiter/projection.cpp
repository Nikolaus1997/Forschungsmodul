#include "projection.h"

Projection::Projection()
{

};

void Projection::project(Array2D &u, const Array2D &x)
{

}

std::pair<Array2D,Array2D> Projection::LR(Array2D &A)
{
    Array2D A_ = Array2D(A.size());
    Array2D L_ = Array2D(A.size());
    Array2D R_ = Array2D(A.size());
    for(int i=0;i<A_.size()[0];i++){
        for(int j=0;j<A_.size()[1];j++){
            A_(i,j) = A(i,j);
        }
    }
    for(int i=0; i<A_.size()[0]-1;i++){
        for(int j=i+1; j<A_.size()[1];j++){
            A_(j,i) = (A_(j,i)/A_(i,i));
            for(int k=i+1;k<A_.size()[0];k++){
                    A_(j,k) = A_(j,k) - A_(i,k)*A_(j,i);
            } 
        }
    }
    for(int i=0; i<A_.size()[0];i++){
        for(int j=0;j<A_.size()[1];j++){
            if(i==j){
                L_(i,j) = 1;
            }else if(i>j){
                L_(i,j) = A_(i,j);
            }else{
                L_(i,j) = 0;
            }
        }
    }

    for(int i=0; i<A_.size()[0];i++){
        for(int j=0;j<A_.size()[1];j++){
            if(i==j){
                R_(i,j) = A_(i,j);
            }else if(i<j){
                R_(i,j) = A_(i,j);
            }else{
                R_(i,j) = 0;
            }
        }
    }
    return std::make_pair(L_,R_);
}

// void Projection::AssA_tilde(Array2D &A)
// {
//     for(int k=0;k<A_tilde.size()[0];k++){
//         for(int l = 0; l<A_tilde.size()[1];l++){
//             for(int i=0;i<A.size()[0];i++){
//                     A_tilde(k,l) += A(i,k)*A(i,l);
//             }
//         }
//     }
// }

void Projection::makeProjection(Array2D &u, const Array2D x, int i, int order)
{
    Array2D f_ = Array2D({u.size()[1],1});

    for(int j=0;j<u.size()[1];j++){
        f_(j,0) = u(i,j);
    }
    Array2D A = MakeMonomBasis(x,i,order);
    Array2D  A_T = MakeTransPosed(A);
    Array2D A_TA = MatMul(A_T,A);
    Array2D A_Tf= MatMul(A_T,f_);
    std::pair<Array2D,Array2D> result = LR(A_TA);
    Array2D y_ = forwardSubstitution(result.first,A_Tf);
    Array2D coeff_ = backwardSubstitution(result.second,y_);
    for(int j=0;j<u.size()[1];j++){
        u(i,j) = 0.0;
        for(int k = 0;k<order;k++){
        u(i,j) += pow(x(i,j),double(k))*coeff_(k,0);
        //std::cout<<" ij"<<i<<" "<<j<<" U "<<u(i,j)<<" pow "<<pow(x(i,j),k)<<" coeff "<< coeff_(k,0)<<std::endl;
        }
    }
    // u.printValues();
    // std::cout<<"U"<<std::endl;

}

Array2D Projection::MatMul(const Array2D &A, const Array2D &B) {
    // Ensure dimensions are valid for multiplication
    if (A.size()[1] != B.size()[0]) {
        throw std::invalid_argument("Matrix dimensions are incompatible for multiplication");
    }

    // Determine output matrix dimensions (A.rows x B.cols)
    int rows = A.size()[0]; // Number of rows in A
    int cols = B.size()[1]; // Number of columns in B
    int common = A.size()[1]; // Shared dimension (A.cols == B.rows)

    Array2D solution({rows, cols}); // Initialize result matrix with correct dimensions

    // Matrix multiplication
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            solution(i, j) = 0; // Ensure initialization
            for (int k = 0; k < common; k++) {
                solution(i, j) += A(i, k) * B(k, j);
            }
        }
    }

    return solution;
}


Array2D Projection::MakeTransPosed(const Array2D &A)
{   
    Array2D solution = Array2D({A.size()[1],A.size()[0]});
    for(int i=0;i<A.size()[0];i++){
        for(int j=0;j<A.size()[1];j++){
            solution(j,i) = A(i,j);
        }
    }
    return solution;
}

Array2D Projection::MakeMonomBasis(const Array2D &u, int i, int order)
{
    Array2D solution = Array2D({u.size()[1],order});

    for(int j=0;j<u.size()[1];j++){
        for(int k=0;k<order;k++){
            solution(j,k) = pow(u(i,j),double(k));
        }
    }

    return solution;
}

Array2D Projection::backwardSubstitution(const Array2D &A, const Array2D &b)
{
    Array2D y = Array2D({A.size()[0],1});
    for(int i=A.size()[0]-1;i>=0;i--){
        y(i,0) = b(i,0);
        for(int j=i+1;j<A.size()[0];j++){
            y(i,0) -= A(i,j)*y(j,0);
        }
        y(i,0) = y(i,0)/A(i,i);
    }
    return y;
}

Array2D Projection::forwardSubstitution(const Array2D &A, const Array2D &b)
{
    Array2D solution = Array2D({A.size()[0],1});

    for(int i=0;i<A.size()[0];i++){
        solution(i,0) = b(i,0);
        for(int j=0;j<i;j++){
            solution(i,0) -= A(i,j)*solution(j,0);
        }
        solution(i,0) = solution(i,0)/A(i,i);
    }

    return solution;
}
