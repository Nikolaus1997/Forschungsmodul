import numpy as np
import matplotlib.pyplot as plt


def LegendrePoly(degree,x):
    if degree ==0:
        return 1.0
    elif degree ==1:
        return x
    else:
        L = (2*degree+1.0)/(degree+1.0)*x*LegendrePoly(degree-1,x)-degree/(degree+1.0)*LegendrePoly(degree-2,x)

    return L


def LegendrePolynomialAndDerivative(N,x):
    if N==0:
        L = 1
        L_prime = 0
        return L,L_prime
    elif N==1:
        L = x
        L_prime= 1
        return L,L_prime
    else:
        L = (2*N-1.0)/N*x*LegendrePolynomialAndDerivative(N-1,x)[0]-(N-1)/N*LegendrePolynomialAndDerivative(N-2,x)[0]
        L_prime = LegendrePolynomialAndDerivative(N-2,x)[1]+(2*N-1)*LegendrePolynomialAndDerivative(N-1,x)[0]
    return L,L_prime

def LegendreGaussNodesAndWeights(N):
    w = np.zeros(N+1)
    x_j= np.copy(w)
    iterations = 4
    if N==0:
        x_j[0] = 0.0
        w[0] = 2.0
    elif N==1:
        x_j[0] = -np.sqrt(1/3)
        x_j[1] = -x_j[0]
        w[0] = 1
        w[1] = w[0] 
    else:
        for j in range(int(np.floor((N+1)/2))):
            x_j[j] = -np.cos((2*j+1)/(2*N+2)*np.pi)
            for k in range(iterations):
                L,L_prime = LegendrePolynomialAndDerivative(N+1,x_j[j])
                dx = -L/L_prime
                x_j[j] = x_j[j]+dx
                k = k+1
                if np.abs(dx) <=1e-8*np.abs(x_j[j]):
                    k = iterations
            L,L_prime = LegendrePolynomialAndDerivative(N+1,x_j[j])
            x_j[N-j] = -x_j[j]
            w[j] = 2.0/((1-x_j[j]**2)*L_prime**2)
            w[N-j] = w[j]
    if (N)%2==0:
        L,L_prime = LegendrePolynomialAndDerivative(N+1,0.0)
        x_j[int(N/2)] = 0.0
        w[int(N/2)]   = 2.0/(L_prime)**2
    return x_j,w

print(LegendreGaussNodesAndWeights(8))
