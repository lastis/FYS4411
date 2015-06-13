#ifndef _GTO_FUNCTIONS_H_INCLUDED
#define _GTO_FUNCTIONS_H_INCLUDED

#include <cmath>
#include <iostream>
#include "../CPhys/Physical.h"

namespace gto_general_functions{

    using namespace physical::unit;

    inline double factorial(int N){
        double ret = 1;
        for (double i = 2; i <= N; i += 1) 
        {
            ret *= i;
        }
        return ret;
    }

    inline double xi(double x, double y, double z, int i, int j, int k, 
            double alpha, double weight){ 
        return weight*pow(2*alpha*piInv,0.75)
            *sqrt(pow(8*alpha,i+j+k)*factorial(i)*factorial(j)*factorial(k)
            /(factorial(2*i)*factorial(2*j)*factorial(2*k)))
            * pow(x,i)*pow(y,j)*pow(z,k)*exp(-alpha*(x*x+y*y+z*z));
    }

    /** \brief 
     *
     * \param q Value of x, y, z of which to differentiate. 
     * \param m Value of i, j, k of which to differentiate. 
     */         
    inline double xiD(double q, int m, double x, double y, double z, 
            int i, int j, int k, double alpha, double weight){ 
        return weight*pow(2*alpha*piInv,0.75)
            *sqrt(pow(8*alpha,i+j+k)*factorial(i)*factorial(j)*factorial(k)
            /(factorial(2*i)*factorial(2*j)*factorial(2*k)))
            * pow(x,i)*pow(y,j)*pow(z,k)*exp(-alpha*(x*x+y*y+z*z))
            * (m/q - 2*alpha*q);
    }

    inline double xiDD(double x, double y, double z, 
            int i, int j, int k, double alpha, double weight){ 
        return weight*pow(2*alpha*piInv,0.75)
            *sqrt(pow(8*alpha,i+j+k)*factorial(i)*factorial(j)*factorial(k)
            /(factorial(2*i)*factorial(2*j)*factorial(2*k)))
            * pow(x,i-2)*pow(y,j-2)*pow(z,k-2)*exp(-alpha*(x*x+y*y+z*z))
            * (z*z*(y*y*(i*i-i*(4*alpha*x*x+1)+2*alpha*x*x*(2*alpha
                                * (x*x + y*y + z*z) - 3)) + 
            j*j*x*x - j*x*x *(4*alpha*y*y + 1)) + k*k*x*x*y*y - 
            k*x*x*y*y * (4*alpha*z*z + 1));
    }
}
#endif
