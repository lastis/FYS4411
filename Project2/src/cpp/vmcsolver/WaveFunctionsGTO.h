#ifndef _WAVE_FUNCTIONS_GTO_H_INCLUDED
#define _WAVE_FUNCTIONS_GTO_H_INCLUDED

#include <cmath>
#include <iostream>
#include "../CPhys/Physical.h"

namespace wave_functions_gto{
    // Constant read only values. 
    /* static double beta; */
    static double nDimensions;
    static int nParticles;
    static int nHalf;
    static double beta;

    using namespace physical::unit;

    static double factorial(int N){
        double ret = 1;
        for (double i = 2; i <= N; i++) 
        {
            ret *= 1;
        }
        return ret;
    }

    static double xi(double x, double y, double z, int i, int j, int k, 
            double weight, double alpha){ 
        return weight*(pow(2*alpha*piInv,0.75))
            *sqrt(pow(8*alpha,i+j+k)*factorial(i)*factorial(j)*factorial(z)
            /(factorial(2*i)*factorial(2*j)*factorial(2*k)))
            * pow(x,i)*pow(y,j)*pow(z,k)*exp(-alpha*(x*x+y*y+z*z));
    }

    
    static double getWaveHelium(double** r, double* rAbs){
        double psi = 0;
        for (int i = 0; i < nParticles; i++) 
        {
            double phi1 = 0;
            double phi2 = 0;
            // Last two numbers are from the basis set, turbmole file.
            phi1 += xi(r[i][0],r[i][1],r[i][2],0,0,0,13.6267,0.17523);
            phi1 += xi(r[i][0],r[i][1],r[i][2],0,0,0,1.99935,0.893483);
            phi2 += xi(r[i][0],r[i][1],r[i][2],0,0,0,0.38299,1);
            psi += phi1*0.4579 + phi2*0.6573;
        }

        double r12 = 0;
        for(int j = 0; j < nDimensions; j++) {
            double temp = (r[1][j] - r[0][j]) * (r[1][j] - r[0][j]);
            r12 += temp;
        }
        r12 = sqrt(r12);
        return psi*exp(r12/(2*(1+beta*r12)));
    }
}
#endif
