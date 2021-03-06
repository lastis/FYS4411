#ifndef _WAVE_FUNCTIONS_H_INCLUDED
#define _WAVE_FUNCTIONS_H_INCLUDED

#include <cmath>
#include <iostream>
#include "../CPhys/Physical.h"
#include "GtoHelium.h"

namespace wave_functions{
    // Constant read only values. 
    static double alpha;
    static double beta;
    static double nDimensions;
    static double h;
    static double hInv;
    static double h2Inv;
    static double charge;
    static int nParticles;
    static int nHalf;
    static double (*getWaveFuncVal)(double** r, double* rAbs);

    using namespace physical::unit;
    using namespace std;

    static double phi1sDAlpha(double r){
        return -r*exp(-alpha*r);
    }

    static double phi2sDAlpha(double r){
        return (r*r*0.25 - r)*exp(-alpha*r*0.5);
    }

    static double phi2pxDAlpha(double x, double r){
        return -x*r*0.5*exp(-alpha*r*0.5);
    }
    
    static double phi2pyDAlpha(double y, double r){
        return -y*r*0.5*exp(-alpha*r*0.5);
    }

    static double phi2pzDAlpha(double z, double r){
        return -z*r*0.5*exp(-alpha*r*0.5);
    }

    static double phiDAlpha(int j, double* r, double rAbs){
        switch (j) {
            case 0 :
                return phi1sDAlpha(rAbs);
            case 1 :
                return phi2sDAlpha(rAbs);
            case 2 :
                return phi2pxDAlpha(r[0],rAbs);
            case 3 : 
                return phi2pyDAlpha(r[1],rAbs);
            case 4 : 
                return phi2pzDAlpha(r[2],rAbs);
            default:
                std::cout << "Index out of bounds in phi()!!!" << std::endl;
                return 0;
        }
    }

    static double f(double r){
        return -0.5*alpha*r;
    }

    static double g(double r){
        return exp(f(r));
    }

    static double phi1s(double r){
        return g(r)*g(r);
    }

    static double phi1sD(double x, double r){
        return -x*alpha*g(r)*g(r)/r;
    }

    static double phi1sDD(double r){
        return alpha*alpha*g(r)*g(r) - 2*alpha*g(r)*g(r)/r;
    }

    static double phi2s(double r){
        return g(r)*(1+f(r));
    }

    static double phi2sD(double x, double r){
        return -x*alpha*g(r)*(1+0.5*f(r))/r;
    }

    static double phi2sDD(double r){
        return 0.75*alpha*alpha*g(r)*(1+f(r)/3) - 2*alpha*g(r)*(1+0.5*f(r))/r;
    }

    static double phi2p(int x, double* r, double rAbs){
        return r[x]*alpha*g(rAbs);
    }

    static double phi2pxD(int x, double* r, double rAbs){
        double tmp = alpha*g(rAbs);
        switch (x) {
            case 0:
                return tmp*(1-alpha*r[0]*r[x]*0.5/rAbs);
            case 1: 
                return -tmp*alpha*r[0]*r[x]*0.5/rAbs;
            case 2:
                return -tmp*alpha*r[0]*r[x]*0.5/rAbs;
            default:
                std::cout << "Out of bound in phi2px!!" << std::endl;
                return 0;
        }
    }

    static double phi2pyD(int x, double* r, double rAbs){
        double tmp = alpha*g(rAbs);
        switch (x) {
            case 0:
                return -tmp*alpha*r[1]*r[x]*0.5/rAbs;
            case 1: 
                return tmp*(1-alpha*r[1]*r[x]*0.5/rAbs);
            case 2:
                return -tmp*alpha*r[1]*r[x]*0.5/rAbs;
            default:
                std::cout << "Out of bound in phi2px!!" << std::endl;
                return 0;
        }
    }

    static double phi2pzD(int x, double* r, double rAbs){
        double tmp = alpha*g(rAbs);
        switch (x) {
            case 0:
                return -tmp*alpha*r[2]*r[x]*0.5/rAbs;
            case 1: 
                return -tmp*alpha*r[2]*r[x]*0.5/rAbs;
            case 2:
                return tmp*(1-alpha*r[2]*r[x]*0.5/rAbs);
            default:
                std::cout << "Out of bound in phi2py!!" << std::endl;
                return 0;
        }
    }

    static double phi2pDD(int x, double* r, double rAbs){
        return phi2p(x,r,rAbs)*alpha*(alpha*rAbs - 8)/(4*rAbs);
    }

    static double phi(int j, double* r, double rAbs){
        switch (j) {
            case 0 :
                return phi1s(rAbs);
            case 1 :
                return phi2s(rAbs);
            case 2 :
                return phi2p(0,r,rAbs);
            case 3 : 
                return phi2p(1,r,rAbs);
            case 4 : 
                return phi2p(2,r,rAbs);
            default:
                std::cout << "Index out of bounds in phi()!!!" << std::endl;
                return 0;
        }
    }

    static double phiD(int j, double* r, double rAbs, int x){
        switch (j) {
            case 0 :
                return phi1sD(r[x],rAbs); 
            case 1 :
                return phi2sD(r[x],rAbs);
            case 2 :
                return phi2pxD(x,r,rAbs); // Cartesian coordinates.
            case 3 : 
                return phi2pyD(x,r,rAbs); // Cartesian coordinates.
            case 4 : 
                return phi2pzD(x,r,rAbs); // Cartesian coordinates.
            default:
                std::cout << "Index out of bounds in phi()!!!" << std::endl;
                return 0;
        }
    }

    static double phiDD(int j, double* r, double rAbs){
        switch (j) {
            case 0 :
                return phi1sDD(rAbs); // Spherical coordinates.
            case 1 :
                return phi2sDD(rAbs);
            case 2 :
                return phi2pDD(0,r,rAbs); // Cartestian coordinates
            case 3 : 
                return phi2pDD(1,r,rAbs);
            case 4 : 
                return phi2pDD(2,r,rAbs);
            default:
                std::cout << "Index out of bounds in phi()!!!" << std::endl;
                return 0;
        }
    }

    static double getLocalEnergyGeneric(double** r, double* rAbs){
        double waveFunctionMinus = 0;
        double waveFunctionPlus = 0;
        double waveFunctionCurrent;
        waveFunctionCurrent = getWaveFuncVal(r, rAbs);

        // Kinetic energy

        double rPlus;
        double rMinus;
        double r0;
        double r0Abs;
        double kineticEnergy = 0;
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rPlus 	= r[i][j] + h;
                rMinus 	= r[i][j] - h;
                r0 		= r[i][j];
                r0Abs = rAbs[i];

                // Update values for rMinus
                r[i][j] = rMinus;
                rAbs[i] = 0;
                for (int x = 0; x < nDimensions; x++) {
                    rAbs[i] += r[i][x]*r[i][x];
                }
                rAbs[i] = sqrt(rAbs[i]);
                // Evaluate at rMins
                waveFunctionMinus = getWaveFuncVal(r,rAbs);
                // Update values for rPlus
                r[i][j] = rPlus;
                rAbs[i] = 0;
                for (int x = 0; x < nDimensions; x++) {
                    rAbs[i] += r[i][x]*r[i][x];
                }
                rAbs[i] = sqrt(rAbs[i]);
                // Evaluate at rPlus
                waveFunctionPlus = getWaveFuncVal(r,rAbs);
                // Reset.
                r[i][j] = r0;
                rAbs[i] = r0Abs;
                kineticEnergy -= (waveFunctionMinus + waveFunctionPlus 
                    - 2 * waveFunctionCurrent);
            }
        }
        kineticEnergy = 0.5 * h2Inv * kineticEnergy / waveFunctionCurrent;

        // Potential energy
        double potentialEnergy = 0;
        double rSingleParticle = 0;
        for(int i = 0; i < nParticles; i++) {
            rSingleParticle = 0;
            for(int j = 0; j < nDimensions; j++) {
                rSingleParticle += r[i][j]*r[i][j];
            }
            potentialEnergy -= charge / sqrt(rSingleParticle);
        }
        // Contribution from electron-electron potential
        double r12 = 0;
        for(int i = 0; i < nParticles; i++) {
            for(int j = i + 1; j < nParticles; j++) {
                r12 = 0;
                for(int k = 0; k < nDimensions; k++) {
                    r12 += (r[i][k] - r[j][k]) * (r[i][k] - r[j][k]);
                }
                potentialEnergy += 1 / sqrt(r12);
            }
        }
        return kineticEnergy + potentialEnergy;
    }

    static double getWaveFuncHeliumNoCor(double** r, double* rAbs){
        double argument = 0;
        using namespace std;
        for(int i = 0; i < nParticles; i++) {
            argument += rAbs[i];
        }
        return exp(-argument * alpha);
    }

    static double getWaveFuncHelium(double** r, double* rAbs){
        double* r1 = r[0];
        double* r2 = r[1];
        double temp = 0;
        double r1Abs = rAbs[0];
        double r2Abs = rAbs[1];
        double r12 = 0;
        for(int j = 0; j < nDimensions; j++) {
            temp = (r2[j] - r1[j]) * (r2[j] - r1[j]);
            r12 += temp;
        }
        r12 = sqrt(r12);
        return exp(-(r1Abs + r2Abs)*alpha)*exp(r12/(2*(1+beta*r12)));
    }

    static double getWaveBerylliumNoCor(double** r, double* rAbs){
        double r1Abs = rAbs[0];
        double r2Abs = rAbs[1];
        double r3Abs = rAbs[2];
        double r4Abs = rAbs[3];
        return (phi1s(r1Abs)*phi2s(r2Abs) - phi1s(r2Abs)*phi2s(r1Abs))
            *(phi1s(r3Abs)*phi2s(r4Abs) - phi1s(r4Abs)*phi2s(r3Abs));
    }

    static double getWaveBeryllium(double** r, double* rAbs){
        double r1Abs = rAbs[0];
        double r2Abs = rAbs[1];
        double r3Abs = rAbs[2];
        double r4Abs = rAbs[3];
        // The value of the slater determinant.
        double phi = (phi1s(r1Abs)*phi2s(r2Abs) -phi1s(r2Abs)*phi2s(r1Abs))
        *(phi1s(r3Abs)*phi2s(r4Abs) -phi1s(r4Abs)*phi2s(r3Abs));

        double cor = 0;
        double rij = 0;
        double a;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (j <= i) continue; 
                // If i < j, calculate something.
                rij = 0;
                for (int k = 0; k < nDimensions; k++) {
                    rij += (r[j][k] - r[i][k])*(r[j][k] - r[i][k]);
                }
                rij = sqrt(rij);
                int spinI = i/nHalf;
                int spinJ = j/nHalf;
                switch (spinI + spinJ){
                    case 0:
                        a = 0.25;
                        break;
                    case 1:
                        a = 0.5;
                        break;
                    case 2:
                        a = 0.25;
                        break;
                }
                double bij = 1/(1 + beta*rij);
                cor += rij*bij*a;
            }
        }
        cor = exp(cor);
        return phi*cor;
    }

    static double getLocalEnergyGenericNoCor(double** r, double* rAbs){
        double waveFunctionMinus = 0;
        double waveFunctionPlus = 0;
        double waveFunctionCurrent;

        waveFunctionCurrent = getWaveFuncVal(r,rAbs);

        // Kinetic energy

        double rPlus;
        double rMinus;
        double r0;
        double r0Abs;
        double kineticEnergy = 0;
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rPlus 	= r[i][j] + h;
                rMinus 	= r[i][j] - h;
                r0 		= r[i][j];
                r0Abs = rAbs[i];


                // Update values for rMinus
                r[i][j] = rMinus;
                rAbs[i] = 0;
                for (int x = 0; x < nDimensions; x++) {
                    rAbs[i] += r[i][x]*r[i][x];
                }
                rAbs[i] = sqrt(rAbs[i]);
                // Evaluate at rMins
                waveFunctionMinus = getWaveFuncVal(r,rAbs);
                // Update values for rPlus
                r[i][j] = rPlus;
                rAbs[i] = 0;
                for (int x = 0; x < nDimensions; x++) {
                    rAbs[i] += r[i][x]*r[i][x];
                }
                rAbs[i] = sqrt(rAbs[i]);
                // Evaluate at rPlus
                waveFunctionPlus = getWaveFuncVal(r,rAbs);
                // Reset
                r[i][j] = r0;
                rAbs[i] = r0Abs;
                kineticEnergy -= (waveFunctionMinus + waveFunctionPlus 
                    - 2 * waveFunctionCurrent);
            }
        }
        kineticEnergy = 0.5 * h2Inv * kineticEnergy / waveFunctionCurrent;

        // Potential energy
        double potentialEnergy = 0;
        double rSingleParticle = 0;
        for(int i = 0; i < nParticles; i++) {
            rSingleParticle = 0;
            for(int j = 0; j < nDimensions; j++) {
                rSingleParticle += r[i][j]*r[i][j];
            }
            potentialEnergy -= charge / sqrt(rSingleParticle);
        }
        return kineticEnergy + potentialEnergy;
    }

    static double getLocalEnergyHydrogen(double** r, double* rAbs){
        return  -1/rAbs[0] - 0.5*alpha*(alpha - 2/rAbs[0]);
    }

    static double getLocalEnergyHeliumNoCor(double** r, double* rAbs){
        double* r1 = r[0];
        double* r2 = r[1];
        double temp = 0;
        double r12Abs = 0;
        double r1Abs = 0;
        double r2Abs = 0;
        double r1r2 = 0; // Dot product.
        for(int j = 0; j < nDimensions; j++) {
            temp = r1[j] * r1[j];
            r1Abs += temp;
            temp = r2[j] * r2[j];
            r2Abs += temp;
            temp = (r1[j] - r2[j]) * (r1[j] - r2[j]);
            r12Abs += temp;
            // Dot product.
            temp = r1[j]*r2[j];
            r1r2 += temp;
        }
        r1Abs = sqrt(r1Abs);
        r2Abs = sqrt(r2Abs);
        r12Abs = sqrt(r12Abs);
        double E1 = (alpha-charge)*(1/r1Abs + 1/r2Abs) - alpha*alpha;
        return E1;
    }

    static double getLocalEnergyHelium(double** r, double* rAbs){
        double* r1 = r[0];
        double* r2 = r[1];
        double temp = 0;
        double r12Abs = 0;
        double r1Abs = 0;
        double r2Abs = 0;
        double r1r2 = 0; // Dot product.
        for(int j = 0; j < nDimensions; j++) {
            temp = r1[j] * r1[j];
            r1Abs += temp;
            temp = r2[j] * r2[j];
            r2Abs += temp;
            temp = (r1[j] - r2[j]) * (r1[j] - r2[j]);
            r12Abs += temp;
            // Dot product.
            temp = r1[j]*r2[j];
            r1r2 += temp;
        }
        r1Abs = sqrt(r1Abs);
        r2Abs = sqrt(r2Abs);
        r12Abs = sqrt(r12Abs);
        double E1 = (alpha-charge)*(1/r1Abs + 1/r2Abs) + 1/r12Abs - alpha*alpha;
        double betaR12 = 1/((1+beta*r12Abs));
        return E1 + betaR12*betaR12/2 * (
            alpha*(r1Abs + r2Abs)/r12Abs*(1-(r1r2/(r1Abs*r2Abs)))
            - betaR12*betaR12/2 - 2/r12Abs + 2*beta*betaR12
            );
    }
}
#endif
