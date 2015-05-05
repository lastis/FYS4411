#ifndef _WAVE_FUNCTIONS_H_INCLUDED
#define _WAVE_FUNCTIONS_H_INCLUDED

#include <cmath>
#include <iostream>

namespace wave_functions{
    // Constant read only values. 
    static double alpha;
    static double beta;
    static double nDimensions;
    static double h;
    static double h2;
    static double charge;
    static int nParticles;
    static int nHalf;
    static double (*getWaveFuncVal)(double** r, double* rAbs);


    static double f(double r){
        return -0.5*alpha*r;
    }

    static double g(double r){
        return exp(f(r));
    }

    static double phi1s(double r){
        return g(r)*g(r);
    }

    static double phi1sD(double r){
        return -alpha*g(r)*g(r);
    }

    static double phi1sDD(double r){
        return alpha*alpha*g(r)*g(r);
    }

    static double phi2s(double r){
        return g(r)*(1+f(r));
    }

    static double phi2sD(double r){
        return -alpha*g(r)*(1+0.5*f(r));
    }

    static double phi2sDD(double r){
        return 0.75*alpha*alpha*g(r)*(1+f(r)/3);
    }

    static double phi2p(double r){
        return alpha*r*exp(-alpha*r/2);
    }

    static double phi(int j, double* r){
        double rAbs = 0;
        for (int i = 0; i < nDimensions; i++) {
            rAbs += r[i]*r[i];
        }
        using namespace std;
        rAbs = sqrt(rAbs);
        switch (j) {
            case 0 :
                return phi1s(rAbs);
            case 1 :
                return phi2s(rAbs);
            case 2 ... 4:
                return phi2p(rAbs);
            default:
                std::cout << "Index out of bounds in phi()!!!" << std::endl;
                return 0;
        }
    }

    static double phiD(int j, double* r){
        double rAbs = 0;
        for (int i = 0; i < nDimensions; i++) {
            rAbs += r[i]*r[i];
        }
        rAbs = sqrt(rAbs);
        switch (j) {
            case 0 :
                return phi1sD(rAbs);
            case 1 :
                return phi2sD(rAbs);
            default:
                std::cout << "Index out of bounds in phi()!!!" << std::endl;
                return 0;
        }
    }

    static double phiDD(int j, double* r){
        double rAbs = 0;
        for (int i = 0; i < nDimensions; i++) {
            rAbs += r[i]*r[i];
        }
        rAbs = sqrt(rAbs);
        switch (j) {
            case 0 :
                return phi1sDD(rAbs) + 2*phi1sD(rAbs)/rAbs;
            case 1 :
                return phi2sDD(rAbs) + 2*phi2sD(rAbs)/rAbs;
            default:
                std::cout << "Index out of bounds in phi()!!!" << std::endl;
                return 0;
        }
    }

    static double berylliumPsiDD(double** r, double* rAbs)
    {
        /* double rAbs[4]; */
        /* for (int i = 0; i < 4; i++) { */
        /*     rAbs[i] = 0; */
        /*     for(int j = 0; j < nDimensions; j++) { */
        /*         rAbs[i] += r[i][j] * r[i][j]; */
        /*     } */
        /*     rAbs[i] = sqrt(rAbs[i]); */
        /* } */
        double tmp = 0;
        double sum = 0;
        tmp += phi2s(rAbs[1])*phi1sDD(rAbs[0]) - phi1s(rAbs[1])*phi2sDD(rAbs[0]);
        tmp -= phi2s(rAbs[0])*phi1sDD(rAbs[1]) - phi1s(rAbs[0])*phi2sDD(rAbs[1]);
        tmp += 2*(phi2s(rAbs[1])*phi1sD(rAbs[0]) 
                - phi1s(rAbs[1])*phi2sD(rAbs[0]))/rAbs[0];
        tmp -= 2*(phi2s(rAbs[0])*phi1sD(rAbs[1]) 
                - phi1s(rAbs[0])*phi2sD(rAbs[1]))/rAbs[1];
        sum += tmp/(phi1s(rAbs[0]) * phi2s(rAbs[1]) - phi2s(rAbs[0]) * phi1s(rAbs[1]));

        tmp = 0;
        tmp += phi2s(rAbs[3])*phi1sDD(rAbs[2]) - phi1s(rAbs[3])*phi2sDD(rAbs[2]);
        tmp -= phi2s(rAbs[2])*phi1sDD(rAbs[3]) - phi1s(rAbs[2])*phi2sDD(rAbs[3]);
        tmp += 2*(phi2s(rAbs[3])*phi1sD(rAbs[2]) 
                - phi1s(rAbs[3])*phi2sD(rAbs[2]))/rAbs[2];
        tmp -= 2*(phi2s(rAbs[2])*phi1sD(rAbs[3]) 
                - phi1s(rAbs[2])*phi2sD(rAbs[3]))/rAbs[3];
        sum += tmp/(phi1s(rAbs[2]) * phi2s(rAbs[3]) - phi2s(rAbs[2]) * phi1s(rAbs[3]));
        return sum;
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
        kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

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
            /* temp = r1[j] * r1[j]; */
            /* r1Abs += temp; */
            /* temp = r2[j] * r2[j]; */
            /* r2Abs += temp; */

            temp = (r2[j] - r1[j]) * (r2[j] - r1[j]);
            r12 += temp;
        }
        /* r1Abs = sqrt(r1Abs); */
        /* r2Abs = sqrt(r2Abs); */
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
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i >= j) continue; 
            // If i < j, calculate something.
            rij = 0;
            for (int k = 0; k < nDimensions; k++) {
                rij += (r[j][k] - r[i][k])*(r[j][k] - r[i][k]);
            }
            rij = sqrt(rij);
            if ((i == 0 || j == 0) && (i == 1 || j == 1))
                cor += 0.25*rij/(1+beta*rij);
            else if ((i == 2 || j == 2) && (i == 3 || j == 3))
                cor += 0.25*rij/(1+beta*rij);
            else
                cor += 0.5*rij/(1+beta*rij);
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
        kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

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
