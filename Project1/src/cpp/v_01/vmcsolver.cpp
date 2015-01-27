#include "vmcsolver.h"
#include <iostream>
#include <stdio.h>
#include <math.h>

using namespace CPhys;
using namespace std;

VMCSolver::VMCSolver() :
    nDimensions(3),
    charge(2),
    stepLength(1.0),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(0),
    alpha(0.5*charge),
    beta(0.5*charge),
    nCycles(1000000)
{
}

void VMCSolver::runMonteCarloIntegration()
{
    rOld = Matrix(nParticles, nDimensions);
    rNew = Matrix(nParticles, nDimensions);

    prOld = rOld.getArrayPointer();
    prNew = rOld.getArrayPointer();

    double waveFunctionOld = 0;
    double waveFunctionNew = 0;

    double energySum = 0;
    double energySquaredSum = 0;

    double deltaE;

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            /* prOld[i][j] = stepLength * (Random::ran2(idum) - 0.5); */
            rOld(i,j) = stepLength * (Random::ran2(idum) - 0.5);
        }
    }
    rNew = rOld;

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++) {

        // Store the current value of the wave function
        waveFunctionOld = waveFunction1(rOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                /* prNew[i][j] = prOld[i][j] + stepLength*(Random::ran2(idum) - 0.5); */
                rNew(i,j) = rOld(i,j) + stepLength*(Random::ran2(idum) - 0.5);
            }

            // Recalculate the value of the wave function
            waveFunctionNew = waveFunction1(rNew);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(Random::ran2(idum) <= (waveFunctionNew*waveFunctionNew) /
			    (waveFunctionOld*waveFunctionOld)) {
                for(int j = 0; j < nDimensions; j++) {
                    /* prOld[i][j] = prNew[i][j]; */
                    rOld(i,j) = rNew(i,j);
                    waveFunctionOld = waveFunctionNew;
                }
            } else {
                for(int j = 0; j < nDimensions; j++) {
                    /* prNew[i][j] = prOld[i][j]; */
                    rNew(i,j) = rOld(i,j);
                }
            }
            // update energies
            deltaE = localEnergy(rNew);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
    }
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
}

double VMCSolver::localEnergy(Matrix &r)
{
    Matrix rPlus = Matrix(nParticles, nDimensions);
    Matrix rMinus = Matrix(nParticles, nDimensions);

    rPlus = rMinus = r;
    double** prPlus = rPlus.getArrayPointer();
    double** prMinus = rMinus.getArrayPointer();
    double** pr = r.getArrayPointer();

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = waveFunction1(r);

    // Kinetic energy

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            /* prPlus[i][j] += h; */
            /* prMinus[i][j] -= h; */
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = waveFunction1(rMinus);
            waveFunctionPlus = waveFunction1(rPlus);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus 
			    - 2 * waveFunctionCurrent);
            /* prPlus[i][j] = pr[i][j]; */
            /* prMinus[i][j] = pr[i][j]; */
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    // Potential energy
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            /* rSingleParticle += pr[i][j]*pr[i][j]; */
            rSingleParticle += r(i,j)*r(i,j);
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                /* r12 += (pr[i][k] - pr[j][k]) * (pr[i][k] - pr[j][k]); */
                r12 += (r(i,j) - r(j,k)) * (r(i,k) - r(j,k));
            }
            potentialEnergy += 1 / sqrt(r12);
        }
    }

    return kineticEnergy + potentialEnergy;
}

double VMCSolver::waveFunction1(Matrix &r)
{
    double argument = 0;
    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }
    return exp(-argument * alpha);
}

double VMCSolver::waveFunction2(Matrix &r)
{
    double argument = 0;
    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }
    return exp(-argument * alpha);
}
