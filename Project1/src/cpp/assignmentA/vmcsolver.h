#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include "CPhys.h"

class VMCSolver
{
public:
    VMCSolver();

    void runMonteCarloIntegration();

private:
    double waveFunction(double** r);
    double localEnergy(Matrix &r);

    int nDimensions;
    int charge;
    double stepLength;
    int nParticles;

    double h;
    double h2;

    long idum;

    double alpha;

    int nCycles;

    Matrix rOld;
    Matrix rNew;
    Matrix rPlus;
    Matrix rMinus;

    double** prOld;
    double** prNew;
    double** prPlus;
    double** prMinus;
};

#endif // VMCSOLVER_H
