#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include "CPhys.h"

class VMCSolver
{
public:
    VMCSolver();

    void runMonteCarloIntegration();

private:
    double waveFunction(Matrix &r);
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

    double** prOld;
    double** prNew;
};

#endif // VMCSOLVER_H
