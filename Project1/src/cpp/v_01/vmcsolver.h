#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include "CPhys.h"

class VMCSolver
{
public:
    VMCSolver();

    void runMonteCarloIntegration();

private:
    double waveFunction1(Matrix &r);
    double waveFunction2(Matrix &r);
    double localEnergy(Matrix &r);

    int nDimensions;
    int charge;
    double stepLength;
    int nParticles;

    double h;
    double h2;

    long idum;

    double alpha;
    double beta;

    int nCycles;

    Matrix rOld;
    Matrix rNew;

    double** prOld;
    double** prNew;
};

#endif // VMCSOLVER_H
