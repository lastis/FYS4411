#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include "CPhys.h"

class VMCSolver
{
public:
    VMCSolver();

    void runMonteCarloIntegration();
    void useWaveType1();
    void useWaveType2();

private:
    int waveFunctionType;
    double waveFunction1(double** r);
    double waveFunction2(double* r1, double* r2);
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
    Matrix rPlus;
    Matrix rMinus;

    double** prOld;
    double** prNew;
    double** prPlus;
    double** prMinus;
};

#endif // VMCSOLVER_H
