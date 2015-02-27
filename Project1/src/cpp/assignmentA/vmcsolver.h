#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include "CPhys.h"

class VMCSolver
{
public:
    VMCSolver();

    bool runMonteCarloIntegration();
    void useWaveType1();
    void useWaveType2();
    bool initFromFile(std::string fName = "main.ini");
    double getStepAcceptance();
    void setStepLength(double stepLength);
    void exportParamters();

private:
    double waveFunction1(double** r);
    double waveFunction2(double* r1, double* r2);
    void reset();
    double localEnergy(Matrix &r);

    bool initialized;
    int waveFunctionType;
    int accepts;
    int rejects;
    int nDimensions;
    int charge;
    int nParticles;
    double stepLength;

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
