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
    bool initFromFile(std::string fName);
    double getStepAcceptance();
    void setStepLength(double stepLength);
    double getStepLength();
    void exportParamters(std::string fName);

private:
    double waveFunction1(double** r);
    double waveFunction2(double** r);
    void reset();
    double localEnergy(Matrix &r);

    bool initialized;
    // Next paramters are gathered from file. 
    int waveFunctionType;
    int accepts;
    int rejects;
    int charge;
    int nDimensions;
    int nCycles;
    int nParticles;
    double stepLength;
    double alpha;
    double beta;
    double h;
    double h2;
    long idum;

    double rSum;
    double mean;



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
