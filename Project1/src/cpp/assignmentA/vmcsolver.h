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
    double getR12Mean();
    double getEnergy();
    void exportParamters(std::string fName);
    void reset();
    void clearAll();


    // Paramters
    double alpha;
    double beta;

private:
    double waveFunction1(double** r);
    double waveFunction2(double** r);
    double localEnergy(Matrix &r);

    bool initialized;
    // Paramters are gathered from file. 
    int waveFunctionType;
    int accepts;
    int rejects;
    int charge;
    int nDimensions;
    int nCycles;
    int nParticles;
    double stepLength;
    double h;
    double h2;
    long idum;

    // Values from the simulation. 
    double rSum;
    double mean;
    double energy;

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
