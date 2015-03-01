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

    bool runIntegration();
    bool initFromFile(std::string fName);
    void exportParamters(std::string fName);
    void useWaveFunction1();
    void useWaveFunction2();
    void useLocalEnergyHelium();
    void useLocalEnergyHydrogen();
    void useLocalEnergyGeneric();
    double getStepAcceptance();
    void setStepLength(double stepLength);
    double getStepLength();
    double getR12Mean();
    double getEnergy();
    void reset();
    void clearAll();
    void supressOutput();


    // Parameters
    double alpha;
    double beta;
    long idum;

private:
    static const int LOCAL_ENERGY_GENERIC = 1;
    static const int LOCAL_ENERGY_HELIUM = 2;
    static const int LOCAL_ENERGY_HYDROGEN = 3;
    static const int WAVE_FUNCTION_1 = 1;
    static const int WAVE_FUNCTION_2 = 2;

    double wave1(double** r);
    double wave2(double* r1, double* r2);
    double localEnergyGeneric(double** r);
    double localEnergyHelium(double** r);
    double localEnergyHydrogen(double* r1);

    bool initialized;
    bool outputSupressed;
    // Paramters are gathered from file. 
    int waveFunction;
    int accepts;
    int rejects;
    int charge;
    int nDimensions;
    int nCycles;
    int nParticles;
    double stepLength;
    double h;
    double h2;
    int localEnergyFunction;

    // Values from the simulation. 
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
