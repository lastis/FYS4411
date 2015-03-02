#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include "CPhys.h"

class VMCSolver
{
    static const int LOCAL_ENERGY_GENERIC = 1;
    static const int LOCAL_ENERGY_HELIUM = 2;
    static const int LOCAL_ENERGY_HYDROGEN = 3;
    static const int WAVE_FUNCTION_1 = 1;
    static const int WAVE_FUNCTION_2 = 2;
public:
    VMCSolver();

    bool runIntegration();
    bool initFromFile(std::string fName);
    void exportParamters(std::string fName);
    void setWaveFunction1();
    void setWaveFunction2();
    void setLocalEnergyHelium();
    void setLocalEnergyHydrogen();
    void setLocalEnergyGeneric();
    void setImportanceSampling();
    double getAcceptanceRatio();
    void setStepLength(double stepLength);
    double getStepLength();
    double getR12Mean();
    double getEnergy();
    void reset();
    void clearAll();
    void supressOutput();


    // Parameters gathered from file TODO.
    double alpha;
    double beta;
    long idum;
    double timeStep;
    double D;

private:

    // Private functions
    double (VMCSolver::*getWaveFuncVal)(double** r);
    double getWaveFunc1Val(double** r);
    double getWaveFunc2Val(double** r);
    double (VMCSolver::*getLocalEnergy)(double** r);
    double getLocalEnergyGeneric(double** r);
    double getLocalEnergyHelium(double** r);
    double getLocalEnergyHydrogen(double** r);
    void   updateQuantumForce(double** r, double ** qForce,double factor);
    bool   initRunVariables();
    inline void runRandomWalk();
    inline void runQuantumWalk();

    bool ready;
    bool outputSupressed;
    bool useImportanceSampling;


    // Paramters are gathered from file. 
    int waveFunction;
    int localEnergyFunction;
    int accepts;
    int rejects;
    int charge;
    int nDimensions;
    int nCycles;
    int nParticles;
    double stepLength;
    double h;
    double h2;

    // Values from the simulation. 
    double mean;
    double energy;
    double energySquared;
    double energySum;
    double energySquaredSum;
    double rAbsSum;
    double deltaE;
    double waveFuncValOld;
    double waveFuncValNew;

    Matrix qForceOld;
    Matrix qForceNew;
    Matrix rOld;
    Matrix rNew;


    double** pqForceOld;
    double** pqForceNew;
    double** prOld;
    double** prNew;
};

#endif // VMCSOLVER_H
