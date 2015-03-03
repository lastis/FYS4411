#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include "../CPhys/CPhys.h"

class VMCSolver
{
    static const int LOCAL_ENERGY_GENERIC = 1;
    static const int LOCAL_ENERGY_HELIUM = 2;
    static const int LOCAL_ENERGY_HYDROGEN = 3;
    static const int WAVE_FUNCTION_1 = 1;
    static const int WAVE_FUNCTION_2 = 2;
    static const int WAVE_FUNCTION_BERYLLIUM = 3;
public:
    VMCSolver();

    bool runIntegration();
    bool initFromFile(std::string fName);
    void exportParamters(std::string fName);
    void exportDensity(std::string fName);
    void exportEnergyArray(std::string fName);
    void setWaveFunction1();
    void setWaveFunction2();
    void setWaveFunctionBeryllium();
    void setLocalEnergyHelium();
    void setLocalEnergyHydrogen();
    void setLocalEnergyGeneric();
    void setImportanceSampling(bool param);
    void setRecordDensity(bool param, int bins = 9, double maxPos = 2);
    void setRecordChargeDensity();
    void setRecordEnergyArray(bool param);
    double getAcceptanceRatio();
    void setStepLength(double stepLength);
    double getStepLength();
    double getR12Mean();
    double getEnergy();
    void clear();
    void supressOutput();


    // Parameters gathered from file. Be careful about changing these directly.
    double alpha;
    double beta;
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
    long idum;
    double timeStep;
    double D;

private:

    // Private functions
    double (VMCSolver::*getWaveFuncVal)(double** r);
    double getWaveFunc1Val(double** r);
    double getWaveFunc2Val(double** r);
    double getWaveBerylliumVal(double** r);
    double phi1s(double r);
    double phi2s(double r);
    double (VMCSolver::*getLocalEnergy)(double** r);
    double getLocalEnergyGeneric(double** r);
    double getLocalEnergyHelium(double** r);
    double getLocalEnergyHydrogen(double** r);
    bool   initRunVariables();
    void   updateQuantumForce(double** r, double ** qForce,double factor);
    void endOfSingleParticleStep(int cycle, int i);
    void endOfCycle(int cycle);
    void runRandomStep(int cycle);
    void runQuantumStep(int cycle);


    // Values from the simulation. 
    double mean;
    double energy;
    double energySquared;

    // Run variables
    double energySum;
    double energySquaredSum;
    double rAbsSum;
    double deltaE;
    double waveFuncValOld;
    double waveFuncValNew;
    int bins;
    double rMax;

    bool ready;
    bool outputSupressed;
    bool useImportanceSampling;
    bool recordDensity;
    bool recordChargeDensity;
    bool recordEnergyArray;

    Matrix qForceOld;
    Matrix qForceNew;
    double** pqForceOld;
    double** pqForceNew;
    Matrix rOld;
    Matrix rNew;
    double** prOld;
    double** prNew;

    Matrix density;
    Matrix densityCharge;
    double** pDensity;
    double** pDensityCharge;

    Vector energyArray;
    double* pEnergyArray;

};

#endif // VMCSOLVER_H
