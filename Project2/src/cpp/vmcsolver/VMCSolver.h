#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include "../CPhys/CPhys.h"
#include "SingleParticleWaveFunctions.h"

class VMCSolver
{
public:
    VMCSolver();
    static const int LOCAL_ENERGY_GENERIC = 1;
    static const int LOCAL_ENERGY_HELIUM_1 = 2;
    static const int LOCAL_ENERGY_HELIUM_2 = 4;
    static const int LOCAL_ENERGY_HYDROGEN = 3;
    static const int LOCAL_ENERGY_GENERIC_NOCOR = 5;
    static const int LOCAL_ENERGY_SLATER = 6;
    static const int LOCAL_ENERGY_SLATER_NOCOR = 7;
    static const int WAVE_FUNCTION_1 = 1;
    static const int WAVE_FUNCTION_2 = 2;
    static const int WAVE_FUNCTION_BERYLLIUM_1 = 3;
    static const int WAVE_FUNCTION_BERYLLIUM_2 = 4;

    bool runIntegration();
    bool initFromFile(std::string fName);
    void exportParamters(std::string fName);
    void exportDensity(std::string fName);
    void exportEnergyArray(std::string fName);
    void exportPositions(std::string fName);
    void setWaveFunction1();
    void setWaveFunction2();
    void setWaveFunctionBeryllium1();
    void setWaveFunctionBeryllium2();
    void setLocalEnergyHelium1();
    void setLocalEnergyHelium2();
    void setLocalEnergyHydrogen();
    void setLocalEnergyGeneric();
    void setLocalEnergyGenericNoCor();
    void setLocalEnergySlater();
    void setLocalEnergySlaterNoCor();
    void setRecordDensity(bool param, int bins = 9, double maxPos = 2);
    void setRecordEnergyArray(bool param);
    void setRecordR12Mean(bool param);
    void setRecordPositions(bool param);
    void setStepLength(double stepLength);
    void useImportanceSampling(bool param);
    void useEfficientSlater(bool param);
    double getAcceptanceRatio();
    double getStepLength();
    double getR12Mean();
    double getEnergy();
    double getEnergySquared();
    void clear();
    void supressOutput();
    bool validateParamters();

    double getLocalEnergyGeneric(double** r);
    double getLocalEnergyGenericNoCor(double** r);
    double getLocalEnergyHelium1(double** r);
    double getLocalEnergyHelium2(double** r);
    double getLocalEnergyHydrogen(double** r);
    double getLocalEnergySlater(double** r);
    double getLocalEnergySlaterNoCor(double** r);
    bool initRunVariables();
    void updateQuantumForce(double** r, double ** qForce,double factor);
    void updateSlater(int i);
    void updateInverse(int i, double ratio, double** mat, double** inv);
    double getWaveFunc1Val(double** r);
    double getWaveFunc2Val(double** r);
    double getWaveBeryllium1Val(double** r);
    double getWaveBeryllium2Val(double** r);


    // Parameters to be set manually or from file.
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
    double (VMCSolver::*getLocalEnergy)(double** r);
    void endOfSingleParticleStep(int cycle, int i);
    void endOfCycle(int cycle);
    void runRandomStep(int cycle);
    void runQuantumStep(int cycle);


public:
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
    bool recordDensity;
    bool recordChargeDensity;
    bool recordEnergyArray;
    bool recordR12Mean;
    bool recordPositions;
    bool importanceSampling;
    bool efficientSlater;

    Matrix slater1;
    Matrix slater1Inv;
    Matrix slater2;
    Matrix slater2Inv;
    double** pslater1;
    double** pslater1Inv;
    double** pslater2;
    double** pslater2Inv;
    Vector vS;
    double* S;

    Matrix qForceOld;
    Matrix qForceNew;
    double** pqForceOld;
    double** pqForceNew;
    Matrix rOld;
    Matrix rNew;
    double** prOld;
    double** prNew;

    Matrix positions;
    double** pPositions;

    Matrix density;
    double** pDensity;

    Vector energyArray;
    double* pEnergyArray;

};

#endif // VMCSOLVER_H
