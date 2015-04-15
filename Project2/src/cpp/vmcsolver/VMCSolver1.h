#ifndef _VMCSOLVER_H_INCLUDED
#define _VMCSOLVER_H_INCLUDED

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include "../CPhys/CPhys.h"
#include "SingleParticleWaveFunctions.h"

class VMCSolver1
{
public:
    VMCSolver1();
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
    bool initRunVariables();
    double getAcceptanceRatio();
    double getStepLength();
    double getR12Mean();
    double getEnergy();
    double getEnergySquared();
    void clear();
    void supressOutput();

    double getLocalEnergyGeneric(double** r);
    double getLocalEnergyGenericNoCor(double** r);
    double getLocalEnergyHelium1(double** r);
    double getLocalEnergyHelium2(double** r);
    double getLocalEnergyHydrogen(double** r);
    double getLocalEnergySlater(double** r);
    double getLocalEnergySlaterNoCor(double** r);
    void updateQuantumForce(double** r, double ** qForce,double factor);
    void updateSlater(int i);
    void updateInverse(int i, double ratio, double** mat, double** inv);
    double getWaveFunc1Val(double** r);
    double getWaveFunc2Val(double** r);
    double getWaveBeryllium1Val(double** r);
    double getWaveBeryllium2Val(double** r);


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
    double (VMCSolver1::*getWaveFuncVal)(double** r);
    double (VMCSolver1::*getLocalEnergy)(double** r);
    void endOfSingleParticleStep(int cycle, int i);
    void endOfCycle(int cycle);
    void runRandomStep(int cycle);
    void runQuantumStep(int cycle);


public:

    double mean;
    double energy;
    double energySquared;
    double energySum;
    double energySquaredSum;
    double rAbsSum;
    double deltaE;
    double waveFuncValOld;
    double waveFuncValNew;
    double rMax;
    int bins;

    bool ready;
    bool outputSupressed;
    bool recordingDensity;
    bool recordingChargeDensity;
    bool recordingEnergyArray;
    bool recordingR12Mean;
    bool recordingPositions;
    bool importanceSampling;
    bool efficientSlater;
    bool parallel;

    // Private variables
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
