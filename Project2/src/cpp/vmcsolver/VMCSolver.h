#ifndef _VMCSOLVER_H_INCLUDED
#define _VMCSOLVER_H_INCLUDED

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include <random>
#include "../CPhys/CPhys.h"
#include "SingleParticleWaveFunctions.h"

class VMCSolver
{
public:
    static const int LOCAL_ENERGY_GENERIC = 1;
    static const int LOCAL_ENERGY_GENERIC_NOCOR = 5;
    static const int LOCAL_ENERGY_HELIUM_1 = 2;
    static const int LOCAL_ENERGY_HELIUM_2 = 4;
    static const int LOCAL_ENERGY_HYDROGEN = 3;
    static const int LOCAL_ENERGY_SLATER = 6;
    static const int LOCAL_ENERGY_SLATER_NOCOR = 7;
    static const int WAVE_FUNCTION_1 = 8;
    static const int WAVE_FUNCTION_2 = 9;
    static const int WAVE_FUNCTION_BERYLLIUM_1 = 10;
    static const int WAVE_FUNCTION_BERYLLIUM_2 = 11;

private:
    void endOfSingleParticleStep(int cycle, int i);
    void endOfCycle(int cycle);
    void updateQuantumForce(double** r, double* rAbs, double** qForce,
                            double factor);
    void updateQuantumForceSlater(double** r, double* rAbs, double** qForce);
    void updateSlater(int i);

public:
    VMCSolver();

    void clear();
    double getAcceptanceRatio();
    double getStepLength();
    double getR12Mean();
    double getEnergy();
    double getEnergySquared();
    void supressOutput();

    bool runIntegration();
    bool initRunVariables();

    void startOfCycle();
    void startOfCycleQuantum();
    void startOfCycleSlaterQuantum();

    void runStep(int cycle);
    void runStepSlater(int cycle);
    void runStepQuantum(int cycle);
    void runStepSlaterQuantum(int cycle);

    void runSingleStep(int i, int cycle);
    void runSingleStepSlater(int i, int cycle);
    void runSingleStepQuantum(int i, int cycle);
    void runSingleStepSlaterQuantum(int i, int cycle);

    double (*getWaveFuncVal)(double** r, double* rAbs);
    double (*getLocalEnergy)(double** r, double* rAbs);
    double getLocalEnergySlater(double** r, double* rAbs);
    double getLocalEnergySlaterNoCor(double** r, double* rAbs);
    double getCorrelationRatio(int i);

    Vector getEnergyArray();

    void setSeed(long seed);
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
    double hInv;
    double h2Inv;
    long idum;
    double timeStep;
    double D;
    double rMax;
    int bins;

    std::mt19937 gen;
    std::uniform_real_distribution<double> dist_uniform;
    std::normal_distribution<double> dist_gauss;

    double mean;
    double energy;
    double energySquared;
    double energySum;
    double energySquaredSum;
    double rAbsSum;
    double deltaE;
    double waveFuncValOld;
    double waveFuncValNew;
    double greensFunction;
    double ratio;
    double potentialEnergy;
    double DD;
    double CC;
    double DC;
    int nHalf;

    bool usingCorrelation;
    bool outputSupressed;
    bool recordingDensity;
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
    Vector rAbsOldVec;
    Vector rAbsNewVec;
    double* rAbsOld;
    double* rAbsNew;
    Matrix positions;
    double** pPositions;
    Matrix density;
    double** pDensity;
    Vector energyArray;
    double* pEnergyArray;

};

#endif  // VMCSOLVER_H
