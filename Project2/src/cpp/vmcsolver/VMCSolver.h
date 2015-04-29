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
    static const int LOCAL_ENERGY_HELIUM_1 = 2;
    static const int LOCAL_ENERGY_HELIUM_2 = 4;
    static const int LOCAL_ENERGY_HYDROGEN = 3;
    static const int LOCAL_ENERGY_GENERIC_NOCOR = 5;
    static const int LOCAL_ENERGY_SLATER = 6;
    static const int LOCAL_ENERGY_SLATER_NOCOR = 7;
    static const int WAVE_FUNCTION_1 = 8;
    static const int WAVE_FUNCTION_2 = 9;
    static const int WAVE_FUNCTION_BERYLLIUM_1 = 10;
    static const int WAVE_FUNCTION_BERYLLIUM_2 = 11;
private:
    bool initRunVariables();
    void endOfSingleParticleStep(int cycle, int i);
    void endOfCycle(int cycle);
    void updateQuantumForce(double** r, double ** qForce,double factor);
    void updateSlater(int i);
public:
    VMCSolver();

    bool runIntegration();
    void clear();
    double getAcceptanceRatio();
    double getStepLength();
    double getR12Mean();
    double getEnergy();
    double getEnergySquared();
    void supressOutput();
    void runRandomStep(int cycle);
    void runQuantumStep(int cycle);

    double (VMCSolver::*getWaveFuncVal)(double** r);
    double (VMCSolver::*getLocalEnergy)(double** r, int i);
    double getLocalEnergyGeneric(double** r, int i);
    double getLocalEnergyGenericNoCor(double** r, int i);
    double getLocalEnergyHelium1(double** r, int i);
    double getLocalEnergyHelium2(double** r, int i);
    double getLocalEnergyHydrogen(double** r, int i);
    double getLocalEnergySlater(double** r, int i);
    double getLocalEnergySlaterNoCor(double** r, int i);
    double getWaveFunc1Val(double** r);
    double getWaveFunc2Val(double** r);
    double getWaveBeryllium1Val(double** r);
    double getWaveBeryllium2Val(double** r);
    Vector getEnergyArray();

    void    setSeed(long seed);
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
    double ratio;
    int nHalf;

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
    Matrix positions;
    double** pPositions;
    Matrix density;
    double** pDensity;
    Vector energyArray;
    double* pEnergyArray;
};

#endif // VMCSOLVER_H

