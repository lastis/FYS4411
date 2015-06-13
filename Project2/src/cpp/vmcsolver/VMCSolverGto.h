#ifndef _VMCSOLVER_GTO_H_INCLUDED
#define _VMCSOLVER_GTO_H_INCLUDED

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include <random>
#include "../CPhys/CPhys.h"
#include "SingleParticleWaveFunctions.h"
#include "GtoHelium.h"
#include "GtoBeryllium.h"
#include "GtoNeon.h"

class VMCSolverGto
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
    static const int WAVE_FUNCTION_HELIUM_GTO = 12;
    static const int WAVE_FUNCTION_BERYLLIUM_GTO = 13;
    static const int WAVE_FUNCTION_NEON_GTO = 14;

private:
    double (*phi)(int i, int j, double** r, double* rAbs);
    double (*phiD)(int i, int j, double** r, double* rAbs, int x);
    double (*phiDD)(int i, int j, double** r, double* rAbs);
    void updateSlater(int i, double** slater1New, double** slater1Old,
                      double** slater2New, double** slater2Old,
                      double** slater1InvNew, double** slater1InvOld,
                      double** slater2InvNew, double** slater2InvOld);

public:
    VMCSolverGto();

    void clear();
    bool initRunVariables();
    void startOfCycle();
    void runSingleStep(int i);
    double getLocalEnergySlater(double** r, double* rAbs);
    double getCorrelationRatio(int i);
    void setSeed(long seed);

    // Initialization variables
    int waveFunction;
    double beta;
    int charge;
    int nDimensions;
    int nParticles;
    double stepLength;
    double h;
    double hInv;
    double h2Inv;
    long idum;

    // Other
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist_uniform;

    // State variables
    int accepts;
    int rejects;
    double deltaE;
    double ratio;
    double potentialEnergy;
    double DD;
    double CC;
    double DC;
    int nHalf;

    Matrix slater1Old;
    Matrix slater1New;
    Matrix slater1InvOld;
    Matrix slater1InvNew;
    Matrix slater2Old;
    Matrix slater2New;
    Matrix slater2InvOld;
    Matrix slater2InvNew;
    double** pslater1Old;
    double** pslater1New;
    double** pslater1InvOld;
    double** pslater1InvNew;
    double** pslater2Old;
    double** pslater2New;
    double** pslater2InvOld;
    double** pslater2InvNew;
    Matrix rOld;
    Matrix rNew;
    double** prOld;
    double** prNew;
    Vector rAbsOldVec;
    Vector rAbsNewVec;
    double* rAbsOld;
    double* rAbsNew;
};

#endif  // VMCSOLVER_H
