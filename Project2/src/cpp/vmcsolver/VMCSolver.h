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
    static const int WAVE_FUNCTION_HELIUM_GTO = 12;

private:
    void endOfSingleParticleStep(int cycle, int i);
    void updateQuantumForce(double** r, double* rAbs, double** qForce,
                            double factor);
    void updateQuantumForceSlater(double** r, double* rAbs, double** qForce,
                                  double** pslater1, double** pslater2,
                                  double** pslater1Inv, double** pslater2Inv);
    void updateSlater(int i, double** slater1New, double** slater1Old,
                      double** slater2New, double** slater2Old,
                      double** slater1InvNew, double** slater1InvOld,
                      double** slater2InvNew, double** slater2InvOld);

public:
    VMCSolver();

    void clear();
    double getAcceptanceRatio();
    double getStepLength();
    double getR12Mean();
    double getEnergy();
    double getEnergySquared();
    void supressOutput();

    bool initRunVariables();

    void startOfCycle();
    void startOfCycleQuantum();
    void startOfCycleSlaterQuantum();

    void runSingleStep(int i, int cycle);
    void runSingleStepSlater(int i, int cycle);
    void runSingleStepQuantum(int i, int cycle);
    void runSingleStepSlaterQuantum(int i, int cycle);

    double calc_dE_dAlpha();
    double calc_dE_dBeta();

    double (*getWaveFuncVal)(double** r, double* rAbs);
    double (*getLocalEnergy)(double** r, double* rAbs);
    double getLocalEnergySlater(double** r, double* rAbs);
    double getLocalEnergySlaterNoCor(double** r, double* rAbs);
    double getCorrelationRatio(int i);

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

    std::mt19937 gen;
    std::uniform_real_distribution<double> dist_uniform;
    std::normal_distribution<double> dist_gauss;

    double deltaE;
    double waveFuncValOld;
    double waveFuncValNew;
    double greensFunction;
    double ratio;
    double potentialEnergy;
    double DD;
    double CC;
    double DC;
    double dE_dAlpha;
    double dE_dBeta;
    int nHalf;

    bool usingCorrelation;
    bool importanceSampling;
    bool efficientSlater;

    // Private variables
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
};

#endif  // VMCSOLVER_H
