#ifndef _VMCWRAPPER_H_INCLUDED_
#define _VMCWRAPPER_H_INCLUDED_

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include "VMCSolver.h"
#include "../CPhys/CPhys.h"

class VMCWrapper
{
public:
    VMCWrapper();
    static const int LOCAL_ENERGY_GENERIC       = VMCSolver::LOCAL_ENERGY_GENERIC;
    static const int LOCAL_ENERGY_HELIUM_1      = VMCSolver::LOCAL_ENERGY_HELIUM_1;
    static const int LOCAL_ENERGY_HELIUM_2      = VMCSolver::LOCAL_ENERGY_HELIUM_2;
    static const int LOCAL_ENERGY_HYDROGEN      = VMCSolver::LOCAL_ENERGY_HYDROGEN;
    static const int LOCAL_ENERGY_GENERIC_NOCOR = VMCSolver::LOCAL_ENERGY_GENERIC_NOCOR;
    static const int LOCAL_ENERGY_SLATER        = VMCSolver::LOCAL_ENERGY_SLATER;
    static const int LOCAL_ENERGY_SLATER_NOCOR  = VMCSolver::LOCAL_ENERGY_SLATER_NOCOR;
    static const int WAVE_FUNCTION_1            = VMCSolver::WAVE_FUNCTION_1;
    static const int WAVE_FUNCTION_2            = VMCSolver::WAVE_FUNCTION_2;
    static const int WAVE_FUNCTION_BERYLLIUM_1  = VMCSolver::WAVE_FUNCTION_BERYLLIUM_1;
    static const int WAVE_FUNCTION_BERYLLIUM_2  = VMCSolver::WAVE_FUNCTION_BERYLLIUM_2;

    bool runIntegration();
    bool initSolver(VMCSolver& solver);
    bool initFromFile(std::string fName);
    void exportParamters(std::string fName);
    void exportDensity(std::string fName);
    void exportEnergyArray(std::string fName);
    void exportPositions(std::string fName);
    void useWaveFunction1();
    void useWaveFunction2();
    void useWaveFunctionBeryllium1();
    void useWaveFunctionBeryllium2();
    void useLocalEnergyHelium1();
    void useLocalEnergyHelium2();
    void useLocalEnergyHydrogen();
    void useLocalEnergyGeneric();
    void useLocalEnergyGenericNoCor();
    void useLocalEnergySlater();
    void useLocalEnergySlaterNoCor();
    void recordDensity(bool param, int bins = 9, double maxPos = 2);
    void recordEnergyArray(bool param);
    void recordR12Mean(bool param);
    void recordPositions(bool param);
    void setStepLength(double stepLength);
    void useParallel(bool param);
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


    int threads;
    // Paramters to the solver.
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
    int bins;
    double rMax;
    bool outputSupressed;
    bool recordingDensity;
    bool recordingEnergyArray;
    bool recordingR12Mean;
    bool recordingPositions;
    bool importanceSampling;
    bool efficientSlater;
    bool parallel;

    Vector energyArray;

    // Results from the solver
    double mean;
    double energy;
    double energySquared;
    double acceptanceRatio;

};

#endif // VMCSOLVER_H
