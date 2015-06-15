#ifndef _VMCWRAPPER_H_INCLUDED_
#define _VMCWRAPPER_H_INCLUDED_

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include "VMCSolver.h"
#include "VMCSolverGto.h"
#include "VMCSolverGtoI.h"
#include "../CPhys/CPhys.h"

class VMCWrapper
{
public:
    VMCWrapper();
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

    bool initSolver(VMCSolver& solver);
    bool initSolver(VMCSolverGto& solver);
    bool initSolver(VMCSolverGtoI& solver);
    bool initFromFile(std::string fName);
    void exportParamters(std::string fName);
    void useWaveFunction1();
    void useWaveFunction2();
    void useWaveFunctionBeryllium1();
    void useWaveFunctionBeryllium2();
    void useWaveFunctionHeliumGTO();
    void useWaveFunctionBerylliumGTO();
    void useWaveFunctionNeonGTO();
    void useLocalEnergyHelium1();
    void useLocalEnergyHelium2();
    void useLocalEnergyHydrogen();
    void useLocalEnergyGeneric();
    void useLocalEnergyGenericNoCor();
    void useLocalEnergySlater();
    void useLocalEnergySlaterNoCor();
    void setStepLength(double stepLength);
    void useImportanceSampling(bool param);
    void useEfficientSlater(bool param);
    void clear();
    bool validateParamters();
    bool validateParamtersGto();
    bool validateParamtersGtoI();
    VMCSolver getInitializedSolver();
    VMCSolverGto getInitializedSolverGto();
    VMCSolverGtoI getInitializedSolverGtoI();

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
    double hInv;
    double h2Inv;
    long idum;
    double timeStep;
    double D;
    bool importanceSampling;
    bool efficientSlater;
};

#endif  // VMCSOLVER_H
