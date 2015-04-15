#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include "../CPhys/CPhys.h"

class VMCWrapper
{
public:
    VMCWrapper();
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
    void setRecordDensity(bool param, int bins = 9, double maxPos = 2);
    void setRecordEnergyArray(bool param);
    void setRecordR12Mean(bool param);
    void setRecordPositions(bool param);
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

    // Shared variables
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
    bool recordDensity;
    bool recordChargeDensity;
    bool recordEnergyArray;
    bool recordR12Mean;
    bool recordPositions;
    bool importanceSampling;
    bool efficientSlater;
    bool parallel;
};

#endif // VMCSOLVER_H
