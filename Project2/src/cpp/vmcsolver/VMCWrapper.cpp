#include "VMCWrapper.h"

using namespace CPhys;
using namespace std;

VMCWrapper::VMCWrapper(){
    clear();
}

VMCSolver VMCWrapper::getInitializedSolver(){
    validateParamters();
    VMCSolver solver = VMCSolver();
    initSolver(solver);
    solver.initRunVariables();
    return solver;
}

bool VMCWrapper::initSolver(VMCSolver& solver){
    solver.alpha = alpha;
    solver.beta = beta;
    solver.waveFunction = waveFunction;
    solver.localEnergyFunction = localEnergyFunction;
    solver.charge = charge;
    solver.nDimensions = nDimensions;
    solver.nCycles = nCycles;
    solver.nParticles = nParticles;
    solver.stepLength = stepLength;
    solver.h = h;
    solver.hInv = hInv;
    solver.h2Inv = h2Inv;
    solver.idum = idum;
    solver.timeStep = timeStep;
    solver.D = D;
    solver.importanceSampling = importanceSampling;
    solver.efficientSlater = efficientSlater;
    return true;
}

void VMCWrapper::setStepLength(double stepLength){
    this->stepLength = stepLength;
}

bool VMCWrapper::initFromFile(std::string fName){
    ifstream myFile;
    string  paramName;
    string  discard;

    string adress = "../../../res/" + fName;
    myFile.open(adress.c_str());
    if (!myFile) {
        cout << fName << " does not exist. Solver could not initialize." << endl;
        return false;
    }
    cout << "Initializing from file : " << fName << endl;
    myFile >> paramName >> discard >> charge;
    myFile >> paramName >> discard >> alpha;
    myFile >> paramName >> discard >> beta;
    myFile >> paramName >> discard >> nDimensions;
    myFile >> paramName >> discard >> nParticles;
    myFile >> paramName >> discard >> stepLength;
    myFile >> paramName >> discard >> nCycles;
    myFile >> paramName >> discard >> waveFunction;
    myFile >> paramName >> discard >> h;
    myFile >> paramName >> discard >> hInv;
    myFile >> paramName >> discard >> h2Inv;
    myFile >> paramName >> discard >> idum;
    myFile >> paramName >> discard >> localEnergyFunction;
    myFile >> paramName >> discard >> timeStep;
    myFile >> paramName >> discard >> D;
    myFile >> paramName >> discard >> importanceSampling;
    myFile >> paramName >> discard >> efficientSlater;

    myFile.close();

    return true;
}

void VMCWrapper::useLocalEnergyHelium1(){
    localEnergyFunction = LOCAL_ENERGY_HELIUM_1;
}

void VMCWrapper::useLocalEnergyHelium2(){
    localEnergyFunction = LOCAL_ENERGY_HELIUM_2;
}

void VMCWrapper::useLocalEnergyHydrogen(){
    localEnergyFunction = LOCAL_ENERGY_HYDROGEN;
}

void VMCWrapper::useLocalEnergySlater(){
    localEnergyFunction = LOCAL_ENERGY_SLATER;
}

void VMCWrapper::useLocalEnergySlaterNoCor(){
    localEnergyFunction = LOCAL_ENERGY_SLATER_NOCOR;
}

void VMCWrapper::useEfficientSlater(bool param){
    efficientSlater = param;
}

void VMCWrapper::useImportanceSampling(bool param){
    importanceSampling = param;
}

void VMCWrapper::useLocalEnergyGeneric(){
    localEnergyFunction = LOCAL_ENERGY_GENERIC;
}

void VMCWrapper::useLocalEnergyGenericNoCor(){
    localEnergyFunction = LOCAL_ENERGY_GENERIC_NOCOR;
}

void VMCWrapper::useWaveFunction1(){
    waveFunction = WAVE_FUNCTION_1;
}

void VMCWrapper::useWaveFunction2(){
    waveFunction = WAVE_FUNCTION_2;
}

void VMCWrapper::useWaveFunctionBeryllium1(){
    waveFunction = WAVE_FUNCTION_BERYLLIUM_1;
}

void VMCWrapper::useWaveFunctionBeryllium2(){
    waveFunction = WAVE_FUNCTION_BERYLLIUM_2;
}

void VMCWrapper::useWaveFunctionHeliumGTO(){
    waveFunction = WAVE_FUNCTION_HELIUM_GTO;
}

void VMCWrapper::clear(){
    waveFunction = 0;
    localEnergyFunction = 0;
    nDimensions = 0;
    charge = 0;
    stepLength = 0;
    nParticles = 0;
    h = 0;
    hInv = 0;
    h2Inv = 0;
    idum = 1;
    alpha = 0;
    beta = 0;
    nCycles = 0;
    timeStep = 0;
    D = 0;

    useImportanceSampling(false);
    useEfficientSlater(false);
}

bool VMCWrapper::validateParamters(){
    bool valid = true;
    if (efficientSlater) {
        if (nParticles > 10) {
            cout << "Error : Slater determinant has not been implemented "
                << "for more than 10 particles!" << endl;
            valid = false;
        }
    }
    if(importanceSampling){
        if (timeStep == 0) {
            cout << "Error : Cannot use importance sampling with timeStep = 0." 
                << endl;
            valid = false;
            useImportanceSampling(false);
        }
        if (D == 0) {
            cout << "Error : Cannot use importance sampling with D = 0." 
                << endl;
            valid = false;
            useImportanceSampling(false);
        }
    }
    if(nDimensions == 0 ) {
        cout << "Cannot simulate with 0 dimensions." << endl;
        valid = false;
    }
    if(localEnergyFunction == LOCAL_ENERGY_HYDROGEN && nParticles != 1) {
        cout << "Cannot use this analytic local energy function " 
            << "for other than 1 particle." << endl;
        valid = false;
    }
    if (localEnergyFunction == LOCAL_ENERGY_GENERIC_NOCOR && waveFunction == 0) {
        cout << "Cannot use generic local energy function "
            << "with no wave function." << endl;
        valid = false;
    }
    if (localEnergyFunction == LOCAL_ENERGY_GENERIC && waveFunction == 0) {
        cout << "Cannot use generic local energy function "
            << "with no wave function." << endl;
        valid = false;
    }
    if(localEnergyFunction == LOCAL_ENERGY_HELIUM_1 && nParticles != 2) {
        cout << "Cannot use this analytic local energy function "
            << "for other than 2 particles." << endl;
        valid = false;
    }
    if(localEnergyFunction == LOCAL_ENERGY_HELIUM_2 && nParticles != 2) {
        cout << "Cannot use this analytic local energy function "
            << "for other than 2 particles." << endl;
        valid = false;
    }

    if (waveFunction == WAVE_FUNCTION_2 && nParticles != 2){
        cout << "Cannot use this wave function  "
            << "for other than 2 particles." << endl;
        valid = false;
    }
    if (waveFunction == WAVE_FUNCTION_BERYLLIUM_1 && nParticles != 4){
        cout << "Cannot use this wave function  "
            << "for other than 4 particles." << endl;
        valid = false;
    }
    if (waveFunction == WAVE_FUNCTION_BERYLLIUM_2 && nParticles != 4){
        cout << "Cannot use this wave function  "
            << "for other than 4 particles." << endl;
        valid = false;
    }
    /* if (recordingR12Mean && nParticles != 2) { */
    /*     cout << "Cannot use record r12 mean   " */
    /*         << "for other than 2 particles." << endl; */
    /*     valid = false; */
    /* } */

    return valid;
}

void VMCWrapper::exportParamters(std::string fName){
    string adress = "../../../res/" + fName;
    ofstream myFile;
    cout << "Dumption paramters to file : " << fName << endl;
    myFile.open(adress.c_str());
    myFile << "charge = " << charge <<  endl;
    myFile << "alpha = " << alpha << endl;
    myFile << "beta = " << beta << endl;
    myFile << "nDimensions = " << nDimensions <<  endl;
    myFile << "nParticles = " << nParticles <<  endl;
    myFile << "stepLength = " << stepLength << endl;
    myFile << "nCycles = " << nCycles << endl;
    myFile << "waveFunction = " << waveFunction << endl;
    myFile << "h = " << h << endl;
    myFile << "hInv = " << hInv << endl;
    myFile << "h2Inv = " << h2Inv << endl;
    myFile << "idum = " << idum << endl;
    myFile << "localEnergyFunction = " << localEnergyFunction << endl;
    myFile << "timeStep = " << timeStep << endl;
    myFile << "D = " << D << endl;
    myFile << "useImportanceSampling = " << importanceSampling << endl;
    /* myFile << "recordDensity = " << recordingDensity <<  endl; */
    /* myFile << "recordEnergyArray = " << recordingEnergyArray <<  endl; */
    /* myFile << "recordR12Mean = " << recordingR12Mean <<  endl; */
    /* myFile << "recordPositions = " << recordingPositions <<  endl; */
    myFile << "useEfficientSlater = " << efficientSlater <<  endl;
    myFile.close();
}

/* void VMCWrapper::exportDensity(std::string fName){ */
/*     /1* string adress = "../../../res/" + fName; *1/ */
/*     /1* ofstream myFile; *1/ */
/*     /1* cout << "Dumption densities to file : " << fName << endl; *1/ */
/*     /1* myFile.open(adress.c_str()); *1/ */
/*     /1* myFile << rMax << endl; *1/ */
/*     /1* for (int i = 0; i < nParticles; i++) { *1/ */
/* 	/1* for (int j = 0; j < bins; j++) { *1/ */
/* 	    /1* myFile << pDensity[i][j] << " "; *1/ */
/* 	/1* } *1/ */
/* 	/1* myFile << endl; *1/ */
/*     /1* } *1/ */
/*     /1* myFile.close(); *1/ */
/* } */

/* void VMCWrapper::exportEnergyArray(std::string fName){ */
/*     /1* string adress = "../../../res/" + fName; *1/ */
/*     /1* ofstream myFile; *1/ */
/*     /1* cout << "Dumption energies to file : " << fName << endl; *1/ */
/*     /1* myFile.open(adress.c_str()); *1/ */
/*     /1* for (int i = 0; i < nCycles; i++) { *1/ */
/* 	/1* myFile << pEnergyArray[i] << " "; *1/ */
/*     /1* } *1/ */
/*     /1* myFile.close(); *1/ */
/* } */

/* void VMCWrapper::exportPositions(std::string fName){ */
/*     /1* string adress = "../../../res/" + fName; *1/ */
/*     /1* ofstream myFile; *1/ */
/*     /1* cout << "Dumption energies to file : " << fName << endl; *1/ */
/*     /1* myFile.open(adress.c_str()); *1/ */
/*     /1* for (int i = 0; i < nParticles; i++) { *1/ */
/*     /1* 	for (int j = 0; j < nCycles; j++) { *1/ */
/* 	    /1* myFile << pPositions[i][j] << " "; *1/ */
/*     /1* 	} *1/ */
/* 	/1* myFile << endl; *1/ */
/*     /1* } *1/ */
/*     /1* myFile.close(); *1/ */
/* } */
