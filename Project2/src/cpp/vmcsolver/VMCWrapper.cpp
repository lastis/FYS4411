#include "VMCWrapper.h"

using namespace CPhys;
using namespace std;

VMCWrapper::VMCWrapper(){
    clear();
}

bool VMCWrapper::runIntegration(){
    bool val;
    if (parallel) {
        VMCSolver1 solver = VMCSolver1();
        #pragma omp parallel private(solver)
        {
            initSolver(solver);
            solver.nCycles = nCycles/omp_get_num_threads();
            val = solver.runIntegration();
            printf("Energy : %f\n", solver.getEnergy());
            solver.nCycles = nCycles;
        }
    }
    else {
        VMCSolver1 solver = VMCSolver1();
        initSolver(solver);
        val = solver.runIntegration();
        mean = solver.getR12Mean();
        energy = solver.getEnergy();
        energySquared = solver.getEnergySquared();
        acceptanceRatio = solver.getAcceptanceRatio();
    }
    return val;
}

bool VMCWrapper::initSolver(VMCSolver1& solver){
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
    solver.h2 = h2;
    solver.idum = idum;
    solver.timeStep = timeStep;
    solver.D = D;
    solver.outputSupressed = outputSupressed;
    solver.recordingDensity = recordingDensity;
    solver.recordingEnergyArray = recordingEnergyArray;
    solver.recordingR12Mean = recordingR12Mean;
    solver.recordingPositions = recordingPositions;
    solver.importanceSampling = importanceSampling;
    solver.efficientSlater = efficientSlater;
    solver.parallel = parallel;
    solver.bins = bins;
    solver.rMax = rMax;
}

void VMCWrapper::supressOutput(){
    outputSupressed = true;
}

void VMCWrapper::setStepLength(double stepLength){
    this->stepLength = stepLength;
}

double VMCWrapper::getStepLength(){
    return stepLength;
}

double VMCWrapper::getEnergy(){
    return energy;
}

double VMCWrapper::getEnergySquared(){
    return energySquared;
}

double VMCWrapper::getR12Mean(){
    return mean;
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
    myFile >> paramName >> discard >> h2;
    myFile >> paramName >> discard >> idum;
    myFile >> paramName >> discard >> localEnergyFunction;
    myFile >> paramName >> discard >> timeStep;
    myFile >> paramName >> discard >> D;
    myFile >> paramName >> discard >> importanceSampling;
    myFile >> paramName >> discard >> recordingDensity;
    myFile >> paramName >> discard >> recordingEnergyArray;
    myFile >> paramName >> discard >> recordingR12Mean;
    myFile >> paramName >> discard >> recordingPositions;
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

void VMCWrapper::useParallel(bool param){
    parallel = param;
}

void VMCWrapper::setRecordEnergyArray(bool param){
    recordingEnergyArray = param;
}

void VMCWrapper::setRecordDensity(bool param, int bins, double maxPos){
    // This is the only place where bins and rMax are set. But 
    // this function is called on clear().
    recordingDensity = param;
    this->bins = bins;
    rMax = maxPos;
}

void VMCWrapper::setRecordR12Mean(bool param){
    recordingR12Mean = param;
}
void VMCWrapper::setRecordPositions(bool param){
    recordingPositions = param;
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

void VMCWrapper::clear(){
    waveFunction = WAVE_FUNCTION_1;
    localEnergyFunction = LOCAL_ENERGY_GENERIC;
    nDimensions = 0;
    charge = 0;
    stepLength = 0;
    nParticles = 0;
    h = 0;
    h2 = 0;
    idum = 1;
    alpha = 0;
    beta = 0;
    nCycles = 0;
    timeStep = 0;
    D = 0;

    outputSupressed = false;
    useImportanceSampling(false);
    useEfficientSlater(false);
    setRecordDensity(false);
    setRecordEnergyArray(false);
    setRecordR12Mean(false);
    setRecordPositions(false);
}

double VMCWrapper::getAcceptanceRatio(){
    return acceptanceRatio;
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
    if(localEnergyFunction == LOCAL_ENERGY_HYDROGEN && nParticles != 1) {
        cout << "Cannot use this analytic local energy function " 
            << "for other than 1 particle." << endl;
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
    if (recordingR12Mean && nParticles != 2) {
        cout << "Cannot use record r12 mean   "
            << "for other than 2 particles." << endl;
        valid = false;
    }

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
    myFile << "h2 = " << h2 << endl;
    myFile << "idum = " << idum << endl;
    myFile << "localEnergyFunction = " << localEnergyFunction << endl;
    myFile << "timeStep = " << timeStep << endl;
    myFile << "D = " << D << endl;
    myFile << "useImportanceSampling = " << importanceSampling << endl;
    myFile << "recordDensity = " << recordingDensity <<  endl;
    myFile << "recordEnergyArray = " << recordingEnergyArray <<  endl;
    myFile << "recordR12Mean = " << recordingR12Mean <<  endl;
    myFile << "recordPositions = " << recordingPositions <<  endl;
    myFile << "useEfficientSlater = " << efficientSlater <<  endl;
    myFile.close();
}

void VMCWrapper::exportDensity(std::string fName){
    /* string adress = "../../../res/" + fName; */
    /* ofstream myFile; */
    /* cout << "Dumption densities to file : " << fName << endl; */
    /* myFile.open(adress.c_str()); */
    /* myFile << rMax << endl; */
    /* for (int i = 0; i < nParticles; i++) { */
	/* for (int j = 0; j < bins; j++) { */
	    /* myFile << pDensity[i][j] << " "; */
	/* } */
	/* myFile << endl; */
    /* } */
    /* myFile.close(); */
}

void VMCWrapper::exportEnergyArray(std::string fName){
    /* string adress = "../../../res/" + fName; */
    /* ofstream myFile; */
    /* cout << "Dumption energies to file : " << fName << endl; */
    /* myFile.open(adress.c_str()); */
    /* for (int i = 0; i < nCycles; i++) { */
	/* myFile << pEnergyArray[i] << " "; */
    /* } */
    /* myFile.close(); */
}

void VMCWrapper::exportPositions(std::string fName){
    /* string adress = "../../../res/" + fName; */
    /* ofstream myFile; */
    /* cout << "Dumption energies to file : " << fName << endl; */
    /* myFile.open(adress.c_str()); */
    /* for (int i = 0; i < nParticles; i++) { */
    /* 	for (int j = 0; j < nCycles; j++) { */
	    /* myFile << pPositions[i][j] << " "; */
    /* 	} */
	/* myFile << endl; */
    /* } */
    /* myFile.close(); */
}
