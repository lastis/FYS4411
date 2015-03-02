#include "vmcsolver.h"

using namespace CPhys;
using namespace std;

VMCSolver::VMCSolver(){
    clearAll();
}


bool VMCSolver::runIntegration(){

    reset();
    ready = initRunVariables();
    if (!ready) {
	cout << "Error: Solver not initialized properly, integration not running."
	    << endl;
	return false;
    }

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++) {
	if (useImportanceSampling) 
	    runQuantumWalk();
	else 
	    runRandomWalk();
    }

    mean = rAbsSum/(nCycles);
    energy = energySum/(nCycles * nParticles);
    energySquared = energySquaredSum/(nCycles * nParticles);

    // Output 
    if (outputSupressed) {
        outputSupressed = false;
        return true;
    }
    cout << "Energy: " << energy << " Energy (squared sum): " 
	<< energySquared << endl;
    cout << "Variance : " << energySquared - energy*energy << endl;
    return true;
}

inline void VMCSolver::runQuantumWalk(){
    double greensFunction;
    // Store the current value of the wave function
    waveFuncValOld = (this->*getWaveFuncVal)(prOld);
    updateQuantumForce(prOld,pqForceOld,waveFuncValOld);

    // New position to test
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            prNew[i][j] = prOld[i][j] + Random::gauss(idum)*sqrt(timeStep) 
		+ pqForceOld[i][j]*timeStep*D;
        }
	//For the other particles we need to set the position to the old
	//position since we move only one particle at the time. 
	for (int k = 0; k < nParticles; k++) {
		if (k != i) {
			for (int j = 0; j < nDimensions; j++) {
				prNew[k][j] = prOld[k][j];
			}
		}
	}

        // Recalculate the value of the wave function
        waveFuncValNew = (this->*getWaveFuncVal)(prNew);
	updateQuantumForce(prNew, pqForceNew, waveFuncValNew);

	// Compute the log ratio of the greens functions to be used in the 
	// Metropolis-Hastings algorithm.
	greensFunction = 0;
	for (int j = 0; j < nDimensions; j++) {
		greensFunction += 0.5*(pqForceOld[i][j]+pqForceNew[i][j])*
		    (D*timeStep*0.5*(pqForceOld[i][j]-pqForceNew[i][j])
		    - prNew[i][j] + prOld[i][j]);
	}
	greensFunction = exp(greensFunction);

        // Check for step acceptance (if yes, 
        // update position, if no, reset position)
	// The metropolis test is performed by moving one particle at the time.
        if(Random::ran2(idum) <= greensFunction*(waveFuncValNew*waveFuncValNew)
                / (waveFuncValOld*waveFuncValOld)) {
            for(int j = 0; j < nDimensions; j++) {
                prOld[i][j] = prNew[i][j];
		pqForceOld[i][j] = pqForceNew[i][j];
                waveFuncValOld = waveFuncValNew;
            }
            accepts++;
        } else {
            for(int j = 0; j < nDimensions; j++) {
                prNew[i][j] = prOld[i][j];
		pqForceNew[i][j] = pqForceOld[i][j];
            }
            rejects++;
        }

        // update energies
        deltaE = (this->*getLocalEnergy)(prNew); 
        
        energySum += deltaE;
        energySquaredSum += deltaE*deltaE;
    }

    double rsq = 0;
    for(int j = 0; j < nDimensions; j++) {
        rsq += (prNew[1][j] - prNew[0][j])*(prNew[1][j] - prNew[0][j]);
    }
    rAbsSum += sqrt(rsq);

}

inline void VMCSolver::runRandomWalk(){
    // Store the current value of the wave function
    waveFuncValOld = (this->*getWaveFuncVal)(prOld);

    // New position to test
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            prNew[i][j] = prOld[i][j] + stepLength*(Random::ran2(idum) - 0.5);
        }

        // Recalculate the value of the wave function
        waveFuncValNew = (this->*getWaveFuncVal)(prNew);
        // Check for step acceptance (if yes, 
        // update position, if no, reset position)
        if(Random::ran2(idum) <= (waveFuncValNew*waveFuncValNew)
                / (waveFuncValOld*waveFuncValOld)) {
            for(int j = 0; j < nDimensions; j++) {
                prOld[i][j] = prNew[i][j];
                waveFuncValOld = waveFuncValNew;
            }
            accepts++;
        } else {
            for(int j = 0; j < nDimensions; j++) {
                prNew[i][j] = prOld[i][j];
            }
            rejects++;
        }

        // update energies
        deltaE = (this->*getLocalEnergy)(prNew); 
        
        energySum += deltaE;
        energySquaredSum += deltaE*deltaE;
    }

    double rsq = 0;
    for(int j = 0; j < nDimensions; j++) {
        rsq += (prNew[1][j] - prNew[0][j])*(prNew[1][j] - prNew[0][j]);
    }
    rAbsSum += sqrt(rsq);
}

void VMCSolver::supressOutput(){
    outputSupressed = true;
}

bool VMCSolver::initRunVariables(){

    // Set the wave function as a function pointer
    if (waveFunction == WAVE_FUNCTION_1)
	getWaveFuncVal = &VMCSolver::getWaveFunc1Val;
    else if (waveFunction == WAVE_FUNCTION_2)
	getWaveFuncVal = &VMCSolver::getWaveFunc2Val;
    else {
	cout << "Error: Wave function not set, integration not running."
	    << endl;
	return false;
    }

    // Set the local energy function as a function pointer
    if (localEnergyFunction == LOCAL_ENERGY_GENERIC)
        getLocalEnergy = &VMCSolver::getLocalEnergyGeneric;
    else if (localEnergyFunction == LOCAL_ENERGY_HELIUM)
        getLocalEnergy = &VMCSolver::getLocalEnergyHelium;
    /* else if (localEnergyFunction == LOCAL_ENERGY_HELIUM); */
	/* double (VMCSolver::*localEnergy)(double** r) */ 
	/*     = &VMCSolver::localEnergyHelium; */
    else {
	cout << "Error: Local energy function not set, integration not running."
	    << endl;
	return false;
    }

    // Initialize arrays

    if (useImportanceSampling) {
	qForceOld = Matrix(nParticles, nDimensions);
	qForceNew = Matrix(nParticles, nDimensions);
	pqForceOld = qForceOld.getArrayPointer();
	pqForceNew = qForceNew.getArrayPointer();
    }
    rOld = Matrix(nParticles, nDimensions);
    rNew = Matrix(nParticles, nDimensions);
    prNew = rNew.getArrayPointer();
    prOld = rOld.getArrayPointer();

    // initial trial positions
    if (useImportanceSampling == true) {
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                prOld[i][j] = Random::gauss(idum)*sqrt(timeStep);
            }
        }
    }
    else {
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                /* prOld[i][j] = stepLength * (Random::ran2(idum) - 0.5); */
                prOld[i][j] = (Random::ran2(idum) - 0.5);
            }
        }
    }

    rNew = rOld;
    prNew = rNew.getArrayPointer();
    prOld = rOld.getArrayPointer();
    // Finished without error (hopefully).
    return true;
}

double VMCSolver::getLocalEnergyHelium(double** r){
    double* r1 = r[0];
    double* r2 = r[1];
    double temp = 0;
    double r12Abs = 0;
    double r1Abs = 0;
    double r2Abs = 0;
    double r1r2 = 0; // Dot product.
    for(int j = 0; j < nDimensions; j++) {
	temp = r1[j] * r1[j];
	r1Abs += temp;
	temp = r2[j] * r2[j];
	r2Abs += temp;
	temp = (r1[j] - r2[j]) * (r1[j] - r2[j]);
	r12Abs += temp;
	temp = r1[j]*r2[j];
	r1r2 += temp;
    }
    r1Abs = sqrt(r1Abs);
    r2Abs = sqrt(r2Abs);
    r12Abs = sqrt(r12Abs);
    double E1 = (alpha-charge)*(1/r1Abs + 1/r2Abs) + 1/r12Abs - alpha*alpha;
    double betaR12 = 1/((1+beta*r12Abs));
    return E1 + betaR12*betaR12/2 * (
	    alpha*(r1Abs + r2Abs)/r12Abs*(1-(r1r2/(r1Abs*r2Abs)))
	    - betaR12*betaR12/2 - 2/r12Abs + 2*beta*betaR12
	    );
}

void VMCSolver::updateQuantumForce(double** r, double ** qForce, double factor){
    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;
    double waveFunctionCurrent;
    waveFunctionCurrent = (this->*getWaveFuncVal)(r);

    double rPlus;
    double rMinus;
    double r0;
    double kineticEnergy = 0;
    // Kinetic energy

    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
	    rPlus 	= r[i][j] + h;
	    rMinus 	= r[i][j] - h;
	    r0 		= r[i][j];

	    r[i][j] = rMinus;
	    waveFunctionMinus = (this->*getWaveFuncVal)(r);
	    r[i][j] = rPlus;
	    waveFunctionPlus = (this->*getWaveFuncVal)(r);
	    r[i][j] = r0;
	    qForce[i][j] = 
		(waveFunctionPlus - waveFunctionMinus)*h/factor;
        }
    }
}

double VMCSolver::getLocalEnergyGeneric(double** r){
    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;
    double waveFunctionCurrent;
    waveFunctionCurrent = (this->*getWaveFuncVal)(r);

    // Kinetic energy

    double rPlus;
    double rMinus;
    double r0;
    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
	    rPlus 	= r[i][j] + h;
	    rMinus 	= r[i][j] - h;
	    r0 		= r[i][j];

	    r[i][j] = rMinus;
	    waveFunctionMinus = (this->*getWaveFuncVal)(r);
	    r[i][j] = rPlus;
	    waveFunctionPlus = (this->*getWaveFuncVal)(r);
	    r[i][j] = r0;
	    kineticEnergy -= (waveFunctionMinus + waveFunctionPlus 
			- 2 * waveFunctionCurrent);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    // Potential energy
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r[i][j]*r[i][j];
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r[i][k] - r[j][k]) * (r[i][k] - r[j][k]);
            }
            potentialEnergy += 1 / sqrt(r12);
        }
    }

    return kineticEnergy + potentialEnergy;
}

void VMCSolver::setStepLength(double stepLength){
    this->stepLength = stepLength;
}

double VMCSolver::getStepLength(){
    return stepLength;
}

double VMCSolver::getEnergy(){
    return energy;
}

double VMCSolver::getR12Mean(){
    return mean;
}

bool VMCSolver::initFromFile(std::string fName){
    ifstream myFile;
    string  paramName;
    string  discard;

    string adress = "../../res/" + fName;
    myFile.open(adress.c_str());
    if (!myFile) {
	cout << fName << "does not exist. Solver could not initialize." << endl;
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

    myFile.close();

    ready = true;
    outputSupressed = false;
    return true;
}

void VMCSolver::setLocalEnergyHelium(){
    if(nParticles != 2) {
	cout << "Cannot use this analytic local energy function "
	    << "for other than 2 particles." << endl;
	cout << "Using generic one." << endl;
	localEnergyFunction = LOCAL_ENERGY_GENERIC;
	return;
    }
    localEnergyFunction = LOCAL_ENERGY_HELIUM;
}

void VMCSolver::setLocalEnergyHydrogen(){
    if(nParticles != 1) {
	cout << "Cannot use this analytic local energy function " 
	    << "for other than 1 particle." << endl;
	cout << "Using generic one." << endl;
	localEnergyFunction = LOCAL_ENERGY_GENERIC;
	return;
    }
    localEnergyFunction = LOCAL_ENERGY_HYDROGEN;
}

void VMCSolver::setImportanceSampling(){
    if (timeStep == 0) {
	cout << "Error : Cannot use importance sampling with timeStep = 0" 
	    << endl;
	return;
    }
    if (D == 0) {
	cout << "Error : Cannot use importance sampling with D = 0" 
	    << endl;
	return;
    }
    useImportanceSampling = true;
}

void VMCSolver::setLocalEnergyGeneric(){
    localEnergyFunction = LOCAL_ENERGY_GENERIC;
}

void VMCSolver::setWaveFunction1(){
    waveFunction = WAVE_FUNCTION_1;
}

void VMCSolver::setWaveFunction2(){
    waveFunction = WAVE_FUNCTION_2;
}

void VMCSolver::reset(){
    mean = 0;
    accepts = 0;
    rejects = 0;

    energy = 0;
    energySquared = 0;
    energySum = 0;
    energySquaredSum = 0;
    rAbsSum = 0;

    rOld.reset();
    rNew.reset();
    qForceOld.reset();
    qForceNew.reset();

    deltaE = 0;
    waveFuncValOld = 0;
    waveFuncValNew = 0;
}

void VMCSolver::clearAll(){
    reset();
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

    ready = false;
    useImportanceSampling = false;
}

double VMCSolver::getAcceptanceRatio(){
    return double(accepts)*100/(rejects+accepts);
}

double VMCSolver::getWaveFunc1Val(double** r)
{
    double argument = 0;
    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r[i][j] * r[i][j];
        }
        argument += sqrt(rSingleParticle);
    }
    return exp(-argument * alpha);
}

double VMCSolver::getWaveFunc2Val(double** r){
    double* r1 = r[0];
    double* r2 = r[1];
    double temp = 0;
    double r1Abs = 0;
    double r2Abs = 0;
    double r12 = 0;
    for(int j = 0; j < nDimensions; j++) {
	temp = r1[j] * r1[j];
	r1Abs += temp;
	temp = r2[j] * r2[j];
	r2Abs += temp;

	temp = (r2[j] - r1[j]) * (r2[j] - r1[j]);
	r12 += temp;
    }
    r1Abs = sqrt(r1Abs);
    r2Abs = sqrt(r2Abs);
    r12 = sqrt(r12);
    return exp(-(r1Abs + r2Abs)*alpha)*exp(r12/(2*(1+beta*r12)));
}

void VMCSolver::exportParamters(std::string fName){
    string adress = "../../res/" + fName;
    ofstream myFile;
    cout << "Dumption to file : " << fName << endl;
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
    myFile.close();
}
