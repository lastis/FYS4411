#include "vmcsolver.h"

using namespace CPhys;
using namespace std;

VMCSolver::VMCSolver(){
    clearAll();
}


bool VMCSolver::runIntegration(){

    if (!initialized) {
	cout << "Error: Solver not initialized, integration not running."
	    << endl;
	return false;
    }
    // Set the wave function as a function pointer
    if (waveFunction == WAVE_FUNCTION_1)
	wave = &VMCSolver::wave1;
    else if (waveFunction == WAVE_FUNCTION_2)
	wave = &VMCSolver::wave2;
    else {
	cout << "Error: Wave function not set, integration not running."
	    << endl;
	return false;
    }
    // Set the local energy function as a function pointer
    if (localEnergyFunction == LOCAL_ENERGY_GENERIC)
	localEnergy = &VMCSolver::localEnergyGeneric;
    else if (localEnergyFunction == LOCAL_ENERGY_HELIUM)
	localEnergy = &VMCSolver::localEnergyHelium;
    /* else if (localEnergyFunction == LOCAL_ENERGY_HELIUM); */
	/* double (VMCSolver::*localEnergy)(double** r) */ 
	/*     = &VMCSolver::localEnergyHelium; */
    else {
	cout << "Error: Local energy function not set, integration not running."
	    << endl;
	return false;
    }


    double waveFunctionOld = 0;
    double waveFunctionNew = 0;

    double energySum = 0;
    double energySquaredSum = 0;

    double deltaE = 0;
    double rSum = 0;

    reset();

    initPositions();

    rNew = rOld;
    prNew = rNew.getArrayPointer();
    prOld = rOld.getArrayPointer();

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++) {

        // Store the current value of the wave function
	waveFunctionOld = (this->*wave)(prOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                prNew[i][j] = prOld[i][j] + stepLength*(Random::ran2(idum) - 0.5);
            }

            // Recalculate the value of the wave function
	    waveFunctionNew = (this->*wave)(prNew);
            // Check for step acceptance (if yes, update position, if no, reset position)
            if(Random::ran2(idum) <= (waveFunctionNew*waveFunctionNew) / 
			    (waveFunctionOld*waveFunctionOld)) {
                for(int j = 0; j < nDimensions; j++) {
                    prOld[i][j] = prNew[i][j];
                    waveFunctionOld = waveFunctionNew;
		    accepts++;
                }
            } else {
                for(int j = 0; j < nDimensions; j++) {
                    prNew[i][j] = prOld[i][j];
		    rejects++;
                }
            }

            // update energies
	    deltaE = (this->*localEnergy)(prNew); 
            
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }

	double rsq = 0;
	for(int j = 0; j < nDimensions; j++) {
	    rsq += (prOld[1][j] - prOld[0][j])*(prOld[1][j] - prOld[0][j]);
	}
	rSum += sqrt(rsq);
    }
    mean = rSum/(nCycles);
    energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);

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

void VMCSolver::supressOutput(){
    outputSupressed = true;
}

void VMCSolver::initPositions(){
    rOld = Matrix(nParticles, nDimensions);
    rNew = Matrix(nParticles, nDimensions);
    prNew = rNew.getArrayPointer();
    prOld = rOld.getArrayPointer();

    rPlus = Matrix(nParticles, nDimensions);
    rMinus = Matrix(nParticles, nDimensions);

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            /* prOld[i][j] = stepLength * (Random::ran2(idum) - 0.5); */
            prOld[i][j] = (Random::ran2(idum) - 0.5);
        }
    }
}

double VMCSolver::localEnergyHelium(double** r){
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

double VMCSolver::localEnergyGeneric(double** r){

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;
    double waveFunctionCurrent;
    waveFunctionCurrent = (this->*wave)(r);

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
	    waveFunctionMinus = (this->*wave)(r);
	    r[i][j] = rPlus;
	    waveFunctionPlus = (this->*wave)(r);
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

    initialized = true;
    outputSupressed = false;
    return true;
}

void VMCSolver::useLocalEnergyHelium(){
    if(nParticles != 2) {
	cout << "Cannot use this analytic local energy function "
	    << "for other than 2 particles." << endl;
	cout << "Using generic one." << endl;
	localEnergyFunction = LOCAL_ENERGY_GENERIC;
	return;
    }
    localEnergyFunction = LOCAL_ENERGY_HELIUM;
}

void VMCSolver::useLocalEnergyHydrogen(){
    if(nParticles != 1) {
	cout << "Cannot use this analytic local energy function " 
	    << "for other than 1 particle." << endl;
	cout << "Using generic one." << endl;
	localEnergyFunction = LOCAL_ENERGY_GENERIC;
	return;
    }
    localEnergyFunction = LOCAL_ENERGY_HYDROGEN;
}

void VMCSolver::useLocalEnergyGeneric(){
    localEnergyFunction = LOCAL_ENERGY_GENERIC;
}

void VMCSolver::useWaveFunction1(){
    waveFunction = WAVE_FUNCTION_1;
}

void VMCSolver::useWaveFunction2(){
    waveFunction = WAVE_FUNCTION_2;
}

void VMCSolver::reset(){
    mean = 0;
    energy = 0;
    accepts = 0;
    rejects = 0;

    rOld.reset();
    rNew.reset();
    rPlus.reset();
    rMinus.reset();

}

void VMCSolver::clearAll(){
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

    accepts = 0;
    rejects = 0;

    energy = 0;
    mean = 0;

    initialized = false;
}

double VMCSolver::getAcceptanceRatio(){
    return double(accepts)*100/(rejects+accepts);
}

double VMCSolver::wave1(double** r)
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

double VMCSolver::wave2(double** r){
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
