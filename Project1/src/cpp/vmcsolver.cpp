#include "vmcsolver.h"

using namespace CPhys;
using namespace std;

VMCSolver::VMCSolver(){
    clearAll();
}


bool VMCSolver::runMonteCarloIntegration(){
    if (!initialized) {
	cout << "Error: Solver not initialized, integration not running."
	    << endl;
	return false;
    }
    rOld = Matrix(nParticles, nDimensions);
    rNew = Matrix(nParticles, nDimensions);
    prNew = rNew.getArrayPointer();
    prOld = rOld.getArrayPointer();

    rPlus = Matrix(nParticles, nDimensions);
    rMinus = Matrix(nParticles, nDimensions);

    double waveFunctionOld = 0;
    double waveFunctionNew = 0;

    double energySum = 0;
    double energySquaredSum = 0;
    
    

    reset();

    double deltaE;

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            prOld[i][j] = stepLength * (Random::ran2(idum) - 0.5);
        }
    }

    rNew = rOld;
    prNew = rNew.getArrayPointer();
    prOld = rOld.getArrayPointer();

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++) {

        // Store the current value of the wave function
	if (waveFunction == WAVE_FUNCTION_1) 
	    waveFunctionOld = wave1(prOld);
	else 
	    waveFunctionOld = wave2(prOld[0],prOld[1]);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                prNew[i][j] = prOld[i][j] + stepLength*(Random::ran2(idum) - 0.5);
            }

            // Recalculate the value of the wave function
            if(waveFunction == WAVE_FUNCTION_1) 
		waveFunctionNew = wave1(prNew);
	    else 
		waveFunctionNew = wave2(prNew[0],prNew[1]);

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
	    if (localEnergyFunction == LOCAL_ENERGY_GENERIC) 
		deltaE = localEnergy(rNew); 
	    else if(localEnergyFunction == LOCAL_ENERGY_HELIUM) 
		deltaE = localEnergyHelium(prNew[0], prNew[1]);
	    else if(localEnergyFunction == LOCAL_ENERGY_HYDROGEN) 
		/* deltaE = localEnergyHydrogen(prNew[0]); */
		cout << "Should go here" << endl;
            
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

    if (outputSupressed) return true;
    cout << "Energy: " << energy << " Energy (squared sum): " 
	<< energySquared << endl;
    cout << "Variance : " << energySquared - energy*energy << endl;
    return true;
}

void VMCSolver::supressOutput(){
    outputSupressed = true;
}

double VMCSolver::localEnergyHelium(double* r1, double* r2){
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

double VMCSolver::localEnergy(Matrix &r)
{

    double** pr = r.getArrayPointer();

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;
    double waveFunctionCurrent;
    if (waveFunction == WAVE_FUNCTION_1) 
	waveFunctionCurrent = wave1(pr);
    else 
	waveFunctionCurrent = wave2(pr[0],pr[1]);

    // Kinetic energy

    double rPlus;
    double rMinus;
    double r0;
    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
	    rPlus 	= pr[i][j] + h;
	    rMinus 	= pr[i][j] - h;
	    r0 		= pr[i][j];

	    pr[i][j] = rMinus;
	    if (waveFunction == WAVE_FUNCTION_1) 
		waveFunctionMinus = wave1(pr);
	    else 
		waveFunctionMinus = wave2(pr[0],pr[1]);
	    pr[i][j] = rPlus;
	    if (waveFunction == WAVE_FUNCTION_1) 
		waveFunctionPlus = wave1(pr);
	    else 
		waveFunctionPlus = wave2(pr[0],pr[1]);
	    pr[i][j] = r0;
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
            rSingleParticle += pr[i][j]*pr[i][j];
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (pr[i][k] - pr[j][k]) * (pr[i][k] - pr[j][k]);
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
    return true;
}

void VMCSolver::useLocalEnergyGeneric(){
    localEnergyFunction = LOCAL_ENERGY_GENERIC;
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
    outputSupressed = false;

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

    rSum = 0;
    energy = 0;
    mean = 0;

    initialized = false;
}

double VMCSolver::getStepAcceptance(){
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

double VMCSolver::wave2(double* r1, double* r2){
    double temp = 0;
    double r = 0;
    double r12 = 0;
    for(int j = 0; j < nDimensions; j++) {
	temp = r1[j] * r1[j];
	r += temp;
	temp = r2[j] * r2[j];
	r += temp;

	temp = (r1[j] - r2[j]) * (r1[j] - r2[j]);
	r12 += temp;
    }
    r = sqrt(r);
    r12 = sqrt(r12);
    return exp(-r*alpha)*exp(r12/(2*(1+beta*r12)));
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
