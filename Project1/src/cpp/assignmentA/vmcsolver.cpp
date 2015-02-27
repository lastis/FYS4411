#include "vmcsolver.h"

using namespace CPhys;
using namespace std;

VMCSolver::VMCSolver(){
    reset();
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
    
    
    // TODO put these "resets" in an own method.
    rSum = 0;
    energy = 0;

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
	if (waveFunctionType == 1) waveFunctionOld = waveFunction1(prOld);
	else waveFunctionOld = waveFunction2(prOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                prNew[i][j] = prOld[i][j] + stepLength*(Random::ran2(idum) - 0.5);
            }

            // Recalculate the value of the wave function
            if(waveFunctionType == 1) waveFunctionNew = waveFunction1(prNew);
	    else waveFunctionNew = waveFunction2(prNew);

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
            deltaE = localEnergy(rNew);
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
    cout << "Energy: " << energy << " Energy (squared sum): " 
	<< energySquared << endl;
    return true;
}

double VMCSolver::localEnergy(Matrix &r)
{

    double** pr = r.getArrayPointer();

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;
    double waveFunctionCurrent;
    if (waveFunctionType == 1) waveFunctionCurrent = waveFunction1(pr);
    else waveFunctionCurrent = waveFunction2(pr);

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
	    if (waveFunctionType == 1) waveFunctionMinus = waveFunction1(pr);
	    else waveFunctionMinus = waveFunction2(pr);
	    pr[i][j] = rPlus;
	    if (waveFunctionType == 1) waveFunctionPlus = waveFunction1(pr);
	    else waveFunctionPlus = waveFunction2(pr);
	    pr[i][j] = r0;
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
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
    string line;

    myFile.open(fName.c_str());
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
    myFile >> paramName >> discard >> waveFunctionType;
    myFile >> paramName >> discard >> h;
    myFile >> paramName >> discard >> h2;
    myFile >> paramName >> discard >> idum;

    myFile.close();

    initialized = true;
    return true;
}

void VMCSolver::useWaveType1(){
    waveFunctionType = 1;
}

void VMCSolver::useWaveType2(){
    waveFunctionType = 2;
}

void VMCSolver::reset(){
    waveFunctionType = 1;
    nDimensions = 0;
    charge = 0;
    stepLength = 0;
    nParticles = 0;
    h = 0;
    h2 = 0;
    idum = -1;
    alpha = 0;
    beta = 0;
    nCycles = 0;
    accepts = 0;
    rejects = 0;

    rSum = 0;

    initialized = false;
}

double VMCSolver::getStepAcceptance(){
    return double(accepts)*100/(rejects+accepts);
}

double VMCSolver::waveFunction1(double** r)
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

double VMCSolver::waveFunction2(double** r){
    double argument = 0;
    double rsq = 0;
    double r12sq = 0;
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
	// Alpha part
        for(int j = 0; j < nDimensions; j++) {
            rsq += r[i][j] * r[i][j];
        }
	argument += sqrt(rsq);
	// Beta part
        for(int j = i + 1; j < nParticles; j++) {
            r12sq = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12sq += (r[i][k] - r[j][k]) * (r[i][k] - r[j][k]);
            }
	    r12 += sqrt(r12sq);
        }
    }
    return exp(-argument*alpha*r12/(2*(1+beta*r12)));
    /* double r1sq = 0; */
    /* double r2sq = 0; */
    /* double r12sq = 0; */
    /* for(int j = 0; j < nDimensions; j++) { */
	/* r1sq += r1[j] * r1[j]; */
	/* r2sq += r2[j] * r2[j]; */
	/* r12sq  += (r2[j] - r1[j])*(r2[j] - r1[j]); */
    /* } */
    /* double argument = sqrt(r1sq) + sqrt(r2sq); */
    /* double r12 = sqrt(r12sq); */
    /* return exp(-argument*alpha*r12/(2*(1+beta*r12))); */
}

void VMCSolver::exportParamters(std::string fName){
    ofstream myFile;
    cout << "Dumption to file : " << fName << endl;
    myFile.open(fName.c_str());
    myFile << "charge = " << charge <<  endl;
    myFile << "alpha = " << alpha << endl;
    myFile << "beta = " << beta << endl;
    myFile << "nDimensions = " << nDimensions <<  endl;
    myFile << "nParticles = " << nParticles <<  endl;
    myFile << "stepLength = " << stepLength << endl;
    myFile << "nCycles = " << nCycles << endl;
    myFile << "waveFunctionType = " << waveFunctionType << endl;
    myFile << "h = " << h << endl;
    myFile << "h2 = " << h2 << endl;
    myFile << "idum = " << idum << endl;
    myFile.close();
}
