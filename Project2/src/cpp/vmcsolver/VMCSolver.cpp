#include "VMCSolver.h"
#include <omp.h>

using namespace CPhys;
using namespace std;
using namespace wave_functions;

VMCSolver::VMCSolver(){
    clear();
}

bool VMCSolver::runIntegration(){
    if (validateParamters() && initRunVariables()) ready = true;
    if (!ready) {
        cout << "Error: Solver not initialized properly, integration not running."
            << endl;
        return false;
    }
    if (parallel) {
        initPositions();
        // Initialize different positions.
        ready = false;
        #pragma omp parallel for private(cycle)
        // Main part of the code. 
        // Loop over Monte Carlo cycles.
        for(int cycle = 0; cycle < nCycles; cycle++) {
            if (importanceSampling) runQuantumStep(cycle);
            else runRandomStep(cycle);
            endOfCycle(cycle);
        }
    }
    else {
        initPositions();
        // Main part of the code. 
        // Loop over Monte Carlo cycles.
        for(int cycle = 0; cycle < nCycles; cycle++) {
            if (importanceSampling) runQuantumStep(cycle);
            else runRandomStep(cycle);
            endOfCycle(cycle);
        }

    }

    // Calculate the density
    if (recordDensity) {
        for (int i = 0; i < nParticles; i++) {
            for (int j = 0; j < bins; j++) {
            pDensity[i][j] /= nCycles;
            }
        }
    }

    // Calculate the mean distance r12 and energy
    mean = rAbsSum/(nCycles);
    energy = energySum/(nCycles * nParticles);
    energySquared = energySquaredSum/(nCycles * nParticles);
    if (recordEnergyArray) {
    	for (int i = 0; i < nCycles; i++) {
    		pEnergyArray[i] /= nParticles;
    	}
    }

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

void VMCSolver::runQuantumStep(int cycle){
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
        } 
        else {
            for(int j = 0; j < nDimensions; j++) {
                prNew[i][j] = prOld[i][j];
                pqForceNew[i][j] = pqForceOld[i][j];
            }
            rejects++;
        }
	endOfSingleParticleStep(cycle, i);
    }
}

void VMCSolver::runRandomStep(int cycle){
    if (efficientSlater){
        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                prNew[i][j] = prOld[i][j] + stepLength*(Random::ran2(idum) - 0.5);
            }
            double ratio = 0;
            for (int j = 0; j < nParticles/2; j++) {
                if (i < nParticles/2) 
                    ratio += phi(j,prNew[i]) * pslater1Inv[j][i];
                else 
                    ratio += phi(j,prNew[i]) * pslater2Inv[j][i - nParticles/2];
            }
            // Check for step acceptance (if yes, 
            // update position, if no, reset position)
            if(Random::ran2(idum) <= ratio) {
                for(int j = 0; j < nDimensions; j++) {
                    prOld[i][j] = prNew[i][j];
                }
                // Update the i'th particle (row) in the slater matrix.
                updateSlater(i);
                // Update the inverse of the slater matrix.
                if (i < nParticles/2) 
                    updateInverse(i, ratio, pslater1, pslater1Inv);
                else 
                    updateInverse(i-nParticles/2, ratio,pslater2,pslater2Inv);
                accepts++;
            } 
            else {
                for(int j = 0; j < nDimensions; j++) {
                    prNew[i][j] = prOld[i][j];
                }
                rejects++;
            }
            endOfSingleParticleStep(cycle, i);
        }
        // All particles moved one step at this point.
    } 
    // Not using efficient slater
    else {
        // Store the current value of the wave function
        waveFuncValOld = (this->*getWaveFuncVal)(prOld);
        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                prNew[i][j] = prOld[i][j] + stepLength*(Random::ran2(idum) - 0.5);
            }
            // Recalculate the value of the wave function
            waveFuncValNew = (this->*getWaveFuncVal)(prNew);
            double ratio = 
                    waveFuncValNew*waveFuncValNew/(waveFuncValOld*waveFuncValOld);
            // Check for step acceptance (if yes, 
            // update position, if no, reset position)
            if(Random::ran2(idum) <= ratio) {
                for(int j = 0; j < nDimensions; j++) {
                    prOld[i][j] = prNew[i][j];
                }
                waveFuncValOld = waveFuncValNew;
                accepts++;
            } else {
                for(int j = 0; j < nDimensions; j++) {
                    prNew[i][j] = prOld[i][j];
                }
                waveFuncValNew = waveFuncValOld;
                rejects++;
            }
            endOfSingleParticleStep(cycle, i);
        }
        // All particles moved one step at this point.
    }
}

void VMCSolver::supressOutput(){
    outputSupressed = true;
}

bool VMCSolver::initRunVariables(){
    // Member variables
    mean = 0;
    accepts = 0;
    rejects = 0;

    energy = 0;
    energySquared = 0;

    energySum = 0;
    energySquaredSum = 0;

    deltaE = 0;
    waveFuncValOld = 0;
    waveFuncValNew = 0;

    wave_functions::alpha = alpha;
    wave_functions::beta = beta;
    wave_functions::nDimensions = nDimensions;

    // Set the wave function as a function pointer
    if (waveFunction == WAVE_FUNCTION_1)
        getWaveFuncVal = &VMCSolver::getWaveFunc1Val;
    else if (waveFunction == WAVE_FUNCTION_2)
        getWaveFuncVal = &VMCSolver::getWaveFunc2Val;
    else if (waveFunction == WAVE_FUNCTION_BERYLLIUM_1)
        getWaveFuncVal = &VMCSolver::getWaveBeryllium1Val;
    else if (waveFunction == WAVE_FUNCTION_BERYLLIUM_2)
        getWaveFuncVal = &VMCSolver::getWaveBeryllium2Val;
    else {
        cout << "Error: Wave function not set, integration not running."
            << endl;
        return false;
    }

    // Set the local energy function as a function pointer
    if (localEnergyFunction == LOCAL_ENERGY_GENERIC)
        getLocalEnergy = &VMCSolver::getLocalEnergyGeneric;
    else if (localEnergyFunction == LOCAL_ENERGY_HELIUM_1)
        getLocalEnergy = &VMCSolver::getLocalEnergyHelium1;
    else if (localEnergyFunction == LOCAL_ENERGY_HELIUM_2)
        getLocalEnergy = &VMCSolver::getLocalEnergyHelium2;
    else if (localEnergyFunction == LOCAL_ENERGY_HYDROGEN)
        getLocalEnergy = &VMCSolver::getLocalEnergyHydrogen;
    else if (localEnergyFunction == LOCAL_ENERGY_GENERIC_NOCOR)
        getLocalEnergy = &VMCSolver::getLocalEnergyGenericNoCor;
    else {
	cout << "Error: Local energy function not set, integration not running."
	    << endl;
        return false;
    }

    // Finished without error (hopefully).
    return true;
}

bool VMCSolver::initPositions(){
    // Initialize arrays
    if (recordPositions) {
        positions = Matrix(nParticles, nCycles);
        positions.reset();
        pPositions = positions.getArrayPointer();
    }
    if (recordEnergyArray) {
    	energyArray = Vector(nCycles);
        energyArray.reset();
        pEnergyArray = energyArray.getArrayPointer();
    }
    if (recordDensity) {
        density = Matrix(nParticles, bins);
        density.reset();
        pDensity = density.getArrayPointer();
    }
    if (importanceSampling) {
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
    if (importanceSampling == true) {
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
    // This is important to get the correct pointer to the new matrix. 
    prNew = rNew.getArrayPointer();
    prOld = rOld.getArrayPointer();

    // Initialize the slater determinant of the initial positions. 
    if (efficientSlater) {
        vS = Vector(nParticles/2);
        S = vS.getArrayPointer();
        if (nParticles == 2 || nParticles == 4 || nParticles == 10){
            slater1 = Matrix(nParticles/2,nParticles/2);
            slater2 = Matrix(nParticles/2,nParticles/2);
            pslater1 = slater1.getArrayPointer();
            pslater2 = slater2.getArrayPointer();
            for (int i = 0; i < nParticles/2; i++) {
                for (int j = 0; j < nParticles/2; j++) {
                    pslater1[i][j] = phi(j,prNew[i]);
                    pslater2[i][j] = phi(j,prNew[i+2]);
                }
            }
            slater1Inv = CPhys::MatOp::getInverse(slater1);
            slater2Inv = CPhys::MatOp::getInverse(slater2);
            pslater1Inv = slater1Inv.getArrayPointer();
            pslater2Inv = slater1Inv.getArrayPointer();
            /* slater1.print(); */
            /* slater1Inv.print(); */
            /* slater2.print(); */
            /* slater2Inv.print(); */
        }
        else {
            slater1 = Matrix(nParticles,nParticles);
            pslater1 = slater1.getArrayPointer();
            for (int i = 0; i < nParticles; i++) {
                for (int j = 0; j < nParticles; j++) {
                    pslater1[i][j] = phi(i,prNew[j]);
                }
            }
            slater1Inv = CPhys::MatOp::getInverse(slater1);
            pslater1Inv = slater1Inv.getArrayPointer();
        }
    }
    return true;
}

double VMCSolver::getLocalEnergyHydrogen(double** r){
    double* r1 = r[0];
    double rAbs = 0;
    for(int j = 0; j < nDimensions; j++) {
	rAbs += r1[j] * r1[j];
    }
    rAbs = sqrt(rAbs);
    return  -1/rAbs - 0.5*alpha*(alpha - 2/rAbs);
}

double VMCSolver::getLocalEnergySlater(double** r){
    double D = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nParticles/2; j++) {
            D += phiDD(j,r[i]);
        }
    }
    
    
}

double VMCSolver::getLocalEnergySlaterNoCor(double** r){
    double sum = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nParticles/2; j++) {
            sum += phiDD(j,r[i]);
        }
    }
    return -0.5*sum;
}

double VMCSolver::getLocalEnergyHelium1(double** r){
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
	// Dot product.
	temp = r1[j]*r2[j];
	r1r2 += temp;
    }
    r1Abs = sqrt(r1Abs);
    r2Abs = sqrt(r2Abs);
    r12Abs = sqrt(r12Abs);
    double E1 = (alpha-charge)*(1/r1Abs + 1/r2Abs) + 1/r12Abs - alpha*alpha;
    return E1;
}

double VMCSolver::getLocalEnergyHelium2(double** r){
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
	// Dot product.
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


void VMCSolver::endOfCycle(int cycle){
    if (!recordR12Mean) return;
    // Calculate the radius of the particle
    double rAbs = 0;
    double rsq = 0;
    for(int j = 0; j < nDimensions; j++) {
        rsq += (prNew[1][j] - prNew[0][j])*(prNew[1][j] - prNew[0][j]);
    }
    rAbs = sqrt(rsq);
    // Add it to a sum so we can calculate the mean.
    rAbsSum += rAbs;

}

void VMCSolver::endOfSingleParticleStep(int cycle, int i){
    // update energies
    deltaE = (this->*getLocalEnergy)(prNew); 
    energySum += deltaE;
    energySquaredSum += deltaE*deltaE;

    // Store in energy array.
    if (recordEnergyArray) {
    	pEnergyArray[cycle] += deltaE;
    }

    if (recordPositions) {
	double rAbs = 0;
	for(int j = 0; j < nDimensions; j++) {
	    rAbs += prNew[i][j]*prNew[i][j];
	}
	rAbs = sqrt(rAbs);
	pPositions[i][cycle] = rAbs;
    }

    // Calculate density
    if (recordDensity) {
        int bin;
        double rsq = 0;
        for(int j = 0; j < nDimensions; j++) {
            rsq += prNew[i][j]*prNew[i][j];
        }
        double rAbs = sqrt(rsq);
        if (rAbs < rMax ) {
            bin = rAbs/rMax*bins;
            pDensity[i][bin] += 1;
        }
    }
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

void VMCSolver::updateSlater(int i){
    for (int j = 0; j < nParticles/2; j++) {
        if (i < nParticles/2) pslater1[i][j] = phi(j,prNew[i]);
        else pslater2[i - nParticles/2][j] = phi(j,prNew[i]);
    }
}

void VMCSolver::updateInverse(int i, double ratio, double** mat, double** inv){
    // Update the inverse matrix for all columns except the i'th.
    for (int j = 0; j < nParticles/2; j++) {
        if (j == i) continue;
        S[j] = 0;
        for (int l = 0; l < nParticles/2; l++) {
            // d_il(new) * dInv_lj(old)
            /* S[j] += phi(l,prNew[i])*inv[l][j]; */
            S[j] += mat[i][l]*inv[l][j];
        }
    }

    for (int j = 0; j < nParticles/2; j++) {
        if (j == i) continue;
        for (int k = 0; k < nParticles/2; k++) {
            inv[k][j] = inv[k][j] - S[j]/ratio*inv[k][i];
        }
    }
    // Update the i'th column.
    for (int k = 0; k < nParticles/2; k++) {
        inv[k][i] = inv[k][i]/ratio;
    }
}

double VMCSolver::getLocalEnergyGenericNoCor(double** r){
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
    return kineticEnergy + potentialEnergy;
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

double VMCSolver::getEnergySquared(){
    return energySquared;
}

double VMCSolver::getR12Mean(){
    return mean;
}

bool VMCSolver::initFromFile(std::string fName){
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
    myFile >> paramName >> discard >> recordDensity;
    myFile >> paramName >> discard >> recordEnergyArray;
    myFile >> paramName >> discard >> recordR12Mean;
    myFile >> paramName >> discard >> recordPositions;
    myFile >> paramName >> discard >> efficientSlater;

    myFile.close();

    return true;
}

void VMCSolver::setLocalEnergyHelium1(){
    localEnergyFunction = LOCAL_ENERGY_HELIUM_1;
}

void VMCSolver::setLocalEnergyHelium2(){
    localEnergyFunction = LOCAL_ENERGY_HELIUM_2;
}

void VMCSolver::setLocalEnergyHydrogen(){
    localEnergyFunction = LOCAL_ENERGY_HYDROGEN;
}

void VMCSolver::setLocalEnergySlater(){
    localEnergyFunction = LOCAL_ENERGY_SLATER;
}

void VMCSolver::setLocalEnergySlaterNoCor(){
    localEnergyFunction = LOCAL_ENERGY_SLATER_NOCOR;
}

void VMCSolver::useEfficientSlater(bool param){
    efficientSlater = param;
}

void VMCSolver::useImportanceSampling(bool param){
    importanceSampling = param;
}

void VMCSolver::useParallel(bool param){
    parallel = param;
}

void VMCSolver::setRecordEnergyArray(bool param){
    recordEnergyArray = param;
}

void VMCSolver::setRecordDensity(bool param, int bins, double maxPos){
    // This is the only place where bins and rMax are set. But 
    // this function is called on clear().
    recordDensity = param;
    this->bins = bins;
    rMax = maxPos;
}

void VMCSolver::setRecordR12Mean(bool param){
    recordR12Mean = param;
}
void VMCSolver::setRecordPositions(bool param){
    recordPositions = param;
}

void VMCSolver::setLocalEnergyGeneric(){
    localEnergyFunction = LOCAL_ENERGY_GENERIC;
}

void VMCSolver::setLocalEnergyGenericNoCor(){
    localEnergyFunction = LOCAL_ENERGY_GENERIC_NOCOR;
}

void VMCSolver::setWaveFunction1(){
    waveFunction = WAVE_FUNCTION_1;
}

void VMCSolver::setWaveFunction2(){
    waveFunction = WAVE_FUNCTION_2;
}

void VMCSolver::setWaveFunctionBeryllium1(){
    waveFunction = WAVE_FUNCTION_BERYLLIUM_1;
}

void VMCSolver::setWaveFunctionBeryllium2(){
    waveFunction = WAVE_FUNCTION_BERYLLIUM_2;
}

void VMCSolver::clear(){
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

double VMCSolver::getAcceptanceRatio(){
    return double(accepts)*100/(rejects+accepts);
}

double VMCSolver::getWaveFunc1Val(double** r){
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

double VMCSolver::getWaveBeryllium1Val(double** r){
    double* r1 = r[0];
    double* r2 = r[1];
    double* r3 = r[2];
    double* r4 = r[3];
    double r1Abs = 0;
    double r2Abs = 0;
    double r3Abs = 0;
    double r4Abs = 0;
    for(int j = 0; j < nDimensions; j++) {
        r1Abs += r1[j] * r1[j];
        r2Abs += r2[j] * r2[j];
        r3Abs += r3[j] * r3[j];
        r4Abs += r4[j] * r4[j];
    }
    r1Abs = sqrt(r1Abs);
    r2Abs = sqrt(r2Abs);
    r3Abs = sqrt(r3Abs);
    r4Abs = sqrt(r4Abs);
    return (phi1s(r1Abs)*phi2s(r2Abs) - phi1s(r2Abs)*phi2s(r1Abs))
	*(phi1s(r3Abs)*phi2s(r4Abs) - phi1s(r4Abs)*phi2s(r3Abs));
}

double VMCSolver::getWaveBeryllium2Val(double** r){
    double* r1 = r[0];
    double* r2 = r[1];
    double* r3 = r[2];
    double* r4 = r[3];
    double r1Abs = 0;
    double r2Abs = 0;
    double r3Abs = 0;
    double r4Abs = 0;
    for(int j = 0; j < nDimensions; j++) {
        r1Abs += r1[j] * r1[j];
        r2Abs += r2[j] * r2[j];
        r3Abs += r3[j] * r3[j];
        r4Abs += r4[j] * r4[j];
    }
    r1Abs = sqrt(r1Abs);
    r2Abs = sqrt(r2Abs);
    r3Abs = sqrt(r3Abs);
    r4Abs = sqrt(r4Abs);
    // The value of the slater determinant.
    double phi = (phi1s(r1Abs)*phi2s(r2Abs) -phi1s(r2Abs)*phi2s(r1Abs))
	*(phi1s(r3Abs)*phi2s(r4Abs) -phi1s(r4Abs)*phi2s(r3Abs));
    double cor = 0;
    double rij = 0;
    for (int i = 0; i < 4; i++) {
    	for (int j = 0; j < 4; j++) {
    	    if (i >= j) continue; 
	    // If i < j, calculate something.
	    rij = 0;
	    for (int k = 0; k < nDimensions; k++) {
	    	rij += (r[j][k] - r[i][k])*(r[j][k] - r[i][k]);
	    }
	    rij = sqrt(rij);
	    if ((i == 0 || i == 0) && (j == 1 || j == 1))
		cor += 0.25*rij/(1+beta*rij);
	    else if ((i == 2 || i == 2) && (j == 3 || j == 3))
		cor += 0.25*rij/(1+beta*rij);
	    else
		cor += 0.5*rij/(1+beta*rij);
    	}
    }
    cor = exp(cor);
    return phi*cor;
    
}

/* double VMCSolver::phi(int j, double* r){ */
/*     double rAbs = 0; */
/*     for (int i = 0; i < nDimensions; i++) { */
/*         rAbs += r[i]*r[i]; */
/*     } */
/*     rAbs = sqrt(rAbs); */
/*     switch (j) { */
/*         case 0 : */
/*             return phi1s(rAbs); */
/*         case 1 : */
/*             return phi2s(rAbs); */
/*         case 2 ... 4: */
/*             return phi2p(rAbs); */
/*         default: */
/*             cout << "Index out of bounds in phi()!!!" << endl; */
/*             return 0; */
/*     } */
/* } */

/* double VMCSolver::phi1s(double r){ */
/*     return exp(-alpha*r); */
/* } */

/* double VMCSolver::phi2s(double r){ */
/*     return (1-alpha*r/2)*exp(-alpha*r/2); */
/* } */

/* double VMCSolver::phi2p(double r){ */
/*     return alpha*r*exp(-alpha*r/2); */
/* } */

bool VMCSolver::validateParamters(){
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
    if (recordR12Mean && nParticles != 2) {
        cout << "Cannot use record r12 mean   "
            << "for other than 2 particles." << endl;
        valid = false;
    }

    return valid;
}

void VMCSolver::exportParamters(std::string fName){
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
    myFile << "recordDensity = " << recordDensity <<  endl;
    myFile << "recordEnergyArray = " << recordEnergyArray <<  endl;
    myFile << "recordR12Mean = " << recordR12Mean <<  endl;
    myFile << "recordPositions = " << recordPositions <<  endl;
    myFile << "useEfficientSlater = " << efficientSlater <<  endl;
    myFile.close();
}

void VMCSolver::exportDensity(std::string fName){
    string adress = "../../../res/" + fName;
    ofstream myFile;
    cout << "Dumption densities to file : " << fName << endl;
    myFile.open(adress.c_str());
    myFile << rMax << endl;
    for (int i = 0; i < nParticles; i++) {
	for (int j = 0; j < bins; j++) {
	    myFile << pDensity[i][j] << " ";
	}
	myFile << endl;
    }
    myFile.close();

}

void VMCSolver::exportEnergyArray(std::string fName){
    string adress = "../../../res/" + fName;
    ofstream myFile;
    cout << "Dumption energies to file : " << fName << endl;
    myFile.open(adress.c_str());
    for (int i = 0; i < nCycles; i++) {
	myFile << pEnergyArray[i] << " ";
    }
    myFile.close();

}

void VMCSolver::exportPositions(std::string fName){
    string adress = "../../../res/" + fName;
    ofstream myFile;
    cout << "Dumption energies to file : " << fName << endl;
    myFile.open(adress.c_str());
    for (int i = 0; i < nParticles; i++) {
    	for (int j = 0; j < nCycles; j++) {
	    myFile << pPositions[i][j] << " ";
    	}
	myFile << endl;
    }
    myFile.close();

}
