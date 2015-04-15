#include "VMCSolver1.h"

using namespace CPhys;
using namespace std;
using namespace wave_functions;

VMCSolver1::VMCSolver1(){
    clear();
}

bool VMCSolver1::runIntegration(){
    bool ready = false;
    if (initRunVariables()) ready = true;
    if (!ready) {
        cout << "Error: Solver not initialized properly, integration not running."
            << endl;
        return false;
    }
    // Main part of the code. 
    // Loop over Monte Carlo cycles.
    for(int cycle = 0; cycle < nCycles; cycle++) {
        if (importanceSampling) runQuantumStep(cycle);
        else runRandomStep(cycle);
        endOfCycle(cycle);
    }
    // Calculate the density
    if (recordingDensity) {
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
    if (recordingEnergyArray) {
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

void VMCSolver1::runQuantumStep(int cycle){
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

void VMCSolver1::runRandomStep(int cycle){
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

void VMCSolver1::supressOutput(){
    outputSupressed = true;
}

bool VMCSolver1::initRunVariables(){
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
        getWaveFuncVal = &VMCSolver1::getWaveFunc1Val;
    else if (waveFunction == WAVE_FUNCTION_2)
        getWaveFuncVal = &VMCSolver1::getWaveFunc2Val;
    else if (waveFunction == WAVE_FUNCTION_BERYLLIUM_1)
        getWaveFuncVal = &VMCSolver1::getWaveBeryllium1Val;
    else if (waveFunction == WAVE_FUNCTION_BERYLLIUM_2)
        getWaveFuncVal = &VMCSolver1::getWaveBeryllium2Val;
    else {
        cout << "Error: Wave function not set, integration not running."
            << endl;
        return false;
    }

    cout << "Local : " << localEnergyFunction << endl;
    // Set the local energy function as a function pointer
    if (localEnergyFunction == LOCAL_ENERGY_GENERIC)
        getLocalEnergy = &VMCSolver1::getLocalEnergyGeneric;
    else if (localEnergyFunction == LOCAL_ENERGY_HELIUM_1)
        getLocalEnergy = &VMCSolver1::getLocalEnergyHelium1;
    else if (localEnergyFunction == LOCAL_ENERGY_HELIUM_2)
        getLocalEnergy = &VMCSolver1::getLocalEnergyHelium2;
    else if (localEnergyFunction == LOCAL_ENERGY_HYDROGEN)
        getLocalEnergy = &VMCSolver1::getLocalEnergyHydrogen;
    else if (localEnergyFunction == LOCAL_ENERGY_GENERIC_NOCOR)
        getLocalEnergy = &VMCSolver1::getLocalEnergyGenericNoCor;
    else if (localEnergyFunction == LOCAL_ENERGY_SLATER)
        getLocalEnergy = &VMCSolver1::getLocalEnergySlater;
    else if (localEnergyFunction == LOCAL_ENERGY_SLATER_NOCOR)
        getLocalEnergy = &VMCSolver1::getLocalEnergySlaterNoCor;
    else {
	cout << "Error: Local energy function not set, integration not running."
	    << endl;
        return false;
    }

    // Initialize arrays
    if (recordingPositions) {
        positions = Matrix(nParticles, nCycles);
        positions.reset();
        pPositions = positions.getArrayPointer();
    }
    if (recordingEnergyArray) {
    	energyArray = Vector(nCycles);
        energyArray.reset();
        pEnergyArray = energyArray.getArrayPointer();
    }
    if (recordingDensity) {
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
    // Finished without error (hopefully).
    return true;
}

double VMCSolver1::getLocalEnergyHydrogen(double** r){
    double* r1 = r[0];
    double rAbs = 0;
    for(int j = 0; j < nDimensions; j++) {
	rAbs += r1[j] * r1[j];
    }
    rAbs = sqrt(rAbs);
    return  -1/rAbs - 0.5*alpha*(alpha - 2/rAbs);
}

double VMCSolver1::getLocalEnergySlater(double** r){
    double D = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nParticles/2; j++) {
            D += phiDD(j,r[i]);
        }
    }
    
    
}

double VMCSolver1::getLocalEnergySlaterNoCor(double** r){
    double sum = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nParticles/2; j++) {
            sum += phiDD(j,r[i]);
        }
    }
    return -0.5*sum;
}

double VMCSolver1::getLocalEnergyHelium1(double** r){
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

double VMCSolver1::getLocalEnergyHelium2(double** r){
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


void VMCSolver1::endOfCycle(int cycle){
    if (!recordingR12Mean) return;
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

void VMCSolver1::endOfSingleParticleStep(int cycle, int i){
    // update energies
    deltaE = (this->*getLocalEnergy)(prNew); 
    energySum += deltaE;
    energySquaredSum += deltaE*deltaE;

    // Store in energy array.
    if (recordingEnergyArray) {
    	pEnergyArray[cycle] += deltaE;
    }

    if (recordingPositions) {
	double rAbs = 0;
	for(int j = 0; j < nDimensions; j++) {
	    rAbs += prNew[i][j]*prNew[i][j];
	}
	rAbs = sqrt(rAbs);
	pPositions[i][cycle] = rAbs;
    }

    // Calculate density
    if (recordingDensity) {
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

void VMCSolver1::updateQuantumForce(double** r, double ** qForce, double factor){
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

void VMCSolver1::updateSlater(int i){
    for (int j = 0; j < nParticles/2; j++) {
        if (i < nParticles/2) pslater1[i][j] = phi(j,prNew[i]);
        else pslater2[i - nParticles/2][j] = phi(j,prNew[i]);
    }
}

void VMCSolver1::updateInverse(int i, double ratio, double** mat, double** inv){
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

double VMCSolver1::getLocalEnergyGenericNoCor(double** r){
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

double VMCSolver1::getLocalEnergyGeneric(double** r){
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

double VMCSolver1::getStepLength(){
    return stepLength;
}

double VMCSolver1::getEnergy(){
    return energy;
}

double VMCSolver1::getEnergySquared(){
    return energySquared;
}

double VMCSolver1::getR12Mean(){
    return mean;
}

void VMCSolver1::clear(){
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
    importanceSampling = false;
    efficientSlater = false;
    parallel = false;
    recordingDensity = false;
    recordingEnergyArray = false;
    recordingR12Mean = false;
    recordingPositions = false;
}

double VMCSolver1::getAcceptanceRatio(){
    return double(accepts)*100/(rejects+accepts);
}

double VMCSolver1::getWaveFunc1Val(double** r){
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

double VMCSolver1::getWaveFunc2Val(double** r){
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

double VMCSolver1::getWaveBeryllium1Val(double** r){
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

double VMCSolver1::getWaveBeryllium2Val(double** r){
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

