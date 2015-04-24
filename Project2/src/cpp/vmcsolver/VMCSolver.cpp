#include "VMCSolver.h"

using namespace CPhys;
using namespace std;
using namespace wave_functions;

VMCSolver::VMCSolver(){
    // Initialize the random generators.
    dist_uniform = std::uniform_real_distribution<double>(0.0,1.0);
    dist_uniform = std::uniform_real_distribution<double>(0.0,sqrt(2));
    clear();
}

bool VMCSolver::runIntegration(){
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
    /* cout << "Energy: " << energy << " Energy (squared sum): " */ 
	/* << energySquared << endl; */
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
            /* prNew[i][j] = prOld[i][j] + Random::gauss(idum)*sqrt(timeStep) */ 
            prNew[i][j] = prOld[i][j] + dist_gauss(gen)*sqrt(timeStep) 
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
        /* if(Random::ran2(idum) <= greensFunction*(waveFuncValNew*waveFuncValNew) */
        if(dist_uniform(gen) <= greensFunction*(waveFuncValNew*waveFuncValNew)
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

void VMCSolver::setSeed(long seed){
    gen.seed(seed);
}

void VMCSolver::runRandomStep(int cycle){
    if (efficientSlater){
        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                prNew[i][j] = prOld[i][j] + stepLength*(dist_uniform(gen) - 0.5);
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
            if(dist_uniform(gen) <= ratio) {
                for(int j = 0; j < nDimensions; j++) {
                    prOld[i][j] = prNew[i][j];
                }
                // Update the i'th particle (row) in the slater matrix.
                updateSlater(i);
                // Update the inverse of the slater matrix.
                /* if (ratio > 100) cout << ratio << endl; */
                if (i < nHalf) 
                    pMatOp::updateInverse(i, ratio, pslater1, pslater1Inv,nHalf);
                else 
                    pMatOp::updateInverse(i-nHalf, ratio,pslater2,pslater2Inv,nHalf);
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
        // ALL PARTICLES MOVED ONE STEP AT THIS POINT.
    } 
    // Not using efficient slater
    else {
        // Store the current value of the wave function
        waveFuncValOld = (this->*getWaveFuncVal)(prOld);
        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                prNew[i][j] = prOld[i][j] + stepLength*(dist_uniform(gen) - 0.5);
            }
            // Recalculate the value of the wave function
            waveFuncValNew = (this->*getWaveFuncVal)(prNew);
            double ratio = 
                    waveFuncValNew*waveFuncValNew/(waveFuncValOld*waveFuncValOld);
            // Check for step acceptance (if yes, 
            // update position, if no, reset position)
            if(dist_uniform(gen) <= ratio) {
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

    nHalf = nParticles/2;

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
    else if (localEnergyFunction == LOCAL_ENERGY_SLATER)
        getLocalEnergy = &VMCSolver::getLocalEnergySlater;
    else if (localEnergyFunction == LOCAL_ENERGY_SLATER_NOCOR)
        getLocalEnergy = &VMCSolver::getLocalEnergySlaterNoCor;
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
                /* prOld[i][j] = Random::gauss(idum)*sqrt(timeStep); */
                /* prOld[i][j] = rng.gauss(idum)*sqrt(timeStep); */
                prOld[i][j] = dist_gauss(gen)*sqrt(timeStep);
            }
        }
    }
    else {
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                prOld[i][j] = stepLength * (dist_uniform(gen) - 0.5);
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

double VMCSolver::getLocalEnergyHydrogen(double** r, int i){
    double* r1 = r[0];
    double rAbs = 0;
    for(int j = 0; j < nDimensions; j++) {
	rAbs += r1[j] * r1[j];
    }
    rAbs = sqrt(rAbs);
    return  -1/rAbs - 0.5*alpha*(alpha - 2/rAbs);
}

double VMCSolver::getLocalEnergySlater(double** r, int i){
    double DD = 0;
    double D = 0;
    double C = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nParticles/2; j++) {
            DD += phiDD(j,r[i]);
            D += phiD(j,r[i]);
        }
    }
    double CC = 0;
    double a1, a2;
    double tmp = 0;
    double rkj, rki;
    double bki, bkj;
    double dot;
    for (int k = 0; k < nParticles; k++) {
        for (int j = 0; j < nParticles; j++) {
            if (j == k) continue;
            rkj = 0;
            for (int x = 0; x < nDimensions; x++) {
                tmp = (r[j][x] - r[k][x]);
                rkj += tmp*tmp;
            }
            rkj = sqrt(rkj);
            switch (k*2 + j){
                case 0:
                    a1 = 0.25;
                case 1:
                    a1 = 0.5;
                case 2:
                    a1 = 0.5;
                case 3:
                    a1 = 0.25;
            }
            bkj = 1/(1 + beta*rkj);
            CC += 2*a1*bkj*bkj/rkj;
            CC -= 2*a1*beta*bkj*bkj*bkj;
            for (int i = 0; i < nParticles; i++) {
                if (i == k) continue;
                rki = 0;
                dot = 0;
                for (int x = 0; x < nDimensions; x++) {
                    tmp = (r[i][x] - r[k][x]);
                    rki += tmp*tmp;
                    dot += r[k][x]*r[k][x] - r[k][x]*r[j][x] 
                        - r[k][x]*r[i][x] + r[i][x]*r[j][x];
                }
                rki = sqrt(rki);
                bki = 1/(1 + beta*rki);
                switch (k*2 + i){
                    case 0:
                        a2 = 0.25;
                    case 1:
                        a2 = 0.5;
                    case 2:
                        a2 = 0.5;
                    case 3:
                        a2 = 0.25;
                }
                CC += dot/(rki*rkj)*a1*a2*bki*bki*bkj*bkj;
            }
        }
    }
    for (int k = 0; k < nParticles; k++) {
        for (int j = 0; j < nParticles; j++) {
            if (j == k) continue;
            rkj = 0;
            for (int x = 0; x < nDimensions; x++) {
                tmp = (r[j][x] - r[k][x]);
                rkj += tmp*tmp;
            }
            rkj = sqrt(rkj);
            switch (k*2 + j){
                case 0:
                    a1 = 0.25;
                case 1:
                    a1 = 0.5;
                case 2:
                    a1 = 0.5;
                case 3:
                    a1 = 0.25;
            }
            bkj = 1/(1 + beta*rkj);
            for (int x = 0; x < nDimensions; x++) {
                // What?
                C += (r[j][x] - r[k][x])*a1*bkj*bkj/rkj;
            }
        }
    }
    

    return DD + CC + 2*D*C;
}

double VMCSolver::getLocalEnergySlaterNoCor(double** r, int i){
    double sum = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nParticles/2; j++) {
            sum += phiDD(j,r[i]);
        }
    }
    return -0.5*sum;
}

double VMCSolver::getLocalEnergyHelium1(double** r, int i){
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

double VMCSolver::getLocalEnergyHelium2(double** r, int i){
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

void VMCSolver::endOfSingleParticleStep(int cycle, int i){
    // update energies
    deltaE = (this->*getLocalEnergy)(prNew, i); 
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


double VMCSolver::getLocalEnergyGenericNoCor(double** r, int i){
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

double VMCSolver::getLocalEnergyGeneric(double** r, int i){
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

Vector VMCSolver::getEnergyArray(){
    return energyArray;
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

    rMax = 0;
    bins = 1;

    mean = 0;
    energy = 0;
    energySquared = 0;
    energySum = 0;
    energySquaredSum = 0;
    rAbsSum = 0;
    deltaE = 0;
    waveFuncValOld = 0;
    waveFuncValNew = 0;


    outputSupressed = false;
    importanceSampling = false;
    efficientSlater = false;
    parallel = false;
    recordingDensity = false;
    recordingEnergyArray = false;
    recordingR12Mean = false;
    recordingPositions = false;

    // Initialize all variables, they are mostly overwritten.
    slater1 = Matrix();
    slater1Inv = Matrix();
    slater2 = Matrix();
    slater2Inv = Matrix();
    qForceOld = Matrix();
    qForceNew = Matrix();
    rOld = Matrix();
    rNew = Matrix();
    positions = Matrix();
    density = Matrix();
    vS = Vector();
    energyArray = Vector();
    pslater1 = slater1.getArrayPointer();
    pslater1Inv = slater1Inv.getArrayPointer();
    pslater2 = slater2.getArrayPointer();
    pslater2Inv = slater2Inv.getArrayPointer();
    pqForceOld = qForceOld.getArrayPointer();
    pqForceNew = qForceNew.getArrayPointer();
    prOld = rOld.getArrayPointer();
    prNew = rNew.getArrayPointer();
    pPositions = positions.getArrayPointer();
    pDensity = density.getArrayPointer();
    S = vS.getArrayPointer();
    pEnergyArray = energyArray.getArrayPointer();
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
	    if ((i == 0 || j == 0) && (i == 1 || j == 1))
		    cor += 0.25*rij/(1+beta*rij);
	    else if ((i == 2 || j == 2) && (i == 3 || j == 3))
		    cor += 0.25*rij/(1+beta*rij);
	    else
		    cor += 0.5*rij/(1+beta*rij);
    	}
    }
    cor = exp(cor);
    return phi*cor;
}

