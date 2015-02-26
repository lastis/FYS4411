#include "vmcsolver.h"

using namespace CPhys;
using namespace std;

VMCSolver::VMCSolver() :
    waveFunctionType(1),
    accepts(0),
    rejects(0),
    nDimensions(3),
    charge(2),
    stepLength(1.0),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(-1),
    alpha(0.5*charge),
    beta(1),
    nCycles(1000000)
{
}

bool VMCSolver::initFromFile(std::string fName){
    ofstream myFile;
    myFile.open(fName.c_str());
    myFile.close();
    return true;
}

void VMCSolver::useWaveType1(){
    waveFunctionType = 1;
}

void VMCSolver::useWaveType2(){
    waveFunctionType = 2;
}

void VMCSolver::reset(){
    accepts = 0;
    rejects = 0;
}

void VMCSolver::runMonteCarloIntegration()
{
    reset();
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
	else waveFunctionOld = waveFunction2(prOld[0],prOld[1]);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                prNew[i][j] = prOld[i][j] + stepLength*(Random::ran2(idum) - 0.5);
            }

            // Recalculate the value of the wave function
            if(waveFunctionType == 1) waveFunctionNew = waveFunction1(prNew);
	    else waveFunctionNew = waveFunction2(prNew[0], prOld[1]);

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
    }
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    cout << "Energy: " << energy << " Energy (squared sum): " 
	<< energySquared << endl;
    cout << "Accepted moves: " << int(double(accepts)*100/(rejects+accepts))
	<< " %" << endl;
    reset();
}

double VMCSolver::localEnergy(Matrix &r)
{

    double** pr = r.getArrayPointer();

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;
    double waveFunctionCurrent;
    if (waveFunctionType == 1) waveFunctionCurrent = waveFunction1(pr);
    else waveFunctionCurrent = waveFunction2(pr[0], pr[1]);

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
	    else waveFunctionMinus = waveFunction2(pr[0], pr[1]);
	    pr[i][j] = rPlus;
	    if (waveFunctionType == 1) waveFunctionPlus = waveFunction1(pr);
	    else waveFunctionPlus = waveFunction2(pr[0], pr[1]);
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

double VMCSolver::waveFunction2(double* r1, double* r2){
    double r1sq = 0;
    double r2sq = 0;
    double r12sq = 0;
    for(int j = 0; j < nDimensions; j++) {
	r1sq += r1[j] * r1[j];
	r2sq += r2[j] * r2[j];
	r12sq  += (r2[j] - r1[j])*(r2[j] - r1[j]);
    }
    double argument = sqrt(r1sq) + sqrt(r2sq);
    double r12 = sqrt(r12sq);
    return exp(-argument*alpha*r12/(2*(1+beta*r12)));
}
