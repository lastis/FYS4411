#include "VMCSolver.h"

using namespace CPhys;
using namespace std;
using namespace wave_functions;

VMCSolver::VMCSolver()
{
    // Initialize the random generators.
    dist_uniform = std::uniform_real_distribution<double>(0.0, 1.0);
    dist_gauss = std::normal_distribution<double>(0.0, sqrt(2));
    clear();
}

bool VMCSolver::runIntegration()
{
    bool ready = false;
    if (initRunVariables()) ready = true;
    if (!ready)
    {
        cout << "Error: Solver not initialized properly, integration not "
                "running." << endl;
        return false;
    }
    // Main part of the code.
    // Loop over Monte Carlo cycles.
    // This implementation might look stupid, but everything breaks all
    // the time and every tiny step requires rigorous testing.
    if (importanceSampling == true && efficientSlater == false)
    {
        for (int cycle = 0; cycle < nCycles; cycle++)
        {
            runStepQuantum(cycle);
        }
    }
    else if (importanceSampling == true && efficientSlater == true)
    {
        for (int cycle = 0; cycle < nCycles; cycle++)
        {
            runStepSlaterQuantum(cycle);
        }
    }
    else if (importanceSampling == false && efficientSlater == true)
    {
        for (int cycle = 0; cycle < nCycles; cycle++)
        {
            runStepSlater(cycle);
        }
    }
    else if (importanceSampling == false && efficientSlater == false)
    {
        for (int cycle = 0; cycle < nCycles; cycle++)
        {
            runStep(cycle);
        }
    }

    // Calculate the density
    if (recordingDensity)
    {
        for (int i = 0; i < nParticles; i++)
        {
            for (int j = 0; j < bins; j++)
            {
                pDensity[i][j] /= nCycles;
            }
        }
    }

    // Calculate the mean distance r12 and energy
    mean = rAbsSum / (nCycles);
    energy = energySum / (nCycles * nParticles);
    energySquared = energySquaredSum / (nCycles * nParticles);
    if (recordingEnergyArray)
    {
        for (int i = 0; i < nCycles; i++)
        {
            pEnergyArray[i] /= nParticles;
        }
    }

    // Output
    if (outputSupressed)
    {
        outputSupressed = false;
        return true;
    }
    /* cout << "Energy: " << energy << " Energy (squared sum): " */
    /* << energySquared << endl; */
    return true;
}

void VMCSolver::runStepSlaterQuantum(int cycle)
{
    startOfCycleSlaterQuantum();
    for (int i = 0; i < nParticles; i++)
    {
        runSingleStepSlaterQuantum(i, cycle);
    }
    // ALL PARTICLES MOVED ONE STEP AT THIS POINT.
    endOfCycle(cycle);
}

void VMCSolver::runSingleStepQuantum(int i, int cycle)
{
    rAbsNew[i] = 0;
    for (int j = 0; j < nDimensions; j++)
    {
        prNew[i][j] = prOld[i][j] + dist_gauss(gen) * sqrt(timeStep) +
                      pqForceOld[i][j] * timeStep * D;
        rAbsNew[i] += prNew[i][j] * prNew[i][j];
    }
    rAbsNew[i] = sqrt(rAbsNew[i]);

    // Recalculate the value of the wave function
    waveFuncValNew = getWaveFuncVal(prNew, rAbsNew);
    updateQuantumForce(prNew, rAbsNew, pqForceNew, waveFuncValNew);

    // Compute the log ratio of the greens functions to be used in the
    // Metropolis-Hastings algorithm.
    greensFunction = 0;
    for (int j = 0; j < nDimensions; j++)
    {
        greensFunction +=
            0.5 * (pqForceOld[i][j] + pqForceNew[i][j]) *
            (D * timeStep * 0.5 * (pqForceOld[i][j] - pqForceNew[i][j]) -
             prNew[i][j] + prOld[i][j]);
    }
    greensFunction = exp(greensFunction);
    ratio = greensFunction * (waveFuncValNew * waveFuncValNew) /
            (waveFuncValOld * waveFuncValOld);

    // Check for step acceptance (if yes,
    // update position, if no, reset position)
    // The metropolis test is performed by moving one particle at the time.
    if (dist_uniform(gen) <= ratio)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            prOld[i][j] = prNew[i][j];
            pqForceOld[i][j] = pqForceNew[i][j];
        }
        rAbsOld[i] = rAbsNew[i];
        waveFuncValOld = waveFuncValNew;
        accepts++;
    }
    else
    {
        for (int j = 0; j < nDimensions; j++)
        {
            prNew[i][j] = prOld[i][j];
            pqForceNew[i][j] = pqForceOld[i][j];
        }
        rAbsNew[i] = rAbsOld[i];
        rejects++;
    }
    // update energies
    endOfSingleParticleStep(cycle, i);
}

void VMCSolver::runStepQuantum(int cycle)
{
    startOfCycleQuantum();
    // New position to test
    for (int i = 0; i < nParticles; i++)
    {
        runSingleStepQuantum(i, cycle);
    }
    endOfCycle(cycle);
}

void VMCSolver::setSeed(long seed)
{
    gen.seed(seed);
}

double VMCSolver::getCorrelationRatio(int i)
{
    double a;
    double rkjNew = 0;
    double rkjOld = 0;
    double valNew = 0;
    double valOld = 0;
    double tmp;
    for (int k = 0; k < nParticles; k++)
    {
        for (int j = k + 1; j < nParticles; j++)
        {
            if (i == j || i == k)
            {
                rkjOld = 0;
                rkjNew = 0;
                for (int x = 0; x < nDimensions; x++)
                {
                    tmp = (prNew[j][x] - prNew[k][x]);
                    rkjNew += tmp * tmp;
                    tmp = (prOld[j][x] - prOld[k][x]);
                    rkjOld += tmp * tmp;
                }
                rkjOld = sqrt(rkjOld);
                rkjNew = sqrt(rkjNew);
                int spinK = k / nHalf;
                int spinJ = j / nHalf;
                switch (spinK + spinJ)
                {
                    case 0:
                        a = 0.25;
                        /* a = 0.5; */
                        break;
                    case 1:
                        a = 0.5;
                        break;
                    case 2:
                        a = 0.25;
                        /* a = 0.5; */
                        break;
                }
                double bkjOld = 1 / (1 + beta * rkjOld);
                double bkjNew = 1 / (1 + beta * rkjNew);
                valOld += a * bkjOld * rkjOld;
                valNew += a * bkjNew * rkjNew;
            }
        }
    }
    return exp(valNew - valOld);
}

void VMCSolver::startOfCycle()
{
    // Store the current value of the wave function
    waveFuncValOld = wave_functions::getWaveFuncVal(prOld, rAbsOld);
}


void VMCSolver::startOfCycleQuantum()
{
    // Store the current value of the wave function
    waveFuncValOld = getWaveFuncVal(prOld, rAbsOld);
    updateQuantumForce(prOld, rAbsOld, pqForceOld, waveFuncValOld);
}

void VMCSolver::startOfCycleSlaterQuantum()
{
    updateQuantumForceSlater(prOld, rAbsOld, pqForceOld);
}

void VMCSolver::runSingleStepSlaterQuantum(int i, int cycle)
{
    // New position to test
    rAbsNew[i] = 0;
    for (int j = 0; j < nDimensions; j++)
    {
        prNew[i][j] = prOld[i][j] + dist_gauss(gen) * sqrt(timeStep) +
                      pqForceOld[i][j] * timeStep * D;
        rAbsNew[i] += prNew[i][j] * prNew[i][j];
    }
    rAbsNew[i] = sqrt(rAbsNew[i]);
    double ratioTmp = 0;
    for (int j = 0; j < nHalf; j++)
    {
        if (i < nHalf)
            ratioTmp += phi(j, prNew[i], rAbsNew[i]) * pslater1Inv[j][i];
        else
            ratioTmp +=
                phi(j, prNew[i], rAbsNew[i]) * pslater2Inv[j][i - nHalf];
    }
    if (usingCorrelation)
    {
        ratio = ratioTmp * getCorrelationRatio(i);
        ratio = ratio * ratio;
    }
    else
        ratio = ratioTmp * ratioTmp;

    // Update the quantum force.
    updateQuantumForceSlater(prNew, rAbsNew, pqForceNew);
    // Compute the log ratio of the greens functions to be used in the
    // Metropolis-Hastings algorithm.
    greensFunction = 0;
    for (int j = 0; j < nDimensions; j++)
    {
        greensFunction +=
            0.5 * (pqForceOld[i][j] + pqForceNew[i][j]) *
            (D * timeStep * 0.5 * (pqForceOld[i][j] - pqForceNew[i][j]) -
             prNew[i][j] + prOld[i][j]);
    }
    greensFunction = exp(greensFunction);

    // Check for step acceptance (if yes,
    // update position, if no, reset position)
    // The metropolis test is performed by moving one particle at the time.
    if (dist_uniform(gen) <= greensFunction * ratio)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            prOld[i][j] = prNew[i][j];
            pqForceOld[i][j] = pqForceNew[i][j];
        }
        rAbsOld[i] = rAbsNew[i];
        // Update the i'th particle (row) in the slater matrix.
        updateSlater(i);
        // Update the inverse of the slater matrix.
        if (i < nHalf)
            pMatOp::updateInverse(i, ratioTmp, pslater1, pslater1Inv, nHalf);
        else
            pMatOp::updateInverse(i - nHalf, ratioTmp, pslater2, pslater2Inv,
                                  nHalf);
        accepts++;
    }
    else
    {
        for (int j = 0; j < nDimensions; j++)
        {
            prNew[i][j] = prOld[i][j];
            pqForceNew[i][j] = pqForceOld[i][j];
        }
        rAbsNew[i] = rAbsOld[i];
        rejects++;
    }
    endOfSingleParticleStep(cycle, i);
}

void VMCSolver::runSingleStepSlater(int i, int cycle)
{
    // New position to test
    rAbsNew[i] = 0;
    for (int j = 0; j < nDimensions; j++)
    {
        prNew[i][j] = prOld[i][j] + stepLength * (dist_uniform(gen) - 0.5);
        rAbsNew[i] += prNew[i][j] * prNew[i][j];
    }
    rAbsNew[i] = sqrt(rAbsNew[i]);
    double ratioTmp = 0;
    for (int j = 0; j < nHalf; j++)
    {
        if (i < nHalf)
            ratioTmp += phi(j, prNew[i], rAbsNew[i]) * pslater1Inv[j][i];
        else
            ratioTmp +=
                phi(j, prNew[i], rAbsNew[i]) * pslater2Inv[j][i - nHalf];
    }
    if (usingCorrelation)
    {
        ratio = ratioTmp * getCorrelationRatio(i);
        ratio = ratio * ratio;
    }
    else
        ratio = ratioTmp * ratioTmp;
    // Check for step acceptance (if yes,
    // update position, if no, reset position)
    if (dist_uniform(gen) <= ratio)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            prOld[i][j] = prNew[i][j];
        }
        rAbsOld[i] = rAbsNew[i];
        // Update the i'th particle (row) in the slater matrix.
        updateSlater(i);
        // Update the inverse of the slater matrix.
        if (i < nHalf)
            pMatOp::updateInverse(i, ratioTmp, pslater1, pslater1Inv, nHalf);
        else
            pMatOp::updateInverse(i - nHalf, ratioTmp, pslater2, pslater2Inv,
                                  nHalf);
        accepts++;
    }
    else
    {
        for (int j = 0; j < nDimensions; j++)
        {
            prNew[i][j] = prOld[i][j];
        }
        rAbsNew[i] = rAbsOld[i];
        rejects++;
    }
    endOfSingleParticleStep(cycle, i);
}

void VMCSolver::runSingleStep(int i, int cycle)
{
    rAbsNew[i] = 0;
    for (int j = 0; j < nDimensions; j++)
    {
        prNew[i][j] = prOld[i][j] + stepLength * (dist_uniform(gen) - 0.5);
        rAbsNew[i] += prNew[i][j] * prNew[i][j];
    }
    rAbsNew[i] = sqrt(rAbsNew[i]);
    // Recalculate the value of the wave function
    waveFuncValNew = getWaveFuncVal(prNew, rAbsNew);
    ratio = waveFuncValNew * waveFuncValNew / (waveFuncValOld * waveFuncValOld);
    // Check for step acceptance (if yes,
    // update position, if no, reset position)
    if (dist_uniform(gen) <= ratio)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            prOld[i][j] = prNew[i][j];
        }
        rAbsOld[i] = rAbsNew[i];
        waveFuncValOld = waveFuncValNew;
        accepts++;
    }
    else
    {
        for (int j = 0; j < nDimensions; j++)
        {
            prNew[i][j] = prOld[i][j];
        }
        rAbsNew[i] = rAbsOld[i];
        rejects++;
    }
    // update energies
    endOfSingleParticleStep(cycle, i);
}

void VMCSolver::runStepSlater(int cycle)
{
    for (int i = 0; i < nParticles; i++)
    {
        runSingleStepSlater(i, cycle);
    }
    // ALL PARTICLES MOVED ONE STEP AT THIS POINT.
    endOfCycle(cycle);
}

void VMCSolver::runStep(int cycle)
{
    // New position to test
    for (int i = 0; i < nParticles; i++)
    {
        runSingleStep(i, cycle);
    }
    // All particles moved one step at this point.
    endOfCycle(cycle);
}

void VMCSolver::supressOutput()
{
    outputSupressed = true;
}

bool VMCSolver::initRunVariables()
{
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
    nHalf = nParticles / 2;

    wave_functions::alpha = alpha;
    wave_functions::beta = beta;
    wave_functions::nDimensions = nDimensions;
    wave_functions::h = h;
    wave_functions::h2 = h2;
    wave_functions::charge = charge;
    wave_functions::nParticles = nParticles;
    wave_functions::nHalf = nHalf;

    // Set the wave function as a function pointer
    if (waveFunction == WAVE_FUNCTION_1)
    {
        getWaveFuncVal = wave_functions::getWaveFuncHeliumNoCor;
    }
    else if (waveFunction == WAVE_FUNCTION_2)
    {
        getWaveFuncVal = wave_functions::getWaveFuncHelium;
        usingCorrelation = true;
    }
    else if (waveFunction == WAVE_FUNCTION_BERYLLIUM_1)
    {
        getWaveFuncVal = wave_functions::getWaveBerylliumNoCor;
    }
    else if (waveFunction == WAVE_FUNCTION_BERYLLIUM_2)
    {
        getWaveFuncVal = wave_functions::getWaveBeryllium;
        usingCorrelation = true;
    }
    wave_functions::getWaveFuncVal = getWaveFuncVal;

    // Set the local energy function as a function pointer
    if (localEnergyFunction == LOCAL_ENERGY_GENERIC)
    {
        getLocalEnergy = wave_functions::getLocalEnergyGeneric;
        usingCorrelation = true;
    }
    else if (localEnergyFunction == LOCAL_ENERGY_HELIUM_1)
    {
        getLocalEnergy = wave_functions::getLocalEnergyHeliumNoCor;
    }
    else if (localEnergyFunction == LOCAL_ENERGY_HELIUM_2)
    {
        getLocalEnergy = wave_functions::getLocalEnergyHelium;
        usingCorrelation = true;
    }
    else if (localEnergyFunction == LOCAL_ENERGY_HYDROGEN)
    {
        getLocalEnergy = wave_functions::getLocalEnergyHydrogen;
        usingCorrelation = true;
    }
    else if (localEnergyFunction == LOCAL_ENERGY_GENERIC_NOCOR)
    {
        getLocalEnergy = wave_functions::getLocalEnergyGenericNoCor;
    }
    else if (localEnergyFunction == LOCAL_ENERGY_SLATER_NOCOR)
    {
    }
    else if (localEnergyFunction == LOCAL_ENERGY_SLATER)
    {
        usingCorrelation = true;
    }

    // Initialize arrays
    if (recordingPositions)
    {
        positions = Matrix(nParticles, nCycles);
        positions.reset();
        pPositions = positions.getArrayPointer();
    }
    if (recordingEnergyArray)
    {
        energyArray = Vector(nCycles);
        energyArray.reset();
        pEnergyArray = energyArray.getArrayPointer();
    }
    if (recordingDensity)
    {
        density = Matrix(nParticles, bins);
        density.reset();
        pDensity = density.getArrayPointer();
    }
    if (importanceSampling)
    {
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
    if (importanceSampling == true)
    {
        for (int i = 0; i < nParticles; i++)
        {
            for (int j = 0; j < nDimensions; j++)
            {
                prOld[i][j] = dist_gauss(gen) * sqrt(timeStep);
            }
        }
    }
    else
    {
        for (int i = 0; i < nParticles; i++)
        {
            for (int j = 0; j < nDimensions; j++)
            {
                prOld[i][j] = stepLength * (dist_uniform(gen) - 0.5);
            }
        }
    }
    rNew = rOld;
    // This is important to get the correct pointer to the new matrix.
    prNew = rNew.getArrayPointer();
    prOld = rOld.getArrayPointer();

    // Initialize the absolute positions.
    rAbsOldVec = Vector(nParticles);
    rAbsNewVec = Vector(nParticles);
    rAbsOld = rAbsOldVec.getArrayPointer();
    rAbsNew = rAbsNewVec.getArrayPointer();
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rAbsOld[i] += prOld[i][j] * prOld[i][j];
            rAbsNew[i] += prNew[i][j] * prNew[i][j];
        }
        rAbsOld[i] = sqrt(rAbsOld[i]);
        rAbsNew[i] = sqrt(rAbsNew[i]);
    }

    // Initialize the slater determinant of the initial positions.
    if (efficientSlater)
    {
        vS = Vector(nParticles / 2);
        S = vS.getArrayPointer();
        slater1 = Matrix(nParticles / 2, nParticles / 2);
        slater2 = Matrix(nParticles / 2, nParticles / 2);
        pslater1 = slater1.getArrayPointer();
        pslater2 = slater2.getArrayPointer();
        for (int i = 0; i < nParticles / 2; i++)
        {
            for (int j = 0; j < nParticles / 2; j++)
            {
                pslater1[i][j] = phi(j, prNew[i], rAbsNew[i]);
                pslater2[i][j] = phi(j, prNew[i + 2], rAbsNew[i + 2]);
            }
        }
        slater1Inv = CPhys::MatOp::getInverse(slater1);
        slater2Inv = CPhys::MatOp::getInverse(slater2);
        pslater1Inv = slater1Inv.getArrayPointer();
        pslater2Inv = slater2Inv.getArrayPointer();
    }
    // Finished without error (hopefully).
    return true;
}

void VMCSolver::endOfCycle(int cycle)
{
    if (!recordingR12Mean) return;
    // Calculate the radius of the particle
    double rAbs = 0;
    double rsq = 0;
    for (int j = 0; j < nDimensions; j++)
    {
        rsq += (prNew[1][j] - prNew[0][j]) * (prNew[1][j] - prNew[0][j]);
    }
    rAbs = sqrt(rsq);
    // Add it to a sum so we can calculate the mean.
    rAbsSum += rAbs;
}

void VMCSolver::endOfSingleParticleStep(int cycle, int i)
{
    // update energies
    if (localEnergyFunction == LOCAL_ENERGY_SLATER)
        deltaE = getLocalEnergySlater(prOld, rAbsOld);
    else if (localEnergyFunction == LOCAL_ENERGY_SLATER_NOCOR)
        deltaE = getLocalEnergySlaterNoCor(prOld, rAbsOld);
    else
        deltaE = getLocalEnergy(prOld, rAbsOld);
    energySum += deltaE;
    energySquaredSum += deltaE * deltaE;

    // Store in energy array.
    if (recordingEnergyArray)
    {
        pEnergyArray[cycle] += deltaE;
    }

    if (recordingPositions)
    {
        double rAbs = 0;
        for (int j = 0; j < nDimensions; j++)
        {
            rAbs += prOld[i][j] * prOld[i][j];
        }
        rAbs = sqrt(rAbs);
        pPositions[i][cycle] = rAbs;
    }

    // Calculate density
    if (recordingDensity)
    {
        int bin;
        double rsq = 0;
        for (int j = 0; j < nDimensions; j++)
        {
            rsq += prOld[i][j] * prOld[i][j];
        }
        double rAbs = sqrt(rsq);
        if (rAbs < rMax)
        {
            bin = rAbs / rMax * bins;
            pDensity[i][bin] += 1;
        }
    }
}

void VMCSolver::updateQuantumForceSlater(double** r, double* rAbs,
                                         double** qForce)
{
    /* double rkVec[nDimensions]; */
    /* double rkjVec[nDimensions]; */
    /* double rkjAbs; */
    /* double rkGrad[nDimensions]; */
    if (usingCorrelation)
    {
        /* for (int k = 0; k < nParticles; k++) */
        /* { */
        /*     // Reset rkVec. */
        /*     for (int x = 0; x < nDimensions; x++) */
        /*     { */
        /*         rkVec[x] = 0; */
        /*     } */
        /*     // Calculate rkVec. */
        /*     for (int j = 0; j < nParticles; j++) */
        /*     { */
        /*         if (j == k) continue; */
        /*         spinK = k / nHalf; */
        /*         spinJ = j / nHalf; */
        /*         switch (spinK + spinJ) */
        /*         { */
        /*             case 0: */
        /*                 a1 = 0.25; */
        /*                 break; */
        /*             case 1: */
        /*                 a1 = 0.5; */
        /*                 break; */
        /*             case 2: */
        /*                 a1 = 0.25; */
        /*                 break; */
        /*         } */
        /*         rkjAbs = 0; */
        /*         for (int x = 0; x < nDimensions; x++) */
        /*         { */
        /*             rkjVec[x] = r[j][x] - r[k][x]; */
        /*             rkjAbs += rkjVec[x] * rkjVec[x]; */
        /*         } */
        /*         rkjAbs = sqrt(rkjAbs); */
        /*         bkj = 1 / (1 + beta * rkjAbs); */
        /*         // Calculate the final vector, the gradient of the correction
         */
        /*         // function. */
        /*         for (int x = 0; x < nDimensions; x++) */
        /*         { */
        /*             rkVec[x] -= rkjVec[x] * a1 * bkj * bkj / rkjAbs; */
        /*         } */
        /*     } */
        /*     tmp = 0; */
        /*     // Reset rkGrad. */
        /*     for (int x = 0; x < nDimensions; x++) */
        /*     { */
        /*         rkGrad[x] = 0; */
        /*     } */
        /*     // Calculate rkGrad. */
        /*     for (int x = 0; x < nDimensions; x++) */
        /*     { */
        /*         /1* rkGrad[x] = r[k][x]*tmp/rAbs[k]; *1/ */
        /*         for (int j = 0; j < nHalf; j++) */
        /*         { */
        /*             if (k < nHalf) */
        /*                 rkGrad[x] += phiD(j, r[k], rAbs[k], x) *
         * pslater1Inv[j][k]; */
        /*             else */
        /*                 rkGrad[x] += */
        /*                     phiD(j, r[k], rAbs[k], x) * pslater2Inv[j][k -
         * nHalf]; */
        /*         } */
        /*     } */
        /*     // Calculate DC. */
        /*     for (int x = 0; x < nDimensions; x++) */
        /*     { */
        /*         DC += rkVec[x] * rkGrad[x]; */
        /*     } */
        /* } */
    }
    // Not using correlation.
    else
    {
        for (int k = 0; k < nParticles; k++)
        {
            // Reset qForce.
            for (int x = 0; x < nDimensions; x++)
            {
                qForce[k][x] = 0;
            }
            // Calculate.
            for (int x = 0; x < nDimensions; x++)
            {
                for (int j = 0; j < nHalf; j++)
                {
                    if (k < nHalf)
                        qForce[k][x] +=
                            phiD(j, r[k], rAbs[k], x) * pslater1Inv[j][k];
                    else
                        qForce[k][x] += phiD(j, r[k], rAbs[k], x) *
                                        pslater2Inv[j][k - nHalf];
                }
            }
        }
    }
}

void VMCSolver::updateQuantumForce(double** r, double* rAbs, double** qForce,
                                   double factor)
{
    double waveFunctionMinus;
    double waveFunctionPlus;
    double r0;
    double rAbs0;
    // Kinetic energy
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            // Save the current position of the particle. 
            r0 = r[i][j];
            rAbs0 = 0;
            for (int x = 0; x < nDimensions; x++) 
            {
                rAbs0 += r[i][x]*r[i][x];
            }
            rAbs0 = sqrt(rAbs0);

            // Evaluate at plus h. 
            r[i][j] = r0 - h;
            rAbs[i] = 0;
            for (int x = 0; x < nDimensions; x++) 
            {
                rAbs[i] += r[i][x]*r[i][x];
            }
            rAbs[i] = sqrt(rAbs[i]);
            waveFunctionMinus = getWaveFuncVal(r, rAbs);

            // Evaluate at minus h.
            r[i][j] = r0 + h;
            rAbs[i] = 0;
            for (int x = 0; x < nDimensions; x++) 
            {
                rAbs[i] += r[i][x]*r[i][x];
            }
            rAbs[i] = sqrt(rAbs[i]);
            waveFunctionPlus = getWaveFuncVal(r, rAbs);

            // Reset position. 
            r[i][j] = r0;
            rAbs[i] = rAbs0;
            qForce[i][j] = (waveFunctionPlus - waveFunctionMinus) * h / factor;
        }
    }
}

void VMCSolver::updateSlater(int i)
{
    for (int j = 0; j < nParticles / 2; j++)
    {
        if (i < nParticles / 2)
            pslater1[i][j] = phi(j, prNew[i], rAbsNew[i]);
        else
            pslater2[i - nHalf][j] = phi(j, prNew[i], rAbsNew[i]);
    }
}

Vector VMCSolver::getEnergyArray()
{
    return energyArray;
}

double VMCSolver::getStepLength()
{
    return stepLength;
}

double VMCSolver::getEnergy()
{
    return energy;
}

double VMCSolver::getEnergySquared()
{
    return energySquared;
}

double VMCSolver::getR12Mean()
{
    return mean;
}

void VMCSolver::clear()
{
    waveFunction = 0;
    localEnergyFunction = 0;
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
    ratio = 0;

    usingCorrelation = false;
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

double VMCSolver::getAcceptanceRatio()
{
    return double(accepts) * 100 / (rejects + accepts);
}

double VMCSolver::getLocalEnergySlaterNoCor(double** r, double* rAbs)
{
    double sum = 0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nHalf; j++)
        {
            if (i < nHalf)
                sum += phiDD(j, r[i], rAbs[i]) * pslater1Inv[j][i];
            else
                sum += phiDD(j, r[i], rAbs[i]) * pslater2Inv[j][i - nHalf];
        }
    }

    // Potential energy.
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for (int i = 0; i < nParticles; i++)
    {
        rSingleParticle = 0;
        for (int j = 0; j < nDimensions; j++)
        {
            rSingleParticle += r[i][j] * r[i][j];
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    return -0.5 * sum + potentialEnergy;
}

double VMCSolver::getLocalEnergySlater(double** r, double* rAbs)
{
    double DD = 0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nHalf; j++)
        {
            if (i < nHalf)
                DD += phiDD(j, r[i], rAbs[i]) * pslater1Inv[j][i];
            else
                DD += phiDD(j, r[i], rAbs[i]) * pslater2Inv[j][i - nHalf];
        }
    }

    // Potential energy.
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for (int i = 0; i < nParticles; i++)
    {
        rSingleParticle = 0;
        for (int j = 0; j < nDimensions; j++)
        {
            rSingleParticle += r[i][j] * r[i][j];
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    double rij = 0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = i + 1; j < nParticles; j++)
        {
            rij = 0;
            for (int k = 0; k < nDimensions; k++)
            {
                rij += (r[i][k] - r[j][k]) * (r[i][k] - r[j][k]);
            }
            potentialEnergy += 1 / sqrt(rij);
        }
    }

    double CC = 0;
    double a1, a2;
    double tmp = 0;
    double rkj, rki;
    double bki, bkj;
    int spinI, spinJ, spinK;
    double dot;
    for (int k = 0; k < nParticles; k++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            if (j == k) continue;
            rkj = 0;
            for (int x = 0; x < nDimensions; x++)
            {
                tmp = (r[j][x] - r[k][x]);
                rkj += tmp * tmp;
            }
            rkj = sqrt(rkj);
            spinK = k / nHalf;
            spinJ = j / nHalf;
            switch (spinK + spinJ)
            {
                case 0:
                    a1 = 0.25;
                    break;
                case 1:
                    a1 = 0.5;
                    break;
                case 2:
                    a1 = 0.25;
                    break;
            }
            bkj = 1 / (1 + beta * rkj);
            CC += 2 * a1 * bkj * bkj / rkj;
            CC -= 2 * a1 * beta * bkj * bkj * bkj;
            for (int i = 0; i < nParticles; i++)
            {
                if (i == k) continue;
                rki = 0;
                dot = 0;
                for (int x = 0; x < nDimensions; x++)
                {
                    tmp = (r[i][x] - r[k][x]);
                    rki += tmp * tmp;
                    dot += r[k][x] * r[k][x] - r[k][x] * r[j][x] -
                           r[k][x] * r[i][x] + r[i][x] * r[j][x];
                }
                rki = sqrt(rki);
                bki = 1 / (1 + beta * rki);
                spinI = i / nHalf;
                switch (spinK + spinI)
                {
                    case 0:
                        a2 = 0.25;
                        break;
                    case 1:
                        a2 = 0.5;
                        break;
                    case 2:
                        a2 = 0.25;
                        break;
                }
                CC += dot / (rki * rkj) * a1 * a2 * bki * bki * bkj * bkj;
            }
        }
    }

    double DC = 0;
    double rkVec[nDimensions];
    double rkjVec[nDimensions];
    double rkjAbs;
    double rkGrad[nDimensions];
    for (int k = 0; k < nParticles; k++)
    {
        // Reset rkVec.
        for (int x = 0; x < nDimensions; x++)
        {
            rkVec[x] = 0;
        }
        // Calculate rkVec.
        for (int j = 0; j < nParticles; j++)
        {
            if (j == k) continue;
            spinK = k / nHalf;
            spinJ = j / nHalf;
            switch (spinK + spinJ)
            {
                case 0:
                    a1 = 0.25;
                    break;
                case 1:
                    a1 = 0.5;
                    break;
                case 2:
                    a1 = 0.25;
                    break;
            }
            rkjAbs = 0;
            for (int x = 0; x < nDimensions; x++)
            {
                rkjVec[x] = r[j][x] - r[k][x];
                rkjAbs += rkjVec[x] * rkjVec[x];
            }
            rkjAbs = sqrt(rkjAbs);
            bkj = 1 / (1 + beta * rkjAbs);
            // Calculate the final vector, the gradient of the correction
            // function.
            for (int x = 0; x < nDimensions; x++)
            {
                rkVec[x] -= rkjVec[x] * a1 * bkj * bkj / rkjAbs;
            }
        }
        tmp = 0;
        // Reset rkGrad.
        for (int x = 0; x < nDimensions; x++)
        {
            rkGrad[x] = 0;
        }
        // Calculate rkGrad.
        for (int x = 0; x < nDimensions; x++)
        {
            for (int j = 0; j < nHalf; j++)
            {
                if (k < nHalf)
                    rkGrad[x] += phiD(j, r[k], rAbs[k], x) * pslater1Inv[j][k];
                else
                    rkGrad[x] +=
                        phiD(j, r[k], rAbs[k], x) * pslater2Inv[j][k - nHalf];
            }
        }
        // Calculate DC.
        for (int x = 0; x < nDimensions; x++)
        {
            DC += rkVec[x] * rkGrad[x];
        }
    }
    return -0.5 * DD - 0.5 * CC - DC + potentialEnergy;
}
