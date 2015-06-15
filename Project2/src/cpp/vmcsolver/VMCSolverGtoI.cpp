#include "VMCSolverGtoI.h"

using namespace CPhys;
using namespace std;

VMCSolverGtoI::VMCSolverGtoI()
{
    // Initialize the random generators.
    dist_uniform = std::uniform_real_distribution<double>(0.0, 1.0);
    dist_gauss = std::normal_distribution<double>(0.0, sqrt(2));
    clear();
}

void VMCSolverGtoI::setSeed(long seed)
{
    gen.seed(seed);
}

inline double VMCSolverGtoI::getCorrelationRatio(int i)
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
                        break;
                    case 1:
                        a = 0.5;
                        break;
                    case 2:
                        a = 0.25;
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

double VMCSolverGtoI::calc_dE_dBeta(){
    double rkjVec[nDimensions];
    double rkjAbs;
    dE_dBeta = 0;
    for (int k = 0; k < nParticles; k++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            // Run the code if j > k.
            if (j <= k) continue;
            int spinK = k / nHalf;
            int spinJ = j / nHalf;
            double a1;
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
                rkjVec[x] = prOld[j][x] - prOld[k][x];
                rkjAbs += rkjVec[x] * rkjVec[x];
            }
            rkjAbs = sqrt(rkjAbs);
            double bkj = 1 / (1 + beta * rkjAbs);

            /* dE_dBeta -= a1 * bkj * bkj * rkjAbs * exp(a1*rkjAbs*bkj); */
            dE_dBeta -= a1 * bkj * bkj * rkjAbs;
        }
    }
}

void VMCSolverGtoI::startOfCycle()
{
    updateQuantumForceSlater(prOld, rAbsOld, pqForceOld, pslater1Old,
                             pslater2Old, pslater1InvOld, pslater2InvOld);
}

/** \brief Calculates the quantum force.
 *
 * The only uknown is the qForce variable. All other paramters must be correct.
 * \param r Particle positions.
 * \param rAbs Particle absolute positions.
 * \param pslater1 Up slater matrix.
 * \param pslater2 Down slater matrix.
 * \param pslater1Inv Up inverse slater matrix.
 * \param pslater2Inv Down inverse slater matrix.
 */
inline void VMCSolverGtoI::updateQuantumForceSlater(double** r, double* rAbs,
                                         double** qForce, double** pslater1,
                                         double** pslater2,
                                         double** pslater1Inv,
                                         double** pslater2Inv)
{
    double rkVec[nDimensions];
    double rkjVec[nDimensions];
    double rkjAbs;
    double rkGrad[nDimensions];
    int spinK;
    int spinJ;
    double a1;
    double bkj;
    double tmp;
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
                    rkGrad[x] +=
                        phiD(k, j, r, rAbs, x) * pslater1Inv[j][k];
                else
                    rkGrad[x] += phiD(k, j, r, rAbs, x) *
                                 pslater2Inv[j][k - nHalf];
            }
        }
        // Calculate DC.
        for (int x = 0; x < nDimensions; x++)
        {
            DC += rkVec[x] * rkGrad[x];
        }
    }
}

/** \brief Copy the slater old variabes to the new variables in an effecient
 * way.
 *
 * \param slater1New Slater matrix that will be overwritten.
 * \param slater1Old Slater matrix that will be copied.
 * \param slater2New Slater matrix that will be overwritten.
 * \param slater2Old Slater matrix that will be copied.
 * \param slater1InvNew Slater matrix that will be overwritten.
 * \param slater1InvOld Slater matrix that will be copied.
 * \param slater2InvNew Slater matrix that will be overwritten.
 * \param slater2InvOld Slater matrix that will be copied.
 */
inline void VMCSolverGtoI::updateSlater(int i, double** slater1New, double** slater1Old,
                             double** slater2New, double** slater2Old,
                             double** slater1InvNew, double** slater1InvOld,
                             double** slater2InvNew, double** slater2InvOld)
{
    if (i < nHalf)
    {
        for (int j = 0; j < nHalf; j++)
        {
            slater1New[i][j] = slater1Old[i][j];
            for (int i = 0; i < nHalf; i++)
            {
                slater1InvNew[i][j] = slater1InvOld[i][j];
            }
        }
    }
    else
    {
        for (int j = 0; j < nHalf; j++)
        {
            slater2New[i - nHalf][j] = slater2Old[i - nHalf][j];
            for (int i = 0; i < nHalf; i++)
            {
                slater2InvNew[i][j] = slater2InvOld[i][j];
            }
        }
    }
}

void VMCSolverGtoI::runSingleStep(int i)
{
    // New position to test
    rAbsNew[i] = 0;
    for (int j = 0; j < nDimensions; j++)
    {
        prNew[i][j] = prOld[i][j] + dist_gauss(gen) * sqrt(timeStep) +
                      2 * pqForceOld[i][j] * timeStep * D;
        rAbsNew[i] += prNew[i][j] * prNew[i][j];
    }
    rAbsNew[i] = sqrt(rAbsNew[i]);

    // Update the slater matrix and calculate the ratio.
    double ratioTmp = 0;
    for (int j = 0; j < nHalf; j++)
    {
        if (i < nHalf)
        {
            pslater1New[i][j] = phi(i, j, prNew, rAbsNew);
            ratioTmp += pslater1New[i][j] * pslater1InvOld[j][i];
        }
        else
        {
            pslater2New[i - nHalf][j] = phi(i, j, prNew, rAbsNew);
            ratioTmp +=
                pslater2New[i - nHalf][j] * pslater2InvOld[j][i - nHalf];
        }
    }

    ratio = ratioTmp * getCorrelationRatio(i);
    ratio = ratio * ratio;

    // Update the inverse slater matrix for the new particle.
    if (i < nHalf)
    {
        pMatOp::updateInverse(i, ratioTmp, pslater1New, pslater1InvNew, nHalf);
    }
    else
    {
        pMatOp::updateInverse(i - nHalf, ratioTmp, pslater2New, pslater2InvNew,
                              nHalf);
    }
    // Update the quantum force.
    updateQuantumForceSlater(prNew, rAbsNew, pqForceNew, pslater1New,
                             pslater2New, pslater1InvNew, pslater2InvNew);
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
        for (int x = 0; x < nDimensions; x++)
        {
            prOld[i][x] = prNew[i][x];
            pqForceOld[i][x] = pqForceNew[i][x];
        }
        rAbsOld[i] = rAbsNew[i];
        // Update the slater matrices.
        updateSlater(i, pslater1Old, pslater1New, pslater2Old, pslater2New,
                     pslater1InvOld, pslater1InvNew, pslater2InvOld,
                     pslater2InvNew);
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
        // Update the slater matrices.
        updateSlater(i, pslater1New, pslater1Old, pslater2New, pslater2Old,
                     pslater1InvNew, pslater1InvOld, pslater2InvNew,
                     pslater2InvOld);
        rejects++;
    }
    // update energies
    deltaE = getLocalEnergySlater(prOld, rAbsOld);
}


bool VMCSolverGtoI::initRunVariables()
{
    accepts = 0;
    rejects = 0;
    deltaE = 0;
    nHalf = nParticles / 2;


    if (waveFunction == WAVE_FUNCTION_HELIUM_GTO){
        phi = gto_helium::phi;
        phiD = gto_helium::phiD;
        phiDD = gto_helium::phiDD;
    } 
    else if (waveFunction == WAVE_FUNCTION_BERYLLIUM_GTO){
        phi = gto_beryllium::phi;
        phiD = gto_beryllium::phiD;
        phiDD = gto_beryllium::phiDD;
    } 
    else if (waveFunction == WAVE_FUNCTION_NEON_GTO){
        phi = gto_neon::phi;
        phiD = gto_neon::phiD;
        phiDD = gto_neon::phiDD;
    } 

    // Initialize arrays that are used in importance sampling.
    qForceOld = Matrix(nParticles, nDimensions);
    qForceNew = Matrix(nParticles, nDimensions);
    pqForceOld = qForceOld.getArrayPointer();
    pqForceNew = qForceNew.getArrayPointer();

    // Initialize position arrays.
    rOld = Matrix(nParticles, nDimensions);
    rNew = Matrix(nParticles, nDimensions);
    prNew = rNew.getArrayPointer();
    prOld = rOld.getArrayPointer();
    // Initial trial positions (different with importance sampling).
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            prOld[i][j] = dist_gauss(gen) * sqrt(timeStep);
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
    {
        // Create the arrays and assign the correct pointers.
        slater1Old = Matrix(nParticles / 2, nParticles / 2);
        slater2Old = Matrix(nParticles / 2, nParticles / 2);
        pslater1Old = slater1Old.getArrayPointer();
        pslater2Old = slater2Old.getArrayPointer();
        // Fill the arrays.
        for (int i = 0; i < nParticles / 2; i++)
        {
            for (int j = 0; j < nParticles / 2; j++)
            {
                pslater1Old[i][j] = phi(i,j, prNew, rAbsNew);
                pslater2Old[i][j] = phi(i+nHalf,j, prNew, rAbsNew);
            }
        }
        // Calculate the inverse and assign the pointers.
        slater1InvOld = CPhys::MatOp::getInverse(slater1Old);
        slater2InvOld = CPhys::MatOp::getInverse(slater2Old);
        pslater1InvOld = slater1InvOld.getArrayPointer();
        pslater2InvOld = slater2InvOld.getArrayPointer();

        // Copy old values to the new values.
        slater1InvNew = slater1InvOld;
        slater2InvNew = slater2InvOld;
        slater1New = slater1Old;
        slater2New = slater2Old;
        // Assign the correct pointers.
        pslater1InvNew = slater1InvNew.getArrayPointer();
        pslater2InvNew = slater2InvNew.getArrayPointer();
        pslater1New = slater1New.getArrayPointer();
        pslater2New = slater2New.getArrayPointer();
    }
    // Finished without error (hopefully).
    return true;
}

void VMCSolverGtoI::clear()
{
    waveFunction = 0;
    nDimensions = 0;
    charge = 0;
    stepLength = 0;
    nParticles = 0;
    h = 0;
    hInv = 0;
    h2Inv = 0;
    idum = 1;
    beta = 0;

    deltaE = 0;
    ratio = 0;

    // Initialize all variables, they are mostly overwritten.
    slater1Old = Matrix();
    slater1New = Matrix();
    slater2Old = Matrix();
    slater2New = Matrix();
    slater1InvOld = Matrix();
    slater1InvNew = Matrix();
    slater2InvOld = Matrix();
    slater2InvNew = Matrix();
    rOld = Matrix();
    rNew = Matrix();
    pslater1Old = slater1Old.getArrayPointer();
    pslater1New = slater1New.getArrayPointer();
    pslater1InvOld = slater1InvOld.getArrayPointer();
    pslater1InvNew = slater1InvNew.getArrayPointer();
    pslater2Old = slater2Old.getArrayPointer();
    pslater2New = slater2New.getArrayPointer();
    pslater2InvOld = slater2InvOld.getArrayPointer();
    pslater2InvNew = slater2InvNew.getArrayPointer();
    prOld = rOld.getArrayPointer();
    prNew = rNew.getArrayPointer();

    qForceOld = Matrix();
    qForceNew = Matrix();
    pqForceOld = qForceOld.getArrayPointer();
    pqForceNew = qForceNew.getArrayPointer();
    vS = Vector();
    S = vS.getArrayPointer();
}

inline double VMCSolverGtoI::getLocalEnergySlater(double** r, double* rAbs)
{
    DD = 0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nHalf; j++)
        {
            if (i < nHalf)
                DD += phiDD(i, j, r, rAbs) * pslater1InvOld[j][i];
            else
                DD += phiDD(i, j, r, rAbs) * pslater2InvOld[j][i - nHalf];
        }
    }

    // Potential energy.
    potentialEnergy = 0;
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

    CC = 0;
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

    double rkVec[nDimensions];
    double rkjVec[nDimensions];
    double rkjAbs;
    double rkGrad[nDimensions];
    DC = 0;
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
                    rkGrad[x] +=
                        phiD(k, j, r, rAbs, x) * pslater1InvOld[j][k];
                else
                    rkGrad[x] += phiD(k, j, r, rAbs, x) *
                                 pslater2InvOld[j][k - nHalf];
            }
        }
        // Calculate DC.
        for (int x = 0; x < nDimensions; x++)
        {
            DC += rkVec[x] * rkGrad[x];
        }
    }

    DD = -0.5 * DD;
    CC = -0.5 * CC;
    DC = -DC;

    return DD + CC + DC + potentialEnergy;
}
