#include <unittest++/UnitTest++.h>
#include <unittest++/TestReporterStdout.h>
#include "../vmcsolver/VMCWrapper.h"
#include <string.h>
#include <stdio.h>

using namespace std;

static void initWaveFunctions(VMCSolver& solver)
{
    wave_functions::alpha = solver.alpha;
    wave_functions::beta = solver.beta;
    wave_functions::nDimensions = solver.nDimensions;
    wave_functions::h = solver.h;
    wave_functions::hInv = solver.hInv;
    wave_functions::h2Inv = solver.h2Inv;
    wave_functions::charge = solver.charge;
    wave_functions::nParticles = solver.nParticles;
    wave_functions::nHalf = solver.nHalf;
    wave_functions::getWaveFuncVal = solver.getWaveFuncVal;

    wave_functions_gto::nDimensions = solver.nDimensions;
    wave_functions_gto::nParticles = solver.nParticles;
    wave_functions_gto::nHalf = solver.nHalf;
    wave_functions_gto::beta = solver.beta;
}

SUITE(CPhys)
{
    TEST(Matrix)
    {
        Matrix mat1;
        mat1 = Matrix(10, 10);
        Matrix mat2 = Matrix(0, 0);
        /* Matrix mat3 = Matrix(3,3); */
        /* mat1 = mat3; */
    }

    TEST(LUdecomposition)
    {
        // LU decomposing mat into LU
        Matrix mat = Matrix(2, 2);
        Matrix L = Matrix(mat);
        Matrix U = Matrix(mat);
        mat(0, 0) = 4;
        mat(0, 1) = 3;
        mat(1, 0) = 6;
        mat(1, 1) = 3;
        CPhys::MatOp::decomposeLU(mat, L, U);
        CHECK_CLOSE(1.5, L(1, 0), 0.0001);

        CHECK_CLOSE(4, U(0, 0), 0.0001);
        CHECK_CLOSE(3, U(0, 1), 0.0001);
        CHECK_CLOSE(-1.5, U(1, 1), 0.0001);
        Matrix res = L * U;

        CHECK_CLOSE(4, res(0, 0), 0.0001);
        CHECK_CLOSE(3, res(0, 1), 0.0001);
        CHECK_CLOSE(6, res(1, 0), 0.0001);
        CHECK_CLOSE(3, res(1, 1), 0.0001);

        Matrix B = Matrix(3, 3);
        B(0, 0) = 7;
        B(0, 1) = 2;
        B(0, 2) = 1;
        B(1, 0) = 0;
        B(1, 1) = 3;
        B(1, 2) = -1;
        B(2, 0) = -3;
        B(2, 1) = 4;
        B(2, 2) = -2;

        CPhys::MatOp::decomposeLU(B, L, U);
        Matrix res1 = L * U;
        CHECK_CLOSE(B(0, 0), res1(0, 0), 0.0001);
        CHECK_CLOSE(B(0, 1), res1(0, 1), 0.0001);
        CHECK_CLOSE(B(0, 2), res1(0, 2), 0.0001);
        CHECK_CLOSE(B(1, 0), res1(1, 0), 0.0001);
        CHECK_CLOSE(B(1, 1), res1(1, 1), 0.0001);
        CHECK_CLOSE(B(1, 2), res1(1, 2), 0.0001);
        CHECK_CLOSE(B(2, 0), res1(2, 0), 0.0001);
        CHECK_CLOSE(B(2, 1), res1(2, 1), 0.0001);
        CHECK_CLOSE(B(2, 2), res1(2, 2), 0.0001);
    }

    TEST(ForwardSubstitution)
    {
        // Solving the equation:
        // LUx = b
        // Ly  = b
        // This means x is unity.
        // With the given paramters, this means x is unity.
        Matrix L = Matrix(2, 2);
        L(0, 0) = 1;
        L(0, 1) = 0;
        L(1, 0) = 1.5;
        L(1, 1) = 1;
        Matrix b = Matrix(2, 2);
        b(0, 0) = 4;
        b(0, 1) = 3;
        b(1, 0) = 6;
        b(1, 1) = 3;
        Matrix y = Matrix(2, 2);
        CPhys::MatOp::substituteForward(L, y, b);
        CHECK_CLOSE(4, y(0, 0), 0.0001);
        CHECK_CLOSE(3, y(0, 1), 0.0001);
        CHECK_CLOSE(0, y(1, 0), 0.0001);
        CHECK_CLOSE(-1.5, y(1, 1), 0.0001);

        // Assuming LU decomposition works, do a decomposition
        // and solve the equation L*Linv = I. And test this.
        Matrix B = Matrix(3, 3);
        B(0, 0) = 7;
        B(0, 1) = 2;
        B(0, 2) = 1;
        B(1, 0) = 0;
        B(1, 1) = 3;
        B(1, 2) = -1;
        B(2, 0) = -3;
        B(2, 1) = 4;
        B(2, 2) = -2;
        Matrix U;
        CPhys::MatOp::decomposeLU(B, L, U);
        Matrix LInv;
        Matrix I = Matrix(L.getN(), L.getN());
        I.eye();
        CPhys::MatOp::substituteForward(L, LInv, I);
        Matrix res = L * LInv;

        CHECK_CLOSE(1, res(0, 0), 0.0001);
        CHECK_CLOSE(0, res(0, 1), 0.0001);
        CHECK_CLOSE(0, res(0, 2), 0.0001);

        CHECK_CLOSE(0, res(1, 0), 0.0001);
        CHECK_CLOSE(1, res(1, 1), 0.0001);
        CHECK_CLOSE(0, res(1, 2), 0.0001);

        CHECK_CLOSE(0, res(2, 0), 0.0001);
        CHECK_CLOSE(0, res(2, 1), 0.0001);
        CHECK_CLOSE(1, res(2, 2), 0.0001);
    }

    TEST(BackwardSubstitution)
    {
        // Solving the equation:
        // LUx = b
        // Ly  = b
        // With the given paramters, this means x is unity.
        Matrix U = Matrix(2, 2);
        U(0, 0) = 4;
        U(0, 1) = 3;
        U(1, 0) = 0;
        U(1, 1) = -1.5;
        Matrix y = Matrix(2, 2);
        y(0, 0) = 4;
        y(0, 1) = 3;
        y(1, 0) = 0;
        y(1, 1) = -1.5;
        Matrix x = Matrix(2, 2);
        CPhys::MatOp::substituteBackward(U, x, y);
        CHECK_CLOSE(1, x(0, 0), 0.0001);
        CHECK_CLOSE(0, x(0, 1), 0.0001);
        CHECK_CLOSE(0, x(1, 0), 0.0001);
        CHECK_CLOSE(1, x(1, 1), 0.0001);

        // Assuming LU decomposition works, do a decomposition
        // and solve the equation U*Uinv = I. And test this.
        Matrix B = Matrix(3, 3);
        B(0, 0) = 7;
        B(0, 1) = 2;
        B(0, 2) = 1;
        B(1, 0) = 0;
        B(1, 1) = 3;
        B(1, 2) = -1;
        B(2, 0) = -3;
        B(2, 1) = 4;
        B(2, 2) = -2;
        Matrix L;
        Matrix UInv;
        Matrix I = Matrix(B.getN(), B.getN());
        I.eye();
        CPhys::MatOp::decomposeLU(B, L, U);
        CPhys::MatOp::substituteBackward(U, UInv, I);
        Matrix res = U * UInv;

        CHECK_CLOSE(1, res(0, 0), 0.0001);
        CHECK_CLOSE(0, res(0, 1), 0.0001);
        CHECK_CLOSE(0, res(0, 2), 0.0001);

        CHECK_CLOSE(0, res(1, 0), 0.0001);
        CHECK_CLOSE(1, res(1, 1), 0.0001);
        CHECK_CLOSE(0, res(1, 2), 0.0001);

        CHECK_CLOSE(0, res(2, 0), 0.0001);
        CHECK_CLOSE(0, res(2, 1), 0.0001);
        CHECK_CLOSE(1, res(2, 2), 0.0001);
    }

    TEST(Multiplication)
    {
        Matrix L = Matrix(2, 2);
        L(0, 0) = 1;
        L(0, 1) = 0;
        L(1, 0) = 1.5;
        L(1, 1) = 1;
        Matrix y = Matrix(2, 2);
        y(0, 0) = 4;
        y(0, 1) = 3;
        y(1, 0) = 0;
        y(1, 1) = -1.5;
        Matrix res = L * y;
        CHECK_CLOSE(4, res(0, 0), 0.0001);
        CHECK_CLOSE(3, res(0, 1), 0.0001);
        CHECK_CLOSE(6, res(1, 0), 0.0001);
        CHECK_CLOSE(3, res(1, 1), 0.0001);
    }

    TEST(Inverse)
    {
        // Check if we get the identity matrix if we multiply back.
        Matrix A = Matrix(2, 2);
        A(0, 0) = 4;
        A(0, 1) = 3;
        A(1, 0) = 6;
        A(1, 1) = 3;
        Matrix AInv = CPhys::MatOp::getInverse(A);

        Matrix res1 = A * AInv;
        CHECK_CLOSE(1, res1(0, 0), 0.0001);
        CHECK_CLOSE(0, res1(0, 1), 0.0001);
        CHECK_CLOSE(0, res1(1, 0), 0.0001);
        CHECK_CLOSE(1, res1(1, 1), 0.0001);

        // Check if we get the identity matrix if we multiply back.
        Matrix B = Matrix(3, 3);
        B(0, 0) = 7;
        B(0, 1) = 2;
        B(0, 2) = 1;
        B(1, 0) = 0;
        B(1, 1) = 3;
        B(1, 2) = -1;
        B(2, 0) = -3;
        B(2, 1) = 4;
        B(2, 2) = -2;
        Matrix BInv = CPhys::MatOp::getInverse(B);

        Matrix res2 = B * BInv;
        CHECK_CLOSE(1, res2(0, 0), 0.0001);
        CHECK_CLOSE(0, res2(0, 1), 0.0001);
        CHECK_CLOSE(0, res2(0, 2), 0.0001);

        CHECK_CLOSE(0, res2(1, 0), 0.0001);
        CHECK_CLOSE(1, res2(1, 1), 0.0001);
        CHECK_CLOSE(0, res2(1, 2), 0.0001);

        CHECK_CLOSE(0, res2(2, 0), 0.0001);
        CHECK_CLOSE(0, res2(2, 1), 0.0001);
        CHECK_CLOSE(1, res2(2, 2), 0.0001);
    }

    TEST(Determinant)
    {
        Matrix AOld = Matrix(2, 2);
        AOld(0, 0) = 4;
        AOld(0, 1) = 3;
        AOld(1, 0) = 3;
        AOld(1, 1) = 2;
        double det = CPhys::MatOp::getDet(AOld);
        CHECK_CLOSE(-1, det, 0.001);
    }

    TEST(UpdateInverse)
    {
// Check if the function updateInverse works (and in parallel)
#pragma omp parallel
        {
            Matrix AOld;
            AOld = Matrix(2, 2);
            AOld(0, 0) = 4;
            AOld(0, 1) = 3;
            AOld(1, 0) = 3;
            AOld(1, 1) = 2;
            Matrix AOldInv = CPhys::MatOp::getInverse(AOld);
            double detAOld = -1;

            Matrix ANew = Matrix(AOld);
            ANew(0, 0) = ANew(0, 0) + 1;
            ANew(0, 1) = ANew(0, 1) + 1;
            Matrix ANewInv = CPhys::MatOp::getInverse(ANew);
            double detANew = ANew(0, 0) * ANew(1, 1) - ANew(1, 0) * ANew(0, 1);
            double ratio = detANew / detAOld;

            Matrix testMat = Matrix(AOldInv);
            double** pANew = ANew.getArrayPointer();
            double** pTestMat = testMat.getArrayPointer();
            CPhys::pMatOp::updateInverse(0, ratio, pANew, pTestMat,
                                         ANew.getN());

            CHECK_CLOSE(ANewInv(0, 0), testMat(0, 0), 0.0001);
            CHECK_CLOSE(ANewInv(0, 1), testMat(0, 1), 0.0001);
            CHECK_CLOSE(ANewInv(1, 0), testMat(1, 0), 0.0001);
            CHECK_CLOSE(ANewInv(1, 1), testMat(1, 1), 0.0001);
        }
    }
}

SUITE(Hydrogen)
{
    VMCWrapper wrapper = VMCWrapper();

    TEST(Instantiate)
    {
        wrapper.charge = 1;
        wrapper.alpha = 1;
        wrapper.beta = 0;
        wrapper.nDimensions = 3;
        wrapper.nParticles = 1;
        wrapper.stepLength = 1.52;
        wrapper.nCycles = 10000;
        wrapper.h = 0.001;
        wrapper.hInv = 1e3;
        wrapper.h2Inv = 1e6;
        wrapper.idum = 1;
    }

    TEST(GroundStateGenergic)
    {
        wrapper.useWaveFunction1();
        wrapper.useLocalEnergyGenericNoCor();
        double energy = 0;
        VMCSolver solver = wrapper.getInitializedSolver();
        // Run simulation.
        for (int cycle = 0; cycle < wrapper.nCycles; cycle++) 
        {
            for (int i = 0; i < wrapper.nParticles; i++) 
            {
                solver.runSingleStep(i,cycle);
                energy += solver.deltaE;
            }
        }
        energy /= (wrapper.nParticles*wrapper.nCycles);
        CHECK_CLOSE(-0.5, energy, 0.01);
    }

    TEST(GroundStateAnalytic)
    {
        wrapper.waveFunction = wrapper.WAVE_FUNCTION_1;
        wrapper.localEnergyFunction = wrapper.LOCAL_ENERGY_HYDROGEN;
        double energy = 0;
        VMCSolver solver = wrapper.getInitializedSolver();
        // Run simulation.
        for (int cycle = 0; cycle < wrapper.nCycles; cycle++) 
        {
            for (int i = 0; i < wrapper.nParticles; i++) 
            {
                solver.runSingleStep(i,cycle);
                energy += solver.deltaE;
            }
        }
        energy /= (wrapper.nParticles*wrapper.nCycles);
        CHECK_CLOSE(-0.5, energy, 0.01);
    }
}

SUITE(VMCWrapper)
{
    TEST(getSolver)
    {
        VMCWrapper solver = VMCWrapper();
        solver.charge = 2;
        solver.alpha = 1.66;
        solver.beta = 0.8;
        solver.nDimensions = 3;
        solver.nParticles = 2;
        solver.stepLength = 1.52;
        solver.nCycles = 10000;
        solver.h = 0.001;
        solver.hInv = 1e3;
        solver.h2Inv = 1e6;
        solver.idum = 1;
        solver.useWaveFunction2();
        solver.useLocalEnergyHelium2();
        VMCSolver solverUnique = solver.getInitializedSolver();
        CHECK_EQUAL(solver.charge, solverUnique.charge);
        CHECK_EQUAL(solver.alpha, solverUnique.alpha);
    }

    TEST(RngPositions)
    {
        VMCWrapper solver = VMCWrapper();
        solver.charge = 4;
        solver.alpha = 3.75;
        solver.beta = 0.8;
        solver.nDimensions = 3;
        solver.nParticles = 4;
        solver.stepLength = 1.52;
        solver.nCycles = 1000;
        solver.h = 0.001;
        solver.hInv = 1e3;
        solver.h2Inv = 1e6;
        solver.idum = 1;
        solver.useEfficientSlater(true);
        solver.useLocalEnergyGenericNoCor();
        solver.useWaveFunctionBeryllium1();
        VMCSolver solverUnique1 = solver.getInitializedSolver();
        solver.useLocalEnergySlater();
        VMCSolver solverUnique2 = solver.getInitializedSolver();
        // Check that the positions are the same for the two different solvers.
        for (int i = 0; i < solverUnique1.nParticles; i++)
        {
            CHECK_EQUAL(solverUnique1.prOld[i][0], solverUnique2.prOld[i][0]);
        }

        // Do the same for importance sampling.
        solver.useImportanceSampling(true);
        solver.timeStep = 0.001;
        solver.D = 0.5;
        solver.useLocalEnergySlater();
        solver.useEfficientSlater(true);
        VMCSolver solverUnique3 = solver.getInitializedSolver();
        solver.useLocalEnergyGeneric();
        solver.useEfficientSlater(false);
        solver.useWaveFunctionBeryllium2();
        VMCSolver solverUnique4 = solver.getInitializedSolver();
        for (int i = 0; i < solverUnique1.nParticles; i++)
        {
            CHECK_EQUAL(solverUnique3.prOld[i][0], solverUnique4.prOld[i][0]);
        }
    }

    TEST(HeliumLocalEnergy)
    {
        VMCWrapper solver = VMCWrapper();
        solver.charge = 2;
        solver.alpha = 2;
        solver.nDimensions = 3;
        solver.nParticles = 2;
        solver.stepLength = 1.52;
        solver.nCycles = 1000;
        solver.h = 0.001;
        solver.hInv = 1e3;
        solver.h2Inv = 1e6;
        solver.idum = 1;
        solver.useWaveFunction1();
        VMCSolver solverUnique1 = solver.getInitializedSolver();
        VMCSolver solverUnique2 = solver.getInitializedSolver();
        // This will initialize the slater matrix and the initial
        // positions.
        solverUnique1.initRunVariables();
        solverUnique2.initRunVariables();
        // These two should be the same.
        double** r1 = solverUnique1.prOld;
        double** r2 = solverUnique2.prOld;
        for (int i = 0; i < solverUnique1.nParticles; i++)
        {
            CHECK_EQUAL(r1[i][0], r2[i][0]);
        }
        double* r1Abs = solverUnique1.rAbsOld;
        double* r2Abs = solverUnique2.rAbsOld;

        initWaveFunctions(solverUnique1);
        double localEnergy1 =
            wave_functions::getLocalEnergyGenericNoCor(r1, r1Abs);
        double localEnergy2 =
            wave_functions::getLocalEnergyGenericNoCor(r2, r2Abs);
        // Just to make sure things are equal.
        CHECK_EQUAL(localEnergy1, localEnergy2);
        localEnergy1 = wave_functions::getLocalEnergyGenericNoCor(r1, r1Abs);
        localEnergy2 = wave_functions::getLocalEnergyHeliumNoCor(r2, r2Abs);
        CHECK_CLOSE(localEnergy1, localEnergy2, 0.0001);
    }

    TEST(QuantumSlaterVsNormalGreensFunctionNoCor)
    {
        VMCWrapper solver = VMCWrapper();
        solver.charge = 4;
        solver.alpha = 4.6;
        solver.beta = 1;
        solver.nDimensions = 3;
        solver.nParticles = 4;
        solver.stepLength = 1.52;
        solver.nCycles = 1000;
        solver.h = 0.001;
        solver.hInv = 1e3;
        solver.h2Inv = 1e6;
        solver.idum = 1;

        solver.useImportanceSampling(true);
        solver.timeStep = 0.001;
        solver.D = 0.5;

        solver.useWaveFunctionBeryllium1();
        solver.useLocalEnergyGenericNoCor();
        VMCSolver solverUnique1 = solver.getInitializedSolver();
        solver.useEfficientSlater(true);
        solver.useLocalEnergySlaterNoCor();
        VMCSolver solverUnique2 = solver.getInitializedSolver();

        // Check that the ratios are the same for the normal step and the
        // efficient slater step.
        double rx1, rx2;
        int cycles = 10;
        int particles = 4;
        double tmp1;
        double tmp2;
        for (int i = 0; i < cycles; i++)
        {
            solverUnique1.startOfCycleQuantum();
            solverUnique2.startOfCycleSlaterQuantum();
            tmp1 = solverUnique1.pqForceOld[0][0];
            tmp2 = solverUnique2.pqForceOld[0][0];
            CHECK_CLOSE(tmp1, tmp2, 0.01);
            for (int j = 0; j < particles; j++)
            {
                // Check that the particles have the same positions.
                rx1 = solverUnique1.prNew[j][0];
                rx2 = solverUnique2.prNew[j][0];
                CHECK_CLOSE(rx1, rx2, 0.001);
                // Check the ratios.
                solverUnique1.runSingleStepQuantum(j, i);
                solverUnique2.runSingleStepSlaterQuantum(j, i);
                tmp1 = solverUnique1.pqForceOld[j][0];
                tmp2 = solverUnique2.pqForceOld[j][0];
                CHECK_CLOSE(tmp1, tmp2, 0.01);
            }
        }
    }

    TEST(SlaterVsNormalRatio)
    {
        VMCWrapper solver = VMCWrapper();
        solver.charge = 4;
        solver.alpha = 4.6;
        solver.beta = 1;
        solver.nDimensions = 3;
        solver.nParticles = 4;
        solver.stepLength = 1.52;
        solver.nCycles = 1000;
        solver.h = 0.001;
        solver.hInv = 1e3;
        solver.h2Inv = 1e6;
        solver.idum = 1;
        solver.useWaveFunctionBeryllium2();
        solver.useLocalEnergyGeneric();
        VMCSolver solverUnique1 = solver.getInitializedSolver();
        solver.useEfficientSlater(true);
        VMCSolver solverUnique2 = solver.getInitializedSolver();

        initWaveFunctions(solverUnique1);
        // Check that the ratios are the same for the normal step and the
        // efficient slater step.
        double ratio1, ratio2;
        int cycles = 10;
        int particles = 4;
        for (int cycle = 0; cycle < cycles; cycle++)
        {
            solverUnique1.startOfCycle();
            for (int i = 0; i < particles; i++)
            {
                solverUnique1.runSingleStep(i, cycle);
                ratio1 = solverUnique1.ratio;

                // check the second solver.
                solverUnique2.runSingleStepSlater(i, cycle);
                ratio2 = solverUnique2.ratio;
                CHECK_CLOSE(ratio1, ratio2, 0.000001);
            }
        }
    }

    TEST(SlaterVsNormalLocalEnergy)
    {
        VMCWrapper solver = VMCWrapper();
        solver.charge = 4;
        solver.alpha = 4.6;
        solver.beta = 1;
        solver.nDimensions = 3;
        solver.nParticles = 4;
        solver.stepLength = 1.52;
        solver.nCycles = 1000;
        solver.h = 0.001;
        solver.hInv = 1e3;
        solver.h2Inv = 1e6;
        solver.idum = 1;
        solver.useWaveFunctionBeryllium2();
        solver.useLocalEnergyGeneric();
        VMCSolver solverUnique1 = solver.getInitializedSolver();
        solver.useEfficientSlater(true);
        solver.useLocalEnergySlater();
        VMCSolver solverUnique2 = solver.getInitializedSolver();
        // This will initialize the slater matrix and the initial
        // positions.
        solverUnique1.initRunVariables();
        solverUnique2.initRunVariables();

        initWaveFunctions(solverUnique1);
        // Check that the ratios are the same for the normal step and the
        // efficient slater step.
        double energy1, energy2;
        int cycles = 1;
        int particles = 1;
        for (int i = 0; i < cycles; i++)
        {
            solverUnique1.startOfCycle();
            for (int j = 0; j < particles; j++)
            {
                solverUnique1.runSingleStep(j, i);
                solverUnique2.runSingleStepSlater(j, i);
                energy1 = solverUnique1.deltaE;
                energy2 = solverUnique2.deltaE;
                CHECK_CLOSE(energy1, energy2, 0.001);
            }
        }
    }

    TEST(SlaterVsNormalRatioNoCor)
    {
        VMCWrapper solver = VMCWrapper();
        solver.charge = 4;
        solver.alpha = 4.6;
        solver.nDimensions = 3;
        solver.nParticles = 4;
        solver.stepLength = 1.52;
        solver.nCycles = 1000;
        solver.h = 0.001;
        solver.hInv = 1e3;
        solver.h2Inv = 1e6;
        solver.idum = 1;
        solver.useWaveFunctionBeryllium1();
        solver.useEfficientSlater(true);
        solver.useLocalEnergyGenericNoCor();
        VMCSolver solverUnique1 = solver.getInitializedSolver();
        VMCSolver solverUnique2 = solver.getInitializedSolver();
        // This will initialize the slater matrix and the initial
        // positions.
        solverUnique1.initRunVariables();
        solverUnique2.initRunVariables();

        double** r = solverUnique2.prOld;
        double rAbs[4];
        for (int i = 0; i < 4; i++)
        {
            rAbs[i] = 0;
            for (int j = 0; j < 3; j++)
            {
                rAbs[i] += r[i][j] * r[i][j];
            }
            rAbs[i] = sqrt(rAbs[i]);
        }
        // Check if the slater matrix is correct.
        using namespace wave_functions;
        initWaveFunctions(solverUnique1);
        CHECK_EQUAL(phi1s(rAbs[0]), solverUnique1.pslater1Old[0][0]);
        CHECK_EQUAL(phi2s(rAbs[0]), solverUnique1.pslater1Old[0][1]);
        CHECK_EQUAL(phi1s(rAbs[1]), solverUnique1.pslater1Old[1][0]);
        CHECK_EQUAL(phi2s(rAbs[1]), solverUnique1.pslater1Old[1][1]);

        // Check that the ratios are the same for the normal step and the
        // efficient slater step.
        double ratio1, ratio2;
        int cycles = 10;
        int particles = 4;
        for (int i = 0; i < cycles; i++)
        {
            // Hack to update the first old wavefunction value.
            double** r = solverUnique1.prOld;
            double* rAbs = solverUnique1.rAbsOld;
            solverUnique1.waveFuncValOld =
                wave_functions::getWaveBerylliumNoCor(r, rAbs);
            for (int j = 0; j < particles; j++)
            {
                solverUnique1.runSingleStep(j, i);
                solverUnique2.runSingleStepSlater(j, i);
                ratio1 = solverUnique1.ratio;
                ratio2 = solverUnique2.ratio;
                CHECK_CLOSE(ratio1, ratio2, 0.000001);
            }
        }
    }

    TEST(SlaterVsNormalLocalEnergyNoCor)
    {
        VMCWrapper solver = VMCWrapper();
        solver.charge = 4;
        solver.alpha = 4.6;
        solver.nDimensions = 3;
        solver.nParticles = 4;
        solver.stepLength = 1.52;
        solver.nCycles = 1000;
        solver.h = 0.001;
        solver.hInv = 1e3;
        solver.h2Inv = 1e6;
        solver.idum = 1;
        solver.useWaveFunctionBeryllium1();
        solver.useLocalEnergyGenericNoCor();
        VMCSolver solverUnique1 = solver.getInitializedSolver();
        solver.useEfficientSlater(true);
        solver.useLocalEnergySlaterNoCor();
        VMCSolver solverUnique2 = solver.getInitializedSolver();
        // This will initialize the slater matrix and the initial
        // positions.
        solverUnique1.initRunVariables();
        solverUnique2.initRunVariables();

        initWaveFunctions(solverUnique1);
        // Check that the ratios are the same for the normal step and the
        // efficient slater step.
        double energy1, energy2;
        int cycles = 10;
        int particles = 4;
        for (int i = 0; i < cycles; i++)
        {
            // Hack to update the first old wavefunction value.
            double** r = solverUnique1.prOld;
            double* rAbs = solverUnique1.rAbsOld;
            solverUnique1.waveFuncValOld =
                wave_functions::getWaveBerylliumNoCor(r, rAbs);
            for (int j = 0; j < particles; j++)
            {
                solverUnique1.runSingleStep(j, i);
                solverUnique2.runSingleStepSlater(j, i);
                energy1 = solverUnique1.deltaE;
                energy2 = solverUnique2.deltaE;
                CHECK_CLOSE(energy1, energy2, 0.001);
            }
        }
    }

    TEST(BerylliumLocalEnergy)
    {
        VMCWrapper solver = VMCWrapper();
        solver.charge = 4;
        solver.alpha = 4.6;
        solver.nDimensions = 3;
        solver.nParticles = 4;
        solver.stepLength = 1.52;
        solver.nCycles = 1000;
        solver.h = 0.001;
        solver.hInv = 1e3;
        solver.h2Inv = 1e6;
        solver.idum = 1;
        solver.useWaveFunctionBeryllium1();
        solver.useEfficientSlater(true);
        solver.useLocalEnergySlaterNoCor();
        VMCSolver solverUnique1 = solver.getInitializedSolver();
        VMCSolver solverUnique2 = solver.getInitializedSolver();
        // This will initialize the slater matrix and the initial
        // positions.
        solverUnique1.initRunVariables();
        solverUnique2.initRunVariables();
        // These two should be the same.
        double** r1 = solverUnique1.prOld;
        double** r2 = solverUnique2.prOld;
        double* r1Abs = solverUnique1.rAbsOld;
        double* r2Abs = solverUnique2.rAbsOld;

        initWaveFunctions(solverUnique1);
        double localEnergy1 =
            wave_functions::getLocalEnergyGenericNoCor(r1, r1Abs);
        double localEnergy2 =
            solverUnique2.getLocalEnergySlaterNoCor(r2, r2Abs);
        CHECK_CLOSE(localEnergy1, localEnergy2, 0.0001);
    }

    TEST(BerylliumSlater)
    {
        VMCWrapper solver = VMCWrapper();
        solver.charge = 4;
        solver.alpha = 4;
        solver.nDimensions = 3;
        solver.nParticles = 4;
        solver.stepLength = 1.52;
        solver.nCycles = 1000;
        solver.h = 0.001;
        solver.hInv = 1e3;
        solver.h2Inv = 1e6;
        solver.idum = 1;
        solver.useWaveFunctionBeryllium1();
        solver.useEfficientSlater(true);
        VMCSolver solverUnique1 = solver.getInitializedSolver();
        // This will initialize the slater matrix and the initial
        // positions.
        solverUnique1.initRunVariables();
        // These two should be the same.
        double* r1 = solverUnique1.prOld[0];
        double* r2 = solverUnique1.prOld[1];
        double* r3 = solverUnique1.prOld[2];
        double* r4 = solverUnique1.prOld[3];

        double rAbs1 = 0;
        double rAbs2 = 0;
        double rAbs3 = 0;
        double rAbs4 = 0;

        for (int i = 0; i < 3; i++)
        {
            rAbs1 += r1[i] * r1[i];
            rAbs2 += r2[i] * r2[i];
            rAbs3 += r3[i] * r3[i];
            rAbs4 += r4[i] * r4[i];
        }
        rAbs1 = sqrt(rAbs1);
        rAbs2 = sqrt(rAbs2);
        rAbs3 = sqrt(rAbs3);
        rAbs4 = sqrt(rAbs4);

        using namespace wave_functions;
        initWaveFunctions(solverUnique1);

        // Check explicit wave functions against the slater matrix.
        CHECK_CLOSE(phi1s(rAbs1), solverUnique1.pslater1Old[0][0], 0.0001);
        CHECK_CLOSE(phi(0, r1, rAbs1), solverUnique1.pslater1Old[0][0], 0.0001);

        CHECK_CLOSE(phi2s(rAbs1), solverUnique1.pslater1Old[0][1], 0.0001);
        CHECK_CLOSE(phi(1, r1, rAbs1), solverUnique1.pslater1Old[0][1], 0.0001);

        CHECK_CLOSE(phi1s(rAbs2), solverUnique1.pslater1Old[1][0], 0.0001);
        CHECK_CLOSE(phi(0, r2, rAbs2), solverUnique1.pslater1Old[1][0], 0.0001);

        CHECK_CLOSE(phi2s(rAbs2), solverUnique1.pslater1Old[1][1], 0.0001);
        CHECK_CLOSE(phi(1, r2, rAbs2), solverUnique1.pslater1Old[1][1], 0.0001);

        // Second slater matrix.
        CHECK_CLOSE(phi1s(rAbs3), solverUnique1.pslater2Old[0][0], 0.0001);
        CHECK_CLOSE(phi(0, r3, rAbs3), solverUnique1.pslater2Old[0][0], 0.0001);

        CHECK_CLOSE(phi2s(rAbs3), solverUnique1.pslater2Old[0][1], 0.0001);
        CHECK_CLOSE(phi(1, r3, rAbs3), solverUnique1.pslater2Old[0][1], 0.0001);

        CHECK_CLOSE(phi1s(rAbs4), solverUnique1.pslater2Old[1][0], 0.0001);
        CHECK_CLOSE(phi(0, r4, rAbs4), solverUnique1.pslater2Old[1][0], 0.0001);

        CHECK_CLOSE(phi2s(rAbs4), solverUnique1.pslater2Old[1][1], 0.0001);
        CHECK_CLOSE(phi(1, r4, rAbs4), solverUnique1.pslater2Old[1][1], 0.0001);
    }

    TEST(Factorial)
    {
        double res;
        res = wave_functions_gto::factorial(2);
        CHECK_CLOSE(2,res, 0.00001);
        res = wave_functions_gto::factorial(3);
        CHECK_CLOSE(6,res, 0.00001);
    }

    TEST(phi)
    {
        VMCWrapper solver = VMCWrapper();
        solver.nDimensions = 3;
        solver.alpha = 3.75;
        wave_functions::alpha = 3.75;
        wave_functions::nDimensions = 3;
        double* r = new double[3];
        r[0] = 0.5;
        r[1] = 0.5;
        r[2] = 0.5;
        double rAbs = 0.8660254;
        // Check that the function phi(j,r_i) works.
        using namespace wave_functions;
        CHECK_CLOSE(phi1s(rAbs), phi(0, r, rAbs), 0.00001);
        CHECK_CLOSE(phi2s(rAbs), phi(1, r, rAbs), 0.00001);

        // Check values of phi.
        CHECK_CLOSE(0.0388675, phi1s(rAbs), 0.00001);
        CHECK_CLOSE(-0.122980, phi2s(rAbs), 0.00001);
        delete[] r;
    }
}

SUITE(Helium)
{
    VMCWrapper wrapper = VMCWrapper();
    double energy;

    TEST(Instantiate)
    {
        wrapper.charge = 2;
        wrapper.alpha = 1.66;
        wrapper.beta = 0.8;
        wrapper.nDimensions = 3;
        wrapper.nParticles = 2;
        wrapper.stepLength = 1.52;
        wrapper.nCycles = 10000;
        wrapper.waveFunction = wrapper.WAVE_FUNCTION_2;
        wrapper.h = 0.001;
        wrapper.hInv = 1e3;
        wrapper.h2Inv = 1e6;
        wrapper.idum = 1;
        wrapper.localEnergyFunction = wrapper.LOCAL_ENERGY_HELIUM_2;
    }

    TEST(h0LocalGenergic)
    {
        wrapper.alpha = 2;
        wrapper.waveFunction = wrapper.WAVE_FUNCTION_1;
        wrapper.localEnergyFunction = wrapper.LOCAL_ENERGY_GENERIC_NOCOR;
        double energy = 0;
        VMCSolver solver = wrapper.getInitializedSolver();
        // Run simulation.
        for (int cycle = 0; cycle < wrapper.nCycles; cycle++) 
        {
            for (int i = 0; i < wrapper.nParticles; i++) 
            {
                solver.runSingleStep(i,cycle);
                energy += solver.deltaE;
            }
        }
        energy /= (wrapper.nParticles*wrapper.nCycles);
        CHECK_CLOSE(-4, energy, 0.1);
    }

    TEST(h0LocalAnalytic)
    {
        wrapper.alpha = 2;
        wrapper.waveFunction = wrapper.WAVE_FUNCTION_1;
        wrapper.useLocalEnergyHelium1();
        double energy = 0;
        VMCSolver solver = wrapper.getInitializedSolver();
        // Run simulation.
        for (int cycle = 0; cycle < wrapper.nCycles; cycle++) 
        {
            for (int i = 0; i < wrapper.nParticles; i++) 
            {
                solver.runSingleStep(i,cycle);
                energy += solver.deltaE;
            }
        }
        energy /= (wrapper.nParticles*wrapper.nCycles);
        CHECK_CLOSE(-4, energy, 0.1);
    }

    TEST(WaveFunction2LocalEnergyGenergic)
    {
        wrapper.alpha = 1.66;
        wrapper.useWaveFunction2();
        wrapper.useLocalEnergyGeneric();
        double energy = 0;
        VMCSolver solver = wrapper.getInitializedSolver();
        // Run simulation.
        for (int cycle = 0; cycle < wrapper.nCycles; cycle++) 
        {
            for (int i = 0; i < wrapper.nParticles; i++) 
            {
                solver.runSingleStep(i,cycle);
                energy += solver.deltaE;
            }
        }
        energy /= (wrapper.nParticles*wrapper.nCycles);
        CHECK_CLOSE(-2.8, energy, 0.2);
    }

    TEST(WaveFunction2LocalEnergyAnalytic)
    {
        wrapper.alpha = 1.66;
        wrapper.useWaveFunction2();
        wrapper.useLocalEnergyHelium2();
        double energy = 0;
        VMCSolver solver = wrapper.getInitializedSolver();
        // Run simulation.
        for (int cycle = 0; cycle < wrapper.nCycles; cycle++) 
        {
            for (int i = 0; i < wrapper.nParticles; i++) 
            {
                solver.runSingleStep(i,cycle);
                energy += solver.deltaE;
            }
        }
        energy /= (wrapper.nParticles*wrapper.nCycles);
        CHECK_CLOSE(-2.8, energy, 0.2);
    }

    TEST(HeliumGTOLocalEnergyGeneric)
    {
        wrapper.useWaveFunctionHeliumGTO();
        wrapper.useLocalEnergyGeneric();
        double energy = 0;
        VMCSolver solver = wrapper.getInitializedSolver();
        // Run simulation.
        for (int cycle = 0; cycle < wrapper.nCycles; cycle++) 
        {
            for (int i = 0; i < wrapper.nParticles; i++) 
            {
                solver.runSingleStep(i,cycle);
                energy += solver.deltaE;
            }
        }
        energy /= (wrapper.nParticles*wrapper.nCycles);
        CHECK_CLOSE(-2.8, energy, 0.2);
    }
}

SUITE(Beryllium)
{
    VMCWrapper wrapper = VMCWrapper();
    TEST(Initialize)
    {
        wrapper.charge = 4;
        wrapper.alpha = 3.75;
        wrapper.beta = 0.8;
        wrapper.nDimensions = 3;
        wrapper.nParticles = 4;
        wrapper.stepLength = 1.52;
        wrapper.nCycles = 1000;
        wrapper.h = 0.001;
        wrapper.hInv = 1e3;
        wrapper.h2Inv = 1e6;
        wrapper.idum = 1;
    }

    TEST(h0LocalGeneric)
    {
        wrapper.waveFunction = wrapper.WAVE_FUNCTION_BERYLLIUM_1;
        wrapper.localEnergyFunction = wrapper.LOCAL_ENERGY_GENERIC_NOCOR;
        double energy = 0;
        VMCSolver solver = wrapper.getInitializedSolver();
        // Run simulation.
        for (int cycle = 0; cycle < wrapper.nCycles; cycle++) 
        {
            for (int i = 0; i < wrapper.nParticles; i++) 
            {
                solver.runSingleStep(i,cycle);
                energy += solver.deltaE;
            }
        }
        energy /= (wrapper.nParticles*wrapper.nCycles);
        CHECK_CLOSE(-20, energy, 0.5);
    }

    TEST(WaveFunction2LocalEnergyGeneric)
    {
        wrapper.useWaveFunctionBeryllium2();
        wrapper.useLocalEnergyGeneric();
        double energy = 0;
        VMCSolver solver = wrapper.getInitializedSolver();
        // Run simulation.
        for (int cycle = 0; cycle < wrapper.nCycles; cycle++) 
        {
            for (int i = 0; i < wrapper.nParticles; i++) 
            {
                solver.runSingleStep(i,cycle);
                energy += solver.deltaE;
            }
        }
        energy /= (wrapper.nParticles*wrapper.nCycles);
        CHECK_CLOSE(-14.3, energy, 0.5);
    }

    TEST(EfficientSlater)
    {
        wrapper.useEfficientSlater(true);
        wrapper.useLocalEnergyGeneric();
        wrapper.useWaveFunctionBeryllium2();
        double energy = 0;
        VMCSolver solver = wrapper.getInitializedSolver();
        // Run simulation.
        for (int cycle = 0; cycle < wrapper.nCycles; cycle++) 
        {
            for (int i = 0; i < wrapper.nParticles; i++) 
            {
                solver.runSingleStep(i,cycle);
                energy += solver.deltaE;
            }
        }
        energy /= (wrapper.nParticles*wrapper.nCycles);
        CHECK_CLOSE(-14.3, energy, 0.5);
    }
}

int main(int argc, const char* argv[])
{
    using namespace UnitTest;

    if (argc > 1)
    {
        // if first arg is "suite", we search for suite names instead of test
        // names
        const bool bSuite = strcmp("suite", argv[1]) == 0;

        // walk list of all tests, add those with a name that
        // matches one of the arguments  to a new TestList
        const TestList& allTests(Test::GetTestList());
        TestList selectedTests;
        Test* p = allTests.GetHead();
        while (p)
        {
            for (int i = 1; i < argc; ++i)
                if (strcmp(
                        bSuite ? p->m_details.suiteName : p->m_details.testName,
                        argv[i]) == 0)
                {
                    selectedTests.Add(p);
                }
            p = p->next;
        }

        // run selected test(s) only
        TestReporterStdout reporter;
        TestRunner runner(reporter);
        return runner.RunTestsIf(selectedTests, 0, True(), 0);
    }
    else
    {
        return RunAllTests();
    }
}
