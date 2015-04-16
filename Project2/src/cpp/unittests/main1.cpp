#include <unittest++/UnitTest++.h>
#include "../vmcsolver/VMCWrapper.h"

using namespace std;

SUITE(CPhys){
    TEST(Matrix){
        Matrix mat1;
        mat1 = Matrix(10,10);
        Matrix mat2 = Matrix(0,0);
        /* Matrix mat3 = Matrix(3,3); */
        /* mat1 = mat3; */
    }

    TEST(LUdecomposition){
        // LU decomposing mat into LU
        Matrix mat = Matrix(2,2);
        Matrix L = Matrix(mat);
        Matrix U = Matrix(mat);
        mat(0,0) = 4;
        mat(0,1) = 3;
        mat(1,0) = 6;
        mat(1,1) = 3;
        CPhys::MatOp::decomposeLU(mat,L,U);
        CHECK_CLOSE(1.5,L(1,0),0.0001);

        CHECK_CLOSE(4,U(0,0),0.0001);
        CHECK_CLOSE(3,U(0,1),0.0001);
        CHECK_CLOSE(-1.5,U(1,1),0.0001);
        Matrix res = L*U;

        CHECK_CLOSE(4,res(0,0),0.0001);
        CHECK_CLOSE(3,res(0,1),0.0001);
        CHECK_CLOSE(6,res(1,0),0.0001);
        CHECK_CLOSE(3,res(1,1),0.0001);

        Matrix B = Matrix(3,3);
        B(0,0) = 7;
        B(0,1) = 2;
        B(0,2) = 1;
        B(1,0) = 0;
        B(1,1) = 3;
        B(1,2) = -1;
        B(2,0) = -3;
        B(2,1) = 4;
        B(2,2) = -2;

        CPhys::MatOp::decomposeLU(B,L,U);
        Matrix res1 = L*U;
        CHECK_CLOSE(B(0,0),res1(0,0), 0.0001);
        CHECK_CLOSE(B(0,1),res1(0,1), 0.0001);
        CHECK_CLOSE(B(0,2),res1(0,2), 0.0001);
        CHECK_CLOSE(B(1,0),res1(1,0), 0.0001);
        CHECK_CLOSE(B(1,1),res1(1,1), 0.0001);
        CHECK_CLOSE(B(1,2),res1(1,2), 0.0001);
        CHECK_CLOSE(B(2,0),res1(2,0), 0.0001);
        CHECK_CLOSE(B(2,1),res1(2,1), 0.0001);
        CHECK_CLOSE(B(2,2),res1(2,2), 0.0001);
    }

    TEST(ForwardSubstitution){
        // Solving the equation:
        // LUx = b
        // Ly  = b
        // This means x is unity. 
        // With the given paramters, this means x is unity. 
        Matrix L = Matrix(2,2);
        L(0,0) = 1;
        L(0,1) = 0;
        L(1,0) = 1.5;
        L(1,1) = 1;
        Matrix b = Matrix(2,2);
        b(0,0) = 4;
        b(0,1) = 3;
        b(1,0) = 6;
        b(1,1) = 3;
        Matrix y = Matrix(2,2);
        CPhys::MatOp::substituteForward(L,y,b);
        CHECK_CLOSE(4,y(0,0),0.0001);
        CHECK_CLOSE(3,y(0,1),0.0001);
        CHECK_CLOSE(0,y(1,0),0.0001);
        CHECK_CLOSE(-1.5,y(1,1),0.0001);


        // Assuming LU decomposition works, do a decomposition
        // and solve the equation L*Linv = I. And test this.
        Matrix B = Matrix(3,3);
        B(0,0) = 7;
        B(0,1) = 2;
        B(0,2) = 1;
        B(1,0) = 0;
        B(1,1) = 3;
        B(1,2) = -1;
        B(2,0) = -3;
        B(2,1) = 4;
        B(2,2) = -2;
        Matrix U;
        CPhys::MatOp::decomposeLU(B,L,U);
        Matrix LInv;
        Matrix I = Matrix(L.getN(),L.getN());
        I.eye();
        CPhys::MatOp::substituteForward(L,LInv,I);
        Matrix res = L*LInv;

        CHECK_CLOSE(1,res(0,0),0.0001);
        CHECK_CLOSE(0,res(0,1),0.0001);
        CHECK_CLOSE(0,res(0,2),0.0001);

        CHECK_CLOSE(0,res(1,0),0.0001);
        CHECK_CLOSE(1,res(1,1),0.0001);
        CHECK_CLOSE(0,res(1,2),0.0001);

        CHECK_CLOSE(0,res(2,0),0.0001);
        CHECK_CLOSE(0,res(2,1),0.0001);
        CHECK_CLOSE(1,res(2,2),0.0001);
    }

    TEST(BackwardSubstitution){
        // Solving the equation:
        // LUx = b
        // Ly  = b
        // With the given paramters, this means x is unity. 
        Matrix U = Matrix(2,2);
        U(0,0) = 4;
        U(0,1) = 3;
        U(1,0) = 0;
        U(1,1) = -1.5;
        Matrix y = Matrix(2,2);
        y(0,0) = 4;
        y(0,1) = 3;
        y(1,0) = 0;
        y(1,1) = -1.5;
        Matrix x = Matrix(2,2);
        CPhys::MatOp::substituteBackward(U,x,y);
        CHECK_CLOSE(1,x(0,0),0.0001);
        CHECK_CLOSE(0,x(0,1),0.0001);
        CHECK_CLOSE(0,x(1,0),0.0001);
        CHECK_CLOSE(1,x(1,1),0.0001);

        // Assuming LU decomposition works, do a decomposition
        // and solve the equation U*Uinv = I. And test this.
        Matrix B = Matrix(3,3);
        B(0,0) = 7;
        B(0,1) = 2;
        B(0,2) = 1;
        B(1,0) = 0;
        B(1,1) = 3;
        B(1,2) = -1;
        B(2,0) = -3;
        B(2,1) = 4;
        B(2,2) = -2;
        Matrix L;
        Matrix UInv;
        Matrix I = Matrix(B.getN(),B.getN());
        I.eye();
        CPhys::MatOp::decomposeLU(B,L,U);
        CPhys::MatOp::substituteBackward(U,UInv,I);
        Matrix res = U*UInv;

        CHECK_CLOSE(1,res(0,0),0.0001);
        CHECK_CLOSE(0,res(0,1),0.0001);
        CHECK_CLOSE(0,res(0,2),0.0001);

        CHECK_CLOSE(0,res(1,0),0.0001);
        CHECK_CLOSE(1,res(1,1),0.0001);
        CHECK_CLOSE(0,res(1,2),0.0001);

        CHECK_CLOSE(0,res(2,0),0.0001);
        CHECK_CLOSE(0,res(2,1),0.0001);
        CHECK_CLOSE(1,res(2,2),0.0001);
    }

    TEST(Multiplication){
        Matrix L = Matrix(2,2);
        L(0,0) = 1;
        L(0,1) = 0;
        L(1,0) = 1.5;
        L(1,1) = 1;
        Matrix y = Matrix(2,2);
        y(0,0) = 4;
        y(0,1) = 3;
        y(1,0) = 0;
        y(1,1) = -1.5;
        Matrix res = L*y;
        CHECK_CLOSE(4,res(0,0),0.0001);
        CHECK_CLOSE(3,res(0,1),0.0001);
        CHECK_CLOSE(6,res(1,0),0.0001);
        CHECK_CLOSE(3,res(1,1),0.0001);
    }

    TEST(Inverse){
        // Check if we get the identity matrix if we multiply back.
        Matrix A = Matrix(2,2);
        A(0,0) = 4;
        A(0,1) = 3;
        A(1,0) = 6;
        A(1,1) = 3;
        Matrix AInv = CPhys::MatOp::getInverse(A);

        Matrix res1 = A*AInv;
        CHECK_CLOSE(1,res1(0,0),0.0001);
        CHECK_CLOSE(0,res1(0,1),0.0001);
        CHECK_CLOSE(0,res1(1,0),0.0001);
        CHECK_CLOSE(1,res1(1,1),0.0001);


        // Check if we get the identity matrix if we multiply back.
        Matrix B = Matrix(3,3);
        B(0,0) = 7;
        B(0,1) = 2;
        B(0,2) = 1;
        B(1,0) = 0;
        B(1,1) = 3;
        B(1,2) = -1;
        B(2,0) = -3;
        B(2,1) = 4;
        B(2,2) = -2;
        Matrix BInv = CPhys::MatOp::getInverse(B);

        Matrix res2 = B*BInv;
        CHECK_CLOSE(1,res2(0,0),0.0001);
        CHECK_CLOSE(0,res2(0,1),0.0001);
        CHECK_CLOSE(0,res2(0,2),0.0001);

        CHECK_CLOSE(0,res2(1,0),0.0001);
        CHECK_CLOSE(1,res2(1,1),0.0001);
        CHECK_CLOSE(0,res2(1,2),0.0001);

        CHECK_CLOSE(0,res2(2,0),0.0001);
        CHECK_CLOSE(0,res2(2,1),0.0001);
        CHECK_CLOSE(1,res2(2,2),0.0001);
    }
}

SUITE(Hydrogen){
    VMCWrapper solver = VMCWrapper();

    TEST(Instantiate){
        solver.charge = 1;
        solver.alpha = 1;
        solver.beta = 0;
        solver.nDimensions = 3;
        solver.nParticles = 1;
        solver.stepLength = 1.52;
        solver.nCycles = 10000;
        solver.waveFunction = solver.WAVE_FUNCTION_1;
        solver.h = 0.001;
        solver.h2 = 1e+06;
        solver.idum = 1;
        //TODO Check all variable names.
    }

    TEST(GroundStateGenergic){
        solver.waveFunction = solver.WAVE_FUNCTION_1;
        solver.localEnergyFunction = solver.LOCAL_ENERGY_GENERIC;
        solver.supressOutput();
        solver.runIntegration();
        double energy = solver.getEnergy();
        CHECK_CLOSE(-0.5,energy,0.01);
    }

    TEST(GroundStateAnalytic){
        solver.waveFunction = solver.WAVE_FUNCTION_1;
        solver.localEnergyFunction = solver.LOCAL_ENERGY_HYDROGEN;
        solver.supressOutput();
        solver.runIntegration();
        double energy = solver.getEnergy();
        CHECK_CLOSE(-0.5,energy,0.01);
    }
}

SUITE(VMCWrapper){
    TEST(phi){
        VMCWrapper solver = VMCWrapper();
        solver.nDimensions = 3;
        solver.alpha = 3.75;
        double* r = new double[3];
        r[0] = 0.5;
        r[1] = 0.5;
        r[2] = 0.5;
        double rAbs = 0.8660254;
        // Check that the function phi(j,r_i) works.
        using namespace wave_functions;
        CHECK_CLOSE(phi1s(rAbs),phi(0,r),0.00001);
        CHECK_CLOSE(phi2s(rAbs),phi(1,r),0.00001);

        // Check values of phi.
        CHECK_CLOSE(0.0388675, phi1s(rAbs), 0.00001);
        CHECK_CLOSE(-0.122980, phi2s(rAbs),0.00001);
    }

    TEST(SlaterDet){
        VMCWrapper solver = VMCWrapper();
        solver.charge = 4; 
        solver.alpha = 3.75;
        solver.beta = 0.8;
        solver.nDimensions = 3;
        solver.nParticles = 4;
        solver.stepLength = 1.52;
        solver.nCycles = 1000;
        solver.h = 0.001;
        solver.h2 = 1e+06;
        solver.idum = 1;
        solver.useWaveFunctionBeryllium2();
        solver.useLocalEnergyGeneric();
        solver.useEfficientSlater(true);

        /* solver.initRunVariables(); */
        /* double r1Abs = 0; */
        /* double r2Abs = 0; */
        /* double r3Abs = 0; */
        /* double r4Abs = 0; */
        /* for (int i = 0; i < solver.nDimensions; i++) { */
        /*     r1Abs += solver.prNew[0][i]*solver.prNew[0][i]; */
        /*     r2Abs += solver.prNew[1][i]*solver.prNew[1][i]; */
        /*     r3Abs += solver.prNew[2][i]*solver.prNew[2][i]; */
        /*     r4Abs += solver.prNew[3][i]*solver.prNew[3][i]; */
        /* } */
        /* r1Abs = sqrt(r1Abs); */
        /* r2Abs = sqrt(r2Abs); */
        /* r3Abs = sqrt(r3Abs); */
        /* r4Abs = sqrt(r4Abs); */
        /* // Check the first "Up" slater det. */
        /* CHECK_CLOSE( 0.33039, solver.phi(0,solver.prNew[0]), 0.0001); */
        /* CHECK_CLOSE( 0.256508, solver.phi(1,solver.prNew[0]), 0.0001); */
        /* CHECK_CLOSE( 0.367256, solver.phi(0,solver.prNew[1]), 0.0001); */
        /* CHECK_CLOSE( 0.302494, solver.phi(1,solver.prNew[1]), 0.0001); */

        /* CHECK_CLOSE( 52.7273, solver.pslater1Inv[0][0], 0.01); */
        /* CHECK_CLOSE(-44.7116, solver.pslater1Inv[0][1], 0.01); */
        /* CHECK_CLOSE(-64.0158, solver.pslater1Inv[1][0], 0.01); */
        /* CHECK_CLOSE( 57.5898, solver.pslater1Inv[1][1], 0.01); */

        // Check if the function updateInverse works. 
        Matrix AOld = Matrix(2,2);
        AOld(0,0) = 4;
        AOld(0,1) = 3;
        AOld(1,0) = 3;
        AOld(1,1) = 2;
        Matrix AOldInv = CPhys::MatOp::getInverse(AOld);
        double detAOld = -1;

        Matrix ANew = Matrix(AOld);
        ANew(0,0) = ANew(0,0) + 1;
        ANew(0,1) = ANew(0,1) + 1;
        Matrix ANewInv = CPhys::MatOp::getInverse(ANew);
        double detANew = ANew(0,0)*ANew(1,1)-ANew(1,0)*ANew(0,1);
        double ratio = detANew/detAOld;

        /* Matrix testMat = Matrix(AOldInv); */
        /* solver.updateInverse(0,ratio,ANew.getArrayPointer(),testMat.getArrayPointer()); */
        /* cout << "Ratio : " << ratio << endl; */

        /* CHECK_CLOSE(ANewInv(0,0),testMat(0,0), 0.0001); */
        /* CHECK_CLOSE(ANewInv(0,1),testMat(0,1), 0.0001); */
        /* CHECK_CLOSE(ANewInv(1,0),testMat(1,0), 0.0001); */
        /* CHECK_CLOSE(ANewInv(1,1),testMat(1,1), 0.0001); */
    }
}

SUITE(Helium){
    VMCWrapper solver = VMCWrapper();
    double energy;

    TEST(Instantiate){
	solver.charge = 2;
	solver.alpha = 1.66;
	solver.beta = 0.8;
	solver.nDimensions = 3;
	solver.nParticles = 2;
	solver.stepLength = 1.52;
	solver.nCycles = 10000;
	solver.waveFunction = solver.WAVE_FUNCTION_2;
	solver.h = 0.001;
	solver.h2 = 1e+06;
	solver.idum = 1;
	solver.localEnergyFunction = solver.LOCAL_ENERGY_HELIUM_2;
	//TODO Check all variable names.
    }

    TEST(h0LocalGenergic){
        solver.alpha = 2;
        solver.waveFunction = solver.WAVE_FUNCTION_1;
        solver.localEnergyFunction = solver.LOCAL_ENERGY_GENERIC_NOCOR;
        solver.supressOutput();
        solver.runIntegration();
        energy = solver.getEnergy();
        CHECK_CLOSE(-4,energy,0.1);
    }

    TEST(h0LocalAnalytic){
        solver.alpha = 2;
        solver.waveFunction = solver.WAVE_FUNCTION_1;
        solver.localEnergyFunction = solver.LOCAL_ENERGY_GENERIC_NOCOR;
        solver.supressOutput();
        solver.runIntegration();
        energy = solver.getEnergy();
        CHECK_CLOSE(-4,energy,0.1);
    }

    TEST(WaveFunc1LocalEnergyGeneric){
        solver.alpha = 1.66;
        solver.useWaveFunction1();
        solver.useLocalEnergyGeneric();
        solver.supressOutput();
        solver.runIntegration();
        energy = solver.getEnergy();
        CHECK_CLOSE(-2.8,energy,0.2);
    }

    TEST(WaveFunc1LocalEnergyAnalytic){
        solver.alpha = 1.66;
        solver.useLocalEnergyHelium1();
        solver.useWaveFunction1();
        solver.supressOutput();
        solver.runIntegration();
        energy = solver.getEnergy();
        CHECK_CLOSE(-2.9,energy,0.1);
    }

    TEST(WaveFunction2LocalEnergyGenergic){
        solver.alpha = 1.66;
        solver.useWaveFunction2();
        solver.useLocalEnergyGeneric();
        solver.supressOutput();
        solver.runIntegration();
        energy = solver.getEnergy();
        CHECK_CLOSE(-2.8,energy,0.2);
    }

    TEST(WaveFunction2LocalEnergyAnalytic){
        solver.alpha = 1.66;
        solver.useWaveFunction2();
        solver.useLocalEnergyHelium2();
        solver.supressOutput();
        solver.runIntegration();
        energy = solver.getEnergy();
        CHECK_CLOSE(-2.8,energy,0.2);
    }
}

SUITE(Beryllium){
    VMCWrapper solver = VMCWrapper();
    TEST(Initialize){
        solver.charge = 4; 
        solver.alpha = 3.75;
        solver.beta = 0.8;
        solver.nDimensions = 3;
        solver.nParticles = 4;
        solver.stepLength = 1.52;
        solver.nCycles = 1000;
        solver.h = 0.001;
        solver.h2 = 1e+06;
        solver.idum = 1;
    }

    TEST(h0LocalGeneric){
        solver.waveFunction = solver.WAVE_FUNCTION_BERYLLIUM_1;
        solver.localEnergyFunction = solver.LOCAL_ENERGY_GENERIC_NOCOR;
        solver.supressOutput();
        solver.runIntegration();
        double energy = solver.getEnergy();
        CHECK_CLOSE(-20,energy,0.5);
    }

    TEST(WaveFunction2LocalEnergyGeneric){
        solver.useWaveFunctionBeryllium2();
        solver.useLocalEnergyGeneric();
        solver.supressOutput();
        solver.runIntegration();
        double energy = solver.getEnergy();
        CHECK_CLOSE(-14.3,energy,0.5);
    }

    TEST(EfficientSlater){
        solver.useWaveFunctionBeryllium2();
        solver.useLocalEnergyGeneric();
        solver.useEfficientSlater(true);
        solver.supressOutput();
        solver.runIntegration();
        double energy = solver.getEnergy();
        CHECK_CLOSE(-14.3,energy,0.5);
    }
}

int main()
{
	return UnitTest::RunAllTests();
	return 0;
}
