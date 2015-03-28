#include <unittest++/UnitTest++.h>
#include "../vmcsolver/VMCSolver.h"
/* #include "../CPhys/CPhys.h" */

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

    /* Test(Inverse){ */
    /*     Matrix AInv = Matrix(2,2); */
    /*     Matrix A = Matrix(2,2); */
    /*     A(0,0) = 4; */
    /*     A(0,1) = 3; */
    /*     A(1,0) = 6; */
    /*     A(1,1) = 3; */

    /*     Matrix res = A*AInv; */

    /*     CHECK_CLOSE(4,res(0,0),0.0001); */
    /*     CHECK_CLOSE(3,res(0,1),0.0001); */
    /*     CHECK_CLOSE(6,res(1,0),0.0001); */
    /*     CHECK_CLOSE(3,res(1,1),0.0001); */
    /* } */
}

SUITE(Hydrogen){
    VMCSolver solver = VMCSolver();

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

SUITE(Helium){
    VMCSolver solver = VMCSolver();
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
        solver.setWaveFunction1();
        solver.setLocalEnergyGeneric();
        solver.supressOutput();
        solver.runIntegration();
        energy = solver.getEnergy();
        CHECK_CLOSE(-2.8,energy,0.2);
    }

    TEST(WaveFunc1LocalEnergyAnalytic){
        solver.alpha = 1.66;
        solver.setLocalEnergyHelium1();
        solver.setWaveFunction1();
        solver.supressOutput();
        solver.runIntegration();
        energy = solver.getEnergy();
        CHECK_CLOSE(-2.9,energy,0.1);
    }

    TEST(WaveFunction2LocalEnergyGenergic){
        solver.alpha = 1.66;
        solver.setWaveFunction2();
        solver.setLocalEnergyGeneric();
        solver.supressOutput();
        solver.runIntegration();
        energy = solver.getEnergy();
        CHECK_CLOSE(-2.8,energy,0.2);
    }

    TEST(WaveFunction2LocalEnergyAnalytic){
        solver.alpha = 1.66;
        solver.setWaveFunction2();
        solver.setLocalEnergyHelium2();
        solver.supressOutput();
        solver.runIntegration();
        energy = solver.getEnergy();
        CHECK_CLOSE(-2.8,energy,0.2);
    }
}

SUITE(BERYLLIUM){
    VMCSolver solver = VMCSolver();
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
}

int main()
{
	return UnitTest::RunAllTests();
	return 0;
}
