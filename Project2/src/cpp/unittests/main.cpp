#include <unittest++/UnitTest++.h>
#include "../vmcsolver/VMCSolver.h"
/* #include "../CPhys/CPhys.h" */


SUITE(CPhys){
    TEST(Matrix){
        Matrix mat1;
        mat1 = Matrix(10,10);
        Matrix mat2 = Matrix(0,0);
        Matrix mat3 = Matrix(3,3);
        mat1 = mat3;
    }
	//TODO Check all variable names.
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
        /* double energy = solver.getEnergy(); */
        /* CHECK_CLOSE(-0.5,energy,0.01); */
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
