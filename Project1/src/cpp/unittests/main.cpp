#include <unittest++/UnitTest++.h>
#include "../vmcsolver/VMCSolver.h"



SUITE(Helium){
    VMCSolver solver;
    double energy;
    TEST(Instantiate){
	solver.initFromFile("unittest/helium2.ini");
	//TODO Check all variable names.
    }

    TEST(WaveFunc1LocalEnergyGeneric){
	solver.setWaveFunction1();
	solver.setLocalEnergyGeneric();
	solver.supressOutput();
	solver.runIntegration();
	energy = solver.getEnergy();
	CHECK_CLOSE(-2.8,energy,0.2);
    }

    TEST(WaveFunc1LocalEnergyAnalytic){
	solver.setLocalEnergyHelium();
	solver.setWaveFunction1();
	solver.supressOutput();
	solver.runIntegration();
	energy = solver.getEnergy();
	CHECK_CLOSE(-2.9,energy,0.1);
    }

    TEST(WaveFunction2LocalEnergyGenergic){
	solver.setWaveFunction2();
	solver.setLocalEnergyGeneric();
	solver.supressOutput();
	solver.runIntegration();
	energy = solver.getEnergy();
	CHECK_CLOSE(-2.8,energy,0.2);
    }

    TEST(WaveFunction2LocalEnergyAnalytic){
	solver.setWaveFunction2();
	solver.setLocalEnergyHelium();
	solver.supressOutput();
	solver.runIntegration();
	energy = solver.getEnergy();
	CHECK_CLOSE(-2.8,energy,0.2);
    }
}

int main()
{
	return UnitTest::RunAllTests();
	return 0;
}
