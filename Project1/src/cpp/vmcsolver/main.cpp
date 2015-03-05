#include "VMCSolver.h"
#include <time.h>
#include <iostream>

using namespace std;

int main()
{
    VMCSolver solver = VMCSolver();
    /* solver.alpha = 3.4; */
    /* solver.D = 0.5; */
    /* solver.timeStep = 0.01; */
    /* solver.setImportanceSampling(true); */
    solver.charge = 2;
    solver.alpha = 1.66;
    solver.beta = 0.8;
    solver.nDimensions = 3;
    solver.nParticles = 2;
    solver.stepLength = 1.52;
    solver.nCycles = 1000000;
    solver.waveFunction = solver.WAVE_FUNCTION_1;
    solver.h = 0.001;
    solver.h2 = 1e+06;
    solver.idum = 1;
    solver.localEnergyFunction = solver.LOCAL_ENERGY_GENERIC_NOCOR;

    clock_t start = clock();
    solver.runIntegration();
    clock_t end = clock();

    double time = double(end - start)/CLOCKS_PER_SEC;
    cout << "Time = " << time << endl;  

    cout << "Accepted moves: " << int(solver.getAcceptanceRatio())
	    << " %" << endl;
    /* cout << "Mean distance: " << solver.getR12Mean() << endl; */

    /* solver.exportEnergyArray("energies.txt"); */
    /* solver.exportDensity("density.txt"); */
    return 0;
}
