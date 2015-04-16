#include "VMCWrapper.h"
#include <time.h>
#include <iostream>

using namespace std;

int main()
{
    VMCWrapper solver = VMCWrapper();
    /* solver.alpha = 3.4; */
    /* solver.D = 0.5; */
    /* solver.timeStep = 0.01; */
    /* solver.setImportanceSampling(true); */
    solver.alpha = 4;
    solver.beta = 0.8;
    solver.nDimensions = 3;
    solver.nParticles = 4;
    solver.charge = 4;
    solver.stepLength = 1.52;
    solver.nCycles = 5e4;
    solver.h = 0.001;
    solver.h2 = 1e+06;
    solver.idum = 2;
    solver.useWaveFunctionBeryllium1();
    solver.useLocalEnergyGenericNoCor();
    /* solver.useEfficientSlater(true); */
    solver.useParallel(true);

    clock_t start = clock();
    solver.runIntegration();
    clock_t end = clock();

    double time = double(end - start)/CLOCKS_PER_SEC;
    cout << "Time = " << time << endl;  

    /* cout << "Mean distance: " << solver.getR12Mean() << endl; */

    /* solver.exportEnergyArray("energies.txt"); */
    /* solver.exportDensity("density.txt"); */
    /* solver.exportParamters("test.ini"); */
    return 0;
}
