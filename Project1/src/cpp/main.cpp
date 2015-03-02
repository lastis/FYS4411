#include "vmcsolver.h"
#include <time.h>
#include <iostream>

using namespace std;

int main()
{
    VMCSolver solver = VMCSolver();
    solver.initFromFile("helium2.ini");

    solver.D = 0.5;
    solver.timeStep = 0.01;
    solver.setImportanceSampling(true);
    solver.setRecordDensity(true,21,2);
    solver.setWaveFunction2();
    solver.setLocalEnergyHelium();
    clock_t start = clock();
    solver.runIntegration();
    clock_t end = clock();

    double time = double(end - start)/CLOCKS_PER_SEC;
    cout << "Time = " << time << endl;  

    cout << "Accepted moves: " << int(solver.getAcceptanceRatio())
	    << " %" << endl;
    cout << "Mean distance: " << solver.getR12Mean() << endl;

    solver.exportDensity("density.txt");
    return 0;
}
