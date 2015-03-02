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
    solver.setImportanceSampling();
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
    return 0;
}
