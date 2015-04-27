#include "VMCWrapper.h"
#include "util.h"
#include <chrono>
#include <iostream>

using namespace std;

int main()
{
    VMCWrapper solver = VMCWrapper();
    /* solver.D = 0.5; */
    /* solver.timeStep = 0.01; */
    /* solver.setImportanceSampling(true); */
    /* solver.alpha = 3.7; */
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

    solver.useWaveFunctionBeryllium1();
    /* solver.localEnergyFunction = solver.LOCAL_ENERGY_GENERIC_NOCOR; */
    solver.useLocalEnergyGenericNoCor();
    solver.supressOutput();
    solver.runIntegration();

    auto start = chrono::high_resolution_clock::now();
    solver.runIntegration();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end-start;
    cout << "Energy : " << solver.getEnergy() <<endl;
    using microfortnights = std::chrono::duration<float, std::ratio<12096,10000>>;
    cout << "Time = " << microfortnights(diff).count() << " micro fortnights." << endl;
    cout << "Acceptance = " << solver.getAcceptanceRatio() << endl;
    return 0;
}
