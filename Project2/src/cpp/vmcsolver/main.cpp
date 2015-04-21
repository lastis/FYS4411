#include "VMCWrapper.h"
#include <chrono>
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
    solver.nCycles = 1e5;
    solver.h = 0.001;
    solver.h2 = 1e+06;
    solver.idum = 2;
    solver.useWaveFunctionBeryllium1();
    solver.useLocalEnergyGenericNoCor();
    solver.useEfficientSlater(true);
    solver.useParallel(true);

    auto start = chrono::high_resolution_clock::now();
    solver.runIntegration();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end-start;
    cout << "Time = " << diff.count() << " seconds." << endl;
    using microfortnights = std::chrono::duration<float, std::ratio<12096,10000>>;
    cout << "Time = " << microfortnights(diff).count() << " micro fortnights." << endl;
    return 0;
}
