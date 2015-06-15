#include "../../vmcsolver/VMCWrapper.h"
#include "../../vmcsolver/util.h"
#include <chrono>
#include <iostream>
#include <omp.h>

using namespace std;

int main(int argc, const char *argv[])
{
    // Take alpha and beta from command line. 
    if (argc != 6) {
        cout << "Not enough arguments." << endl;
        return -1;
    }
    double alpha = atof(argv[1]);
    double beta = atof(argv[2]);
    double timeStep = atof(argv[3]);
    string fileNameBins = argv[4];
    string fileNameVariance = argv[5];

    int nParticles = 2;
    double nCycles = 1e4;
    int idum = 1;
    int threads = 4;

    string adress = "../../../../res/plot/blocking_01/";

    VMCWrapper wrapper = VMCWrapper();
    wrapper.alpha = alpha;
    wrapper.beta = beta;
    wrapper.nDimensions = 3;
    wrapper.nParticles = nParticles;
    wrapper.charge = nParticles;
    wrapper.stepLength = 1.52;
    wrapper.nCycles = nCycles;
    wrapper.h = 0.001;
    wrapper.hInv = 1000;
    wrapper.h2Inv = 1e+06;
    wrapper.idum = idum;
    wrapper.useWaveFunction2();
    wrapper.useLocalEnergyHelium2();
    wrapper.useImportanceSampling(true);
    wrapper.timeStep = timeStep;
    wrapper.D = 0.5;

    Vector vEnergyArray = Vector(nCycles);
    double* energyArray = vEnergyArray.getArrayPointer();

    // Set the max number of threads that can be run.
    omp_set_num_threads(threads);

    auto start = chrono::high_resolution_clock::now();
    #pragma omp parallel 
    {
        VMCSolver solver = wrapper.getInitializedSolver();
        solver.setSeed(idum + omp_get_thread_num());

        // Run simulation.
        for (int cycle = 0; cycle < nCycles; cycle++) 
        {
            solver.startOfCycleQuantum();
            for (int i = 0; i < nParticles; i++) 
            {
                solver.runSingleStepQuantum(i,cycle);
                energyArray[cycle] += solver.deltaE;
            }
        }
    }
    vEnergyArray /= (nParticles*threads);

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end-start;
    cout << "Time = " << diff.count() << " seconds." << endl;

    // Manipulate data. 
    Vector binSizes = Vector();
    Vector energyVariance = Vector();
    util::blockingVar(1,vEnergyArray,energyVariance,binSizes);
    util::appendToFile(adress,fileNameBins,binSizes);
    util::appendToFile(adress,fileNameVariance,energyVariance);

    return 0;
}
