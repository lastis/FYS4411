#include "../../vmcsolver/VMCWrapper.h"
#include "../../vmcsolver/util.h"
#include <chrono>
#include <iostream>
#include <omp.h>

using namespace std;

using microfortnights = std::chrono::duration<float, std::ratio<12096,10000>>;
int main(int argc, const char *argv[])
{
    // Take alpha and beta from command line. 
    if (argc != 3) {
        cout << "Not enough arguments." << endl;
        return -1;
    }
    double alpha = atof(argv[1]);
    double beta = atof(argv[2]);

    int nParticles = 10;
    double nCycles = 1e4;
    int binSize = nCycles;
    int idum = 1;

    string adress = "../../../../res/blocking/neon_prep/";

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
    wrapper.useEfficientSlater(true);
    wrapper.useLocalEnergySlater();
    wrapper.useImportanceSampling(true);
    wrapper.timeStep = 0.001;
    wrapper.D = 0.5;

    VMCSolver solver = wrapper.getInitializedSolver();

    Vector vEnergyArray = Vector(nCycles);
    double* energyArray = vEnergyArray.getArrayPointer();

    int threads = 4;
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
            for (int i = 0; i < nParticles; i++) 
            {
                solver.runSingleStepSlater(i,cycle);
                energyArray[cycle] += solver.deltaE;
            }
        }
    }
    vEnergyArray /= (nParticles*threads);

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end-start;
    /* cout << "Time = " << microfortnights(diff).count() << " micro fortnights." << endl; */
    cout << "Time = " << diff.count() << " seconds." << endl;

    // Manipulate data. 
    Vector blocking = Vector();
    Vector blockingStd = Vector();
    util::blockingVar(1,vEnergyArray,blockingStd,blocking);
    util::appendToFile(adress,"energies_blocking.txt",blocking);
    util::appendToFile(adress,"energies_blocking_std.txt",blockingStd);

    return 0;
}
