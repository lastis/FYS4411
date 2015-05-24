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

    int nParticles = 4;
    double nCycles = 2e5;
    int idum = 1;

    string adress = "../../../../res/variableOptimazation/beryllium/";

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
    wrapper.useImportanceSampling(true);
    wrapper.timeStep = 0.001;
    wrapper.D = 0.5;
    wrapper.useEfficientSlater(true);
    wrapper.useLocalEnergySlater();

    VMCSolver solver = wrapper.getInitializedSolver();

    int threads = 4;
    // Set the max number of threads that can be run.
    omp_set_num_threads(threads);

    /* auto start = chrono::high_resolution_clock::now(); */
    double energy = 0;
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
            }
        }
        energy += solver.energySum / (nCycles * nParticles);
    }
    energy /= threads;
    cout << energy << endl;

    /* auto end = chrono::high_resolution_clock::now(); */
    /* chrono::duration<double> diff = end-start; */
    /* cout << "Time = " << microfortnights(diff).count() << " micro fortnights." << endl; */

    return 0;
}
