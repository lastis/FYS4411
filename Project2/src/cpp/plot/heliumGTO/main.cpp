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
    if (argc != 2) {
        cout << "Not enough arguments." << endl;
        return -1;
    }
    double beta = atof(argv[1]);
    int nCycles = 1e4;
    int nParticles = 2;
    int threads = 4;
    int trials = 4;
    int totalTrials = threads*trials;
    int idum = 1000;
    Vector vEnergiesMean = Vector(totalTrials);
    double* energiesMean = vEnergiesMean.getArrayPointer();

    string adress = "../../../../res/plot/heliumGTO/";
    VMCWrapper wrapper = VMCWrapper();
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
    wrapper.useWaveFunctionHeliumGTO();
    wrapper.useLocalEnergyGeneric();

    // Set the max number of threads that can be run.
    omp_set_num_threads(threads);
    auto start = chrono::high_resolution_clock::now();
    #pragma omp parallel 
    {
        for (int trial = 0; trial < trials; trial++) 
        {
            double energy = 0;
            VMCSolver solver = wrapper.getInitializedSolver();
            // Give an unique seed to the solver.
            int thread = omp_get_thread_num();
            solver.setSeed(wrapper.idum + trials*thread + trial);
            // Run simulation.
            for (int cycle = 0; cycle < nCycles; cycle++) 
            {
                for (int i = 0; i < nParticles; i++) 
                {
                    solver.runSingleStep(i,cycle);
                    energy += solver.deltaE;
                }
            }
            energy /= (nParticles*nCycles);
            energiesMean[trials*thread + trial] = energy;
        }
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end-start;
    cout << totalTrials << " trials completed." << endl;
    cout << "Time = " << diff.count() << " seconds." << endl;

    util::appendToFile(adress,"energies_mean.txt",vEnergiesMean);

    return 0;
}
