#include "../../vmcsolver/VMCWrapper.h"
#include "../../vmcsolver/util.h"
#include <chrono>
#include <iostream>
#include <omp.h>

using namespace std;

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
    double nCycles = 1000;
    int threads = 4;
    int trials = 200;
    int totalTrials = threads*trials;
    int idum = 1000;

    string adress = "../../../../res/plot/beryllium_04/";
    string energyFileName = "energies_01.txt";

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

    Matrix mEnergies = Matrix(totalTrials, nCycles);
    double** energies = mEnergies.getArrayPointer();

    // Set the max number of threads that can be run.
    omp_set_num_threads(threads);

    auto start = chrono::high_resolution_clock::now();
    #pragma omp parallel 
    {
        for (int trial = 0; trial < trials; trial++) 
        {
            VMCSolver solver = wrapper.getInitializedSolver();
            // Give an unique seed to the solver.
            int thread = omp_get_thread_num();
            int uniqueID = trials*thread + trial;
            solver.setSeed(wrapper.idum + uniqueID);
            // Run simulation.
            for (int cycle = 0; cycle < nCycles; cycle++) 
            {
                solver.startOfCycleSlaterQuantum();
                double energy = 0;
                for (int i = 0; i < nParticles; i++) 
                {
                    solver.runSingleStepSlaterQuantum(i,cycle);
                    energy += solver.deltaE;
                }
                energy /= nParticles;
                energies[uniqueID][cycle] = energy;
            }
        }
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end-start;
    cout << totalTrials << " trials completed." << endl;
    cout << "Time = " << diff.count() << " seconds." << endl;
    util::writeToFile(adress,energyFileName,mEnergies);
    return 0;
}
