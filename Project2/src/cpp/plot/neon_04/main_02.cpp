#include "../../vmcsolver/VMCWrapper.h"
#include "../../vmcsolver/util.h"
#include <chrono>
#include <iostream>
#include <omp.h>

using namespace std;

int main(int argc, const char *argv[])
{
    int nParticles = 10;
    double nCycles = 1000;
    int threads = 4;
    int trials = 100;
    int totalTrials = threads*trials;
    int idum = 1000;

    string adress = "../../../../res/plot/neon_04/";
    string energyFileName = "energies_02.txt";

    VMCWrapper wrapper = VMCWrapper();
    wrapper.beta = 0.02;
    wrapper.nDimensions = 3;
    wrapper.nParticles = nParticles;
    wrapper.charge = nParticles;
    wrapper.stepLength = 1.52;
    wrapper.nCycles = nCycles;
    wrapper.h = 0.001;
    wrapper.hInv = 1000;
    wrapper.h2Inv = 1e+06;
    wrapper.idum = idum;
    wrapper.useWaveFunctionNeonGTO();
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
            VMCSolverGtoI solver = wrapper.getInitializedSolverGtoI();
            // Give an unique seed to the solver.
            int thread = omp_get_thread_num();
            int uniqueID = trials*thread + trial;
            solver.setSeed(wrapper.idum + uniqueID);
            // Run simulation.
            for (int cycle = 0; cycle < nCycles; cycle++) 
            {
                solver.startOfCycle();
                double energy = 0;
                for (int i = 0; i < nParticles; i++) 
                {
                    solver.runSingleStep(i);
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
