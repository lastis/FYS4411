#include "../../vmcsolver/VMCWrapper.h"
#include "../../vmcsolver/util.h"
#include <chrono>
#include <iostream>
#include <omp.h>
#include <stdio.h>

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

    int nCycles = 1e4;
    int nParticles = 2;
    int threads = 4;
    int trials = 10;
    int totalTrials = threads*trials;
    int idum = 1000;

    string adress = "../../../../res/parameter_optimization/helium_01/";

    VMCWrapper wrapper = VMCWrapper();
    wrapper.charge = nParticles;
    wrapper.alpha = alpha;
    wrapper.beta = beta;
    wrapper.nDimensions = 3;
    wrapper.nParticles = nParticles;
    wrapper.stepLength = 1.52;
    wrapper.h = 0.001;
    wrapper.hInv = 1e3;
    wrapper.h2Inv = 1e6;
    wrapper.idum = 1;
    wrapper.useEfficientSlater(true);
    wrapper.useLocalEnergySlater();
    wrapper.useImportanceSampling(true);
    wrapper.timeStep = 0.001;
    wrapper.D = 0.5;

    // Set the max number of threads that can be run.
    omp_set_num_threads(threads);

    double energy = 0;
    double dE_dAlpha = 0;
    double dE_dBeta = 0;
    double factored1 = 0;
    double factored2 = 0;
    auto start = chrono::high_resolution_clock::now();
    #pragma omp parallel 
    {
        for (int trial = 0; trial < trials; trial++) 
        {
            int thread = omp_get_thread_num();

            VMCSolver solver = wrapper.getInitializedSolver();
            // Give an unique seed to the solver.
            solver.setSeed(wrapper.idum + trials*thread + trial);

            // Run simulation.
            for (int cycle = 0; cycle < nCycles; cycle++) 
            {
                solver.startOfCycleSlaterQuantum();
                for (int i = 0; i < nParticles; i++) 
                {
                    solver.runSingleStepSlaterQuantum(i,cycle);
                    solver.calc_dE_dAlpha();
                    solver.calc_dE_dBeta();
                    dE_dAlpha += solver.dE_dAlpha;
                    dE_dBeta += solver.dE_dBeta;
                    energy += solver.deltaE;
                    factored1 += solver.dE_dAlpha*solver.deltaE;
                    factored2 += solver.dE_dBeta*solver.deltaE;
                }
            }
        }
    }
    energy /= (nParticles*nCycles*totalTrials);
    dE_dAlpha /= (nParticles*nCycles*totalTrials);
    dE_dBeta /= (nParticles*nCycles*totalTrials);
    factored1 /= (nParticles*nCycles*totalTrials);
    factored2 /= (nParticles*nCycles*totalTrials);

    double retAlpha = 2*(factored1 - dE_dAlpha*energy);
    double retBeta = 2*(factored2 - dE_dBeta*energy);

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end-start;
    cout << totalTrials << " trials completed." << endl;
    cout << "Time = " << diff.count() << " seconds." << endl;
    cout << "Energy = " << energy 
        << " dE_dAlpha = " << retAlpha 
        << " dE_dBeta = " << retBeta  << endl;
    printf("%f %f %f", energy, retAlpha, retBeta);

    return 0;
}
