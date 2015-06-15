#include "../../vmcsolver/VMCWrapper.h"
#include "../../vmcsolver/util.h"
#include <chrono>
#include <iostream>
#include <omp.h>

using namespace std;

using microfortnights = std::chrono::duration<float, std::ratio<12096,10000>>;
int main(int argc, const char *argv[])
{

    int nParticles = 2;
    double nCycles = 1e6;
    int threads = 4;
    int trials = 1;
    int totalTrials = threads*trials;
    int idum = 1000;
    int bins = 100;
    double rMax = 3.0;

    string adress = "../../../../res/plot/helium_03/";
    string fileNameDensity = "density_01.txt";

    VMCWrapper wrapper = VMCWrapper();
    wrapper.alpha = 1.66;
    wrapper.beta = 0.5;
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
    wrapper.timeStep = 0.001;
    wrapper.D = 0.5;

    // Make a matrix with a density histogram for each trials.
    // This way we can get the variance of the density plot. 
    Matrix mDensity1D = Matrix(totalTrials, bins);
    double** density1D = mDensity1D.getArrayPointer();

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
            int uniqueID = trials*thread + trial;
            solver.setSeed(wrapper.idum + uniqueID);
            // Run simulation.
            for (int cycle = 0; cycle < nCycles; cycle++) 
            {
                solver.startOfCycleQuantum();
                for (int i = 0; i < nParticles; i++) 
                {
                    solver.runSingleStepQuantum(i,cycle);
                    int rBin = solver.rAbsOld[i]*bins/rMax;
                    if (rBin < 0 || rBin >= bins) continue;
                    density1D[uniqueID][rBin] += 1;
                }
            }
        }
    }

    // Normalize the histograms.
    for (int i = 0; i < totalTrials; i++) 
    {
        int sum = 0;
        for (int j = 0; j < bins; j++) 
        {
            sum += density1D[i][j];
        }
        for (int j = 0; j < bins; j++) 
        {
            density1D[i][j] /= sum;
        }
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end-start;
    cout << totalTrials << " trials completed." << endl;
    cout << "Time = " << diff.count() << " seconds." << endl;
    util::writeToFile(adress,fileNameDensity,mDensity1D);
    return 0;
}
