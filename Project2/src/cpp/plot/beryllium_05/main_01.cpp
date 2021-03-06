#include "../../vmcsolver/VMCWrapper.h"
#include "../../vmcsolver/util.h"
#include <chrono>
#include <iostream>
#include <omp.h>

using namespace std;

int main(int argc, const char *argv[])
{

    int nParticles = 4;
    double nCycles = 1e5;
    int threads = 4;
    int trials = 10;
    int totalTrials = threads*trials;
    int idum = 1000;
    int bins = 100;
    double rMax = 5.0;

    string adress = "../../../../res/plot/beryllium_05/";
    string fileNameDensity = "density_02.txt";

    VMCWrapper wrapper = VMCWrapper();
    wrapper.alpha = 3.9190;
    wrapper.beta = 0.1101;
    wrapper.nDimensions = 3;
    wrapper.nParticles = nParticles;
    wrapper.charge = nParticles;
    wrapper.stepLength = 1.52;
    wrapper.nCycles = nCycles;
    wrapper.h = 0.001;
    wrapper.hInv = 1000;
    wrapper.h2Inv = 1e+06;
    wrapper.idum = idum;
    wrapper.useWaveFunctionBerylliumGTO();
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
            VMCSolverGtoI solver = wrapper.getInitializedSolverGtoI();
            // Give an unique seed to the solver.
            int thread = omp_get_thread_num();
            int uniqueID = trials*thread + trial;
            solver.setSeed(wrapper.idum + uniqueID);
            // Run simulation.
            for (int cycle = 0; cycle < nCycles; cycle++) 
            {
                solver.startOfCycle();
                for (int i = 0; i < nParticles; i++) 
                {
                    solver.runSingleStep(i);
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
