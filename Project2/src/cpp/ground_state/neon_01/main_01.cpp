#include "../../vmcsolver/VMCWrapper.h"
#include "../../vmcsolver/util.h"
#include <chrono>
#include <iostream>
#include <omp.h>

using namespace std;

using microfortnights = std::chrono::duration<float, std::ratio<12096,10000>>;
int main(int argc, const char *argv[])
{
    int nParticles = 10;
    int nCycles = 1e5;
    int threads = 4;
    int idum = 1000;
    int skipps = 1000;
    int n_b = 1e3;
    int nBlocks = nCycles/n_b;

    string adress = "../../../../res/ground_state/neon_01/";
    string energyFileName = "energies_01.txt";

    VMCWrapper wrapper = VMCWrapper();
    wrapper.alpha = 10.2592;
    wrapper.beta = 0.11;
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

    Vector* vEnergiesMean = new Vector[threads];

    // Set the max number of threads that can be run.
    omp_set_num_threads(threads);

    auto start = chrono::high_resolution_clock::now();
    #pragma omp parallel 
    {
        int thread = omp_get_thread_num();
        vEnergiesMean[thread] = Vector(nBlocks);
        double* energies = vEnergiesMean[thread].getArrayPointer();

        VMCSolver solver = wrapper.getInitializedSolver();
        // Give an unique seed to the solver.
        solver.setSeed(wrapper.idum + thread);
        int block = 0;
        // Run simulation.
        for (int cycle = 0; cycle < nCycles; cycle++) 
        {
            solver.startOfCycleSlaterQuantum();
            for (int i = 0; i < nParticles; i++) 
            {
                solver.runSingleStepSlaterQuantum(i,cycle);
                if (cycle < skipps) continue;
                energies[block] += solver.deltaE;
            }
            if ((cycle+1) % n_b == 0)
            {
                block++;
            }
        }
    }
    Vector mergedVector = Vector(nBlocks*threads);
    double* tmp = mergedVector.getArrayPointer();
    for (int i = 0; i < threads; i++) 
    {
        for (int j = 0; j < nBlocks; j++) 
        {
            double* vector = vEnergiesMean[i].getArrayPointer();
            tmp[i*nBlocks + j] = vector[j]/(nParticles*n_b);
        }
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end-start;
    cout << "Time = " << diff.count() << " seconds." << endl;
    util::writeToFile(adress,energyFileName,mergedVector);
    /* util::writeToFile(adress,energyFileName,vEnergiesMean[0]); */
    delete vEnergiesMean;
    return 0;
}
