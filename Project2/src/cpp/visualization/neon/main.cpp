#include "../../vmcsolver/VMCWrapper.h"
#include "../../vmcsolver/util.h"
#include <chrono>
#include <iostream>
#include <omp.h>

using namespace std;

using microfortnights = std::chrono::duration<float, std::ratio<12096, 10000>>;
int main(int argc, const char* argv[])
{
    // Take alpha and beta from command line.
    if (argc != 3)
    {
        cout << "Not enough arguments." << endl;
        return -1;
    }
    double alpha = atof(argv[1]);
    double beta = atof(argv[2]);

    int nParticles = 10;
    double nCycles = 5e3;
    int threads = 4;
    int trials = 16;
    int totalTrials = threads * trials;
    int idum = 1000;
    int X = 100;
    int Y = 100;
    int Z = 100;

    string adress = "../../../../res/visualization/neon/";

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

    // Density plot in 3D. Along the columns are x, y, z, density
    Matrix mDensity3D = Matrix(X * Y * Z, 4);
    double** density3D = mDensity3D.getArrayPointer();

    // Set the max number of threads that can be run.
    omp_set_num_threads(threads);

    int xBin;
    int yBin;
    int zBin;
    double xMin = -2;
    double yMin = -2;
    double zMin = -2;
    double xMax = -xMin;
    double yMax = -yMin;
    double zMax = -zMin;
    double dx = (xMax-xMin)/X;
    double dy = (yMax-yMin)/Y;
    double dz = (zMax-zMin)/Z;
    auto start = chrono::high_resolution_clock::now();
    #pragma omp parallel
    {
        for (int trial = 0; trial < trials; trial++)
        {
            double energy = 0;

            VMCSolver solver = wrapper.getInitializedSolver();
            // Give an unique seed to the solver.
            int thread = omp_get_thread_num();
            solver.setSeed(wrapper.idum + trials * thread + trial);
            // Run simulation.
            for (int cycle = 0; cycle < nCycles; cycle++)
            {
                for (int i = 0; i < nParticles; i++)
                {
                    solver.runSingleStepSlater(i, cycle);
                    xBin = (solver.prOld[i][0] - xMin)*X / (xMax - xMin);
                    yBin = (solver.prOld[i][1] - yMin)*Y / (yMax - yMin);
                    zBin = (solver.prOld[i][2] - zMin)*Z / (zMax - zMin);
                    if (xBin < 0 || xBin >= X) continue;
                    if (yBin < 0 || yBin >= Y) continue;
                    if (zBin < 0 || zBin >= Z) continue;
                    density3D[xBin * Y * Z + Z * yBin + zBin][3] += 1;
                }
            }
        }
    }
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end - start;
    cout << totalTrials << " trials completed." << endl;
    cout << "Time = " << diff.count() << " seconds." << endl;

    for (int x = 0; x < X; x++)
    {
        for (int y = 0; y < Y; y++)
        {
            for (int z = 0; z < Z; z++)
            {
                density3D[x * Y * Z + y * Z + z][0] = x*dx + xMin;
                density3D[x * Y * Z + y * Z + z][1] = y*dy + yMin;
                density3D[x * Y * Z + y * Z + z][2] = z*dz + zMin;
            }
        }
    }

    util::writeToFile(adress, "density3D.txt", mDensity3D);

    return 0;
}
