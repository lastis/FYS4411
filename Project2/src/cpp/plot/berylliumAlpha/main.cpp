#include "../../vmcsolver/VMCWrapper.h"
#include "../../vmcsolver/util.h"
#include <chrono>
#include <iostream>

using namespace std;
int main(int argc, const char *argv[])
{
    // Take alpha and beta from command line. 
    if (argc != 5) return -1;
    string fName = string(argv[1]);
    double alpha = atof(argv[2]);
    double beta = atof(argv[3]);
    int binSize = atof(argv[4]);

    ofstream myFile;
    /* string fName; */
    string adress;

    // Dump variance
    /* fName = "energies_mean.dat"; */
    adress = "../../../../res/plot/berylliumAlpha/" + fName;

    VMCWrapper solver = VMCWrapper();
    solver.alpha = alpha;
    solver.beta = beta;
    solver.nDimensions = 3;
    solver.nParticles = 4;
    solver.charge = 4;
    solver.stepLength = 1.52;
    solver.nCycles = 1e4;
    solver.h = 0.001;
    solver.h2 = 1e+06;
    solver.idum = 2;
    solver.useWaveFunctionBeryllium2();
    solver.useLocalEnergyGeneric();
    /* solver.useEfficientSlater(true); */
    solver.useParallel(true);
    solver.recordEnergyArray(true);

    // Run simulation.
    auto start = chrono::high_resolution_clock::now();
    solver.runIntegration();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end-start;
    /* cout << "Time = " << diff.count() << " seconds." << endl; */
    using microfortnights = std::chrono::duration<float, std::ratio<12096,10000>>;
    cout << "Time = " << microfortnights(diff).count() << " micro fortnights." << endl;

    // Manipulate data. 
    Vector energyArray = solver.getEnergyArray();
    Vector meanArray = util::getMeanArray(binSize,energyArray);
    double* mean = meanArray.getArrayPointer();
    // Dump results to the end of the file. 
    myFile.open(adress.c_str(),std::ios_base::app);
    for (int i = 0; i < meanArray.getLength(); i++) {
        myFile << mean[i] << " ";
    }
    myFile << endl;
    myFile.close();

    return 0;
}
