#include "../../vmcsolver/VMCWrapper.h"
#include "../../vmcsolver/util.h"
#include <chrono>
#include <iostream>

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
    double nCycles = 1e5;
    int binSize = nCycles/10;

    string adress = "../../../../res/plot/testplot/";

    VMCWrapper wrapper = VMCWrapper();
    wrapper.alpha = alpha;
    wrapper.beta = beta;
    wrapper.nDimensions = 3;
    wrapper.nParticles = nParticles;
    wrapper.charge = 4;
    wrapper.stepLength = 1.52;
    wrapper.nCycles = nCycles;
    wrapper.h = 0.001;
    wrapper.hInv = 1000;
    wrapper.h2Inv = 1e+06;
    wrapper.idum = 1;
    wrapper.useEfficientSlater(true);
    wrapper.useLocalEnergySlater();

    VMCSolver solver = wrapper.getInitializedSolver();

    Vector vEnergyArray = Vector(nCycles);
    Vector vDD = Vector(nCycles);
    Vector vCC = Vector(nCycles);
    Vector vDC = Vector(nCycles);
    Vector vEPot = Vector(nCycles);
    double* energyArray = vEnergyArray.getArrayPointer();
    double* DD = vDD.getArrayPointer();
    double* CC = vCC.getArrayPointer();
    double* DC = vDC.getArrayPointer();
    double* EPot = vEPot.getArrayPointer();


    // Run simulation.
    auto start = chrono::high_resolution_clock::now();
    for (int cycle = 0; cycle < nCycles; cycle++) 
    {
        for (int i = 0; i < nParticles; i++) 
        {
            solver.runSingleStepSlater(i,cycle);
            energyArray[cycle] += solver.deltaE;
            DD[cycle] += solver.DD;
            CC[cycle] += solver.CC;
            DC[cycle] += solver.DC;
            EPot[cycle] += solver.potentialEnergy;
        }
    }
    vDD /= nParticles;
    vCC /= nParticles;
    vDC /= nParticles;
    vEnergyArray /= nParticles;
    vEPot /= nParticles;
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end-start;
    cout << "Time = " << microfortnights(diff).count() << " micro fortnights." << endl;

    // Manipulate data. 
    Vector meanDD = util::getMeanArray(binSize,vDD);
    Vector meanCC = util::getMeanArray(binSize,vCC);
    Vector meanDC = util::getMeanArray(binSize,vDC);
    Vector meanArray = util::getMeanArray(binSize,vEnergyArray);
    Vector meanEPot = util::getMeanArray(binSize,vEPot);

    util::appendToFile(adress,"DD.txt",meanDD);
    util::appendToFile(adress,"CC.txt",meanCC);
    util::appendToFile(adress,"DC.txt",meanDC);
    util::appendToFile(adress,"energies_mean.txt",meanArray);
    util::appendToFile(adress,"energies_mean_potential.txt",meanEPot);

    return 0;
}
