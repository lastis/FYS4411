#include "../vmcsolver/VMCWrapper.h"
#include "../vmcsolver/util.h"
#include <chrono>
#include <iostream>

using namespace std;
int main()
{
    VMCWrapper solver = VMCWrapper();
    solver.alpha = 4;
    solver.beta = 0.8;
    solver.nDimensions = 3;
    solver.nParticles = 4;
    solver.charge = 4;
    solver.stepLength = 1.52;
    solver.nCycles = 1e5;
    solver.h = 0.001;
    solver.h2 = 1e+06;
    solver.idum = 2;
    solver.useWaveFunctionBeryllium2();
    solver.useLocalEnergyGeneric();
    solver.useEfficientSlater(true);
    solver.useParallel(true);
    solver.recordEnergyArray(true);

    // Run simulation.
    auto start = chrono::high_resolution_clock::now();
    solver.runIntegration();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end-start;
    cout << "Time = " << diff.count() << " seconds." << endl;
    using microfortnights = std::chrono::duration<float, std::ratio<12096,10000>>;
    cout << "Time = " << microfortnights(diff).count() << " micro fortnights." << endl;

    // Manipulate data. 
    int bins = 10;
    Vector energyArray = solver.getEnergyArray();
    Vector blockSizes = Vector();
    Vector stdArray = Vector();
    util::blockingVar(bins,energyArray,stdArray,blockSizes);
    double* pblockSizes = blockSizes.getArrayPointer();
    double* pstdArray = stdArray.getArrayPointer();

    string fName = "energyStd.dat";
    string dir = "berylliumAlpha/";
    string adress = "../../../res/plot/" + dir + fName;
    ofstream myFile;
    cout << "Dumption energies to file : " << adress << endl;
    myFile.open(adress.c_str());
    for (int i = 0; i < blockSizes.getLength(); i++) {
        myFile << pblockSizes[i] << " ";
    }
    myFile << endl;
    for (int i = 0; i < stdArray.getLength(); i++) {
        myFile << pstdArray[i] << " ";
    }
    myFile.close();

    return 0;
}
