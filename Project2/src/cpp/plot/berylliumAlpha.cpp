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
    solver.nCycles = 1e4;
    solver.h = 0.001;
    solver.h2 = 1e+06;
    solver.idum = 2;
    solver.useWaveFunctionBeryllium2();
    solver.useLocalEnergyGeneric();
    /* solver.useEfficientSlater(true); */
    /* solver.useParallel(true); */
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
    Vector energyArray = solver.getEnergyArray();
    Vector blockSizes = Vector();
    Vector stdArray = Vector();
    util::blockingVar(1,energyArray,stdArray,blockSizes);
    double* pEnergyArray = energyArray.getArrayPointer();
    double* pblockSizes = blockSizes.getArrayPointer();
    double* pstdArray = stdArray.getArrayPointer();



    ofstream myFile;
    string dir;
    string fName;
    string adress;
    dir = "berylliumAlpha/";



    // Dump raw energy file
    fName = "energy.dat";
    adress = "../../../res/plot/" + dir + fName;
    cout << "Dumption energies to file : " << adress << endl;
    myFile.open(adress.c_str());
    for (int i = 0; i < energyArray.getLength(); i++) {
        myFile << pEnergyArray[i] << " ";
    }
    myFile.close();
    


    // Dump variance
    fName = "energyStd.dat";
    adress = "../../../res/plot/" + dir + fName;
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
