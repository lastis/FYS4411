#include "../../vmcsolver/VMCWrapper.h"
#include "../../vmcsolver/util.h"
#include <chrono>
#include <iostream>

using namespace std;
int main(int argc, const char *argv[])
{
    // Arguments are, 1.output_file_name, 2.alpha, 3.beta 
    if (argc != 4) return -1;
    string fName = string(argv[1]);
    double alpha = atof(argv[2]);
    double beta = atof(argv[3]);

    // Output variables.
    string adress;
    adress = "../../../../res/plot/berylliumBlocking/"  + fName;

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
    /* solver.useParallel(true); */
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
    Vector blockSizes = Vector();
    Vector stdArray = Vector();
    // Do the blocking analysis.
    util::blockingVar(1,energyArray,stdArray,blockSizes);
    /* double* pEnergyArray = energyArray.getArrayPointer(); */
    double* pblockSizes = blockSizes.getArrayPointer();
    double* pstdArray = stdArray.getArrayPointer();

    // Dump variance
    ofstream myFile;
    myFile.open(adress.c_str());
    cout << "dumping to adress" << adress << endl;
    for (int i = 0; i < stdArray.getLength(); i++) {
        myFile << pstdArray[i] << " ";
    }
    myFile << endl;
    for (int i = 0; i < blockSizes.getLength(); i++) {
        myFile << pblockSizes[i] << " ";
    }
    myFile.close();

    return 0;
}
