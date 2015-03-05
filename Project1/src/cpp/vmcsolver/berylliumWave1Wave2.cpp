#include "VMCSolver.h"
#include <time.h>
#include <iostream>

using namespace std;

int main()
{
    // Initialize the system 1.
    VMCSolver solver1 = VMCSolver();
    solver1.charge = 4;
    solver1.alpha = 3.75;
    solver1.beta = 0.5;
    solver1.nDimensions = 3;
    solver1.nParticles = 4;
    solver1.stepLength = 1.52;
    solver1.nCycles = 1000000;
    solver1.waveFunction = solver1.WAVE_FUNCTION_BERYLLIUM_1;
    solver1.h = 0.001;
    solver1.h2 = 1e+06;
    solver1.idum = 1;
    solver1.localEnergyFunction = solver1.LOCAL_ENERGY_GENERIC_NOCOR;

    solver1.setRecordEnergyArray(true);
    solver1.setRecordR12Mean(true);


    // Initialize the system 2.
    VMCSolver solver2 = VMCSolver();
    solver2.charge = 4; 
    solver2.alpha = 3.75;
    solver2.beta = 0.8;
    solver2.nDimensions = 3;
    solver2.nParticles = 4;
    solver2.stepLength = 1.52;
    solver2.nCycles = 1000000;
    solver2.waveFunction = solver2.WAVE_FUNCTION_BERYLLIUM_2;
    solver2.h = 0.001;
    solver2.h2 = 1e+06;
    solver2.idum = 2;
    solver2.localEnergyFunction = solver2.LOCAL_ENERGY_GENERIC;

    solver2.setRecordEnergyArray(true);
    solver2.setRecordR12Mean(true);

    // Data collection
    ofstream myFile;
    string base = "../../../res/";
    string directory = "berylliumWave1Wave2/";
    string fName;
    string adress;

    // Integration
    clock_t start = clock();
    solver1.runIntegration();
    clock_t end = clock();
    double time = double(end - start)/CLOCKS_PER_SEC;
    cout << "Time = " << time << endl;  
    cout << "Accepted moves: " << int(solver1.getAcceptanceRatio()) 
	<< " %" << endl;
    cout << "Mean distance: " << solver1.getR12Mean() << endl;

    // Store data.
    fName = "wave1Energies.txt";
    solver1.exportEnergyArray(directory + fName);
    fName = "wave1TimeAccpetanceRatio.txt";
    adress = base + directory + fName;
    myFile.open(adress.c_str());
    cout << "Creating other data file : " << "/res/" + directory + fName << endl;
    myFile.open(adress.c_str());
    myFile << time << " ";
    myFile << endl;
    myFile << solver1.getAcceptanceRatio() << " ";
    myFile << endl;
    myFile.close();

    // Store paramters.
    fName = "wave1Paramters.ini";
    solver1.exportParamters(directory + fName);

    // Again with system 2.
    // Integration
    start = clock();
    solver2.runIntegration();
    end = clock();
    time = double(end - start)/CLOCKS_PER_SEC;
    cout << "Time = " << time << endl;  
    cout << "Accepted moves: " << int(solver2.getAcceptanceRatio()) 
	<< " %" << endl;
    cout << "Mean distance: " << solver2.getR12Mean() << endl;

    // Store data.
    fName = "wave2Energies.txt";
    solver2.exportEnergyArray(directory + fName);
    fName = "wave2TimeAccpetanceRatio.txt";
    adress = base + directory + fName;
    myFile.open(adress.c_str());
    cout << "Creating other data file : " << "/res/" + directory + fName << endl;
    myFile.open(adress.c_str());
    myFile << time << " ";
    myFile << endl;
    myFile << solver2.getAcceptanceRatio() << " ";
    myFile << endl;
    myFile.close();

    // Store paramters.
    fName = "wave2Paramters.ini";
    solver2.exportParamters(directory + fName);

    return 0;
}
