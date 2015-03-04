#include "VMCSolver.h"
#include <time.h>
#include <iostream>

using namespace std;

int main()
{
    // Initialize the system 1.
    VMCSolver solver1 = VMCSolver();
    solver1.charge = 2;
    solver1.alpha = 1.7;
    solver1.beta = 0.75;
    solver1.nDimensions = 3;
    solver1.nParticles = 2;
    solver1.stepLength = 1.52;
    solver1.nCycles = 1000000;
    solver1.waveFunction = solver1.WAVE_FUNCTION_1;
    solver1.h = 0.001;
    solver1.h2 = 1e+06;
    solver1.idum = 1;
    solver1.localEnergyFunction = solver1.LOCAL_ENERGY_HELIUM;

    solver1.setImportanceSampling(false);
    solver1.setRecordEnergyArray(true);
    solver1.setRecordR12Mean(true);


    // Initialize the system 2.
    VMCSolver solver2 = VMCSolver();
    solver2.charge = 2;
    solver2.alpha = 1.7;
    solver2.beta = 0.75;
    solver2.nDimensions = 3;
    solver2.nParticles = 2;
    solver2.stepLength = 1.52;
    solver2.nCycles = 1000000;
    solver2.waveFunction = solver2.WAVE_FUNCTION_2;
    solver2.h = 0.001;
    solver2.h2 = 1e+06;
    solver2.idum = 1;
    solver2.localEnergyFunction = solver2.LOCAL_ENERGY_HELIUM;

    solver2.setImportanceSampling(true);
    solver2.D = 0.5;
    solver2.timeStep = 0.01;
    solver2.setRecordEnergyArray(true);
    solver2.setRecordR12Mean(true);

    // Initialize the system 3.
    VMCSolver solver3 = VMCSolver();
    solver3.charge = 2;
    solver3.alpha = 1.7;
    solver3.beta = 0.75;
    solver3.nDimensions = 3;
    solver3.nParticles = 2;
    solver3.stepLength = 1.52;
    solver3.nCycles = 1000000;
    solver3.waveFunction = solver3.WAVE_FUNCTION_2;
    solver3.h = 0.001; // This paramter is changed
    solver3.h2 = 1e+06;
    solver3.idum = 1;
    solver3.localEnergyFunction = solver3.LOCAL_ENERGY_HELIUM;

    solver3.setImportanceSampling(true);
    solver3.D = 0.5;
    solver3.timeStep = 0.001;
    solver3.setRecordEnergyArray(true);
    solver3.setRecordR12Mean(true);
    // Data collection
    ofstream myFile;
    string base = "../../../res/";
    string directory = "heliumWave1Wave2Importance/";
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
    fName = "noImportanceEnergies.txt";
    solver1.exportEnergyArray(directory + fName);
    fName = "noImportanceTimeAccpetanceRatio.txt";
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
    fName = "noImportanceParamters.ini";
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
    fName = "importance1Energies.txt";
    solver2.exportEnergyArray(directory + fName);
    fName = "importance1TimeAccpetanceRatio.txt";
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
    fName = "importance1Paramters.ini";
    solver2.exportParamters(directory + fName);

    // Again with system 3.
    // Integration
    start = clock();
    solver3.runIntegration();
    end = clock();
    time = double(end - start)/CLOCKS_PER_SEC;
    cout << "Time = " << time << endl;  
    cout << "Accepted moves: " << int(solver3.getAcceptanceRatio()) 
	<< " %" << endl;
    cout << "Mean distance: " << solver3.getR12Mean() << endl;

    // Store data.
    fName = "importance2Energies.txt";
    solver3.exportEnergyArray(directory + fName);
    fName = "importance2TimeAccpetanceRatio.txt";
    adress = base + directory + fName;
    myFile.open(adress.c_str());
    cout << "Creating other data file : " << "/res/" + directory + fName << endl;
    myFile.open(adress.c_str());
    myFile << time << " ";
    myFile << endl;
    myFile << solver3.getAcceptanceRatio() << " ";
    myFile << endl;
    myFile.close();

    // Store paramters.
    fName = "importance2Paramters.ini";
    solver3.exportParamters(directory + fName);
    return 0;
}
