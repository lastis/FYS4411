#include "VMCSolver.h"
#include <time.h>
#include <iostream>

using namespace std;

// Descripttion: Find good values for alpha when using wave function 1. 
int main()
{
    // Initialize the system.
    VMCSolver solver = VMCSolver();
    solver.charge = 2;
    solver.alpha = 0; // To be set
    solver.beta = 0.75; 
    solver.nDimensions = 3;
    solver.nParticles = 2;
    solver.stepLength = 1.52;
    solver.nCycles = 1000000;
    solver.waveFunction = solver.WAVE_FUNCTION_1;
    solver.h = 0.001;
    solver.h2 = 1e+06;
    solver.idum = 1;
    solver.localEnergyFunction = solver.LOCAL_ENERGY_GENERIC;

    // Run with different alpha
    int N = 51;
    Vector alpha = Vector(N);
    Vector energy = Vector(N);
    Vector variance = Vector(N);
    Vector acceptedMoves = Vector(N);
    Vector times = Vector(N);
    alpha.linspace(0.1,3.1);

    for (int i = 0; i < N; i++) {
	solver.alpha = alpha(i);
	solver.supressOutput();
	clock_t start = clock();
	solver.runIntegration();
	clock_t end = clock();
	times(i) = double(end - start)/CLOCKS_PER_SEC;
	energy(i) = solver.getEnergy();
	variance(i) = solver.getEnergySquared() - energy(i)*energy(i);
	acceptedMoves(i) = solver.getAcceptanceRatio();
	cout << "Energy : " << energy(i) << endl;
	cout << "Alpha  : " << alpha(i) << endl;
    }


    ofstream myFile;
    string base = "../../../res/";
    string directory = "heliumWave2Alpha/";
    string fName;
    string adress;
    // Store the data.
    fName = "alpha_energy_variance.txt";
    adress = base + directory + fName;
    cout << "Creating plot file : " << "/res/" + fName << endl;
    myFile.open(adress.c_str());
    for (int i = 0; i < N; i++) {
    	myFile << alpha(i) << " ";
    }
    myFile << endl;
    for (int i = 0; i < N; i++) {
    	myFile << energy(i) << " ";
    }
    myFile << endl;
    for (int i = 0; i < N; i++) {
    	myFile << variance(i) << " ";
    }
    myFile.close();

    // Store other data.
    solver.exportParamters(directory + "paramters.ini");
    fName = "runtimes_acceptedMoveRatios.txt";
    adress = base + directory + fName;
    cout << "Creating other data file : " << "/res/" + directory + fName << endl;
    myFile.open(adress.c_str());
    for (int i = 0; i < N; i++) {
    	myFile << times(i) << " ";
    }
    myFile << endl;
    for (int i = 0; i < N; i++) {
    	myFile << acceptedMoves(i) << " ";
    }
    myFile.close();
    return 0;
}
