#include "VMCSolver.h"
#include <time.h>
#include <iostream>

using namespace std;

// Descripttion: Find good values for beta when using wave function 1. 
ofstream myFile;
VMCSolver solver;
int N;
void run();
Vector beta;
Vector energy;
Vector variance;
Vector acceptedMoves;
Vector times;

string base = "../../../res/";
string directory = "heliumWave2Beta/";
string fName;
string adress;

int main()
{
    // Initialize the system.
    solver = VMCSolver();
    solver.charge = 2;
    solver.alpha = 1.65; 
    solver.beta = 0; // To be set
    solver.nDimensions = 3;
    solver.nParticles = 2;
    solver.stepLength = 1.52;
    solver.nCycles = 1000000;
    solver.waveFunction = solver.WAVE_FUNCTION_2;
    solver.h = 0.001;
    solver.h2 = 1e+06;
    solver.idum = 1; // Will be varied.
    solver.localEnergyFunction = solver.LOCAL_ENERGY_GENERIC;

    // Run with different beta
    N = 51; // Sample points
    beta = Vector(N);
    energy = Vector(N);
    variance = Vector(N);
    acceptedMoves = Vector(N);
    times = Vector(N);
    beta.linspace(0.1,2.1);
    for (int i = 0; i < 100; i++) {
    	run();
    }
    return 0;
}

void run(){
    for (int i = 0; i < N; i++) {
	solver.beta = beta(i);
	solver.supressOutput();
	clock_t start = clock();
	solver.runIntegration();
	clock_t end = clock();
	times(i) = double(end - start)/CLOCKS_PER_SEC;
	energy(i) = solver.getEnergy();
	variance(i) = solver.getEnergySquared() - energy(i)*energy(i);
	acceptedMoves(i) = solver.getAcceptanceRatio();
	cout << "Energy : " << energy(i) << endl;
	cout << "Beta  : " << beta(i) << endl;
    }


    // Store the data.
    fName = "beta.txt";
    adress = base + directory + fName;
    myFile.open(adress.c_str());
    for (int i = 0; i < N; i++) {
    	myFile << beta(i) << " ";
    }
    myFile << endl;
    myFile.close();

    // Energy store.
    fName = "energy.txt";
    adress = base + directory + fName;
    myFile.open(adress.c_str(), ios::out | ios::app);
    for (int i = 0; i < N; i++) {
    	myFile << energy(i) << " ";
    }
    myFile << endl;
    myFile.close();

    // Variance store.
    fName = "variance.txt";
    adress = base + directory + fName;
    myFile.open(adress.c_str(), ios::out | ios::app);
    for (int i = 0; i < N; i++) {
    	myFile << variance(i) << " ";
    }
    myFile << endl;
    myFile.close();

    // Store other data.
    solver.exportParamters(directory + "paramters.ini");
    fName = "runtimes.txt";
    adress = base + directory + fName;
    myFile.open(adress.c_str(), ios::out | ios::app);
    for (int i = 0; i < N; i++) {
    	myFile << times(i) << " ";
    }
    myFile << endl;
    myFile.close();

    fName = "acceptedMoveRatios.txt";
    adress = base + directory + fName;
    myFile.open(adress.c_str(), ios::out | ios::app);
    for (int i = 0; i < N; i++) {
    	myFile << acceptedMoves(i) << " ";
    }
    myFile << endl;
    myFile.close();
}

