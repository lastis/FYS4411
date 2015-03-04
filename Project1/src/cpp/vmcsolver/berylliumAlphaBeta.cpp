#include "VMCSolver.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main()
{
    VMCSolver solver = VMCSolver();
    solver.charge = 4;
    solver.alpha = 0; // To be set
    solver.beta = 0; // To be set
    solver.nDimensions = 3;
    solver.nParticles = 4;
    solver.stepLength = 1.40;
    solver.nCycles = 10000;
    solver.waveFunction = solver.WAVE_FUNCTION_BERYLLIUM_2;
    solver.h = 0.001;
    solver.h2 = 1e+06;
    solver.idum = 1;
    solver.localEnergyFunction = solver.LOCAL_ENERGY_GENERIC;

    int N = 11;

    Vector alpha = Vector(N);
    Vector beta = Vector(N);
    Matrix energy = Matrix(N,N);

    ofstream myFile;
    string base = "../../../res/";
    string directory = "berylliumAlphaBeta/";
    string fName;
    string adress;

    alpha.linspace(2,5);
    beta.linspace(0.1,2.1);
    double* pAlpha = alpha.getArrayPointer();
    double* pBeta = beta.getArrayPointer();
    double** pEnergy = energy.getArrayPointer();

    // Run simulation
    for (int i = 0; i < N; i++) {
    	for (int j = 0; j < N; j++) {
	    solver.alpha = pAlpha[i];
	    solver.beta = pBeta[j];
	    solver.runIntegration();
	    pEnergy[i][j] = solver.getEnergy();
	    cout << "Acceptance ratio : " << solver.getAcceptanceRatio() 
		<< endl;
	    cout << "Alpha : " << pAlpha[i] << endl;
	    cout << "Beta  : "<< pBeta[j] << endl;
    	}
    }

    // Store the data.
    fName = "alpha.txt";
    cout << "Creating plot file : " << directory + fName << endl;
    adress = base + directory + fName;
    myFile.open(adress.c_str());
    for (int i = 0; i < N; i++) {
    	myFile << pAlpha[i] << " ";
    }
    myFile.close();

    fName = "beta.txt";
    adress = base + directory + fName;
    myFile.open(adress.c_str());
    cout << "Creating plot file : " << directory + fName << endl;
    for (int i = 0; i < N; i++) {
    	myFile << pBeta[i] << " ";
    }
    myFile.close();

    fName = "energy.txt";
    cout << "Creating plot file : " << directory + fName << endl;
    adress = base + directory + fName;
    myFile.open(adress.c_str());
    for (int i = 0; i < N; i++) {
	for (int j = 0; j < N; j++) {
	    myFile << pEnergy[i][j] << " ";
	}
	myFile << endl;
    }
    myFile.close();
    return 0;
}

