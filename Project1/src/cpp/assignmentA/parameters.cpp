#include "vmcsolver.h"
#include <iostream>
#include <fstream>

using namespace std;

// Functions
void exportParamters();
void adjustStepLength();

// Paramters 
// NB! These paramters can be overwritten during the program.
int nDimensions = 3;
int charge 	= 2;
int nParticles = 2;
double stepLength = 1;
double h = 0.001;
double h2 = 1000000;
double alpha = 0.5*charge;
double beta = 1;
int nCycles = 1000000;
int waveFunctionType = 1;

int main()
{

    adjustStepLength();
    exportParamters();

    return 0;
}

void adjustStepLength(){
    VMCSolver solver = VMCSolver();
    exportParamters();
    solver.initFromFile();
    solver.runMonteCarloIntegration();
    double ratio = solver.getStepAcceptance();
    while (ratio < 48 || ratio > 52) {
        if (ratio < 48) stepLength -= 0.5;
        else            stepLength += 0.5;
        solver.setStepLength(stepLength);
        solver.runMonteCarloIntegration();
        ratio = solver.getStepAcceptance();
        cout << "Adjust step length: " << endl;
        cout << "\tRatio \t= " << ratio << endl;
        cout << "\tStep length \t= " <<  stepLength << endl;
    }
    solver.runMonteCarloIntegration();
}


void exportParamters(){
    ofstream myFile;
    myFile.open("main.ini");
    myFile << "charge = " << charge <<  endl;
    myFile << "alpha = " << alpha << endl;
    myFile << "beta = " << beta << endl;
    myFile << "nDimensions = " << nDimensions <<  endl;
    myFile << "nParticles = " << nParticles <<  endl;
    myFile << "stepLength = " << stepLength << endl;
    myFile << "nCycles = " << nCycles << endl;
    myFile << "waveFunctionType = " << waveFunctionType << endl;
    myFile << "h = " << h << endl;
    myFile << "h2 = " << h2 << endl;
    myFile.close();

}
