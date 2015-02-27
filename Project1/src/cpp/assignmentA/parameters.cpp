#include "vmcsolver.h"
#include <iostream>
#include <fstream>

using namespace std;

// Functions
void adjustStepLength(double deltaL = 0.5, 
        double epsilon = 2, double targetRatio = 50);


// Paramters 
// NB! These paramters can be overwritten during the program.

VMCSolver solver;

int nDimensions = 3;
int charge 	= 2;
int nParticles = 2;
double stepLength = 2.5;
double h = 0.001;
double h2 = 1000000;
double alpha = 0.5*charge;
double beta = 1;
int nCycles = 1000000;
int waveFunctionType = 2;

int main()
{

    solver = VMCSolver();
    adjustStepLength(0.01, 1);
    solver.exportParamters();

    return 0;
}

void adjustStepLength(double deltaL, double epsilon, double targetRatio){
    // Initialize the solver from paramters in this file
    /* exportParamters(); */
    solver.initFromFile();

    // Find better value for stepLength
    solver.runMonteCarloIntegration();
    double ratio = solver.getStepAcceptance();
    cout << "\tRatio \t\t= " << ratio << " %" << endl;
    cout << "\tStep length \t= " <<  stepLength << endl;
    while (ratio < targetRatio - epsilon || ratio > targetRatio + epsilon) {
        cout << "Adjust step length: " << endl;
        if (ratio < targetRatio - epsilon)  stepLength -= deltaL;
        else                                stepLength += deltaL;
        solver.setStepLength(stepLength);
        solver.runMonteCarloIntegration();
        ratio = solver.getStepAcceptance();
        cout << "\tRatio \t\t= " << ratio << " %" << endl;
        cout << "\tStep length \t= " <<  stepLength << endl;
    }
}

