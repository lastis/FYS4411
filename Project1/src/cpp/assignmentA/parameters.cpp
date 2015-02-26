#include "vmcsolver.h"
#include <iostream>
#include <fstream>

using namespace std;

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
    // TODO Varify that the values are written "exactly" so they can be read
    // correctly. 
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

    VMCSolver solver = VMCSolver();
    solver.runMonteCarloIntegration();
    return 0;
}
