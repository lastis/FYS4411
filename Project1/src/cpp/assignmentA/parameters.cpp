#include "vmcsolver.h"
#include <iostream>
#include <fstream>

using namespace std;

int nDimensions = 3;
int charge 	= 2;

int main()
{
    ofstream myFile;
    myFile.open("main.ini");
    myFile << "nDimensions = " << nDimensions <<  endl;
    myFile << "charge = " << charge <<  endl;
    myFile << "nParticles = " << 2 <<  endl;
    myFile.close();
    VMCSolver solver = VMCSolver();
    solver.useWaveType1();
    solver.runMonteCarloIntegration();
    return 0;
}
