#include "vmcsolver.h"
#include <iostream>
#include <fstream>

using namespace std;

// Functions
void adjustStepLength(double deltaL = 0.5, 
        double epsilon = 2, double targetRatio = 50);


VMCSolver solver;

int main()
{
    solver = VMCSolver();
    solver.initFromFile("trial1.ini");
    adjustStepLength(0.01, 1);
    solver.exportParamters("trial1.ini");

    return 0;
}

void adjustStepLength(double deltaL, double epsilon, double targetRatio){
    // Find better value for stepLength
    double stepLength = solver.getStepLength();
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

