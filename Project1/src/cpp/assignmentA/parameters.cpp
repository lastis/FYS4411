#include "vmcsolver.h"
#include <iostream>
#include <fstream>

using namespace std;

// Functions
void adjustStepLength(double deltaL = 0.5, 
        double epsilon = 2, double targetRatio = 50);
void adjustAlpha(double deltaAlpha);
void adjustBeta(double deltaBeta);


VMCSolver solver;

int main()
{
    solver = VMCSolver();
    solver.initFromFile("trial2.ini");
    adjustBeta(1);
    solver.exportParamters("trial2.ini");

    return 0;
}

void adjustBeta(double deltaBeta){
    double deltaEnergy;

    // Set previous beta value
    double betaOld = solver.beta;
    solver.runMonteCarloIntegration();
    double energyOld = solver.getEnergy();
    cout << "Beta : " << betaOld << endl;

    // Set new beta value
    double betaNew = betaOld + deltaBeta;
    bool goUp = true;
    solver.runMonteCarloIntegration();
    double energyNew = solver.getEnergy();
    cout << "Beta : " << betaNew << endl;

    // Itterate one time to find which way we must adjust beta. 
    // New beta value should give lower energy
    if (energyNew > energyOld) {
	goUp = false;
    	betaNew = betaOld - deltaBeta;
	solver.runMonteCarloIntegration();
	energyNew = solver.getEnergy();
	cout << "Beta : " << betaNew << endl;
    }
    // If it doesn't, we are at a minima (With the current deltaAlpha).
    if (energyNew > energyOld) {
	cout << "Could not adjust beta with the current adjustment step." 
	    << endl;
	return;
    }

    // Adjust as long as we have a lower energy.
    while (energyNew < energyOld) {
	// Copy new to old values
	betaOld = betaNew;
	energyOld = energyNew;
	// Update new values
	if(goUp) betaNew = betaOld + deltaBeta;
	else betaNew = betaOld - deltaBeta;
	solver.runMonteCarloIntegration();
	energyNew = solver.getEnergy();
	cout << "Beta : " << betaNew << endl;
    }
    
    solver.beta = betaOld;
}

void adjustAlpha(double deltaAlpha){
    double deltaEnergy;

    // Set previous alpha value
    double alphaOld = solver.alpha;
    solver.runMonteCarloIntegration();
    double energyOld = solver.getEnergy();
    cout << "Alpha : " << alphaOld << endl;

    // Set new alpha value
    double alphaNew = alphaOld + deltaAlpha;
    bool goUp = true;
    solver.runMonteCarloIntegration();
    double energyNew = solver.getEnergy();

    // Itterate one time to find which way we must adjust alpha. 
    // New alpha value should give lower energy
    if (energyNew > energyOld) {
	goUp = false;
    	alphaNew = alphaOld - deltaAlpha;
	solver.runMonteCarloIntegration();
	energyNew = solver.getEnergy();
    }
    // If it doesn't, we are at a minima (With the current deltaAlpha).
    if (energyNew > energyOld) {
	cout << "Could not adjust alpha with the current adjustment step." 
	    << endl;
	return;
    }
    cout << "Alpha : " << alphaNew << endl;

    // Adjust as long as we have a lower energy.
    while (energyNew < energyOld) {
	// Copy new to old values
	alphaOld = alphaNew;
	energyOld = energyNew;
	// Update new values
	if(goUp) alphaNew = alphaOld + deltaAlpha;
	else alphaNew = alphaOld - deltaAlpha;
	solver.runMonteCarloIntegration();
	energyNew = solver.getEnergy();
	if (energyNew > energyOld) cout << "Alpha : " << alphaNew << endl;
    }
    
    solver.alpha = alphaOld;
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

