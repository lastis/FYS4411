#include "vmcsolver.h"
#include <iostream>
#include <fstream>

using namespace std;

// Functions
void adjustStepLength(double deltaL = 0.05, 
        double epsilon = 1, double targetRatio = 50);
void adjustAlpha(double deltaAlpha);
void adjustBeta(double deltaBeta);
void createAlphaBetaData(int N, double alphaStart, double alphaEnd, 
	double betaStart, double betaEnd);
void createAlphaData(int N, double alphaStart, double alphaEnd);
void createBetaData(int N, double betaStart, double betaEnd);
void calculateSolverVariance(int N);


VMCSolver solver;

int main()
{

    solver = VMCSolver();
    solver.initFromFile("helium2.ini");
    solver.useLocalEnergyHelium();
    createBetaData(21, 0.001, 0.321);
    return 0;
}

void createBetaData(int N, double betaStart, double betaEnd){
    Vector beta = Vector(N);
    Vector energy = Vector(N);

    beta.linspace(betaStart,betaEnd);
    double* pBeta = beta.getArrayPointer();
    double* pEnergy = energy.getArrayPointer();
    for (int i = 0; i < N; i++) {
	solver.beta = pBeta[i];
	solver.supressOutput();
	solver.runMonteCarloIntegration();
	pEnergy[i] = solver.getEnergy();
	cout << "Energy : " << pEnergy[i] << endl;
	cout << "Beta  : " << pBeta[i] << endl;
    }
    ofstream myFile;
    cout << "Creating plot file : " << "/res/betaPlot.txt" << endl;
    myFile.open("../../res/betaPlot.txt");
    for (int i = 0; i < N; i++) {
    	myFile << pBeta[i] << " ";
    }
    myFile << endl;
    for (int i = 0; i < N; i++) {
    	myFile << pEnergy[i] << " ";
    }
    myFile.close();
}

void createAlphaData(int N, double alphaStart, double alphaEnd){
    Vector alpha = Vector(N);
    Vector energy = Vector(N);

    alpha.linspace(alphaStart,alphaEnd);
    double* pAlpha = alpha.getArrayPointer();
    double* pEnergy = energy.getArrayPointer();
    for (int i = 0; i < N; i++) {
	solver.alpha = pAlpha[i];
	solver.supressOutput();
	solver.runMonteCarloIntegration();
	pEnergy[i] = solver.getEnergy();
	cout << "Energy : " << pEnergy[i] << endl;
	cout << "Alpha  : " << pAlpha[i] << endl;
    }
    ofstream myFile;
    cout << "Creating plot file : " << "/res/alphaPlot.txt" << endl;
    myFile.open("../../res/alphaPlot.txt");
    for (int i = 0; i < N; i++) {
    	myFile << pAlpha[i] << " ";
    }
    myFile << endl;
    for (int i = 0; i < N; i++) {
    	myFile << pEnergy[i] << " ";
    }
    myFile.close();
}

void calculateSolverVariance(int N){
    double tmp;
    double energy;
    double energysq;
    for (int i = 0; i < N; i++) {
	solver.idum++;
	solver.supressOutput();
	solver.runMonteCarloIntegration();
	tmp = solver.getEnergy();
	energy += tmp;
	energysq += tmp*tmp;

	cout << "Energy : " << tmp << endl;
    }
    energy /= N;
    energysq /= N;
    cout << "Variance : " << sqrt(energysq - energy*energy) << endl;
}

void createAlphaBetaData(int N, double alphaStart, double alphaEnd, 
	double betaStart, double betaEnd){
    Vector alpha = Vector(N);
    Vector beta = Vector(N);
    Matrix energy = Matrix(N,N);

    alpha.linspace(alphaStart,alphaEnd);
    beta.linspace(betaStart,betaEnd);
    double* pAlpha = alpha.getArrayPointer();
    double* pBeta = beta.getArrayPointer();
    double** pEnergy = energy.getArrayPointer();
    for (int i = 0; i < N; i++) {
    	for (int j = 0; j < N; j++) {
	    solver.alpha = pAlpha[i];
	    solver.beta = pBeta[j];
	    solver.runMonteCarloIntegration();
	    pEnergy[i][j] = solver.getEnergy();
    	}
    }
    ofstream myFile;
    cout << "Creating plot file : " << "/res/alphadata.txt" << endl;
    myFile.open("../../res/alphadata.txt");
    for (int i = 0; i < N; i++) {
    	myFile << pAlpha[i] << " ";
    }
    myFile.close();

    myFile.open("../../res/betadata.txt");
    cout << "Creating plot file : " << "/res/betadata.txt" << endl;
    for (int i = 0; i < N; i++) {
    	myFile << pBeta[i] << " ";
    }
    myFile.close();

    cout << "Creating plot file : " << "/res/energydata.txt" << endl;
    myFile.open("../../res/energydata.txt");
    for (int i = 0; i < N; i++) {
	for (int j = 0; j < N; j++) {
	    myFile << pEnergy[i][j] << " ";
	}
	myFile << endl;
    }
    myFile.close();
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
    cout << "Alpha : " << alphaNew << endl;

    // Itterate one time to find which way we must adjust alpha. 
    // New alpha value should give lower energy
    if (energyNew > energyOld) {
	goUp = false;
    	alphaNew = alphaOld - deltaAlpha;
	solver.runMonteCarloIntegration();
	energyNew = solver.getEnergy();
	cout << "Alpha : " << alphaNew << endl;
    }
    // If it doesn't, we are at a minima (With the current deltaAlpha).
    if (energyNew > energyOld) {
	cout << "Could not adjust alpha with the current adjustment step." 
	    << endl;
	return;
    }

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

