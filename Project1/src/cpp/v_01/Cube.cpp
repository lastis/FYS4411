#include <iostream>
#include <cstdio>
#include <cmath>
#include <math.h>
#include "Cube.h"

using namespace std;

Cube::Cube(){
	mN = 0;
	mM = 0;
	mO = 0;
	allocateMemory(mN,mM,mO);
}


Cube::Cube(int N, int M, int O){
	allocateMemory(N,M,O);
	mN = N;
	mM = M;
	mO = O;
}

Cube::Cube(const Cube& mat){
	// Note: copy constructor is needed for the assignment 
	// operator (operator=(Cube other)) "The rule of three"
	mN = mat.mN;
	mM = mat.mM;
	mO = mat.mO;
	allocateMemory(mN,mM,mO);
	for (int i = 0; i < mN; i++) {
		for (int j = 0; j < mM; j++) {
			for (int k = 0; k < mO; k++) {
				mMat[i][j][k] = mat.mMat[i][j][k];
			}
		}
	}
}


Cube&	Cube::operator =(double num){
	for (int i = 0; i < mN; i++) {
		for (int j = 0; j < mM; j++) {
			for (int k = 0; k < mO; k++) {
				mMat[i][j][k] = num;
			}
		}
	}
	return *this;
}

Cube& Cube::operator =(Cube other){
	// We need to do a deep copy so we don't lose
	// our dynamic array (pointer) in the switch;
	swap(*this, other);
	
	return *this;
}

void	Cube::swap(Cube& c1, Cube& c2){
	std::swap(c1.mN  , c2.mN);
	std::swap(c1.mM  , c2.mM);
	std::swap(c1.mO  , c2.mO);
	std::swap(c1.mMat, c2.mMat);
}

double&	Cube::operator()(int i, int j, int k){
	// Make the array circular by taking the modulus
	/*
	cout << "Boxes = " << mN << endl;
	cout << "i = " << i << endl;
	cout << "j = " << j << endl;
	cout << "k = " << k << endl;
	*/
	i -= floor(double(i)/mN)*mN;
	j -= floor(double(j)/mM)*mM;
	k -= floor(double(k)/mO)*mO;
	/*
	cout << "i = " << i << endl;
	cout << "j = " << j << endl;
	cout << "k = " << k << endl;
	*/
	return mMat[i][j][k];
}


double*** Cube::getArrayPointer(){
	return mMat;
}

int	Cube::getN(){
	return mN;
}	

int 	Cube::getM(){
	return mM;
}	

int 	Cube::getO(){
	return mO;
}	

void 	Cube::reset(){
	// Set every element to 0
	for (int i = 0; i < mN; i++) {
		for (int j = 0; j < mM; j++) {
			for (int k = 0; k < mO; k++) {
				mMat[i][j][k] = 0;
			}
		}
	}
}

void 	Cube::freeMemory(){
	//delete[] mMat[0][0];
	//delete[] mMat[0];
	//delete[] mMat;
}

void 	Cube::allocateMemory(int N, int M, int O){
	// Allocate memory for pointers
	double ***ptr1 	= new double**[N];
	double **ptr2 	= new double*[N*M];
	double *ptr3	= new double[N*M*O];
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M ; j++) {
			ptr2[j] = ptr3;
			ptr3 += O;
		}
		ptr1[i] = ptr2;
		ptr2 += M;
	}
	mMat = ptr1;
}

Cube::~Cube(){
	freeMemory();
}
