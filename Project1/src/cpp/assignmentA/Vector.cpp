#include <algorithm>
#include <cstdio>
#include <iostream>
#include "Vector.h"

using namespace std;

Vector::Vector(){
	mN = 0;
	allocateMemory(mN);
}

Vector::Vector(int N){
	mN = N;
	allocateMemory(N);
}

Vector::Vector(const Vector& vec){
	// Note: copy constructor is needed for the assignment 
	// operator (operator=(Vector other)) "The rule of three"
	mN = vec.mN;
	allocateMemory(vec.mN);
	for (int i = 0; i < mN; i++) {
		mVec[i] = vec.mVec[i];
	}
}

void	Vector::sort(){
	std::sort(mVec, mVec + mN);
}

Vector&	Vector::operator=(Vector other){
	// We need to do a deep copy so we don't lose
	// our dynamic array (pointer) in the switch;
	copy(other);
	return *this;
}

Vector&	Vector::operator=(double num){
	for (int i = 0; i < mN; i++) {
		mVec[i] = num;
	}
	return *this;
}

Vector&	Vector::operator+=(double num){
	for (int i = 0; i < mN; i++) {
		mVec[i] += num;
	}
	return *this;
}

double&	Vector::operator()(int i){
	// We want to do a check so the user is "protected"
	// against errors
	if(i > mN || i < 0){
		cout << "Index out of bounds i = " << i;
		cout << " N = " << mN << "." << endl;
		return mVec[0];
	}
	double& num = mVec[i];
	return num;
}

double 	Vector::linspace(double start, double end){
	int N = this->getLength();
	double delta = end - start;
	double h  = delta/(N-1);       // step defined here
	for (int i = 0; i < N; i++) {
		mVec[i] = start + i*h;
	}
	return h;
}

void	Vector::copy(Vector& vec){
	// Copy length
	mN = vec.mN;
	// Delte this array
	freeMemory();
	// Make new array and copy
	allocateMemory(mN);
	for (int i = 0; i < mN; i++) {
		mVec[i] = vec.mVec[i];
	}
}
  
double*	Vector::getArrayPointer(){
	// This function is the only one that let's the user
	// play with the array. This class is created so the 
	// user will have as little contact with this variable
	// as possible
	return mVec;
}

int	Vector::getLength(){
	return mN;
}

void	Vector::reset(){
	for (int i = 0; i < mN; i++) {
		mVec[i] = 0;
	}
}

void	Vector::print(){
	double num = 0;
	for (int i = 0; i < mN; i++) {
		num = mVec[i];
		if(num < 0) printf("%.4f\t", num);
		else 	    printf(" %.4f\t", num);
		cout << std::endl;
	}
}

void Vector::allocateMemory(int length){
	// Allocate memory for pointers
	double *arr = new double [length];
	if(!arr) { 
		cout << "Segmentation fault in vector creation" << endl;
		mVec = NULL;
		return;
	}
	mVec = arr;
}
void	Vector::freeMemory(){
	if(mVec) delete[] mVec;
}

Vector::~Vector(){
	freeMemory();
}
