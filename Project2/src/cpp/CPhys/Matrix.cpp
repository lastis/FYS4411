#include <iostream>
#include <cstdio>
#include "Matrix.h"

using namespace std;

Matrix::Matrix(){
	mN = 1;
	mM = 1;
	allocateMemory(mN,mM);
    reset();
}

Matrix::Matrix(int N, int M){
	mN = N;
	mM = M;
    if(mN < 1) mN = 1;
    if(mM < 1) mM = 1;
	allocateMemory(N,M);
}

Matrix::Matrix(const Matrix &mat){
    copy(mat);
}

Matrix	Matrix::operator+(double num){
	Matrix A = Matrix(*this);
	for (int i = 0; i < mN; i++) {
		for (int j = 0; j < mM; j++) {
			A.mMat[i][j] += num;
		}
	}
	return A;
}

Matrix&	Matrix::operator+=(double num){
	for (int i = 0; i < mN; i++) {
		for (int j = 0; j < mM; j++) {
			mMat[i][j] += num;
		}
	}
	return *this;
}

Matrix&	Matrix::operator =(double num){
	for (int i = 0; i < mN; i++) {
		for (int j = 0; j < mM; j++) {
			mMat[i][j] = num;
		}
	}
	return *this;
}

Matrix& Matrix::operator=(const Matrix& other){
	// We need to do a deep copy so we don't lose
	// our dynamic array (pointer) in the switch;
    copy(other);
	return *this;
}

/* void    Matrix::copy(Matrix& obj){ */
/*     freeMemory(); */
/*     mN = obj.mN; */
/*     mM = obj.mM; */
/*     allocateMemory(mN,mM); */
/*     for (int i = 0; i < mN; i++) { */
/*         for (int j = 0; j < mM; j++) { */
/*             mMat[i][j] = obj.mMat[i][j]; */
/*         } */
/*     } */
/* } */

void    Matrix::copy(const Matrix& obj){
    freeMemory();
    mN = obj.mN;
    mM = obj.mM;
    allocateMemory(mN,mM);
    for (int i = 0; i < mN; i++) {
        for (int j = 0; j < mM; j++) {
            mMat[i][j] = obj.mMat[i][j];
        }
    }
}

double&	Matrix::operator()(int i, int j){
	if(i > mN || i < 0){
		cout << "Index out of bounds i = " << i;
		cout << " N = " << mN << "." << endl;
		return mMat[0][0];
	}
	if(j > mM || j < 0){
		cout << "Index out of bounds j = " << j;
		cout << " M = " << mM << "." << endl;
		return mMat[0][0];
	}
	return mMat[i][j];
}

void 	Matrix::diag(Vector& vec, int k){
	// Tod0 Put in check for square matrix
	int a, b;
	int d;
	double *arr = vec.getArrayPointer();
	// Decide the length of the diagonal
	d = mN;
	// Decide wether we will insert the vector above
	// or beneath the diagonal
	if (k > 0){	  // Above
		a = 0;
		b = k;
		d = d - k;
	}
	else if (k == 0){ // Center
		a = 0;
		b = 0;
	}
	else{		  // Beneath
		a = -k;
		b = 0;
		d = d + k;
	}
	// Copy vector to our diagonal
	for (int i = 0; i < d; i++) {
		// k is offset from the diagonal
		mMat[i+a][i+b] = arr[i];
	}
}

Vector	Matrix::diag(int k){
	// Tod0 Put in check for square matrix
	int d;
	int a, b;
	// Decide the length of the diagonal
	if(mN < mM) d = mN;
	else	    d = mM;
	// Decide wether we will aquire the vector above
	// or beneath the diagonal
	if (k > 0){	  // Above
		a = 0;
		b = k;
		d = d - k;
	}
	else if (k == 0){ // Center
		a = 0;
		b = 0;
	}
	else{		  // Beneath
		a = -k;
		b = 0;
		d = d + k ;
	}
	// Create a vector and return the diagonal
	Vector vec = Vector(d);
	double *arr = vec.getArrayPointer();
	for (int i = 0; i < d; i++) {
		// k is offset from the diagonal
		arr[i] = mMat[i+a][i+b];
	}
	return vec;
}

double** Matrix::getArrayPointer(){
	return mMat;
}
Vector	Matrix::getRow(int i){
	Vector vec   = Vector(mM);
	double* pVec = vec.getArrayPointer();
	for (int j = 0; j < mM; j++) {
		pVec[j] = mMat[i][j];
	}
	return vec;
}
void	Matrix::setRow(int i, Vector& v){
	double* pV = v.getArrayPointer();
	for (int j = 0; j < mM; j++) {
		mMat[i][j] = pV[j];
	}
}

void  	Matrix::setCol(int j, Vector& v){
	double* pV = v.getArrayPointer();
	for (int i = 0; i < mN; i++) {
		mMat[i][j] = pV[i];
	}
}

Vector  Matrix::getCol(int j){
	Vector vec   = Vector(mN);
	double* pVec = vec.getArrayPointer();
	for (int i = 0; i < mN; i++) {
		pVec[i] = mMat[i][j];
	}
	return vec;
}

int	Matrix::getN(){
	return mN;
}	

int 	Matrix::getM(){
	return mM;
}	

void 	Matrix::print(){
	double num = 0;
	for (int i = 0; i < mN; i++) {
		for (int j = 0; j < mM; j++) {
			num = mMat[i][j];
			if(num < 0) printf("%.4f\t", num);
			else 	    printf(" %.4f\t", num);
		}
		std::cout << std::endl;
	}
}

void 	Matrix::reset(){
	for (int i = 0; i < mN; i++) {
		for (int j = 0; j < mM; j++) {
			mMat[i][j] = 0;
		}
	}
}

void	Matrix::t(){
	// Transpose
	// Inverted size
	Matrix matrix = Matrix(mM,mN);
	double** nMat = matrix.getArrayPointer();
	for (int i = 0; i < mN; i++) {
		for (int j = 0; j < mM; j++) {
			nMat[j][i] = mMat[i][j];
		}
	}
	swap(*this,matrix);
}

void 	Matrix::eye(){
	// Create the idintity matrix
	for (int i = 0; i < mN; i++) {
		for (int j = 0; j < mM; j++) {
			if(i == j) mMat[i][j] = 1;
			else 	   mMat[i][j] = 0;
		}
	}
}

void 	Matrix::freeMemory(){
	delete[] ref1;
	delete[] ref2;
}

// 	Create matrix in row major order. 
void Matrix::allocateMemory(int row, int col){
    /* double* ptr1 = new (nothrow) double[row*col]; */
    /* double** ptr2 = new (nothrow) double*[row]; */
    double* ptr1 = new double[row*col];
    double** ptr2 = new double*[row];
    ref1 = ptr1;
    ref2 = ptr2;
    mMatFlat = ref1;
	mMat = ref2;
	if(ptr1 == NULL) { 
		cout << "Could not allocate memory in matrix creation" << endl;
		mMat = NULL;
		return;
	}
	if(ptr2 == NULL) { 
		cout << "Could not allocate memory in matrix creation" << endl;
		mMat = NULL;
		return;
	}
	// Assign pointers to the correct rows
	for (int i = 0; i < row; i++) {
		ptr2[i] = ptr1; 
		ptr1   += col;
	}
}

Matrix::~Matrix(){
	freeMemory();
}
