#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include "CPhys.h"
#include <cstdio>
#include <assert.h>

using namespace std;
using namespace CPhys;


Matrix	Lattice::getFCC(int Nc, double b){
	// NOTES:
	// Nc cannot be less than 1
	// Counting is done by adding all atoms from all unit cells
	// then adding the sides except the edges, then add the 
	// edges (they overlap on the last atom)
	int 	atoms 	= 4*Nc*Nc*Nc;
	Matrix 	r = Matrix(atoms,3);
	// Initiate r vector
	int index = 0;
	for (int k = 0; k < Nc; k++) {
		for (int j = 0; j < Nc; j++) {
			for (int i = 0; i < Nc; i++) {
				r(index,0) = i*b;
				r(index,1) = j*b;
				r(index,2) = k*b;
				index++;

				r(index,0) = i*b + 0.5*b;
				r(index,1) = j*b + 0.5*b;
				r(index,2) = k*b;
				index++;

				r(index,0) = i*b;
				r(index,1) = j*b + 0.5*b;
				r(index,2) = k*b + 0.5*b;
				index++;

				r(index,0) = i*b + 0.5*b;
				r(index,1) = j*b;
				r(index,2) = k*b + 0.5*b;
				index++;
			}
		}
	}
	return r;
}

Vector	LinAlg::tridiagSolve(double a, double b, double c, Vector y){
	int N = y.getLength();
	Vector x = Vector(N);
	double* xVec  = x.getArrayPointer();
	double* yVec  = y.getArrayPointer();
	pLinAlg::tridiagSolve(a,b,c,xVec,yVec,N);
	// return the vector
	return x;
}

void EigVal::jacobiMethod(Matrix& A, Matrix& R, int N){
	using namespace pEigVal;
	double** ppR = R.getArrayPointer();
	double** ppA = A.getArrayPointer();
	int 	 k, l;
	int 	 iterations  	= 0;
	double	 eps 		= 10e-8;
	double 	 maxIterations 	= (double) N * N * N;
	// Prepare eigenvector matrix
	R.eye();
	// Decide which elements to rotate
	double 	maxoff		= maxoffdiag(ppA, &k, &l, N);
	// Do the rotations until we converge at a solution
	while (		maxoff > eps 	
		 && iterations < maxIterations ) {
		
		int num = iterations;
		if(num%20000 == 0){
			cout << "Computing - "<< iterations << " iterations ";
			cout << " max off diag = " << maxoff << endl;
		}
		rotate(ppA, ppR, k, l, N );
		maxoff = maxoffdiag(ppA, &k, &l, N );
		iterations++;
	}
	cout << "Number of iterations: " << iterations << "\n";
}



double 	VecOp::normalize(Vector& v, double dx){
	int     N   = v.getLength();
	double* pV  = v.getArrayPointer();
	double  sum = 0;
	double  Z   = 0;
	for (int i = 0; i < N; i++) {
		sum += pV[i]*pV[i]*dx;
	}
	Z = 1/sqrt(sum);
	for (int i = 0; i < N; i++) {
		pV[i] *= Z;
	}
	return Z;

}

void	MatOp::sortCol(Matrix& A, Vector& v){
	// We want to sort the colums in A by the elements in v.
	int rows = A.getN();
	int cols = A.getM();
	// Make a new matrix wit dim N = A.N+1 and M = A.M
	Matrix tmp = Matrix(rows+1,cols);
	// Insert v in the first row
	tmp.setRow(0,v);
	// Copy A to tmp after the first row
	double** ppTmp = tmp.getArrayPointer();
	double** ppA   =    A.getArrayPointer();
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ppTmp[i+1][j] = ppA[i][j];
		}
	}
	// Our matrix is row major so we need to transpose
	// to be able to sort by columns
	// Now **arr[0,1,2, ... , i] points to whole arrays 
	// of the columns
	tmp.t();
	// Sort by first elements of A (sorting by vector v)
	ppTmp = tmp.getArrayPointer();
	qsort((void*)ppTmp,cols,sizeof(double),&pMatOp::compareTwoRows);
	// Transpose back
	tmp.t();
	ppTmp = tmp.getArrayPointer();
	// Now copy the sorted elements to A
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			 ppA[i][j] = ppTmp[i+1][j];
		}
	}
	// Then sort v
	v.sort();
}

double  MatOp::getDet(Matrix& A)
{
    int n = A.getN();
    Matrix L = Matrix(n,n);
    Matrix U = Matrix(n,n);
    double** a = A.getArrayPointer();
    double** l = L.getArrayPointer();
    double** u = U.getArrayPointer();
    pMatOp::decomposeLU(a,l,u,n);
    double res = 1;
    for (int i = 0; i < n; i++) 
    {
        res *= u[i][i];
    }
    return res;
}

void MatOp::decomposeLU(Matrix& mat, Matrix& L, Matrix& U){
    int n = mat.getN();
    L = Matrix(n,n);
    U = Matrix(n,n);
    double** l = L.getArrayPointer();
    double** u = U.getArrayPointer();
    double** a = mat.getArrayPointer();
    pMatOp::decomposeLU(a,l,u,n);
}

/** \brief Update the inverse matrix when only the i'th row has been changed.
 *
 * \param i The i'th row that has been changed.
 * \param ratio The ratio between the new matrix and the old one. 
 * \param mat The matrix to update the inverse for. 
 * \param inv The old inverse matrix.
 * \param N Size of the matrix.
 */         
void pMatOp::updateInverse(int i, double ratio, double** mat, double** inv, int N){
    // Update the inverse matrix for all columns except the i'th.
    for (int j = 0; j < N; j++) {
        if (j == i) continue;
        double Sj = 0;
        for (int l = 0; l < N; l++) {
            // d_il(new) * dInv_lj(old)
            Sj += mat[i][l]*inv[l][j];
        }
        for (int k = 0; k < N; k++) {
            inv[k][j] = inv[k][j] - Sj/ratio*inv[k][i];
        }
    }
    // Update the i'th column.
    for (int k = 0; k < N; k++) {
        inv[k][i] = inv[k][i]/ratio;
    }
}

void pMatOp::decomposeLU(double** a, double** l, double** u, int n){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i<=j){
                u[i][j]=a[i][j];
                for(int k=0; k<i; k++)
                    u[i][j]-=l[i][k]*u[k][j];
                if(i==j)
                    l[i][j]=1;
                else
                    l[i][j]=0;
            }
            else{
                l[i][j]=a[i][j];
                for(int k=0; k<=j-1; k++)
                    l[i][j]-=l[i][k]*u[k][j];
                l[i][j]/=u[j][j];
                u[i][j]=0;
            }
        }
    }
}

Matrix  MatOp::getInverse(Matrix mA){
    int N = mA.getN();
    Matrix mAInv = Matrix(N,N);
    Matrix mU = Matrix(N,N);
    Matrix mL = Matrix(N,N);
    Matrix mLInv = Matrix(N,N);
    Matrix mI = Matrix(N,N);
    mI.eye();
    double** AInv = mAInv.getArrayPointer();
    double** A = mA.getArrayPointer();
    double** U = mU.getArrayPointer();
    double** L = mL.getArrayPointer();
    double** LInv = mLInv.getArrayPointer();
    double** I = mI.getArrayPointer();
    pMatOp::decomposeLU(A,L,U,N);
    pMatOp::substituteForward(L,LInv,I,N);
    pMatOp::substituteBackward(U,AInv,LInv,N);
    return mAInv;
}

void MatOp::substituteBackward(Matrix& U, Matrix& x, Matrix& y){
    int n = U.getN();
    x = Matrix(n,n);
    pMatOp::substituteBackward(U.getArrayPointer(),x.getArrayPointer(), 
            y.getArrayPointer(),n);
}

void pMatOp::substituteBackward(double** U, double** x, double** y, int N){
    for (int k = N-1; k >= 0; k--) {
        for (int i = N-1; i >= 0; i--) {
            double a = y[i][k];
            for (int j = i+1; j < N; j++) {
                a = a - x[j][k]*U[i][j];
            }
            x[i][k] = a/U[i][i];
        }
    }
}

void MatOp::substituteForward(Matrix& L, Matrix& y, Matrix& b){
    int n = L.getN();
    y = Matrix(n,n);
    pMatOp::substituteForward(L.getArrayPointer(),y.getArrayPointer(), 
            b.getArrayPointer(),n);
}

void pMatOp::substituteForward(double** L, double** y, double** b, int N){
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < N; i++) {
            double a = b[i][k];
            int j = 0;
            while (j < i) {
                a = a - y[j][k]*L[i][j];
                j++;
            }
            y[i][k] = a/L[i][i];
        }
    }
}

void	pLinAlg::tridiagSolve(double  a, double  b, double c, 
			      double* x, double* y, int    N){
	double* temp  = new double[N];
	double* aVec  = new double[N];
	double* bVec  = new double[N];
	double* cVec  = new double[N];
	double  bTemp = 0;
	// Init arrays
	for (int i = 0; i < N; i++) {
		temp[i] = 0;
		aVec[i] = a;
		bVec[i] = b;
		cVec[i] = c;
	}
	// Forward substitution
	bTemp = bVec[1];
	x[1]  = y[1]/bTemp;
	for (int i = 2; i < N; i++) {
		temp[i] =  cVec[i-1]/bTemp;
		bTemp   =  bVec[i] - aVec[i]*temp[i];
		x[i] = (y[i] - aVec[i]*x[i-1])/bTemp;
	}
	// Backward substitution
	for (int i = N-1; i >= 1; i--) {
		x[i] -= temp[i+1]*x[i+1];
	}
}

void	pDiffusion::stepJacobi2D(double** u, double** uNext, double alpha,
				 int	      N, int       M, double	T){

	double 	factor 	= 1/(1+4*alpha);
	double	eps 	= 0.00001;
	// Make sure diff is bigger than eps
	double 	diff 	= eps + 1;
	Matrix	mat	= Matrix(N,M);
	double** temp	= mat.getArrayPointer();
	// Our guess of the next step is the previous step, we therefor
	// copy u to a temp array
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			temp[i][j] = u[i][j];
		}
	}
	// We want the two arrays uNext and temp to converge. When the
	// difference between them is below a given threshold we stop
	// the loop
	while (diff > eps){
		// Make suggestions
		for (int i = 1; i < N-1; i++) {
			for (int j = 1; j < M-1; j++) {
				uNext[i][j]  = factor*(
		alpha*(temp[i+1][j]+temp[i-1][j]+temp[i][j+1]+temp[i][j-1])
					    +u[i][j]);
			}
		}
		// Insert boundary conditions
		for (int k = 0; k < N; k++) {
			double xyz = float(k)/N;
			uNext[0][k] 	= (1-xyz)*exp(T);
			uNext[N-1][k] 	= (1-xyz)*exp(1+T);
			uNext[k][0] 	= exp(xyz+T);
			uNext[k][N-1] 	= 0;
		}
		// Calculate difference
		diff = 0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < M; j++) {
				diff += abs(temp[i][j] - uNext[i][j]);
			}
		}
		diff /= N*M;
		// Copy uNew to an temp array
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < M; j++) {
				temp[i][j] = uNext[i][j];
			}
		}
	}
}

void	pDiffusion::stepEuler2D(double** u, double** un, double alpha,
				int      N, int       M ){
	for (int i = 1; i < N-1; i++) {
		for (int j = 1; j < M-1; j++) {
			un[i][j]  = u[i][j] + alpha*(u[i+1][j] + u[i-1][j] + 
				    u[i][j+1] + u[i][j-1] - 4*u[i][j]);
		}
	}
}

void	pDiffusion::stepEuler(double* u, double* un, 
			      double  a, double   b, int N){
	// Forward Euler scheme for the diffusion equation
	for (int i = 1; i < N-1; i++) {
		un[i] = a*u[i-1] + b*u[i] + a*u[i+1];
	}
	un[0] = u[0];
}

void	pDiffusion::stepBackwardEuler(double* u, double* un, 
				      double  a, double   b, int N){
	// Backward Euler scheme for the diffusion equation
	pLinAlg::tridiagSolve(a,b,a,un,u,N);
	un[N-1] = u[N-1];
}

int 	pMatOp::compareTwoRows(const void* rowA, const void* rowB){
	// Compare the first element in row A (this is what ** does) and 
	// compare it to the first element in row B
	return (**(double**) rowA - **(double**) rowB);
}

double 	pEigVal::maxoffdiag(double** A, int* k, int* l, int N){
	double max  = 0.0;
	double a_ij = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			a_ij = fabs(A[i][j]);
			if (a_ij > max){
				max	= a_ij;
				*l 	= i;
				*k 	= j;
			}
		}
	}
	return max;
}

void 	pEigVal::rotate(double** A, double** R,int k, int l, int N){
	double s, c;
	double tau, t;
	tau = (A[l][l] -A[k][k])/(2*A[k][l]);
	if(tau > 0) t = -tau + sqrt(1+tau*tau);
	else 	    t = -tau - sqrt(1+tau*tau);
	c = 1 / sqrt(1+t*t);
	s = t * c;
	double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
	a_kk = A[k][k];
	a_ll = A[l][l];
	// changing the matrix elements with indices k and l
	A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
	A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
	A[k][l] = 0.0; // hard-coding of the zeros
	A[l][k] = 0.0;
	// and then we change the remaining elements
	for ( int i = 0; i < N; i++ ) {
		if ( i != k && i != l ) {
			a_ik 	= A[i][k];
			a_il 	= A[i][l];
			A[i][k]	= c*a_ik - s*a_il;
			A[k][i]	= A[i][k];
			A[i][l]	= c*a_il + s*a_ik;
			A[l][i]	= A[i][l];
		}
		// Finally, we compute the new eigenvectors
		r_ik 	= R[i][k];
		r_il 	= R[i][l];
		R[i][k]	= c*r_ik - s*r_il;
		R[i][l]	= c*r_il + s*r_ik;
	}
}
