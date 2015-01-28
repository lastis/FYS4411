#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include "CPhys.h"
#include <assert.h>

using namespace std;
using namespace CPhys;

double** CPhys::matrix(int N, int M){
	/*
    double** mat = new double*[N]; 
	for(int i = 0; i < N; i++){
		mat[i] = new double[M]; 
	}
    */
    double **ptr1 	= new double*[N];
	double *ptr2	= new double[N*M];
    for(int i = 0; i < N; i++){
		ptr1[i] = ptr2;
        ptr2 += M;
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			ptr1[i][j] = 0;
		}
	}
	return ptr1;
}

void	PeriodicBounds::correctPos(double& x, double& L, double& move){
	if(x > L){
		x -= L*(int(x/L));
        move++;
	}
	else if(x < 0){
		x -= L*(int(x/L)-1);
        move--;
	}
}

double	PeriodicBounds::getClosestDist(double& x1, double& x2, double& L){
	double rij = x2 - x1;
	double lHalf = L/2;
	if(rij > lHalf){
		rij -= L;
	}
	else if(rij < -lHalf){
		rij += L;
	}
	return rij;
}

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

double Random::ran0(long &seed)
{

	#define IA 16807
	#define IM 2147483647
	#define AM (1.0/IM)
	#define IQ 127773
	#define IR 2836
	#define MASK 123459876

	long* idum = &seed;
	long     k;
	double   ans;
	*idum ^= MASK;
	k = (*idum)/IQ;
	*idum = IA*(*idum - k*IQ) - IR*k;
	if(*idum < 0) *idum += IM;
	ans=AM*(*idum);
	*idum ^= MASK;
	return ans;

	#undef IA
	#undef IM
	#undef AM
	#undef IQ
	#undef IR
	#undef MASK

}



double Random::ran2(long &seed) {
 	#define IM1 2147483563
  	#define IM2 2147483399
  	#define AM (1.0/IM1)
  	#define IMM1 (IM1-1)
  	#define IA1 40014
  	#define IA2 40692
	#define IQ1 53668
	#define IQ2 52774
	#define IR1 12211
	#define IR2 3791
	#define NTAB 32
	#define NDIV (1+IMM1/NTAB)
	#define EPS 1.2e-7
	#define RNMX (1.0-EPS)
	int            j;
	long           k;
	long* 	 idum = &seed;
	static long    idum2 = 123456789;
	static long    iy=0;
	static long    iv[NTAB];
	double         temp;

	if(*idum <= 0) {
		if(-(*idum) < 1) *idum = 1;
		else             *idum = -(*idum);
		idum2 = (*idum);
		for(j = NTAB + 7; j >= 0; j--) {
			k     = (*idum)/IQ1;
			*idum = IA1*(*idum - k*IQ1) - k*IR1;
			if(*idum < 0) *idum +=  IM1;
			if(j < NTAB)  iv[j]  = *idum;
		}
		iy=iv[0];
	}
	k     = (*idum)/IQ1;
	*idum = IA1*(*idum - k*IQ1) - k*IR1;
	if(*idum < 0) *idum += IM1;
	k     = idum2/IQ2;
	idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
	if(idum2 < 0) idum2 += IM2;
	j     = iy/NDIV;
	iy    = iv[j] - idum2;
	iv[j] = *idum;
	if(iy < 1) iy += IMM1;
	if((temp = AM*iy) > RNMX) return RNMX;
	else return temp;
	#undef IM1
	#undef IM2
	#undef AM
	#undef IMM1
	#undef IA1
	#undef IA2
	#undef IQ1
	#undef IQ2
	#undef IR1
	#undef IR2
	#undef NTAB
	#undef NDIV
	#undef EPS
	#undef RNMX
}


// random numbers with gaussian distribution
double Random::gauss(long &seed) {
	using namespace Random;
	static int iset = 0;
	static double gset;
	long*	 idum = &seed;
	double fac, rsq, v1, v2;

	if ( idum < 0) iset =0;
	if (iset == 0) {
		do {
			v1 = 2.*ran2(*idum) -1.0;
			v2 = 2.*ran2(*idum) -1.0;
			rsq = v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.);
		fac = sqrt(-2.*log(rsq)/rsq);
		gset = v1*fac;
		iset = 1;
		return v2*fac;
	} else {
		iset = 0;
		return gset;
	}
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
