#include "Matrix.h"
#include "Cube.h"
namespace CPhys{

	namespace Lattice{
		Matrix	getFCC(int Nc, double dist);
	}
	namespace LinAlg{
		Vector	tridiagSolve(double a, double b, double c, Vector y);
	}
	namespace EigVal{
		void 	jacobiMethod(Matrix& A, Matrix& R, int N);
	}
	namespace MatOp{
		void	sortCol(Matrix& A, Vector& v);
        void    decomposeLU(Matrix& mat, Matrix& L, Matrix& U);
        void    substituteForward(Matrix& L, Matrix& y, Matrix& b);
        void    substituteBackward(Matrix& U, Matrix& x, Matrix& y);
        Matrix  getInverse(Matrix A);
	}
	namespace VecOp{
		double	normalize(Vector& v, double dx);
	}
	namespace Random{
		double  ran0(long& seed);
		double	ran2(long& seed);
		double	gauss(long& seed);
	}


	// prefix p (abbr. for 'pointer') is used to avoid name collision
	namespace pMatOp{
        void    decomposeLU(double** mat, double** L, double** U, int n);
        void    substituteForward(double** L, double** y, double** b, int N);
        void    substituteBackward(double** U, double** x, double** y, int N);
        void    getInverse(double** A, double** B);
		int 	compareTwoRows(const void* rowA, const void* rowB);
	}
	namespace pLinAlg{
		void	tridiagSolve(double  a, double  b, double c, 
				     double* x, double* y, int    N);
	}
	namespace pEigVal{
		double 	maxoffdiag        (double** A, int* k, int* l, int N);
		void	rotate(double** A, double** R, int  k, int  l, int N);
	}
	namespace pDiffusion{
		void	stepEuler(double* u, double* un, 
				  double  a, double   b, int N);
		void	stepBackwardEuler(double* u, double* un,
				          double  a, double   b, int N);
		void	stepEuler2D(double** u, double** un, double alpha,
				    int      N, int       M);
		void	stepJacobi2D(double** u, double **un, double alpha,
				     int      N, int      M, double T);
	}
}
