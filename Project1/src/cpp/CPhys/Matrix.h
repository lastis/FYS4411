#include "Vector.h"
class Matrix {
public:
	Matrix();
	Matrix(int N, int M);
	Matrix(const Matrix& mat);
	~Matrix();

	void	eye  ();
	void	t   ();
	//void	inv ();
	Vector	diag(int k = 0);
	void	diag(Vector& vec, int k = 0);
	int	getN();
	int	getM();
	Vector  getCol(int j           );
	void  	setCol(int j, Vector& v);
	Vector  getRow(int i           );
	void  	setRow(int i, Vector& v);

	double	**getArrayPointer();

	void	reset();
	void	print();

	// Operator overrides
	//Matrix&	operator *(double num);
	//Matrix& operator -(Matrix other);
	Matrix  operator +(double num);
	//Matrix&	operator-=(double num);
	Matrix&	operator+=(double num);
	Matrix&	operator =(double num);
	Matrix& operator =(Matrix other);
	// Return by refrence so the values can be changed
	double& operator()(int i, int j);
private:
	void 	allocateMemory(int     row, int col);
	void	freeMemory    ();
	void	swap	      (Matrix& m1 , Matrix& m2);

	// Variables
	double	** mMat;
	int	mN;
	int 	mM;

};
