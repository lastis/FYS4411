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
	void	reset();
	void	print();
	double	**getArrayPointer();

	// Operator overrides
	//Matrix&	operator *(double num);
	//Matrix& operator -(Matrix other);
	Matrix  operator +(double num);
	//Matrix&	operator-=(double num);
	Matrix&	operator+=(double num);
	Matrix&	operator =(double num);
	Matrix& operator =(const Matrix& other);
	// Return by refrence so the values can be changed
	double& operator()(int i, int j);
private:
	void 	allocateMemory(int     row, int col);
	void	freeMemory    ();
    void    copy(const Matrix& obj);
    /* void    copy(Matrix& obj); */

	// Variables
    double*  mMatFlat;
	double** mMat;
    double*  ref1;
	double** ref2;
	int	    mN;
	int 	mM;

};
