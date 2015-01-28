class Cube {
public:
	Cube();
	Cube(int N, int M, int O);
	Cube(const Cube& mat);
	~Cube();

	int	getN();
	int	getM();
	int	getO();

	double	***getArrayPointer();

	void	reset();

	// Operator overrides
	//Cube&	operator *(double num);
	//Cube& operator -(Matrix other);
	//Cube  operator +(double num);
	//Cube&	operator-=(double num);
	//Cube&	operator+=(double num);
	Cube&	operator =(double num);
	Cube& 	operator =(Cube other);
	// Return by refrence so the values can be changed
	double& operator()(int i, int j, int k);
private:
	void 	allocateMemory(int i, int j, int k);
	void	freeMemory    ();
	void	swap	      (Cube& m1 , Cube& m2);

	// Variables
	double	*** mMat;
	int	mN;
	int 	mM;
	int	mO;

};
