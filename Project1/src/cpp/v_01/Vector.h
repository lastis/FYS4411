class Vector {
public:
	Vector();
	Vector(int  N);
	Vector(const Vector& vec);
	~Vector();
	void 	sort();
	void	reset();
	void	print();
	void	copy(Vector& vec);
	int	getLength();
	double	linspace (double  start, double end);
	double*	getArrayPointer();
	Vector&	operator+=(double num);
	Vector&	operator =(double num);
	Vector& operator =(Vector other);
	// Return by refrence so the values can be changed
	double&	operator()(int i);
private:
	void 	allocateMemory   (int length);
	void	freeMemory       ();

	// Variables
	double*	mVec;
	int	mN;
};
