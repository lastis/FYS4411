main:
	g++ -fopenmp VMCSolver.cpp main.cpp ../CPhys/CPhys.cpp  \
		../CPhys/Vector.cpp ../CPhys/Matrix.cpp ../CPhys/Cube.cpp 
param:
	g++ -fopenmp VMCSolver.cpp parameters.cpp ../CPhys/CPhys.cpp  \
		../CPhys/Vector.cpp ../CPhys/Matrix.cpp ../CPhys/Cube.cpp 
alt:
	g++ -std=c++0x -fopenmp VMCWrapper.cpp VMCSolver1.cpp main1.cpp ../CPhys/CPhys.cpp  \
		../CPhys/Vector.cpp ../CPhys/Matrix.cpp ../CPhys/Cube.cpp 
