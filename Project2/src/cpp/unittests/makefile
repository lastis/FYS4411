compile:
	g++ ../vmcsolver/VMCSolver.cpp main.cpp ../CPhys/CPhys.cpp  \
		../CPhys/Vector.cpp ../CPhys/Matrix.cpp ../CPhys/Cube.cpp -lunittest++
alt:
	g++ --std=c++0x -fopenmp    \
    main1.cpp                   \
    ../vmcsolver/VMCWrapper.cpp \
    ../vmcsolver/VMCSolver1.cpp \
    ../CPhys/CPhys.cpp          \
    ../CPhys/Vector.cpp         \
    ../CPhys/Matrix.cpp         \
    ../CPhys/Cube.cpp           \
    -lunittest++
