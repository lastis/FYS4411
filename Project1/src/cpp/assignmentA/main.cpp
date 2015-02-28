#include "vmcsolver.h"
#include <time.h>
#include <iostream>

using namespace std;

int main()
{
    VMCSolver solver = VMCSolver();
    solver.initFromFile("trial2.ini");

    clock_t start = clock();
    solver.runMonteCarloIntegration();
    clock_t end = clock();

    double time = double(end - start)/CLOCKS_PER_SEC;
    cout << "Time = " << time << endl;  

    cout << "Accepted moves: " << int(solver.getStepAcceptance())
	    << " %" << endl;
    /* cout << "Mean distance: " << solver.getR12Mean() << endl; */
    return 0;
}
