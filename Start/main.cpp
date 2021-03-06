#include "vmcsolver.h"
#include <time.h>
#include <iostream>

using namespace std;

int main()
{
    VMCSolver *solver = new VMCSolver();
    clock_t start = clock();
    solver->runMonteCarloIntegration();
    clock_t end = clock();
    double time = double(end - start)/CLOCKS_PER_SEC;
    cout << "Time = " << time << endl;  
    return 0;
}
