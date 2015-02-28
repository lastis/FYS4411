#include "vmcsolver.h"
#include <iostream>

using namespace std;

double localEnergy1(){
    return (alpha - charge)*(1/r)
}
int main(){
    VMCSolver solver = VMCSolver();
    solver.initFromFile("helium1.ini");

    return 0;
}
