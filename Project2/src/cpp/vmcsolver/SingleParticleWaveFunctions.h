#include <cmath>
#include <iostream>

namespace wave_functions{
    static double alpha;
    static double beta;
    static double nDimensions;

    static double f(double r){
        return -0.5*alpha*r;
    }

    static double g(double r){
        return exp(f(r));
    }

    static double phi1s(double r){
        return g(r)*g(r);
    }

    static double phi1sD(double r){
        return -alpha*g(r)*g(r);
    }

    static double phi1sDD(double r){
        return alpha*alpha*g(r)*g(r);
    }

    static double phi2s(double r){
        return g(r)*(1+f(r));
    }

    static double phi2sD(double r){
        return -alpha*g(r)*(1+0.5*f(r));
    }

    static double phi2sDD(double r){
        return 0.75*alpha*alpha*g(r)*(1+f(r)/3);
    }

    static double phi2p(double r){
        return alpha*r*exp(-alpha*r/2);
    }

    static double phi(int j, double* r){
        double rAbs = 0;
        for (int i = 0; i < nDimensions; i++) {
            rAbs += r[i]*r[i];
        }
        using namespace std;
        rAbs = sqrt(rAbs);
        switch (j) {
            case 0 :
                return phi1s(rAbs);
            case 1 :
                return phi2s(rAbs);
            case 2 ... 4:
                return phi2p(rAbs);
            default:
                std::cout << "Index out of bounds in phi()!!!" << std::endl;
                return 0;
        }
    }

    static double phiD(int j, double* r){
        double rAbs = 0;
        for (int i = 0; i < nDimensions; i++) {
            rAbs += r[i]*r[i];
        }
        rAbs = sqrt(rAbs);
        switch (j) {
            case 0 :
                return phi1sD(rAbs);
            case 1 :
                return phi2sD(rAbs);
            default:
                std::cout << "Index out of bounds in phi()!!!" << std::endl;
                return 0;
        }
    }

    static double phiDD(int j, double* r){
        double rAbs = 0;
        for (int i = 0; i < nDimensions; i++) {
            rAbs += r[i]*r[i];
        }
        rAbs = sqrt(rAbs);
        switch (j) {
            case 0 :
                return phi1sDD(rAbs) + 2*phi1sD(rAbs)/rAbs;
            case 1 :
                return phi2sDD(rAbs) + 2*phi2sD(rAbs)/rAbs;
            default:
                std::cout << "Index out of bounds in phi()!!!" << std::endl;
                return 0;
        }
    }
}
