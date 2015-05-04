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

    static double berylliumPsiDD(double** r)
    {
        double rAbs[4];
        for (int i = 0; i < 4; i++) {
            rAbs[i] = 0;
            for(int j = 0; j < nDimensions; j++) {
                rAbs[i] += r[i][j] * r[i][j];
            }
            rAbs[i] = sqrt(rAbs[i]);
        }
        double tmp = 0;
        double sum = 0;
        tmp += phi2s(rAbs[1])*phi1sDD(rAbs[0]) - phi1s(rAbs[1])*phi2sDD(rAbs[0]);
        tmp -= phi2s(rAbs[0])*phi1sDD(rAbs[1]) - phi1s(rAbs[0])*phi2sDD(rAbs[1]);
        tmp += 2*(phi2s(rAbs[1])*phi1sD(rAbs[0]) 
                - phi1s(rAbs[1])*phi2sD(rAbs[0]))/rAbs[0];
        tmp -= 2*(phi2s(rAbs[0])*phi1sD(rAbs[1]) 
                - phi1s(rAbs[0])*phi2sD(rAbs[1]))/rAbs[1];
        sum += tmp/(phi1s(rAbs[0]) * phi2s(rAbs[1]) - phi2s(rAbs[0]) * phi1s(rAbs[1]));

        tmp = 0;
        tmp += phi2s(rAbs[3])*phi1sDD(rAbs[2]) - phi1s(rAbs[3])*phi2sDD(rAbs[2]);
        tmp -= phi2s(rAbs[2])*phi1sDD(rAbs[3]) - phi1s(rAbs[2])*phi2sDD(rAbs[3]);
        tmp += 2*(phi2s(rAbs[3])*phi1sD(rAbs[2]) 
                - phi1s(rAbs[3])*phi2sD(rAbs[2]))/rAbs[2];
        tmp -= 2*(phi2s(rAbs[2])*phi1sD(rAbs[3]) 
                - phi1s(rAbs[2])*phi2sD(rAbs[3]))/rAbs[3];
        sum += tmp/(phi1s(rAbs[2]) * phi2s(rAbs[3]) - phi2s(rAbs[2]) * phi1s(rAbs[3]));
        return sum;
    }
}
