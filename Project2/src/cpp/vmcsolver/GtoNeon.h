#ifndef _GTO_NEON_H_INCLUDED
#define _GTO_NEON_H_INCLUDED

#include <cmath>
#include <iostream>
#include "../CPhys/Physical.h"
#include "GtoGeneralFunctions.h"

using namespace gto_general_functions;
using namespace std;

namespace gto_neon{
    static const double C(int q, int j) 
    {
        static const double array[][5] = {
            { -9.8077e-01 ,-2.6062e-01 ,1.1596e-16 ,-8.3716e-18 ,-1.9554e-17 },
            { -9.3714e-02 ,2.5858e-01 ,-2.0106e-16 ,-9.7173e-17 ,-7.3738e-17 },
            { 2.2863e-02, 8.1619e-01 ,-3.2361e-16, 1.3237e-16 ,1.5789e-16 },
            { -9.9519e-19, -5.6186e-18 ,2.7155e-02, -4.0320e-01 ,3.9171e-01 },
            { -1.2125e-18, -2.8615e-16 ,-5.6207e-01, -2.5833e-02 ,1.2375e-02 },
            { -4.1800e-19, 4.6199e-17, 9.1139e-03 ,-3.9180e-01 ,-4.0392e-01 },
            { -1.6696e-19, -4.2405e-18, 2.8890e-02, -4.2895e-01 ,4.1673e-01 },
            { 1.2125e-18, -2.9426e-16, -5.9797e-01 ,-2.7482e-02, 1.3166e-02 },
            { 3.8779e-19, 5.0519e-17, 9.6959e-03 ,-4.1683e-01 ,-4.2972e-01 }
        };
        return array[q][j];
    }
    
    inline double phi(int i, int j, double** r, double* rAbs){
        // Last two numbers are from the basis set, turbmole file.
        double chi[9] = {};
        chi[0] += xi(r[i][0],r[i][1],r[i][2],0,0,0,515.724,0.058143);
        chi[0] += xi(r[i][0],r[i][1],r[i][2],0,0,0,77.6538,0.347951);
        chi[0] += xi(r[i][0],r[i][1],r[i][2],0,0,0,16.8136,0.710714);

        chi[1] += xi(r[i][0],r[i][1],r[i][2],0,0,0,12.483,-0.409922);
        chi[1] += xi(r[i][0],r[i][1],r[i][2],0,0,0,2.66451,1.22431);

        chi[2] += xi(r[i][0],r[i][1],r[i][2],0,0,0,0.60625,1.0);

        chi[3] += xi(r[i][0],r[i][1],r[i][2],1,0,0,12.483,0.24746);
        chi[3] += xi(r[i][0],r[i][1],r[i][2],1,0,0,2.66451,0.851743);
        chi[4] += xi(r[i][0],r[i][1],r[i][2],1,0,0,0.60625,1.0);

        chi[5] += xi(r[i][0],r[i][1],r[i][2],0,1,0,12.483,0.24746);
        chi[5] += xi(r[i][0],r[i][1],r[i][2],0,1,0,2.66451,0.851743);
        chi[6] += xi(r[i][0],r[i][1],r[i][2],0,1,0,0.60625,1.0);

        chi[7] += xi(r[i][0],r[i][1],r[i][2],0,0,1,12.483,0.24746);
        chi[7] += xi(r[i][0],r[i][1],r[i][2],0,0,1,2.66451,0.851743);
        chi[8] += xi(r[i][0],r[i][1],r[i][2],0,0,1,0.60625,1.0);

        double psi = 0;
        for (int q = 0; q < 9; q++) 
        {
            psi += chi[q]*C(q,j);
        }
        return psi;
    }

    inline double phiD(int i, int j, double** r, double* rAbs, int x){
        // Last two numbers are from the basis set, turbmole file.
        double chi[9] = {};
        chi[0] += xiD(r[i][x], 0, r[i][0],r[i][1],r[i][2],0,0,0,515.724,0.058143);
        chi[0] += xiD(r[i][x], 0, r[i][0],r[i][1],r[i][2],0,0,0,77.6538,0.347951);
        chi[0] += xiD(r[i][x], 0, r[i][0],r[i][1],r[i][2],0,0,0,16.8136,0.710714);

        chi[1] += xiD(r[i][x], 0, r[i][0],r[i][1],r[i][2],0,0,0,12.483,-0.409922);
        chi[1] += xiD(r[i][x], 0, r[i][0],r[i][1],r[i][2],0,0,0,2.66451,1.22431);

        chi[2] += xiD(r[i][x], 0, r[i][0],r[i][1],r[i][2],0,0,0,0.60625,1.0);

        int qNum[3];
        qNum[0] = 1; qNum[1] = 0; qNum[2] = 0;
        chi[3] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],1,0,0,12.483,0.24746);
        chi[3] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],1,0,0,2.66451,0.851743);
        chi[4] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],1,0,0,0.60625,1.0);

        qNum[0] = 0; qNum[1] = 1; qNum[2] = 0;
        chi[5] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],0,1,0,12.483,0.24746);
        chi[5] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],0,1,0,2.66451,0.851743);
        chi[6] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],0,1,0,0.60625,1.0);

        qNum[0] = 0; qNum[1] = 0; qNum[2] = 1;
        chi[7] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],0,0,1,12.483,0.24746);
        chi[7] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],0,0,1,2.66451,0.851743);
        chi[8] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],0,0,1,0.60625,1.0);

        double psi = 0;
        for (int q = 0; q < 9; q++) 
        {
            psi += chi[q]*C(q,j);
        }
        return psi;
    }

    inline double phiDD(int i, int j, double** r, double* rAbs){
        // Last two numbers are from the basis set, turbmole file.
        double chi[9] = {};
        chi[0] += xiDD(r[i][0],r[i][1],r[i][2],0,0,0,515.724,0.058143);
        chi[0] += xiDD(r[i][0],r[i][1],r[i][2],0,0,0,77.6538,0.347951);
        chi[0] += xiDD(r[i][0],r[i][1],r[i][2],0,0,0,16.8136,0.710714);

        chi[1] += xiDD(r[i][0],r[i][1],r[i][2],0,0,0,12.483,-0.409922);
        chi[1] += xiDD(r[i][0],r[i][1],r[i][2],0,0,0,2.66451,1.22431);

        chi[2] += xiDD(r[i][0],r[i][1],r[i][2],0,0,0,0.60625,1.0);

        chi[3] += xiDD(r[i][0],r[i][1],r[i][2],1,0,0,12.483,0.24746);
        chi[3] += xiDD(r[i][0],r[i][1],r[i][2],1,0,0,2.66451,0.851743);
        chi[4] += xiDD(r[i][0],r[i][1],r[i][2],1,0,0,0.60625,1.0);

        chi[5] += xiDD(r[i][0],r[i][1],r[i][2],0,1,0,12.483,0.24746);
        chi[5] += xiDD(r[i][0],r[i][1],r[i][2],0,1,0,2.66451,0.851743);
        chi[6] += xiDD(r[i][0],r[i][1],r[i][2],0,1,0,0.60625,1.0);

        chi[7] += xiDD(r[i][0],r[i][1],r[i][2],0,0,1,12.483,0.24746);
        chi[7] += xiDD(r[i][0],r[i][1],r[i][2],0,0,1,2.66451,0.851743);
        chi[8] += xiDD(r[i][0],r[i][1],r[i][2],0,0,1,0.60625,1.0);

        double psi = 0;
        for (int q = 0; q < 9; q++) 
        {
            psi += chi[q]*C(q,j);
        }
        return psi;
    }

}
#endif
