#ifndef _GTO_BERYLLIUM_H_INCLUDED
#define _GTO_BERYLLIUM_H_INCLUDED

#include <cmath>
#include <iostream>
#include "../CPhys/Physical.h"
#include "GtoGeneralFunctions.h"

using namespace gto_general_functions;
using namespace std;

namespace gto_beryllium{

    static const double C(int q, int j) 
    {
        static const double array[][2] = 
        {{-9.9281e-01,-2.1571e-01},
            {-7.6425e-02,2.2934e-01},
            {2.8727e-02,8.2235e-01},
            {1.2898e-16, 5.1721e-16},
            {-2.3257e-19 ,4.5670e-18},
            {5.6097e-19 ,-1.1040e-17},
            {1.2016e-16 ,8.5306e-16},
            {-4.6874e-19 ,7.0721e-18},
            {1.1319e-18 ,-1.7060e-17}};
        /* static const double array[][2] = */ 
        /* {{-9.9281e-01,-2.1571e-01}, */
        /*     {-7.6425e-02,2.2934e-01}, */
        /*     {2.8727e-02,8.2235e-01}, */
        /*     {1.2898e-16, 5.1721e-16}, */
        /*     {5.6097e-19 ,-1.1040e-17}, */
        /*     {-4.6874e-19 ,7.0721e-18}, */
        /*     {-2.3257e-19 ,4.5670e-18}, */
        /*     {1.2016e-16 ,8.5306e-16}, */
        /*     {1.1319e-18 ,-1.7060e-17}}; */
        return array[q][j];
    }
    
    inline double phi(int i, int j, double** r, double* rAbs){
        // Last two numbers are from the basis set, turbmole file.
        double chi[9] = {};
        chi[0] += xi(r[i][0],r[i][1],r[i][2],0,0,0,71.8876,0.0644263);
        chi[0] += xi(r[i][0],r[i][1],r[i][2],0,0,0,10.7289,0.366096);
        chi[0] += xi(r[i][0],r[i][1],r[i][2],0,0,0,2.22205,0.695934);

        chi[1] += xi(r[i][0],r[i][1],r[i][2],0,0,0,1.29548,-0.421064);
        chi[1] += xi(r[i][0],r[i][1],r[i][2],0,0,0,0.268881,1.22407);

        chi[2] += xi(r[i][0],r[i][1],r[i][2],0,0,0,0.07735,1.00);

        chi[3] += xi(r[i][0],r[i][1],r[i][2],1,0,0,1.29548,0.205132);
        chi[3] += xi(r[i][0],r[i][1],r[i][2],1,0,0,0.2688810,0.8825280);
        chi[4] += xi(r[i][0],r[i][1],r[i][2],1,0,0,0.07735,1.0);

        chi[5] += xi(r[i][0],r[i][1],r[i][2],0,1,0,1.29548,0.205132);
        chi[5] += xi(r[i][0],r[i][1],r[i][2],0,1,0,0.2688810,0.8825280);
        chi[6] += xi(r[i][0],r[i][1],r[i][2],0,1,0,0.07735,1.0);

        chi[7] += xi(r[i][0],r[i][1],r[i][2],0,0,1,1.29548,0.205132);
        chi[7] += xi(r[i][0],r[i][1],r[i][2],0,0,1,0.2688810,0.8825280);
        chi[8] += xi(r[i][0],r[i][1],r[i][2],0,0,1,0.07735,1.0);

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
        chi[0] += xiD(r[i][x], 0, r[i][0],r[i][1],r[i][2],0,0,0,71.8876,0.0644263);
        chi[0] += xiD(r[i][x], 0, r[i][0],r[i][1],r[i][2],0,0,0,10.7289,0.366096);
        chi[0] += xiD(r[i][x], 0, r[i][0],r[i][1],r[i][2],0,0,0,2.22205,0.695934);

        chi[1] += xiD(r[i][x], 0, r[i][0],r[i][1],r[i][2],0,0,0,1.29548,-0.421064);
        chi[1] += xiD(r[i][x], 0, r[i][0],r[i][1],r[i][2],0,0,0,0.268881,1.22407);

        chi[2] += xiD(r[i][x], 0, r[i][0],r[i][1],r[i][2] ,0,0,0,0.07735,1.00);

        int qNum[3];
        qNum[0] = 1; qNum[1] = 0; qNum[2] = 0;
        chi[3] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],1,0,0,1.29548,0.205132);
        chi[3] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],1,0,0,0.2688810,0.8825280);
        chi[4] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],1,0,0,0.07735,1.0);

        qNum[0] = 0; qNum[1] = 1; qNum[2] = 0;
        chi[5] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],0,1,0,1.29548,0.205132);
        chi[5] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],0,1,0,0.2688810,0.8825280);
        chi[6] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],0,1,0,0.07735,1.0);

        qNum[0] = 0; qNum[1] = 0; qNum[2] = 1;
        chi[7] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],0,0,1,1.29548,0.205132);
        chi[7] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],0,0,1,0.2688810,0.8825280);
        chi[8] += xiD(r[i][x], qNum[x], r[i][0],r[i][1],r[i][2],0,0,1,0.07735,1.0);

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
        chi[0] += xiDD(r[i][0],r[i][1],r[i][2],0,0,0,71.8876,0.0644263);
        chi[0] += xiDD(r[i][0],r[i][1],r[i][2],0,0,0,10.7289,0.366096);
        chi[0] += xiDD(r[i][0],r[i][1],r[i][2],0,0,0,2.22205,0.695934);

        chi[1] += xiDD(r[i][0],r[i][1],r[i][2],0,0,0,1.29548,-0.421064);
        chi[1] += xiDD(r[i][0],r[i][1],r[i][2],0,0,0,0.268881,1.22407);

        chi[2] += xiDD(r[i][0],r[i][1],r[i][2],0,0,0,0.07735,1.00);

        chi[3] += xiDD(r[i][0],r[i][1],r[i][2],1,0,0,1.29548,0.205132);
        chi[3] += xiDD(r[i][0],r[i][1],r[i][2],1,0,0,0.2688810,0.8825280);
        chi[4] += xiDD(r[i][0],r[i][1],r[i][2],1,0,0,0.07735,1.0);

        chi[5] += xiDD(r[i][0],r[i][1],r[i][2],0,1,0,1.29548,0.205132);
        chi[5] += xiDD(r[i][0],r[i][1],r[i][2],0,1,0,0.2688810,0.8825280);
        chi[6] += xiDD(r[i][0],r[i][1],r[i][2],0,1,0,0.07735,1.0);

        chi[7] += xiDD(r[i][0],r[i][1],r[i][2],0,0,1,1.29548,0.205132);
        chi[7] += xiDD(r[i][0],r[i][1],r[i][2],0,0,1,0.2688810,0.8825280);
        chi[8] += xiDD(r[i][0],r[i][1],r[i][2],0,0,1,0.07735,1.0);

        double psi = 0;
        for (int q = 0; q < 9; q++) 
        {
            psi += chi[q]*C(q,j);
        }
        return psi;
    }

}
#endif
