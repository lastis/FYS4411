#ifndef _GTO_HELIUM_H_INCLUDED
#define _GTO_HELIUM_H_INCLUDED

#include <cmath>
#include <iostream>
#include "../CPhys/Physical.h"
#include "GtoGeneralFunctions.h"

using namespace gto_general_functions;

namespace gto_helium{
    inline double phi1s(double* r, double rAbs){
        double psi = 0;
        double phi1 = 0;
        double phi2 = 0;
        // Last two numbers are from the basis set, turbmole file.
        phi1 += xi(r[0],r[1],r[2],0,0,0,13.6267,0.17523);
        phi1 += xi(r[0],r[1],r[2],0,0,0,1.99935,0.893483);
        phi2 += xi(r[0],r[1],r[2],0,0,0,0.38299,1);
        psi += phi1*0.4579 + phi2*0.6573;
        return psi;
    }

    inline double phi1sD(int x, double* r, double rAbs){
        double psi = 0;
        double phi1 = 0;
        double phi2 = 0;
        // Last two numbers are from the basis set, turbmole file.
        phi1 += xiD(r[x],0,r[0],r[1],r[2],0,0,0,13.6267,0.17523);
        phi1 += xiD(r[x],0,r[0],r[1],r[2],0,0,0,1.99935,0.893483);
        phi2 += xiD(r[x],0,r[0],r[1],r[2],0,0,0,0.38299,1);
        psi += phi1*0.4579 + phi2*0.6573;
        return psi;
    }

    inline double phi1sDD(double* r, double rAbs){
        double psi = 0;
        double phi1 = 0;
        double phi2 = 0;
        // Last two numbers are from the basis set, turbmole file.
        phi1 += xiDD(r[0],r[1],r[2],0,0,0,13.6267,0.17523);
        phi1 += xiDD(r[0],r[1],r[2],0,0,0,1.99935,0.893483);
        phi2 += xiDD(r[0],r[1],r[2],0,0,0,0.38299,1);
        psi += phi1*0.4579 + phi2*0.6573;
        return psi;
    }
}
#endif
