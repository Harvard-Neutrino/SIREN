/*
  CRunDec.cpp, version 3
  alphas_cms.cpp
*/

#include <iostream>
#include <cmath>
#include "CRunDec.h"

using namespace std;

// Numerical input values for alpha_s^(6)(mu) and 'mu' from arXiv:1304:7498:
const double muin  = 896;
const double as6mu = 0.0889;

// Numerical input values for alpha_s^(5)(Mz) and 'Q' from arXiv:1609.05331:
const double asmz = 0.1162;
const double asmzplus = asmz + 0.007;
const double asmzminus = asmz - 0.0062;
const double Q = 1508.04;

// vary decoupling scale between mu/f and f*mu:
const double f = 4.0;

void FindExtrema(double& max, double& min, CRunDec* crd, double asmz, double f) {
    min = 1.0;
    max = 0.0;
    TriplenfMmu dec[4];
    for(int i = 1; i < 4; i++) {
        dec[i].nf = 0;
    }
    for(int i = 0; i < 10000*(f-1/f); i++) {
        dec[0].nf = 6;
        dec[0].Mth = NumDef.Mt;
        dec[0].muth = NumDef.Mt/f + 0.0001*i*NumDef.Mt;
        double as = crd->AlL2AlH(asmz, NumDef.Mz, dec, Q, 5);
        if(as > max) {
            max = as;
        }
        if(as < min) {
            min = as;
        }
    }
}


int main() {
    CRunDec* crd = new CRunDec(5);

    TriplenfMmu dec[4];
    dec[0].nf = 6;
    dec[0].Mth = NumDef.Mt;
    dec[0].muth = 2*NumDef.Mt;
    for(int i = 1; i < 4; i++) {
        dec[i].nf = 0;
    }


    double as5Mz = crd->AlH2AlL(as6mu,muin, dec, NumDef.Mz, 5);
    cout << "\\alpha_s^(5)(M_Z) interpreting \\alpha_s(896 GeV) as \\alpha_s^{(6)}: " << as5Mz << endl;
    as5Mz = crd->AlphasExact(as6mu, muin, NumDef.Mz, 5, 5);
    cout << "\\alpha_s^(5)(M_Z) interpreting \\alpha_s(896 GeV) as \\alpha_s^{(5)}: " << as5Mz << endl;


    int nloops = 2;

    double asQ = crd->AlphasExact(asmz, NumDef.Mz, Q, 5, nloops);
    double asQplus = abs(crd->AlphasExact(asmzplus, NumDef.Mz, Q, 5, nloops)-asQ);
    double asQminus = abs(crd->AlphasExact(asmzminus, NumDef.Mz, Q, 5, nloops)-asQ);
    double truncuncert = abs(crd->AlphasExact(asmz, NumDef.Mz, Q, 5, nloops-1)-asQ);
    double totuncertplus = sqrt(asQplus*asQplus + truncuncert*truncuncert);
    double totuncertminus = sqrt(asQminus*asQminus + truncuncert*truncuncert);
    cout << "\\alpha_s^(5)(1508.04 GeV) using 2-loop evolution: " << asQ << endl;
    cout << "uncertainty from \\delta\\alpha_s:                + " 
	 << asQplus << " - " << asQminus << endl;
    cout << "uncertainty due to truncation of pert. series: +/-" << truncuncert << endl;
    cout << "total uncertainty:                              + " 
	 << totuncertplus << " - " << totuncertminus << endl;

    cout << endl;

    nloops = 5;
    asQ = crd->AlphasExact(asmz, NumDef.Mz, Q, 5, nloops);
    asQplus = abs(crd->AlphasExact(asmzplus, NumDef.Mz, Q, 5, nloops)-asQ);
    asQminus = abs(crd->AlphasExact(asmzminus, NumDef.Mz, Q, 5, nloops)-asQ);
    truncuncert = abs(crd->AlphasExact(asmz, NumDef.Mz, Q, 5, nloops-1)-asQ);
    totuncertplus = sqrt(asQplus*asQplus + truncuncert*truncuncert);
    totuncertminus = sqrt(asQminus*asQminus + truncuncert*truncuncert);
    cout << "\\alpha_s^(5)(1508.04 GeV) using 5-loop evolution: " << asQ << endl;
    cout << "uncertainty from \\delta\\alpha_s:                + " 
	 << asQplus << " - " << asQminus << endl;
    cout << "uncertainty due to truncation of pert. series: +/-" << truncuncert << endl;
    cout << "total uncertainty:                              + " 
	 << totuncertplus << " - " << totuncertminus << endl << endl;

    dec[0].nf = 6;
    dec[0].Mth = NumDef.Mt;
    dec[0].muth = 2*NumDef.Mt;
    for(int i = 1; i < 4; i++) {
        dec[i].nf = 0;
    }
    asQ = crd->AlL2AlH(asmz, NumDef.Mz, dec, Q, nloops);
    dec[0].nf = 6;
    dec[0].Mth = NumDef.Mt;
    dec[0].muth = 2*NumDef.Mt;
    asQplus = abs(asQ-crd->AlL2AlH(asmzplus, NumDef.Mz,dec, Q, nloops));
    dec[0].nf = 6;
    dec[0].Mth = NumDef.Mt;
    dec[0].muth = 2*NumDef.Mt;
    asQminus = abs(asQ-crd->AlL2AlH(asmzminus, NumDef.Mz,dec, Q, nloops));
    dec[0].nf = 6;
    dec[0].Mth = NumDef.Mt;
    dec[0].muth = 2*NumDef.Mt;
    truncuncert = abs(asQ - crd->AlL2AlH(asmz, NumDef.Mz,dec, Q, nloops - 1));
    double max, min;
    FindExtrema(max,min,crd,asmz,f);
    double scaleuncert = max - min;
    totuncertplus  = sqrt(truncuncert*truncuncert+asQplus*asQplus  +scaleuncert*scaleuncert);
    totuncertminus = sqrt(truncuncert*truncuncert+asQminus*asQminus+scaleuncert*scaleuncert);

    cout << "\\alpha_s^(6)(1508.04 GeV) using 5-loop evolution: " << asQ << endl;
    cout << "uncertainty from \\delta\\alpha_s:                + " 
	 << asQplus << " - " << asQminus << endl;
    cout << "uncertainty due to truncation of pert. series: +/-" << truncuncert << endl;

    cout << "minimal and maximal values of alpha_s in interval" << endl;
    cout << "       " << NumDef.Mt/f << " <= mu <= " << f*NumDef.Mt << "." << endl;
    cout << "max: " << max << " min: " << min << endl;
    cout << "uncertainty due to variation of decoupling scale: +/-" << scaleuncert << endl;
    cout << "total uncertainty,                                 + " 
	 << totuncertplus << " - " << totuncertminus << endl << endl;


    dec[0].nf = 6;
    dec[0].Mth = NumDef.Mt;
    dec[0].muth = NumDef.Mz;
    asQ = crd->AlL2AlH(asmz, NumDef.Mz,dec, Q, nloops);
    cout << "alpha_s^(6)(Q) for decoupling scale M_Z:   " << asQ << endl;
    dec[0].nf = 6;
    dec[0].Mth = NumDef.Mt;
    dec[0].muth = NumDef.Mt;
    asQ = crd->AlL2AlH(asmz, NumDef.Mz,dec, Q, nloops);
    cout << "alpha_s^(6)(Q) for decoupling scale M_t:   " << asQ << endl;
    dec[0].nf = 6;
    dec[0].Mth = NumDef.Mt;
    dec[0].muth = 2*NumDef.Mt;
    asQ = crd->AlL2AlH(asmz, NumDef.Mz,dec, Q, nloops);
    cout << "alpha_s^(6)(Q) for decoupling scale 2*M_t: " << asQ << endl;
    dec[0].nf = 6;
    dec[0].Mth = NumDef.Mt;
    dec[0].muth = 5*NumDef.Mt;
    asQ = crd->AlL2AlH(asmz, NumDef.Mz,dec, Q, nloops);
    cout << "alpha_s^(6)(Q) for decoupling scale 5*M_t: " << asQ << endl;
    dec[0].nf = 6;
    dec[0].Mth = NumDef.Mt;
    dec[0].muth = Q;
    asQ = crd->AlL2AlH(asmz, NumDef.Mz,dec, Q, nloops);
    cout << "alpha_s^(6)(Q) for decoupling scale Q:     " << asQ << endl;

    return 0;
}
