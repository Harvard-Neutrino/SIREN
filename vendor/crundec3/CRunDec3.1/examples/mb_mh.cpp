/*
  CRunDec.cpp, version 3
  mb_mh.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>
#include "CRunDec.h"

using namespace std;

struct Mass {
    double CentralVal;
    double Err;
};

// Module to calculate mb^{(5)}(M_h) using alpha_s^{(5)}(M_Z) and its uncertainty as input.
Mass mbMh(double as, double aserr, int nloop) {
    CRunDec* crd = new CRunDec(5);
    Mass ret;
    // Calculate mb(10 GeV) and its uncertainties.
    double mb10           = 3.610 - 12./1000.*(as - 0.1189)/0.002;
    double mb10mberrplus  = mb10 + 11./1000.;
    double mb10mberrminus = mb10 - 11./1000.;
    double mb10aserrplus  = mb10 + 12./1000.*aserr/0.002;
    double mb10aserrminus = mb10 - 12./1000.*aserr/0.002;

    // Calculate alpha_s at 10 GeV and Mh with its uncertainties.
    double as10      = crd->AlphasExact(as, NumDef.Mz, 10.0, 5, nloop);
    double asMh      = crd->AlphasExact(as, NumDef.Mz, NumDef.Mh, 5, nloop);
    double as10plus  = crd->AlphasExact(as + aserr, NumDef.Mz, 10.0, 5, nloop);
    double asMhplus  = crd->AlphasExact(as + aserr, NumDef.Mz, NumDef.Mh, 5, nloop);
    double as10minus = crd->AlphasExact(as - aserr, NumDef.Mz, 10.0, 5, nloop);
    double asMhminus = crd->AlphasExact(as - aserr, NumDef.Mz, NumDef.Mh, 5, nloop);

    // Calculate mb(Mh) and its uncertainty.
    ret.CentralVal     = crd->mMS2mMS(mb10, as10, asMh, 5, nloop);
    double mbmhplus    = abs(crd->mMS2mMS(mb10mberrplus, as10, asMh, 5, nloop) - ret.CentralVal);
    double mbmhminus   = abs(crd->mMS2mMS(mb10mberrminus, as10, asMh, 5, nloop) - ret.CentralVal);
    double mbmhasplus  = abs(crd->mMS2mMS(mb10aserrplus, as10plus, asMhplus, 5, nloop) - ret.CentralVal);
    double mbmhasminus = abs(crd->mMS2mMS(mb10aserrminus, as10minus, asMhminus, 5, nloop) - ret.CentralVal);
    ret.Err            = sqrt(max(mbmhplus, mbmhminus)*max(mbmhplus, mbmhminus) +
                              max(mbmhasplus, mbmhasminus)*max(mbmhasplus, mbmhasminus));
    delete crd;
    return ret;
}

int main() {
    CRunDec* crd = new CRunDec(5);
    cout << "Reproduce mb(Mh) from arXiv:1502.00509:" << endl;
    double asmzuncert = 0.0011;
    double asold = 0.1189;
    double asolderr = 0.002;

    Mass mbmh = mbMh(asold, asolderr, 5);
    cout << "mb(Mh) = " << mbmh.CentralVal << " \\pm " << mbmh.Err << endl;
    cout << "Use input values from RunDec.m:" << endl;
    mbmh = mbMh(NumDef.asMz, asmzuncert, 5);
    cout << "mb(Mh) = " << mbmh.CentralVal << " \\pm " << mbmh.Err << endl << endl;

    // Example from the paper
    cout << "Compute mb^(6)(mt) from mb^(5)(10 GeV):" << endl;
    double mb10 = 3.610 - 12./1000.*((NumDef.asMz) - 0.1189)/0.002;
    double mu10 = 10.0;
    int nloops  = 5;
    double as10 = crd->AlphasExact(NumDef.asMz, NumDef.Mz, mu10, 5, nloops);
    TriplenfMmu dec[4];
    dec[0].nf   = 6;
    dec[0].Mth  = NumDef.Mt;
    dec[0].muth = 2*NumDef.Mt;
    for(int i = 1; i < 4; i++) {
        dec[i].nf = 0;
    }
    double mbmt = crd->mL2mH(mb10, as10, mu10, dec, NumDef.Mt, nloops);
    cout << "mb(Mt) = " << mbmt << endl;
    return 0;
}
