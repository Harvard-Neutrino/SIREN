#include <iostream>
#include <cmath>
#include "CRunDec.h"

using namespace std;

// Numerical input values from arXiv:0801.1821:
const double mudecc = 3.0;
const double mudecb = NumDef.mub;
const double asmtauerror = 0.016;

// vary decoupling scale between mu/f and f*mu: 
const double f = 3.0;

void FindExtrema1(double& max, double& min, CRunDec* crd, double asin, double mu, double thr1, double thr2, int n, double f) {
    min = 1.0;
    max = 0.0;
    for(int i = 0; i < 10000*(f-1/f); i++) {
        double as3c = crd->AlphasExact(asin, mu, thr1/f + 0.0001*i*thr1, 3, n);
        double as4c = crd->DecAsUpSI(as3c, NumDef.muc, thr1/f + 0.0001*i*thr1, 3, n);
        double as4b = crd->AlphasExact(as4c, thr1/f + 0.0001*i*thr1, thr2, 4, n);
        double as5b = crd->DecAsUpSI(as4b, NumDef.mub, thr2, 4, n);
        double as   = crd->AlphasExact(as5b, thr2, NumDef.Mz, 5, n);
        if(as > max) {
            max = as;
        }
        if(as < min) {
            min = as;
        }
    }
}

void FindExtrema2(double& max, double& min, CRunDec* crd, double asin, double mu, double thr1, double thr2, int n, double f) {
    min = 1.0;
    max = 0.0;
    for(int i = 0; i < 10000*(f-1/f); i++) {
        double as3c = crd->AlphasExact(asin, mu, thr1, 3, n);
        double as4c = crd->DecAsUpSI(as3c, NumDef.muc, thr1, 3, n);
        double as4b = crd->AlphasExact(as4c, thr1, thr2/f + 0.0001*i*thr2, 4, n);
        double as5b = crd->DecAsUpSI(as4b, NumDef.mub, thr2/f + 0.0001*i*thr2, 4, n);
        double as   = crd->AlphasExact(as5b, thr2/f + 0.0001*i*thr2, NumDef.Mz, 5, n);
        if(as > max) {
            max = as;
        }
        if(as < min) {
            min = as;
        }
    }
}

// Calculate alpha_s(M_Z) from alpha_s(mu) using n-loop running and n-loop decoupling.
double As5MZ(double asin, double mu, double thr1, double thr2, int n) {
    CRunDec* crd = new CRunDec(5);
    double as3c  = crd->AlphasExact(asin, mu, thr1, 3, n);
    double as4c  = crd->DecAsUpSI(as3c, NumDef.muc, thr1, 3, n);
    double as4b  = crd->AlphasExact(as4c, thr1, thr2, 4, n);
    double as5b  = crd->DecAsUpSI(as4b, NumDef.mub, thr2, 4, n);
    double as5mz = crd->AlphasExact(as5b, thr2, NumDef.Mz, 5, n);
    delete crd;
    return as5mz;
}


int main() {
    CRunDec* crd = new CRunDec(5);
    cout << "Calculate \\alpha_s(M_Z) for two- to five-loop running." << endl << endl;
    for(int i = 2; i <= 5; i++) {
        double as5      = As5MZ(NumDef.asMtau, NumDef.Mtau, mudecc, mudecb, i);
        double truncerr = abs(as5 - As5MZ(NumDef.asMtau, NumDef.Mtau, mudecc, mudecb, i-1));
        double max, min;
        FindExtrema1(max, min, crd, NumDef.asMtau, NumDef.Mtau, mudecc, mudecb, i, f);
        double scaleuncertc = max - min;
        FindExtrema2(max, min, crd, NumDef.asMtau, NumDef.Mtau, mudecc, mudecb, i, f);
        double scaleuncertb = max - min;
        double totaluncert  = sqrt(truncerr*truncerr + scaleuncertc*scaleuncertc + scaleuncertb*scaleuncertb);

        cout << "\\alpha_s(M_Z) with " << i << "-loop running and " << i-1 << "-loop decoupling: "
             << as5 << endl;
        cout << "take into account error on \\alpha_s(M_\\tau):   ";
        cout << "         +" << abs(as5 - As5MZ(NumDef.asMtau + asmtauerror, NumDef.Mtau, mudecc, mudecb, i)) << endl
             << "                                                        -"
             << abs(as5 - As5MZ(NumDef.asMtau - asmtauerror, NumDef.Mtau, mudecc, mudecb, i)) << endl;
        cout << "uncertainty due to truncation of the perturbation series:      " << truncerr << endl;
        cout << "decoupling scale uncertainty for \\mu_dec,c:                    " << scaleuncertc << endl;
        cout << "decoupling scale uncertainty for \\mu_dec,b:                    " << scaleuncertb << endl;
        cout << "total uncertainty (truncation and dec. scale):        +/-" << totaluncert << endl;
        cout << endl;
    }
    delete crd;
    return 0;
}
