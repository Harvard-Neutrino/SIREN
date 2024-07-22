/*
  CRunDec.cpp, version 3
  lambda.cpp
*/

#include <iostream>
#include <cmath>
#include <algorithm>
#include "CRunDec.h"

using namespace std;

const double mudecc = 3.0;
const double mudecb = NumDef.mub;

// vary decoupling scale between mu/f and f*mu:
const double f = 3.0;

struct Lambda {
    double mean;
    double expuncert;
    double truncuncert;
    double diffuncert;
};

void FindExtremaMt(double& max, double& min, CRunDec* crd, double asmz, int n, double f) {
    min = 1.0;
    max = 0.0;
    TriplenfMmu dec[4];
    for(int i = 1; i < 4; i++) {
        dec[i].nf = 0;
    }
    for(int i = 0; i < 10000*(f-1/f); i++) {
        dec[0].nf   = 5;
        dec[0].Mth  = NumDef.Mt;
        dec[0].muth = NumDef.Mt/f + 0.0001*i*NumDef.Mt;
        double as   = crd->AlL2AlH(asmz, NumDef.Mz, dec, NumDef.Mt, n);
        if(as > max) {
            max = as;
        }
        if(as < min) {
            min = as;
        }
    }
}

void FindExtremamub(double& max, double& min, CRunDec* crd, double asmz, int n, double f) {
    min = 1.0;
    max = 0.0;

    for(int i = 0; i < 10000*(f-1/f); i++) {
        double as5thr = crd->AlphasExact(asmz, NumDef.Mz, NumDef.mub/f + 0.0001*i*NumDef.mub, 5, n);
        double as4thr = crd->DecAsDownSI(as5thr, NumDef.mub, NumDef.mub/f + 0.0001*i*NumDef.mub, 4, n);
        double as     = crd->AlphasExact(as4thr, NumDef.mub/f + 0.0001*i*NumDef.mub, NumDef.mub, 4, n);
        if(as > max) {
            max = as;
        }
        if(as < min) {
            min = as;
        }
    }
}

void FindExtremamuc(double& max, double& min, CRunDec* crd, double asmz, int n, double f) {
    min = 1.0;
    max = 0.0;

    for(int i = 0; i < 10000*(f-1/f); i++) {
        double as5thr1 = crd->AlphasExact(asmz, NumDef.Mz, NumDef.mub, 5, n);
        double as4thr1 = crd->DecAsDownSI(as5thr1, NumDef.mub, NumDef.mub, 4, n);
        double as4thr2 = crd->AlphasExact(as4thr1, NumDef.mub, mudecc/f + 0.0001*i*mudecc, 4, n);
        double as3thr2 = crd->DecAsDownSI(as4thr2, NumDef.muc, mudecc/f + 0.0001*i*mudecc, 3, n);
        double as      = crd->AlphasExact(as3thr2, mudecc/f + 0.0001*i*mudecc, 3.0, 3, n);
        if(as > max) {
            max = as;
        }
        if(as < min) {
            min = as;
        }
    }
}

// Calculate alpha_s^{(3)}(mu) from Lambda^{(3)}.
double As3Lambda(double lamin, double scale1, double scale2, int n) {
    CRunDec* crdtmp = new CRunDec(5);
    double as31 = crdtmp->AlphasLam(lamin, scale1, 3, n);
    double as32 = crdtmp->AlphasExact(as31, scale1, scale2, 3, n);
    delete crdtmp;
    return as32;
}

// Calculate alpha_s(M_Z) from alpha_s(mu) using n-loop running and n-loop decoupling.
double As5MZ(double asin, double mu, double thr1, double thr2, int n) {
    CRunDec* crdtmp = new CRunDec(5);
    double as3c  = crdtmp->AlphasExact(asin, mu, thr1, 3, n);
    double as4c  = crdtmp->DecAsUpSI(as3c, NumDef.muc, thr1, 3, n);
    double as4b  = crdtmp->AlphasExact(as4c, thr1, thr2, 4, n);
    double as5b  = crdtmp->DecAsUpSI(as4b, NumDef.mub, thr2, 4, n);
    double as5mz = crdtmp->AlphasExact(as5b, thr2, NumDef.Mz, 5, n);
    delete crdtmp;
    return as5mz;
}


// Calculate alpha_s^{(5)}(M_Z) from Lambda^{(3)}.
double As5MZfromLambda(double lamin, double scale1, double thr1, double thr2, int n) {
    double as3 = As3Lambda(lamin, scale1, thr1, n);
    double as5 = As5MZ(as3, thr1, thr1, thr2, n);
    return as5;
}

void FindExtremaLam(double& max, double& min, CRunDec* crd, double lam3, int n, double f) {
    min = 1.0;
    max = 0.0;

    for(int i = 0; i < 10000*(f-1/f); i++) {
        double as = As5MZfromLambda(lam3, mudecc/f + 0.0001*i*mudecc, mudecc, mudecb, n);
        if(as > max) {
            max = as;
        }
        if(as < min) {
            min = as;
        }
    }
}

void FindExtremamuc2(double& max, double& min, CRunDec* crd, double lam3, int n, double f) {
    min = 1.0;
    max = 0.0;

    for(int i = 0; i < 10000*(f-1/f); i++) {
        double as = As5MZfromLambda(lam3, mudecc, mudecc/f + 0.0001*i*mudecc, mudecb, n);
        if(as > max) {
            max = as;
        }
        if(as < min) {
            min = as;
        }
    }
}

void FindExtremamub2(double& max, double& min, CRunDec* crd, double lam3, int n, double f) {
    min = 1.0;
    max = 0.0;

    for(int i = 0; i < 10000*(f-1/f); i++) {
        double as = As5MZfromLambda(lam3, mudecc, mudecc, mudecb/f + 0.0001*i*mudecb, n);
        if(as > max) {
            max = as;
        }
        if(as < min) {
            min = as;
        }
    }
}


// Calculate alpha_s^{(6)}(M_t) from alpha_s^{(5)}(M_Z).
double As6Mt(double asin, double thr, int n) {
    double ret = 0.0;
    CRunDec* crdtmp = new CRunDec(5);
    double as5thr = crdtmp->AlphasExact(asin, NumDef.Mz, thr, 5, n);
    double as6thr = crdtmp->DecAsUpOS(as5thr, NumDef.Mt, thr, 5, n);
    ret = crdtmp->AlphasExact(as6thr, thr, NumDef.Mt, 6, n);
    delete crdtmp;
    return ret;
}

// Calculate alpha_s^{(4)}(mu_b) from alpha_s^{(5)}(M_Z).
double As4mub(double asin, double thr, int n) {
    CRunDec* crdtmp = new CRunDec(5);
    double as5thr = crdtmp->AlphasExact(asin, NumDef.Mz, thr, 5, n);
    double as4thr = crdtmp->DecAsDownSI(as5thr, NumDef.mub, thr, 4, n);
    double as4mb  = crdtmp->AlphasExact(as4thr, thr, NumDef.mub, 4, n);
    delete crdtmp;
    return as4mb;
}

// Calculate alpha_s^{(3)}(3 GeV) from alpha_s^{(5)}(M_Z).
double As33gev(double asin, double thr1, double thr2, int n) {
    CRunDec* crdtmp = new CRunDec(5);
    double as5thr1 = crdtmp->AlphasExact(asin, NumDef.Mz, thr1, 5, n);
    double as4thr1 = crdtmp->DecAsDownSI(as5thr1, NumDef.mub, thr1, 4, n);
    double as4thr2 = crdtmp->AlphasExact(as4thr1, thr1, thr2, 4, n);
    double as3thr2 = crdtmp->DecAsDownSI(as4thr2, NumDef.muc, thr2, 3, n);
    double as33gev = crdtmp->AlphasExact(as3thr2, thr2, 3.0, 3, n);
    delete crdtmp;
    return as33gev;
}

// Calculate Lambda^{(nf)} from alpha_s^{(nf)}(mu) and its uncertainties.
Lambda CalcLam(double as, double aserr, double mu, int nf, int nloops) {
    CRunDec* crdtmp = new CRunDec(5);
    Lambda ret;
    double expl4     = crdtmp->LamExpl(as, mu, nf, nloops-1);
    double impl4     = crdtmp->LamImpl(as, mu, nf, nloops-1);
    double expl      = crdtmp->LamExpl(as, mu, nf, nloops);
    double impl      = crdtmp->LamImpl(as, mu, nf, nloops);
    double explplus  = abs(crdtmp->LamExpl(as + aserr, mu, nf, nloops) - expl);
    double implplus  = abs(crdtmp->LamImpl(as + aserr, mu, nf, nloops) - impl);
    double explminus = abs(crdtmp->LamExpl(as - aserr, mu, nf, nloops) - expl);
    double implminus = abs(crdtmp->LamImpl(as - aserr, mu, nf, nloops) - impl);
    ret.mean         = (expl + impl)/2.;
    double mean4     = (expl4 + impl4)/2.;
    ret.truncuncert  = abs(ret.mean - mean4);
    ret.diffuncert   = abs(expl - impl)/2.;
    ret.expuncert    = max(max(explplus, implplus), max(explminus, implminus));
    delete crdtmp;
    return ret;
}

int main() {
    CRunDec* crd = new CRunDec(5);

    cout << "Comparison between LamImpl, LamExpl for nf = 5:" << endl << endl;
    cout << "#loops\tLamExpl\tLamImpl\tdifference (MeV)" << endl << endl;

    double expl[5], impl[5], diff[5];
    for(int i = 1; i <= 5; i++) {
        expl[i-1] = crd->LamExpl(NumDef.asMz, NumDef.Mz, 5, i);
        impl[i-1] = crd->LamImpl(NumDef.asMz, NumDef.Mz, 5, i);
        diff[i-1] = abs(expl[i-1] - impl[i-1]);
        cout << i << "\t" << expl[i-1]*1000. << "\t" << impl[i-1]*1000. << "\t"
             << diff[i-1]*1000. << endl;
    }
    cout << endl;
    cout << "Calculate values in table 3." << endl;
    cout << "In the last column in addition the uncertainty from the variation" << endl 
         << "of the matching scale is shown." << endl << endl;
    cout << "nf\t\\Lambda^(nf)\texp. uncert.\t\\delta_trunc\t\\delta_diff\t\\delta_scale" << endl;
    Lambda lam5 = CalcLam(NumDef.asMz, 0.0011, NumDef.Mz, 5, 5);
    cout << "5\t" << lam5.mean << "\t" << lam5.expuncert
         << "\t" << lam5.truncuncert << "\t" << lam5.diffuncert << "\t0.0" << endl;

    double as6      = As6Mt(NumDef.asMz, NumDef.Mt, 5);
    double as6l4    = As6Mt(NumDef.asMz, NumDef.Mt, 4);
    double as6plus  = abs(As6Mt((NumDef.asMz) + 0.0011, NumDef.Mt, 5) - as6);
    double as6minus = abs(As6Mt((NumDef.asMz) - 0.0011, NumDef.Mt, 5) - as6);
    Lambda lam6     = CalcLam(as6, max(as6plus, as6minus), NumDef.Mt, 6, 5);
    Lambda lam6l4   = CalcLam(as6l4, 0, NumDef.Mt, 6, 4);
    double Max, Min;
    FindExtremaMt(Max, Min, crd, NumDef.asMz, 5, f);
    double scaleuncert = Max - Min;
    cout << "6\t" << lam6.mean << "\t" << lam6.expuncert
         << "\t" << abs(lam6.mean - lam6l4.mean) << "\t" << lam6.diffuncert << "\t" << scaleuncert << endl;

    double as4      = As4mub(NumDef.asMz, NumDef.mub, 5);
    double as4l4    = As4mub(NumDef.asMz, NumDef.mub, 4);
    double as4plus  = abs(As4mub((NumDef.asMz) + 0.0011, NumDef.mub, 5) - as4);
    double as4minus = abs(As4mub((NumDef.asMz) - 0.0011, NumDef.mub, 5) - as4);
    Lambda lam4     = CalcLam(as4, max(as4plus, as4minus), NumDef.mub, 4, 5);
    Lambda lam4l4   = CalcLam(as4l4, 0, NumDef.mub, 4, 4);
    FindExtremamub(Max, Min, crd, NumDef.asMz, 5, f);
    scaleuncert = Max - Min;
    cout << "4\t" << lam4.mean << "\t" << lam4.expuncert
         << "\t" << abs(lam4.mean - lam4l4.mean) << "\t" << lam4.diffuncert << "\t" << scaleuncert << endl;

    double as3      = As33gev(NumDef.asMz, NumDef.mub, 3.0, 5);
    double as3l4    = As33gev(NumDef.asMz, NumDef.mub, 3.0, 4);
    double as3plus  = abs(As33gev((NumDef.asMz) + 0.0011, NumDef.mub, 3.0, 5) - as3);
    double as3minus = abs(As33gev((NumDef.asMz) - 0.0011, NumDef.mub, 3.0, 5) - as3);
    Lambda lam3     = CalcLam(as3, max(as3plus, as3minus), 3.0, 3, 5);
    Lambda lam3l4   = CalcLam(as3l4, 0, 3.0, 3, 4);
    FindExtremamuc(Max, Min, crd, NumDef.asMz, 5, f);
    scaleuncert = Max - Min;
    cout << "3\t" << lam3.mean << "\t" << lam3.expuncert
         << "\t" << abs(lam3.mean - lam3l4.mean) << "\t" << lam3.diffuncert << "\t" << scaleuncert << endl;
    cout << endl;

    cout << endl;
    cout << "Compute \\alpha_s similar to arXiv:1701.03075:" << endl << endl;
    double Lambda3      = 0.332;
    double Lambda3plus  = Lambda3 + 0.014;
    double Lambda3minus = Lambda3 - 0.014;
    int nloops = 5;

    double Lambda4 = crd->DecLambdaUp(Lambda3, NumDef.muc, 3, nloops);
    double Lambda5 = crd->DecLambdaUp(Lambda4, NumDef.mub, 4, nloops);
    double aslam   = crd->AlphasLam(Lambda5, NumDef.Mz, 5, nloops);

    double Lambda4l4 = crd->DecLambdaUp(Lambda3, NumDef.muc, 3, nloops-1);
    double Lambda5l4 = crd->DecLambdaUp(Lambda4l4, NumDef.mub, 4, nloops-1);
    double aslaml4   = crd->AlphasLam(Lambda5l4, NumDef.Mz, 5, nloops-1);

    double Lambda4l3 = crd->DecLambdaUp(Lambda3, NumDef.muc, 3, nloops-2);
    double Lambda5l3 = crd->DecLambdaUp(Lambda4l3, NumDef.mub, 4, nloops-2);
    double aslaml3   = crd->AlphasLam(Lambda5l3, NumDef.Mz, 5, nloops-2);

    double Lambda4plus = crd->DecLambdaUp(Lambda3plus, NumDef.muc, 3, nloops);
    double Lambda5plus = crd->DecLambdaUp(Lambda4plus, NumDef.mub, 4, nloops);
    double aslamplus   = crd->AlphasLam(Lambda5plus, NumDef.Mz, 5, nloops);

    double Lambda4minus = crd->DecLambdaUp(Lambda3minus, NumDef.muc, 3, nloops);
    double Lambda5minus = crd->DecLambdaUp(Lambda4minus, NumDef.mub, 4, nloops);
    double aslamminus   = crd->AlphasLam(Lambda5minus, NumDef.Mz, 5, nloops);

    cout << "\\Lambda^{(4)} = " << Lambda4 << "+" << abs(Lambda4 - Lambda4plus)
         << "-" << abs(Lambda4 - Lambda4minus) << endl;
    cout << "\\Lambda^{(5)} = " << Lambda5 << "+" << abs(Lambda5 - Lambda5plus)
         << "-" << abs(Lambda5 - Lambda5minus) << endl;
    cout << "\\alpha_s^{(5)}(M_Z) = " << aslam << endl << " + "
         << abs(aslam - aslamplus) << " - " << abs(aslam - aslamminus) << " (\\delta\\Lambda^(3))"
         << endl << " +/- " << abs(aslam - aslaml4) + abs(aslaml3 - aslaml4) 
         << " (truncation)" << endl << endl;

    cout << "Compute \\alpha_s as proposed in the paper:" << endl << endl;
    double as5  = As5MZfromLambda(Lambda3, mudecc, mudecc, mudecb, nloops);
    double truncuncert = abs(as5 - As5MZfromLambda(Lambda3, mudecc, mudecc, mudecb, nloops-1));
    double expuncert   = max(abs(as5 - As5MZfromLambda(Lambda3plus, mudecc, mudecc, mudecb, nloops)),
                      abs(as5 - As5MZfromLambda(Lambda3minus, mudecc, mudecc, mudecb, nloops)));
    cout << "\\alpha_s(M_Z) with " << 5 << "-loop running and " << 4 << "-loop decoupling:  " << as5 << endl;
    cout << "uncertainty due to \\delta\\Lambda^{(3)}:                   +/-" << expuncert << endl;
    cout << "uncertainty due to truncation of the perturbation series: +/-" << truncuncert << endl;
    FindExtremaLam(Max, Min, crd, Lambda3, nloops, f);
    double scaleuncertlam = Max - Min;
    cout << "scale uncertainty for \\mu_lam,c:                              " << scaleuncertlam << endl;
    FindExtremamuc2(Max, Min, crd, Lambda3, nloops, f);
    double scaleuncertc = Max - Min;
    cout << "scale uncertainty for \\mu_dec,c:                              " << scaleuncertc << endl;
    FindExtremamub2(Max, Min, crd, Lambda3, nloops, f);
    double scaleuncertb = Max - Min;
    double totalscaleuncert = sqrt(scaleuncertc*scaleuncertc + scaleuncertb*scaleuncertb
                                + scaleuncertlam*scaleuncertlam);
    cout << "scale uncertainty for \\mu_dec,b:                              " << scaleuncertb << endl;
    cout << "uncertainty introduced through evolution:                 +/-" << totalscaleuncert << endl << endl;


    cout << "Compute \\alpha_s^{(3)}(M_\\tau) using \\Lambda^{(3)} from arXiv:1701.03075:" << endl << endl;
    double astau       = crd->AlphasLam(Lambda3, NumDef.Mtau, 3, 5);
    double astau4      = crd->AlphasLam(Lambda3, NumDef.Mtau, 3, 4);
    double astauplus   = abs(crd->AlphasLam(Lambda3plus, NumDef.Mtau, 3, 5) - astau);
    double astauminus  = abs(crd->AlphasLam(Lambda3minus, NumDef.Mtau, 3, 5) - astau);
    truncuncert = abs(astau - astau4);
    expuncert   = max(astauplus, astauminus);
    cout << "\\Lambda^{(3)} -> \\alpha_s^{(3)}(M_\\tau)" << endl;
    cout << "\\alpha_s^{(3)}(M_\\tau) =            " << astau << endl;
    cout << "uncertainty due to \\delta\\Lambda^{(3)}: +/-" << expuncert << endl;
    cout << "uncertainty due to truncation:          +/-" << truncuncert << endl << endl;
    double as3Gev      = crd->AlphasLam(Lambda3, 3.0, 3, 5);
    double as3Gev4     = crd->AlphasLam(Lambda3, 3.0, 3, 4);
    double as3Gevplus  = crd->AlphasLam(Lambda3plus, 3.0, 3, 5);
    double as3Gevminus = crd->AlphasLam(Lambda3minus, 3.0, 3, 5);
    astau              = crd->AlphasExact(as3Gev, 3.0, NumDef.Mtau, 3, 5);
    astau4             = crd->AlphasExact(as3Gev4, 3.0, NumDef.Mtau, 3, 4);
    astauplus          = abs(crd->AlphasExact(as3Gevplus, 3.0, NumDef.Mtau, 3, 5) - astau);
    astauminus         = abs(crd->AlphasExact(as3Gevplus, 3.0, NumDef.Mtau, 3, 5) - astau);
    truncuncert        = abs(astau - astau4);
    expuncert          = max(astauplus, astauminus);
    cout << "\\Lambda^{(3)} -> \\alpha_s^{(3)}(3 GeV) -> \\alpha_s^{(3)}(M_\\tau)" << endl << endl;
    cout << "\\alpha_s^{(3)}(M_\\tau) =            " << astau << endl;
    cout << "uncertainty due to \\delta\\Lambda^{(3)}: +/-" << expuncert << endl;
    cout << "uncertainty due to truncation:          +/-" << truncuncert << endl << endl;
    return 0;
}
