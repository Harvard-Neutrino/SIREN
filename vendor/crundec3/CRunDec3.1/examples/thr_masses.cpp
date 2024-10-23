/*
  CRunDec.cpp, version 3
  thr_masses.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>
#include <string>
#include "CRunDec.h"


// Compare to arXiv:1606.06754v2
using namespace std;

const double xerr = 0.002;
const double mbSI = NumDef.mub;
const double mtOS = NumDef.Mt;
const double mcSI = NumDef.muc;

double asmzglobal;

struct Mass {
    double CentralVal;
    double UpperErr;
    double LowerErr;
    double AsErr;
    double MSOSErr;
};

int ToInt(string scheme) {
    if(scheme == "1S")
        return 1;
    else if(scheme == "PS")
        return 2;
    else if(scheme == "RS")
        return 3;
    else if(scheme == "RSp")
        return 4;
}

double asnl5(double mu) {
    double ret = 0.0;
    CRunDec* crdtmp = new CRunDec(5);
    ret = crdtmp->AlphasExact(asmzglobal, NumDef.Mz, mu, 5, 4);
    delete crdtmp;
    return ret;
}

double asnl4(double mu) {
    double ret = 0.0;
    CRunDec* crdtmp = new CRunDec(5);
    double astmp = crdtmp->DecAsDownSI(crdtmp->AlphasExact(asmzglobal, NumDef.Mz, 2*mbSI, 5, 4), mbSI, 2*mbSI, 4, 4);
    ret = crdtmp->AlphasExact(astmp, 2*mbSI, mu, 4, 4);
    delete crdtmp;
    return ret;
}

double asnl3(double mu) {
    double ret = 0.0;
    CRunDec* crdtmp = new CRunDec(5);
    double astmp = crdtmp->DecAsDownSI(crdtmp->AlphasExact(asmzglobal, NumDef.Mz, 2*mbSI, 5, 4), mbSI, 2*mbSI, 4, 4);
    double astmp2 = crdtmp->DecAsDownSI(crdtmp->AlphasExact(astmp, 2*mbSI, 3.0, 4, 4), mcSI, 3.0, 3, 4);
    ret = crdtmp->AlphasExact(astmp2, 3.0, mu, 3, 4);
    delete crdtmp;
    return ret;
}


/* Routine to calculate the  scaleinvariant mass for given
alpha_s(M_Z) and 1S, PS, RS or RS' mass including their uncertainties. */
Mass mThr2mSI(double mThr, pair<double,double> merr, double asmz, double aserr, int nl, int nloop, string scheme) {
    double (*asfct)(double) = 0;
    CRunDec* crdtmp = new CRunDec(5);  
    Mass ret;
    double muf;

    switch(nl) {
        case 5:
            asfct = &asnl5;
            muf = 80.0;
            break;
        case 4:
            asfct = &asnl4;
            muf = 2.0;
            break;
        case 3:
            asfct = &asnl3;
            muf = 2.0;
            break;
    }

    switch(ToInt(scheme)) {
        case 1:
            {
            asmzglobal = asmz;           
            ret.CentralVal = crdtmp->m1S2mSI(mThr, crdtmp->mq, asfct, nl, nloop, 1.0);
            ret.MSOSErr = abs(crdtmp->m1S2mSI(mThr, crdtmp->mq, asfct, nl, nloop, 1.0 + xerr) - ret.CentralVal);
            ret.UpperErr = abs(crdtmp->m1S2mSI(mThr + merr.first, crdtmp->mq, asfct, nl, nloop, 1.0) - ret.CentralVal);
            ret.LowerErr = abs(crdtmp->m1S2mSI(mThr - merr.second, crdtmp->mq, asfct, nl, nloop, 1.0) - ret.CentralVal);
            asmzglobal = asmz + aserr;
            double delasplus = abs(crdtmp->m1S2mSI(mThr, crdtmp->mq, asfct, nl, nloop, 1.0) - ret.CentralVal); 
            asmzglobal = asmz - aserr;
            double delasminus = abs(crdtmp->m1S2mSI(mThr, crdtmp->mq, asfct, nl, nloop, 1.0) - ret.CentralVal);
            ret.AsErr = (delasplus > delasminus) ? delasplus : delasminus;
            }
            break;
        case 2:
            asmzglobal = asmz;
            {
            ret.CentralVal = crdtmp->mPS2mSI(mThr, crdtmp->mq, asfct, muf, nl, nloop, 1.0);
            ret.MSOSErr = abs(crdtmp->mPS2mSI(mThr, crdtmp->mq, asfct, muf, nl, nloop, 1.0 + xerr) - ret.CentralVal);
            ret.UpperErr = abs(crdtmp->mPS2mSI(mThr + merr.first, crdtmp->mq, asfct, muf, nl, nloop, 1.0) - ret.CentralVal);
            ret.LowerErr = abs(crdtmp->mPS2mSI(mThr - merr.second, crdtmp->mq, asfct, muf, nl, nloop, 1.0) - ret.CentralVal);
            asmzglobal = asmz + aserr;
            double delasplus = abs(crdtmp->mPS2mSI(mThr, crdtmp->mq, asfct, muf, nl, nloop, 1.0) - ret.CentralVal); 
            asmzglobal = asmz - aserr;
            double delasminus = abs(crdtmp->mPS2mSI(mThr, crdtmp->mq, asfct, muf, nl, nloop, 1.0) - ret.CentralVal);
            ret.AsErr = (delasplus > delasminus) ? delasplus : delasminus;
            }
            break;
        case 3:
            asmzglobal = asmz;         
            {  
            ret.CentralVal = crdtmp->mRS2mSI(mThr, crdtmp->mq, asfct, muf, nl, nloop, 1.0);
            ret.MSOSErr = abs(crdtmp->mRS2mSI(mThr, crdtmp->mq, asfct, muf, nl, nloop, 1.0 + xerr) - ret.CentralVal);
            ret.UpperErr = abs(crdtmp->mRS2mSI(mThr + merr.first, crdtmp->mq, asfct, muf, nl, nloop, 1.0) - ret.CentralVal);
            ret.LowerErr = abs(crdtmp->mRS2mSI(mThr - merr.second, crdtmp->mq, asfct, muf, nl, nloop, 1.0) - ret.CentralVal);
            asmzglobal = asmz + aserr;
            double delasplus = abs(crdtmp->mRS2mSI(mThr, crdtmp->mq, asfct, muf, nl, nloop, 1.0) - ret.CentralVal); 
            asmzglobal = asmz - aserr;
            double delasminus = abs(crdtmp->mRS2mSI(mThr, crdtmp->mq, asfct, muf, nl, nloop, 1.0) - ret.CentralVal);
            ret.AsErr = (delasplus > delasminus) ? delasplus : delasminus;
            }
            break; 
        case 4:
            asmzglobal = asmz; 
            {          
            ret.CentralVal = crdtmp->mRSp2mSI(mThr, crdtmp->mq, asfct, muf, nl, nloop, 1.0);
            ret.MSOSErr = abs(crdtmp->mRSp2mSI(mThr, crdtmp->mq, asfct, muf, nl, nloop, 1.0 + xerr) - ret.CentralVal);
            ret.UpperErr = abs(crdtmp->mRSp2mSI(mThr + merr.first, crdtmp->mq, asfct, muf, nl, nloop, 1.0) - ret.CentralVal);
            ret.LowerErr = abs(crdtmp->mRSp2mSI(mThr - merr.second, crdtmp->mq, asfct, muf, nl, nloop, 1.0) - ret.CentralVal);
            asmzglobal = asmz + aserr;
            double delasplus = abs(crdtmp->mRSp2mSI(mThr, crdtmp->mq, asfct, muf, nl, nloop, 1.0) - ret.CentralVal); 
            asmzglobal = asmz - aserr;
            double delasminus = abs(crdtmp->mRSp2mSI(mThr, crdtmp->mq, asfct, muf, nl, nloop, 1.0) - ret.CentralVal);
            ret.AsErr = (delasplus > delasminus) ? delasplus : delasminus;
            }
            break;  
    }
    delete crdtmp;
    return ret;
}

/* Routine to calculate the MS-mass at a scale mu for given
alpha_s(M_Z) and 1S, PS, RS or RS' mass including their uncertainties. */
Mass mThr2mMS(double mThr, pair<double,double> merr, double asmz, double aserr, double scale, int nl, int nloop, string scheme) {
    double as, asplus, asminus;
    CRunDec* crdtmp = new CRunDec(5);  
    Mass ret;
    double muf;

    switch(nl) {
        case 5:
            asmzglobal = asmz;
            as = asnl5(scale);
            asmzglobal = asmz + aserr;
            asplus = asnl5(scale);
            asmzglobal = asmz - aserr;
            asminus = asnl5(scale);
            muf = 80.0;
            break;
        case 4:
            asmzglobal = asmz;
            as = asnl4(scale);
            asmzglobal = asmz + aserr;
            asplus = asnl4(scale);
            asmzglobal = asmz - aserr;
            asminus = asnl4(scale);
            muf = 2.0;
            break;
        case 3:
            asmzglobal = asmz;
            as = asnl3(scale);
            asmzglobal = asmz + aserr;
            asplus = asnl3(scale);
            asmzglobal = asmz - aserr;
            asminus = asnl3(scale);
            muf = 2.0;
            break;
    }

    switch(ToInt(scheme)) {
        case 1:
            {          
            ret.CentralVal = crdtmp->m1S2mMS(mThr, crdtmp->mq, as, scale, nl, nloop, 1.0);
            ret.MSOSErr = abs(crdtmp->m1S2mMS(mThr, crdtmp->mq, as, scale, nl, nloop, 1.0 + xerr) - ret.CentralVal);
            ret.UpperErr = abs(crdtmp->m1S2mMS(mThr + merr.first, crdtmp->mq, as, scale, nl, nloop, 1.0) - ret.CentralVal);
            ret.LowerErr = abs(crdtmp->m1S2mMS(mThr - merr.second, crdtmp->mq, as, scale, nl, nloop, 1.0) - ret.CentralVal);
            double delasplus = abs(crdtmp->m1S2mMS(mThr, crdtmp->mq, asplus, scale, nl, nloop, 1.0) - ret.CentralVal); 
            double delasminus = abs(crdtmp->m1S2mMS(mThr, crdtmp->mq, asminus, scale, nl, nloop, 1.0) - ret.CentralVal);
            ret.AsErr = (delasplus > delasminus) ? delasplus : delasminus;
            }
            break;
        case 2:
            {
            ret.CentralVal = crdtmp->mPS2mMS(mThr, crdtmp->mq, as, scale, muf, nl, nloop, 1.0);
            ret.MSOSErr = abs(crdtmp->mPS2mMS(mThr, crdtmp->mq, as, scale, muf, nl, nloop, 1.0 + xerr) - ret.CentralVal);
            ret.UpperErr = abs(crdtmp->mPS2mMS(mThr + merr.first, crdtmp->mq, as, scale, muf, nl, nloop, 1.0) - ret.CentralVal);
            ret.LowerErr = abs(crdtmp->mPS2mMS(mThr - merr.second, crdtmp->mq, as, scale, muf, nl, nloop, 1.0) - ret.CentralVal);
            double delasplus = abs(crdtmp->mPS2mMS(mThr, crdtmp->mq, asplus, scale, muf, nl, nloop, 1.0) - ret.CentralVal); 
            double delasminus = abs(crdtmp->mPS2mMS(mThr, crdtmp->mq, asminus, scale, muf, nl, nloop, 1.0) - ret.CentralVal);
            ret.AsErr = (delasplus > delasminus) ? delasplus : delasminus;
            }
            break;
        case 3:       
            {  
            ret.CentralVal = crdtmp->mRS2mMS(mThr, crdtmp->mq, as, scale, muf, nl, nloop, 1.0);
            ret.MSOSErr = abs(crdtmp->mRS2mMS(mThr, crdtmp->mq, as, scale, muf, nl, nloop, 1.0 + xerr) - ret.CentralVal);
            ret.UpperErr = abs(crdtmp->mRS2mMS(mThr + merr.first, crdtmp->mq, as, scale, muf, nl, nloop, 1.0) - ret.CentralVal);
            ret.LowerErr = abs(crdtmp->mRS2mMS(mThr - merr.second, crdtmp->mq, as, scale, muf, nl, nloop, 1.0) - ret.CentralVal);
            double delasplus = abs(crdtmp->mRS2mMS(mThr, crdtmp->mq, asplus, scale, muf, nl, nloop, 1.0) - ret.CentralVal); 
            double delasminus = abs(crdtmp->mRS2mMS(mThr, crdtmp->mq, asminus, scale, muf, nl, nloop, 1.0) - ret.CentralVal);
            ret.AsErr = (delasplus > delasminus) ? delasplus : delasminus;
            }
            break; 
        case 4:
            {          
            ret.CentralVal = crdtmp->mRSp2mMS(mThr, crdtmp->mq, as, scale, muf, nl, nloop, 1.0);
            ret.MSOSErr = abs(crdtmp->mRSp2mMS(mThr, crdtmp->mq, as, scale, muf, nl, nloop, 1.0 + xerr) - ret.CentralVal);
            ret.UpperErr = abs(crdtmp->mRSp2mMS(mThr + merr.first, crdtmp->mq, as, scale, muf, nl, nloop, 1.0) - ret.CentralVal);
            ret.LowerErr = abs(crdtmp->mRSp2mMS(mThr - merr.second, crdtmp->mq, as, scale, muf, nl, nloop, 1.0) - ret.CentralVal);
            double delasplus = abs(crdtmp->mRSp2mMS(mThr, crdtmp->mq, asplus, scale, muf, nl, nloop, 1.0) - ret.CentralVal); 
            double delasminus = abs(crdtmp->mRSp2mMS(mThr, crdtmp->mq, asminus, scale, muf, nl, nloop, 1.0) - ret.CentralVal);
            ret.AsErr = (delasplus > delasminus) ? delasplus : delasminus;
            }
            break;  
    }
    delete crdtmp;
    return ret;
}

int main() {
    CRunDec* crd = new CRunDec(5);
    
    Mass mfromPS, mfrom1S, mfromRS, mfromRSp;
    pair<double,double> zp;
    zp.first  = 0.0;
    zp.second = 0.0;

    cout << "Reproduce table 3 of arXiv:1606.06754:" << endl
         << "(Deviations in the last digit might occur due to rounding" << endl
         << "effects in the input.)" << endl;
    cout << endl << "mt(mt)" << endl;
    cout << "input    mtPS =  mt1S =  mtRS =  mtRSp =" << endl;
    cout << "#loops   168.049 172.060 166.290 171.785" << endl;
    for(int i = 1; i <= 4; i++) {
        mfromPS  = mThr2mSI(168.049, zp, NumDef.asMz, 0.0, 5, i, "PS");
        mfrom1S  = mThr2mSI(172.060, zp, NumDef.asMz, 0.0, 5, i, "1S");
        mfromRS  = mThr2mSI(166.290, zp, NumDef.asMz, 0.0, 5, i, "RS");
        mfromRSp = mThr2mSI(171.785, zp, NumDef.asMz, 0.0, 5, i, "RSp");
        cout << i << "       " << mfromPS.CentralVal  << "  " << mfrom1S.CentralVal
             << "  " << mfromRS.CentralVal  << "  " << mfromRSp.CentralVal << endl;
    }
    cout << "4 (x1.002) " << mfromPS.CentralVal - mfromPS.MSOSErr << "  " <<  mfrom1S.CentralVal - mfrom1S.MSOSErr
         << "  " << mfromRS.CentralVal - mfromRS.MSOSErr << "  " << mfromRSp.CentralVal - mfromRSp.MSOSErr << endl;

    cout << endl << "mb(mb)" << endl;
    cout << "input    mbPS =  mb1S =  mbRS =  mbRSp =" << endl;
    cout << "#loops   4.481   4.668   4.364   4.692"   << endl;
    for(int i = 1; i <= 4; i++) {
        mfromPS  = mThr2mSI(4.481, zp, NumDef.asMz, 0.0, 4, i, "PS");
        mfrom1S  = mThr2mSI(4.668, zp, NumDef.asMz, 0.0, 4, i, "1S");
        mfromRS  = mThr2mSI(4.364, zp, NumDef.asMz, 0.0, 4, i, "RS");
        mfromRSp = mThr2mSI(4.692, zp, NumDef.asMz, 0.0, 4, i, "RSp");
        cout << i << "       " << mfromPS.CentralVal  << "  " << mfrom1S.CentralVal
             << "  " << mfromRS.CentralVal  << "  " << mfromRSp.CentralVal << endl;
    }
    cout << "4 (x1.002) " << mfromPS.CentralVal - mfromPS.MSOSErr << "  " <<  mfrom1S.CentralVal - mfrom1S.MSOSErr
         << "  " << mfromRS.CentralVal - mfromRS.MSOSErr << "  " << mfromRSp.CentralVal - mfromRSp.MSOSErr << endl;

    cout << endl << "mc(mc)" << endl;
    cout << "input    mcPS =  mc1S =  mcRS =  mcRSp =" << endl;
    cout << "#loops   1.130   1.513   1.035   1.351" << endl;
    for(int i = 1; i <= 4; i++) {
        mfromPS  = mThr2mSI(1.130, zp, NumDef.asMz, 0.0, 3, i, "PS");
        mfrom1S  = mThr2mSI(1.513, zp, NumDef.asMz, 0.0, 3, i, "1S");
        mfromRS  = mThr2mSI(1.035, zp, NumDef.asMz, 0.0, 3, i, "RS");
        mfromRSp = mThr2mSI(1.351, zp, NumDef.asMz, 0.0, 3, i, "RSp");
        cout << i << "       " << mfromPS.CentralVal  << "  " << mfrom1S.CentralVal
             << "  " << mfromRS.CentralVal  << "  " << mfromRSp.CentralVal << endl;
    }
    cout << "4 (x1.002) " << mfromPS.CentralVal - mfromPS.MSOSErr << "  " << mfrom1S.CentralVal - mfrom1S.MSOSErr
         << "  " << mfromRS.CentralVal - mfromRS.MSOSErr << "  " << mfromRSp.CentralVal - mfromRSp.MSOSErr << endl;

    cout << endl << "mc(3 GeV)" << endl;
    cout << "input    mcPS =  mc1S =  mcRS =  mcRSp =" << endl;
    cout << "#loops   1.153   1.5145  1.043   1.357"   << endl;
    for(int i = 1; i <= 4; i++) {
        mfromPS  = mThr2mMS(1.153, zp, NumDef.asMz, 0.0, 3.0, 3, i, "PS");
        mfrom1S  = mThr2mMS(1.5145, zp, NumDef.asMz, 0.0, 3.0, 3, i, "1S");
        mfromRS  = mThr2mMS(1.043, zp, NumDef.asMz, 0.0, 3.0, 3, i, "RS");
        mfromRSp = mThr2mMS(1.357, zp, NumDef.asMz, 0.0, 3.0, 3, i, "RSp");
        cout << i << "       " << mfromPS.CentralVal  << "  " << mfrom1S.CentralVal
             << "  " << mfromRS.CentralVal  << "  " << mfromRSp.CentralVal << endl;
    }
    cout << "4 (x1.002) " << mfromPS.CentralVal - mfromPS.MSOSErr << "  " <<  mfrom1S.CentralVal - mfrom1S.MSOSErr
         << "  " << mfromRS.CentralVal - mfromRS.MSOSErr << "  " << mfromRSp.CentralVal - mfromRSp.MSOSErr << endl;


    cout << endl << "Calculate m_t(m_t) for a given PS-mass:" << endl;
    pair<double,double> mtPSerr;
    mtPSerr.first  = 0.1;
    mtPSerr.second = 0.1;
    Mass mtMS      = mThr2mSI(168.049, mtPSerr, NumDef.asMz, 0.0011, 5, 4, "PS");
    cout << "mt_MS = " << mtMS.CentralVal << endl;
    cout << " \\delta_mPS:      +/- " << mtMS.UpperErr << endl;
    cout << " \\delta_\\alpha_s: +/- " << mtMS.AsErr << endl;

    cout << endl << "Calculate m_b(m_b) from PS mass provided in arXiv:1411.3132:" << endl;
    double mbPS      = 4.532;
    pair<double,double> mbPSerr;
    mbPSerr.first    = 0.013;
    mbPSerr.second   = 0.039;
    double alsuncert = 0.0011;
    Mass mbmb        = mThr2mSI(mbPS, mbPSerr, NumDef.asMz, alsuncert, 4, 4, "PS");
    cout << mbmb.CentralVal 
	 << "+" << mbmb.UpperErr << "-" << mbmb.LowerErr << "(\\delta m_PS)"
	 << "+/-" << mbmb.AsErr << "(\\delta\\alpha_s)+/-" << mbmb.MSOSErr << "(4l-MS-OS)" 
	 << "\n=" << mbmb.CentralVal 
	 << "+" << sqrt(mbmb.UpperErr*mbmb.UpperErr + mbmb.AsErr*mbmb.AsErr + mbmb.MSOSErr*mbmb.MSOSErr)
	 << "-" << sqrt(mbmb.LowerErr*mbmb.LowerErr + mbmb.AsErr*mbmb.AsErr + mbmb.MSOSErr*mbmb.MSOSErr)
	 << endl;

    return 0;
}
