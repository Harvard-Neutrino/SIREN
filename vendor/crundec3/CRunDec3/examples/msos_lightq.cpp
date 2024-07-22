/*
  CRunDec.cpp, version 3
  msos_lightq.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include "CRunDec.h"


using namespace std;
double asnl6(double mu) {
    double ret = 0.0;
    TriplenfMmu dec[4];
    dec[0].nf   = 6;
    dec[0].Mth  = NumDef.Mt;
    dec[0].muth = 2*NumDef.Mt;
    for(int i = 1; i < 4; i++) {
        dec[i].nf = 0;
    }
    CRunDec* crdtmp = new CRunDec(5);
    ret = crdtmp->AlL2AlH(NumDef.asMz, NumDef.Mz, dec, mu, 5);
    delete crdtmp;
    return ret;
}

double asnl5(double mu) {
    double ret = 0.0;
    CRunDec* crdtmp = new CRunDec(5);
    ret = crdtmp->AlphasExact(NumDef.asMz,NumDef.Mz,mu,5,5);
    delete crdtmp;
    return ret;
}


int main() {
    CRunDec* crd = new CRunDec(5);

    double mut = 163.0; // mt(mt) in GeV
    std::pair<double,double>* mbc = new pair<double,double>[4];
    std::pair<double,double>* mb  = new pair<double,double>[4];
    std::pair<double,double>* mc  = new pair<double,double>[4];
    mbc[0].first  = NumDef.mub;
    mbc[0].second = NumDef.mub;
    mbc[1].first  = NumDef.mc3;
    mbc[1].second = 3.0;
    mbc[2].first  = 0.0;
    mbc[2].second = 0.0;
    mbc[3].first  = 0.0;
    mbc[3].second = 0.0;
    mb[0].first   = NumDef.mub;
    mb[0].second  = NumDef.mub;
    mb[1].first   = 0.0;
    mb[1].second  = 0.0;
    mb[2].first   = 0.0;
    mb[2].second  = 0.0;
    mb[3].first   = 0.0;
    mb[3].second  = 0.0;
    mc[0].first   = NumDef.mc3;
    mc[0].second  = 3.0;
    mc[1].first   = 0.0;
    mc[1].second  = 0.0;
    mc[2].first   = 0.0;
    mc[2].second  = 0.0;
    mc[3].first   = 0.0;
    mc[3].second  = 0.0;

    cout << "Calculate on-shell top quark mass from MS-bar top quark mass." << endl
         << "Separate effects of finite bottom and charm masses." << endl;
    double Mt4  = crd->mMS2mOS(mut, 0, asnl6(mut), mut, 6, 4);
    double Mtb  = crd->mMS2mOS(mut,mb,asnl6(mut),mut, 6, 4);
    double Mtbc = crd->mMS2mOS(mut,mbc,asnl6(mut),mut, 6, 4);
    cout << "M_t = " << Mtbc << endl;
    cout << "without light quark mass effects: " << Mt4 << endl;
    cout << "influence of bottom quark: " << Mtb - Mt4 << endl;
    cout << "influence of charm quark: " << Mtbc - Mtb << endl;
    cout << endl;

    cout << "Calculate on-shell bottom quark mass from MS-bar bottom quark mass." << endl
         << "Separate finite charm mass effects." << endl;
    double Mb4 = crd->mMS2mOS(NumDef.mub, 0, asnl5(NumDef.mub), NumDef.mub, 5, 4);
    double Mbc = crd->mMS2mOS(NumDef.mub, mc, asnl5(NumDef.mub), NumDef.mub, 5, 4);

    cout << "M_b = " << Mbc << endl;
    cout << "without charm quark mass effects: " << Mb4 << endl;
    cout << "influence of charm quark: " << (Mbc - Mb4) << endl;
    return 0;
}
