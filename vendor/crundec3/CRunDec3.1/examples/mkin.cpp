/*
  CRunDec.cpp, version 3
  mkin.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include "CRunDec.h"


using namespace std;

double as5(double mu) {
    double ret = 0.0;
    CRunDec* crdtmp = new CRunDec(5);
    ret = crdtmp->AlphasExact(0.1179,91.1876,mu,5,4);
    delete crdtmp;
    return ret;
}

double as4(double mu) {
    double ret = 0.0;
    CRunDec* crdtmp = new CRunDec(5);
    ret = crdtmp->AlphasExact(crdtmp->DecAsDownSI(as5(2*4.163), 4.163, 2*4.163, 4, 4),2*4.163, mu, 4, 4);
    delete crdtmp;
    return ret;
}

double as3(double mu) {
    double ret = 0.0;
    CRunDec* crdtmp = new CRunDec(5);
    ret = crdtmp->AlphasExact(crdtmp->DecAsDownMS(as4(3.0), 0.993, 3.0, 3, 4), 3.0, mu, 3, 4);
    delete crdtmp;
    return ret;
}


int main() {
    CRunDec* crd = new CRunDec(5);

    double mbkinin = 4.550;
    double musmcact = 3.0;
    double mufin = 0.5;
    double mus = mbkinin;
    double mc3MSin  = 0.993;
    double mc2MSin  = 1.0986725567495788;
    double mcmcMSin = 1.2785614182835914;
    double mcMSmbMSin = crd->mMS2mMS(0.993,as4(3.0),as4(4.163),4,5);
    double musact = 4.163;
    double mbMSact = crd->mMS2mMS(4.163, as5(4.163), as5(musact), 5, 5);
    std::pair<double,double>* mcMS  = new pair<double,double>[1];
    
    cout << "Convert m_c from MS mass to kinetic mass:" << endl;
    mcMS[0].first = 0.0;
    mcMS[0].second = 0.0;
    
    double murenin = 3.0;    
    cout << "mu_s = " << murenin << ", m_c(mu_s) = " << mc3MSin << " -> m_c^{kin} = " << crd->mMS2mKIN(mc3MSin, mcMS, as3(murenin), murenin, mufin, 3, 3, 3, "") << endl;
    murenin = 2.0;
    cout << "mu_s = " << murenin << ", m_c(mu_s) = " << mc2MSin << " -> m_c^{kin} = " << crd->mMS2mKIN(mc2MSin, mcMS, as3(murenin), murenin, mufin, 3, 3, 3, "") << endl;    
    murenin = mcmcMSin;
    cout << "mu_s = " << murenin << ", m_c(mu_s) = " << mcmcMSin << " -> m_c^{kin} = " << crd->mMS2mKIN(mcmcMSin, mcMS, as3(murenin), murenin, mufin, 3, 3, 3, "") << endl;
    
    cout << endl;

    cout << "Convert m_b from MS mass to kinetic mass:" << endl;
    mcMS[0].first = crd->mMS2mMS(mcMSmbMSin,as4(4.163),as4(musmcact),4,5);
    mcMS[0].second = musmcact;
    mufin = 1.0;
    cout << "m_b(m_b) = " << mbMSact << endl;
    cout << "Case A -> m_b^{kin} = " << crd->mMS2mKIN(mbMSact, mcMS, as3(musact), musact, mufin, 3, "A") << endl;
    cout << "Case B -> m_b^{kin} = " << crd->mMS2mKIN(mbMSact, mcMS, as4(musact), musact, mufin, 3, "B") << endl;
    cout << "Case C -> m_b^{kin} = " << crd->mMS2mKIN(mbMSact, mcMS, as4(musact), musact, mufin, 3, "C") << endl;
    cout << "Case D -> m_b^{kin} = " << crd->mMS2mKIN(mbMSact, mcMS, as3(musact), musact, mufin, 3, "D") << endl;
    
    cout << endl;

    cout << "Convert m_b from kinetic mass to MS mass:" << endl;
    musmcact = 2.0;
    mcMS[0].first = crd->mMS2mMS(mcMSmbMSin,as4(4.163),as4(musmcact),4,5);
    mcMS[0].second = musmcact;
    cout << "m_b^{kin} = " << mbkinin << endl;
    cout << "Case A -> m_b(" << mus << ") = " << crd->mKIN2mMS(mbkinin, mcMS, as3(mus), mus, mufin, 3, "A") << endl;
    cout << "Case B -> m_b(" << mus << ") = " << crd->mKIN2mMS(mbkinin, mcMS, as4(mus), mus, mufin, 3, "B") << endl;
    cout << "Case C -> m_b(" << mus << ") = " << crd->mKIN2mMS(mbkinin, mcMS, as4(mus), mus, mufin, 3, "C") << endl;
    cout << "Case D -> m_b(" << mus << ") = " << crd->mKIN2mMS(mbkinin, mcMS, as3(mus), mus, mufin, 3, "D") << endl;
    return 0;
}
