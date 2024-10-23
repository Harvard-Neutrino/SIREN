/*
  CRunDec.cpp, v3.1

  Authors: Florian Herren and Matthias Steinhauser (May 2021)

  Changes since v3.0:
  May 2021: fix array lengths in AlphasLam and fSetAsL. Thanks to Florian Bernlochner for pointing out this bug.
  Jul 2020: implemented mMS2mKIN and mKIN2mMS (arXiv:2005:06478)
            Coefficients ctil[3][{3,4,5}] updated. Thanks to Christopher Lepenik for providing the new values.
  Jan 2019: Changed license from GPLv3 to MIT
*/

/*
  CRunDec.cpp, v3.0

  Authors: Florian Herren and Matthias Steinhauser (Feb 2017)

  F. Herren, M. Steinhauser,
  "Version 3 of RunDec and CRunDec"
  Comput.Phys.Commun. 224 (2018) 333
  arXiv:1703.03751, TTP17-011

  Changes since v2.1:
  Dec 2016: added five-loop beta function (see arXiv:1606.08659);
            extended LamExpl, LamImpl, AlphasLam to five loops;
            ported DecAsDownMS, DecMqDownMS and DecMq{Up,Down}SI from RunDec;
            extended Dec{As,Mq}{Up,Down}{MS,SI} to four loops;
  Jan 2017: fixed Segmentation fault in case null pointer is passed as "mq" in
            mass conversion routines, minor bugfixes in threshold mass routines
            added DecLambdaUp/Down to 5 loops
            extended mRGI2mMS/mMS2mRGI/mMS2mRGImod to 5 loops
            included light quark mass effects in MS-OS relation up to three loops (arXiv:0708.1729)
*/

/*
  CRunDec.cpp, v2.1

  Authors: Florian Herren and Matthias Steinhauser (Jun 2016)

  Changes since v2.0:
  Jun 2016: improved 4-loop coefficients in mMS2mOS, mOS2mSI, mOS2mMS
*/

/*
  Based on:

  CRunDec.cpp, v2.0

  Authors: Florian Herren and Matthias Steinhauser (Jan 2016)

  Changes since v1.1:
  Sep 2015: extended mOS2mMS, mMS2mOS to four loops (see arXiv:1502.01030);  
            mMS2mMS to five loops (see arXiv:1402.6611); 
            DecMqUp/DownOS to four loops (see hep-ph/0512058, hep-ph/0512060
            and arXiv:1502.04719); 
            prepared AsmMSrunexact, AlphasExact for five loops
  Oct 2015: implemented PS, 1S and RS masses, now requires C++11 for tgamma, 
            extended mOS2mSI to four loops (see arXiv:1502.01030)
  Dec 2015: bugfix: nl -> nlq in fas6to5ms, ported DecAsDownSI/DecAsUpSI from RunDec
  Jan 2016: Add RS' mass relations (see arXiv:1502.01030)
*/

/*
  Based on:

  CRunDec.cpp, v1.1

  Author: Barbara Schmidt

  For documentation see

  CRunDec: a C++ package for running and decoupling of the
  strong coupling and quark masses

  by

  Barbara Schmidt, Matthias Steinhauser

  Comput.Phys.Commun. 183 (2012) 1845-1848
  arXiv:1201.6149 
  SFB/CPP-12-03
  TTP12-02

  See also:
  K.G. Chetyrkin, J.H. Kuhn and M. Steinhauser,
  "RunDec: A Mathematica package for running and decoupling of the strong
  coupling and quark masses"
  Comput. Phys. Commun.  133 (2000) 43,
  arXiv:hep-ph/0004189

  May 2013: minor modifications to avoid "-Wunused-parameter"
            during compilation with g++ -Wall -Wextra
  Sep 2013: bug fix in fRungeKuttaImpl (thanks to Stephen Jones)
  Oct 2014: bug fix in AsmMSrunexact
*/

/*
License:

Copyright (c) 2019 Florian Herren & Matthias Steinhauser

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "CRunDec.h"

// Some constants:
#define cf 4./3.
#define ca 3.
#define tr 1./2.
#define B4 -1.762800087073770864061897634679818807215137274389016762629478603776
#define A4 0.5174790616738993863307581618988629456223774751413792582443193479770
#define A5 0.5084005792422687074591088492585899413195411256648216487244977963526
#define Pi M_PI
#define Zeta2 (Pi*Pi)/6.
#define Zeta3 1.20205690315959428539973816151144999076498629234049888179227155534
#define Zeta4 (Pi*Pi*Pi*Pi)/90.
#define Zeta5 1.03692775514336992633136548645703416805708091950191281197419267790
#define Zeta6 (Pi*Pi*Pi*Pi*Pi*Pi)/945.
#define Zeta7 1.00834927738192282683979754984979675959986356056523870641728313657160147831735573534609696891385132
#define ln2 log(2)

using namespace std;

// Default constructor:
CRunDec::CRunDec(){
     for(int i=0; i<4; i++){
       mq[i].first = 0.;
       mq[i].second = 0.;
       nfMmu[i].Mth=0.;
       nfMmu[i].muth=0.;
       nfMmu[i].nf=0;
     }
     AM.Asexact=0.;
     AM.mMSexact=0.;
}

// Constructor called with number for active flavours nf
CRunDec::CRunDec(int n){
     this->SetConstants(n);
     for(int i=0; i<4; i++) {
       mq[i].first = 0.;
       mq[i].second = 0.;
       nfMmu[i].Mth=0.;
       nfMmu[i].muth=0.;
       nfMmu[i].nf=0;
     } 
     AM.Asexact=0.;
     AM.mMSexact=0.;        
}

// Define constants (if not already done in constructor)
void CRunDec::SetConstants(int n){
     double nf=(double) n;
     Nf=nf;
     double nl = nf-1;
     Beta[0]= (double)0.25*(11.0 - 2.*nf/3.0);
     Beta[1]= (double)(102.0 - 38.*nf/3.0)/16.0;
     Beta[2]= (double)(0.5*2857. - 5033.*nf/18.0 + 325.*nf*nf/54.0)/64.0;
     Beta[3]= (double)(149753./6.0 + 3564.*Zeta3 + 
              (-1078361./162.0 - 6508.*Zeta3/27.0)*nf + 
              (50065./162.0 + 6472.*Zeta3/81.0)*nf*nf + 
              1093.*nf*nf*nf/729.0)/256.0;
     Beta[4]= (double)(8157455.0/16.0 + 621885.0/2.0*Zeta3  - 88209.0/2.0*Zeta4 - 288090.0*Zeta5 +
              (-336460813./1944. - 4811164./81.*Zeta3 + 33935./6.*Zeta4 + 1358995./27.*Zeta5)*nf +
              (25960913./1944. + 698531./81.*Zeta3 - 10526./9.*Zeta4 - 381760./81.*Zeta5)*nf*nf +
              (-630559./5832. - 48722./243.*Zeta3 + 1618./27.*Zeta4 + 460./9.*Zeta5)*nf*nf*nf +
              (1205./2916. - 152./81.*Zeta3)*nf*nf*nf*nf)/1024.0;
     Betap[0]= (double)0.25*(11.0 - 2.*nl/3.0);
     Betap[1]= (double)(102.0 - 38.*nl/3.0)/16.0;
     Betap[2]= (double)(0.5*2857. - 5033.*nl/18.0 + 325.*nl*nl/54.0)/64.0;
     Betap[3]= (double)(149753./6.0 + 3564.*Zeta3 + 
              (-1078361./162.0 - 6508.*Zeta3/27.0)*nl + 
              (50065./162.0 + 6472.*Zeta3/81.0)*nl*nl + 
              1093.*nl*nl*nl/729.0)/256.0;
     Betap[4]= (double)(8157455.0/16.0 + 621885.0/2.0*Zeta3  - 88209.0/2.0*Zeta4 - 288090.0*Zeta5 +
              (-336460813./1944. - 4811164./81.*Zeta3 + 33935./6.*Zeta4 + 1358995./27.*Zeta5)*nl +
              (25960913./1944. + 698531./81.*Zeta3 - 10526./9.*Zeta4 - 381760./81.*Zeta5)*nl*nl +
              (-630559./5832. - 48722./243.*Zeta3 + 1618./27.*Zeta4 + 460./9.*Zeta5)*nl*nl*nl +
              (1205./2916. - 152./81.*Zeta3)*nl*nl*nl*nl)/1024.0;
     
     Gamma[0]=(double)1.0;
     Gamma[1]=(double)(202./3.-20.*nf/9.)/16.; 
     Gamma[2]=(double)(1249. + (-2216./27. - 160.*Zeta3/3.)*nf-
              140.*nf*nf/81.)/64.;
     Gamma[3]=(double)(4603055./162. + 135680.*Zeta3/27. - 8800.*Zeta5 +
              (-91723./27. - 34192.*Zeta3/9. + 
              880.*Zeta4 + 18400.*Zeta5/9.)*nf +
              (5242./243. + 800.*Zeta3/9. - 160.*Zeta4/3.)*nf*nf +
              (-332./243. + 64.*Zeta3/27.)*nf*nf*nf)/256.;
     Gamma[4]=(double)(99512327./162. + 46402466.*Zeta3/243. + 96800.*Zeta3*Zeta3 - 698126.*Zeta4/9.
              -231757160.*Zeta5/243. + 242000.*Zeta6 + 412720.*Zeta7
              +nf*(-150736283./1458. - 12538016.*Zeta3/81. - 75680.*Zeta3*Zeta3/9. + 2038742.*Zeta4/27.
              + 49876180.*Zeta5/243. - 638000.*Zeta6/9. - 1820000.*Zeta7/27.)
              +nf*nf*(1320742./729. + 2010824.*Zeta3/243. + 46400.*Zeta3*Zeta3/27. - 166300.*Zeta4/27. - 264040.*Zeta5/81. + 92000.*Zeta6/27.)
              +nf*nf*nf*(91865./1458. + 12848.*Zeta3/81. + 448.*Zeta4/9. - 5120.*Zeta5/27.)
              +nf*nf*nf*nf*(-260./243 - 320.*Zeta3/243. + 64.*Zeta4/27.)
              )/(4*4*4*4*4);

       
     for(int i=0; i<5; i++) {
       B[i]=Beta[i]/Beta[0];
       Bp[i]=Betap[i]/Betap[0];
       C[i]=Gamma[i]/Beta[0];
     }                          
}


// Function int CRunDec::GetNf()
// Returns number of active flavours.
int CRunDec::GetNf(){
     return (int)Nf;
}

// Function void CRunDec::SetNf(int nf)
// Set the private component Nf to the number of active flavours.
void CRunDec::SetNf(int nf){
     this->SetConstants(nf);
}

// Aux. function to exit function.
int CRunDec::Abbruch(void){
     RETURN
}

// PolyLog[n,x] for -1 < x < 1 and positive integer n
double CRunDec::PolyLog(unsigned int n, double x) {
     if(x >= 1 || x <= -1 || n == 0){
       cout << "PolyLog NOT IMPLEMENTED FOR n = " << n << " and x = " << x << endl;
       RETURN
     }
     double m = n;
     double ret = 0.0;
     double sav1 = 1.0;
     double sav2 = 0.0;
     double diff = 1.0;
     double k = 1.0;
     while(abs(diff) > 1e-16) {
       sav2 = sav1;
       sav1 = pow(x,k)*pow(k,-m);
       ret += sav1;
       diff = sav2 - sav1;
       k += 1.0;
     }
     return ret;
}

// Function double CRunDec::LamExpl(double AlphaS, double Mu, int nl)
// Compute \Lambda using eq.(4) of [RunDec].
double CRunDec::LamExpl(double AlphaS, double Mu, int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN
     }
     double A=AlphaS/Pi;
     double sum[5];
     sum[0]= 1./(A*Beta[0]);
     sum[1]= (B[1]*log(A))/Beta[0]
              + (B[1]/Beta[0])*log(Beta[0]);
     sum[2]= (B[2]*A - B[1]*B[1]*A )/Beta[0];
     sum[3]= (0.5*B[3]*A*A-B[1]*B[2]*A*A+ 0.5*B[1]*B[1]*B[1]*A*A)/Beta[0];
     sum[4]= (B[4]/3.-2.*B[3]*B[1]/3. - B[2]*B[2]/3. + B[2]*B[1]*B[1] - B[1]*B[1]*B[1]*B[1]/3.)*A*A*A/Beta[0];

     double LogM2L2=0.0;
     for(int i=1; i<=nl; i++){
       LogM2L2+=sum[i-1];
     }
               
     double Lambda= Mu*exp(-0.5*LogM2L2);
     return Lambda;
}

// Function double CRunDec::AlphasLam(double Lambda, double Mu, int nl)
// Compute \alpha_s using eq.(5) of [RunDec]. 
double CRunDec::AlphasLam(double Lambda, double Mu, int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(Mu/Lambda<1.5){
       cout<<"WARNING: the ratio \\mu/\\lambda = "<< Mu/Lambda
             <<" is very small!"<<endl; 
       RETURN  
     }
     double L=log(Mu*Mu/(Lambda*Lambda));
     double h=1/(L*Beta[0]);
     double c=log(L);
     double sum[5];
     sum[0]= h;
     sum[1]= -h*h*B[1]*c;
     sum[2]= + h*h*h*(B[1]*B[1]*(c*c-c-1)+B[2]);
     sum[3]= h*h*h*h*(B[1]*B[1]*B[1]*(-c*c*c+2.5*c*c+2*c-0.5) 
                   -3*B[1]*B[2]*c +0.5*B[3]);
     sum[4]= h*h*h*h*h*(B[4]/3. + B[3]*B[1]*(-1./6.-2.*c) + 5.*B[2]*B[2]/3. +
                        B[2]*B[1]*B[1]*(6.*c*c-3.*c-3.)+B[1]*B[1]*B[1]*B[1]*(c*c*c*c-13.*c*c*c/3.-3.*c*c/2.+4.*c+7./6.));
     double a=0.0;
     for(int i=1; i<=nl; i++){
       a+=sum[i-1];
     }
     return a*Pi;
}

// Eq.(5) rewritten in a form suitable to determine zero.
double CRunDec::fSetAsL(double Lambda, double Mu, int nl, double AlphaS){
     double L=log(Mu*Mu/(Lambda*Lambda));
     double h=1/(L*Beta[0]);
     double c=log(L);
     double sum[5];
     sum[0]= h;
     sum[1]= -h*h*B[1]*c;
     sum[2]= + h*h*h*(B[1]*B[1]*(c*c-c-1)+B[2]);
     sum[3]= h*h*h*h*(B[1]*B[1]*B[1]*(-c*c*c+2.5*c*c+2*c-0.5) 
                   -3*B[1]*B[2]*c +0.5*B[3]);
     sum[4]= h*h*h*h*h*(B[4]/3. + B[3]*B[1]*(-1./6.-2.*c) + 5.*B[2]*B[2]/3. +
                        B[2]*B[1]*B[1]*(6.*c*c-3.*c-3.)+B[1]*B[1]*B[1]*B[1]*(c*c*c*c-13.*c*c*c/3.-3.*c*c/2.+4.*c+7./6.));
     double Add=0.0;
     for(int i=1; i<=nl; i++) Add+=sum[i-1];
     return (Add-(AlphaS/Pi));
}

// Function double CRunDec::LamImpl(double AlphaS, double Mu,int nl)
// Compute \Lambda using eq.(5) of [RunDec]. 
double CRunDec::LamImpl(double AlphaS, double Mu,int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double epsilonX= 1e-8;
     double Lambda0=LamExpl(AlphaS,Mu,nl);
     double x0=Lambda0 - 0.2*Lambda0;
     double x1=Lambda0 + 0.2*Lambda0;
     double f0= this->fSetAsL(x0,Mu,nl,AlphaS);
     double f1= this->fSetAsL(x1,Mu,nl,AlphaS);
     if(f0*f1>0){
       cout<<"WARNING: No root can be calculatet!"<<endl;
       RETURN
     }
     double xTest;
     double fTest;
     do{
       xTest= (x0+x1)/2;
       fTest= fSetAsL(xTest,Mu,nl,AlphaS);
       if(f0*fTest<0){x1= xTest;}
       else {x0= xTest;}
     }
     while(abs(x1-x0)>= epsilonX);
     double Lambda=xTest;
     return Lambda;
}

// Right-hand side of differential equation for \alpha_s times 1/\mu^2,
// see eq.(1) of [1].
double fSetdydx(CRunDec S, double A,int nl){ 
     double f=0.0;
     double sum[5];
     double B=A*A;
     sum[0]=-S.Beta[0]*B;    
     sum[1]=-S.Beta[1]*B*A;
     sum[2]=-S.Beta[2]*B*B;
     sum[3]=-S.Beta[3]*B*B*A; 
     sum[4]=-S.Beta[4]*B*B*B;   
     for(int i=1; i<= nl; i++) {
       f+=sum[i-1];
     }
     return (f*2);
}

// Implicit Runge-Kutte step (4th order)
// Call with x=\mu^2, y=\alpha_s(\mu)/\pi, step size h, number of loops nl 
// New y value is returned at x+h
double CRunDec::fRungeKuttaImpl(double &x, double y, double &htry, int nl, 
                              double (*f)(CRunDec, double, int)){
     // Precision
     double eps=1e-10;
     double yerr,ytemp,htemp, hnext;
     double h=htry;
     double xnew; // new variable   
     double k1,k2,k3,k4,k5,k6;
     for(;;){
       k1=h*f(*this,y,nl);
       k2=h*f(*this,y+b21*k1,nl);
       k3=h*f(*this,y+b31*k1+b32*k2,nl);
       k4=h*f(*this,y+b41*k1+b42*k2+b43*k3,nl);
       k5=h*f(*this,y+b51*k1+b52*k2+b53*k3+b54*k4,nl);
       k6=h*f(*this,y+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5,nl);
       // y value at x+h as a sum of the previous value and the
       // correspondingly weighted function evaluations
       ytemp= y+ c1*k1+ c2*k2+ c3*k3+ c4*k4+ c5*k5+ c6*k6;
       // Estimate of uncertainty
       yerr=dc1*k1 + dc2*k2 + dc3*k3 + dc4*k4 + dc5*k5 + dc6*k6;
       double err=0.;
       err=fmax(err,fabs(yerr/eps));
       
       // Uncertainty too big? -> Discard result and reduce step size
       if(err>1.){      
         htemp=0.9*h*pow(err,-0.25);
         if(h>=0.){h=fmax(htemp,0.1*h);}
         else{h=fmin(htemp,0.1*h);}
         xnew=x+h;  // modification to previous code
         //decide whether reduced stepsize is still big enough 
         //(in order to prevent a closed loop)
         if(xnew==x){cout<<"stepsize too small"<<endl; RETURN} 
         continue;          
       }
       else{
         if(err>1.89e-4){
         hnext=0.9*h*pow(err,-0.2);
         }
         // Uncertainty OK? -> take y value, increase h
         else{
           hnext=5.*h;
         }
         x+=h;
         
         y=ytemp;
         htry=hnext;
         break;
       }       
     }
     return y; 
}

// Function: double CRunDec::AlphasExact(double AlphaS0, double Mu0, 
//                           double MuEnd, int nl)
// Compute \alpha_s using eq.(1) of [RunDec]
double CRunDec::AlphasExact(double AlphaS0, double Mu0, double MuEnd, int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN
     }
     double Lambda=LamExpl(AlphaS0,Mu0,nl);
     if(MuEnd/Lambda<1.5){
       cout<<"WARNING: the ratio \\mu/\\lambda = "<< MuEnd/Lambda
         <<" is very small!"<<endl;
       RETURN
     }
     double x,y;
     x=log(Mu0);
     y=AlphaS0/Pi;
     double h;
     
     if(Mu0<MuEnd){
       h=1e-4;
       while(x<(log(MuEnd))){
         y=this->fRungeKuttaImpl(x,y,h,nl,fSetdydx);
         if(x+h>=log(MuEnd)){
           h=log(MuEnd)-x;
         }
       }
     return(y*Pi);
     }
     else{h=-1e-4;}
     while(x>(log(MuEnd))){
       y=this->fRungeKuttaImpl(x,y,h,nl,fSetdydx);
       if(x+h<=log(MuEnd)){
         h=log(MuEnd)-x;
       }
     }
     return(y*Pi);
} 

// Eq.(10) of [RunDec]
double CRunDec::fSetcx(double x, int nl){
     double sum[5];
     sum[0]=1;
     sum[1]=(C[1]-B[1]*C[0])*x;
     sum[2]=0.5*((C[1]-B[1]*C[0])*(C[1]-B[1]*C[0]) + C[2] - B[1]*C[1]+
           B[1]*B[1]*C[0] - B[2]*C[0])*x*x;
     sum[3]=((C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])/6. +
           0.5*(C[1]-B[1]*C[0])*(C[2]-B[1]*C[1]+B[1]*B[1]*C[0]-B[2]*C[0])+
           (C[3]-B[1]*C[2]+B[1]*B[1]*C[1]-B[2]*C[1]-B[1]*B[1]*B[1]*C[0] +
            2.*B[1]*B[2]*C[0] - B[3]*C[0])/3.)*x*x*x;
     sum[4]=((C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])/24.
             +(C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])*(C[2]/2. - B[1]*C[1]/2.+ B[1]*B[1]*C[0]/2. - B[2]*C[0]/2.)/2.
             +(C[2]/2. - B[1]*C[1]/2.+ B[1]*B[1]*C[0]/2. - B[2]*C[0]/2.)*(C[2]/2. - B[1]*C[1]/2.+ B[1]*B[1]*C[0]/2. - B[2]*C[0]/2.)/2.
             +((C[1]-B[1]*C[0])*(C[3]-B[1]*C[2]+B[1]*B[1]*C[1]-B[2]*C[1]-B[1]*B[1]*B[1]*C[0] + 2.*B[1]*B[2]*C[0] - B[3]*C[0])/3.)
             +B[1]*B[1]*B[1]*B[1]*C[0]/4. - 3.*B[1]*B[1]*B[2]*C[0]/4. + B[2]*B[2]*C[0]/4. + B[1]*B[3]*C[0]/2. - B[4]*C[0]/4.
             -B[1]*B[1]*B[1]*C[1]/4. + B[1]*B[2]*C[1]/2. - B[3]*C[1]/4. + B[1]*B[1]*C[2]/4. - B[2]*C[2]/4. - B[1]*C[3]/4. + C[4]/4.)*x*x*x*x;
     double erg=0.0;
     for(int i=1; i<=nl; i++){
       erg+=sum[i-1];
     }
     return (pow(x,C[0])*erg);             
}

// Function double CRunDec::mMS2mMS(double mmu0, double AlphaS0, 
//                          double AlphaSEnd, int nl)
// Compute m_q(\mu) using eqs.(9) and (10) of [RunDec]
double CRunDec::mMS2mMS(double mmu0, double AlphaS0, double AlphaSEnd, int nl){
     if(nl<0||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl==0){
       return mmu0;
     }
     double cAlphaS0= this->fSetcx(AlphaS0/Pi, nl);
     double cAlphaSEnd= this->fSetcx(AlphaSEnd/Pi, nl);
     return mmu0*cAlphaSEnd/cAlphaS0;        
}

// Aux. functions (r.h.s of diff. eqs. for alpha_s and m_q)
double fSetdydxM1(CRunDec S, double A, double M){
     return (M*(S.Gamma[0])/(S.Beta[0]*A));
}

double fSetdydxa1(CRunDec S, double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A));      
}

double fSetdydxM2(CRunDec S, double A, double M){
     return (M*(S.Gamma[0]+S.Gamma[1]*A)/(S.Beta[0]*A+S.Beta[1]*A*A));   
}

double fSetdydxa2(CRunDec S, double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A+S.Beta[1]*A*A*A));      
}

double fSetdydxM3(CRunDec S, double A, double M){
     return (M*(S.Gamma[0]+S.Gamma[1]*A+S.Gamma[2]*A*A)/
           (S.Beta[0]*A+S.Beta[1]*A*A+S.Beta[2]*A*A*A));     
}

double fSetdydxa3(CRunDec S, double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A+S.Beta[1]*A*A*A+S.Beta[2]*A*A*A*A));      
}

double fSetdydxM4(CRunDec S, double A, double M){
     return (M*(S.Gamma[0]+S.Gamma[1]*A+S.Gamma[2]*A*A+S.Gamma[3]*A*A*A)/
           (S.Beta[0]*A+S.Beta[1]*A*A+S.Beta[2]*A*A*A+S.Beta[3]*A*A*A*A)); 
}

double fSetdydxa4(CRunDec S, double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A+S.Beta[1]*A*A*A+S.Beta[2]*A*A*A*A+
                S.Beta[3]*A*A*A*A*A));      
}

double fSetdydxM5(CRunDec S, double A, double M){
     return (M*(S.Gamma[0]+S.Gamma[1]*A+S.Gamma[2]*A*A+S.Gamma[3]*A*A*A + S.Gamma[4]*A*A*A*A)/
           (S.Beta[0]*A+S.Beta[1]*A*A+S.Beta[2]*A*A*A+S.Beta[3]*A*A*A*A+S.Beta[4]*A*A*A*A*A)); 
}

double fSetdydxa5(CRunDec S, double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A+S.Beta[1]*A*A*A+S.Beta[2]*A*A*A*A+
                S.Beta[3]*A*A*A*A*A+S.Beta[4]*A*A*A*A*A*A));      
}  

// Runge-Kutta step for implicit procedure
double CRunDec::fRKSchritt(double x, double y, double h, double &yerr,
                    double (*f)(CRunDec, double , double)){
     double k1,k2,k3,k4,k5,k6; 
     k1=h*f(*this,x,y);
     k2=h*f(*this,x+a2*h,y+b21*k1);
     k3=h*f(*this,x+a3*h,y+b31*k1+b32*k2);
     k4=h*f(*this,x+a4*h,y+b41*k1+b42*k2+b43*k3);
     k5=h*f(*this,x+a5*h,y+b51*k1+b52*k2+b53*k3+b54*k4);
     k6=h*f(*this,x+a6*h,y+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5);

     yerr=(dc1*k1 + dc2*k2 + dc3*k3 + dc4*k4 + dc5*k5 + dc6*k6);
     return (y+ c1*k1+ c2*k2+ c3*k3+ c4*k4+ c5*k5+ c6*k6);      
}  

// Function AsmMS CRunDec::AsmMSrunexact(double mMu, double AlphaS0, double Mu0, 
//                         double MuEnd, int nl)
// Compute \alpha_s and m_q solving diff. eqs. simultaneously.
AsmMS CRunDec::AsmMSrunexact(double mMu, double AlphaS0, double Mu0, 
                             double MuEnd, int nl){
     AsmMS Erg;
     if(nl<0||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       this->Abbruch();
     }
     if(nl==0){
       Erg.Asexact=AlphaS0;
       Erg.mMSexact=mMu;
       return Erg;
     }
     double yerr0,ytemp0,xnew,h=1e-3;
     float eps=1e-15;
     float errmax;
     double x0=log(Mu0);
     double y0=AlphaS0/Pi;
     double xEnd=log(MuEnd);
     double yscal0=abs(x0)+abs(h*y0);
     
     double (*falpha)(CRunDec, double , double);
     double (*fmMS)(CRunDec, double , double);
     
     if(nl==1){
       falpha=fSetdydxa1;
       fmMS=fSetdydxM1;
     }
     if(nl==2){
       falpha=fSetdydxa2;
       fmMS=fSetdydxM2;
     } 
     if(nl==3){
       falpha=fSetdydxa3;
       fmMS=fSetdydxM3;
     }
     if(nl==4){
       falpha=fSetdydxa4;
       fmMS=fSetdydxM4;
     }
     if(nl==5){
       falpha=fSetdydxa5;
       fmMS=fSetdydxM5;
     }
     
     if(Mu0<MuEnd){
       h=1e-2;
       while(x0<xEnd){
         for(;;){
           ytemp0=fRKSchritt(x0,y0,h,yerr0,falpha);
           errmax=0.;
           errmax=fmax(errmax,fabs((float)yerr0/yscal0));
           errmax/=eps;
           if(errmax>1){
             h*=0.9;
             xnew=x0+h;
             if(xnew==x0){cout<<"stepsize too small!"<<endl;}
             continue;
           } //if
           else{
             x0+=h;
             y0=ytemp0;  
             if(errmax>1.89e-4){h=0.9*h*pow((double)errmax,-0.2);}
             else{h=5.*h;}
             break;
           }//else
         }//for
         if(x0+h>=xEnd){
           h=xEnd-x0;
         }//if
       }//while
       Erg.Asexact=y0*Pi;   
       
       x0=AlphaS0/Pi;
       xEnd=y0;
       y0=mMu;
       yscal0=abs(x0)+abs(h*y0);
       eps=1e-10;
       h=-1e-3;

       while(x0>xEnd){
         for(;;){
           ytemp0=fRKSchritt(x0,y0,h,yerr0,fmMS);
           errmax=0.;
           errmax=fmax(errmax,fabs((float)yerr0/yscal0));
           errmax/=eps;
           if(errmax>1){
             h*=0.9;
             xnew=x0+h;
             if(xnew==x0){cout<<"stepsize too small!"<<endl;}
             continue;
           } //if
           else{
             x0+=h;
             y0=ytemp0;
             if(errmax>1.89e-4){h=0.9*h*pow((double)errmax,-0.2);}
             else{h=5.*h;}
             break;
           }//else
         }//for
         if(x0+h<=xEnd){
           h=xEnd-x0;
         }//if
       }//while  
       Erg.mMSexact=y0;
       return Erg;
     }//if
     else{
       h=-1e-2;
       while(x0>xEnd){
         for(;;){
           ytemp0=fRKSchritt(x0,y0,h,yerr0,falpha);
           errmax=0.;
           errmax=fmax(errmax,fabs((float)yerr0/yscal0));
           errmax/=eps;
           if(errmax>1){
             h*=0.9;
             xnew=x0+h;
             if(xnew==x0){cout<<"stepsize too small!"<<endl;}
             continue;
           } //if
           else{
             x0+=h;
             y0=ytemp0;
             if(errmax>1.89e-4){h=0.9*h*pow((double)errmax,-0.2);}
             else{h=5.*h;}
             break;
           }//else
         }//for
         if(x0+h<=xEnd){
           h=xEnd-x0;
         }//if
       }//while
       Erg.Asexact=y0*Pi; 

       x0=AlphaS0/Pi;     
       xEnd=y0;
       y0=mMu;
       yscal0=abs(x0)+abs(h*y0);
       eps=1e-10;
       h=1e-3;
       while(x0<xEnd){
         for(;;){
           ytemp0=fRKSchritt(x0,y0,h,yerr0,fmMS);
           errmax=0.;
           errmax=fmax(errmax,fabs((float)yerr0/yscal0));
           errmax/=eps;
           if(errmax>1){
             h*=0.9;
             xnew=x0+h;
             if(xnew==x0){cout<<"stepsize too small!"<<endl;}
             continue;
           } //if
           else{
             x0+=h;
             y0=ytemp0;
             if(errmax>1.89e-4){h=0.9*h*pow((double)errmax,-0.2);}
             else{h=5.*h;}
             break;
           }//else
         }//for
         if(x0+h>=xEnd){
           h=xEnd-x0;
         }//if
       }//while
       Erg.mMSexact=y0;
       return Erg;
     }//else
}   

// Coefficients of eq.(13) of [RunDec]
double CRunDec::fMsFromOs1(double mu, double M){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) (-cf - (3.*cf*lmM)/4.);
     return erg;
}

double CRunDec::fMsFromOs2(double mu, double M, double nl){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) ((-1111.*ca*cf)/384. + (7.*cf*cf)/128. -
     (185.*ca*cf*lmM)/96. + (21.*cf*cf*lmM)/32. - (11.*ca*cf*lmM*lmM)/32. +
     (9.*cf*cf*lmM*lmM)/32. + (143.*cf*tr)/96. + (13.*cf*lmM*tr)/24. + 
     (cf*lmM*lmM*tr)/8. +
     (71.*cf*nl*tr)/96. + (13.*cf*lmM*nl*tr)/24. + (cf*lmM*lmM*nl*tr)/8. +
     (ca*cf*Zeta2)/2 - (15.*cf*cf*Zeta2)/8. - (3.*ca*cf*log(2)*Zeta2)/2. + 
     3.*cf*cf*log(2)*Zeta2 -
     cf*tr*Zeta2 + (cf*nl*tr*Zeta2)/2. + (3.*ca*cf*Zeta3)/8. - 
     (3.*cf*cf*Zeta3)/4.);
     return erg;
}

double CRunDec::fMsFromOs3(double mu, double M, double nl){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) (lmM*lmM*(-2341.*ca*ca*cf + 1962.*ca*cf*cf - 243.*cf*cf*cf 
     + 1492.*ca*cf*tr -
      468.*cf*cf*tr + 1492.*ca*cf*nl*tr - 468.*cf*cf*nl*tr - 208.*cf*tr*tr -
      416.*cf*nl*tr*tr - 208.*cf*nl*nl*tr*tr))/1152. +
     (lmM*lmM*lmM*(-242.*ca*ca*cf + 297.*ca*cf*cf - 81.*cf*cf*cf + 
     176.*ca*cf*tr - 108.*cf*cf*tr + 176.*ca*cf*nl*tr - 108.*cf*cf*nl*tr - 
     32.*cf*tr*tr - 64.*cf*nl*tr*tr - 32.*cf*nl*nl*tr*tr))/1152. +
     (lmM*(-105944.*ca*ca*cf + 52317.*ca*cf*cf - 13203.*cf*cf*cf + 
     74624.*ca*cf*tr -
     5436.*cf*cf*tr + 55616.*ca*cf*nl*tr + 2340.*cf*cf*nl*tr - 
     12608.*cf*tr*tr -
     18304.*cf*nl*tr*tr - 5696.*cf*nl*nl*tr*tr + 12672.*ca*ca*cf*Zeta2 -
     52704.*ca*cf*cf*Zeta2 + 19440.*cf*cf*cf*Zeta2 - 
     38016.*ca*ca*cf*log(2)*Zeta2 +
     91584.*ca*cf*cf*log(2)*Zeta2 - 31104.*cf*cf*cf*log(2)*Zeta2 - 
     29952.*ca*cf*tr*Zeta2 +
     27648.*cf*cf*tr*Zeta2 + 13824.*ca*cf*log(2)*tr*Zeta2 - 
     27648.*cf*cf*log(2)*tr*Zeta2 +
     8064.*ca*cf*nl*tr*Zeta2 + 12096.*cf*cf*nl*tr*Zeta2 + 
     13824.*ca*cf*log(2)*nl*tr*Zeta2 -
     27648.*cf*cf*log(2)*nl*tr*Zeta2 + 9216.*cf*tr*tr*Zeta2 + 
     4608.*cf*nl*tr*tr*Zeta2 -
     4608.*cf*nl*nl*tr*tr*Zeta2 + 9504.*ca*ca*cf*Zeta3 - 
     22896.*ca*cf*cf*Zeta3 +
     7776.*cf*cf*cf*Zeta3 + 6912.*ca*cf*tr*Zeta3 - 3456.*cf*cf*tr*Zeta3 +
     6912.*ca*cf*nl*tr*Zeta3 - 3456.*cf*cf*nl*tr*Zeta3))/13824.;
     return erg;
}

//z[3,m](M) according to eq.(15)
double CRunDec::fZmM(double nl){
     double erg;
     erg= (double) -9478333./93312. + 55.*log(2)*log(2)*log(2)*log(2)/162. +
            (-644201./6480. + 587.*log(2)/27. + 44.*log(2)*log(2)/27.)*Zeta2 -
            61.*Zeta3/27. + 3475*Zeta4/432. + 1439.*Zeta2*Zeta3/72. -
            1975.*Zeta5/216. + 220.*A4/27. + nl*(246643./23328. - 
            log(2)*log(2)*log(2)*log(2)/81. +(967./108. + 22.*log(2)/27. -
            4.*log(2)*log(2)/27.)*Zeta2 + 241.*Zeta3/72. - 305.*Zeta4/108. -
            8.*A4/27.) + nl*nl*(-2353./23328. - 13.*Zeta2/54 - 7.*Zeta3/54.);
     return erg;
}

double CRunDec::fMsFromOs4(double mu, double M, double nl, double err){
     double lmM=log((mu*mu)/(M*M));
     double erg;

     erg =  - 3654.15040757339*err - 1524.2292266911543*lmM - 288.778291935394*lmM*lmM - 32.54735725308642*lmM*lmM*lmM - 
            1.85546875*lmM*lmM*lmM*lmM +
            nl*nl*nl*(0. + 0.678141025604516*err + 0.3205521521864135*lmM + 0.0800290327210927*lmM*lmM +      
            0.010030864197530864*lmM*lmM*lmM) + 
            nl*nl*(0. - 43.48241924867489*err - 19.82672048099557*lmM - 4.482957520194182*lmM*lmM - 
            0.5270061728395061*lmM*lmM*lmM - 0.04108796296296297*lmM*lmM*lmM*lmM) +
            nl*(0. + 756.9421565599532*err + 330.1770776731065*lmM + 67.99849534415492*lmM*lmM + 
            7.595293209876542*lmM*lmM*lmM + 0.48119212962962954*lmM*lmM*lmM*lmM);
     return erg;
}

// Function: double CRunDec::mOS2mMS(double mOS, std::pair<double,double>* mq, double asmu,
//                           double Mu, int nl, double fdelm)
double CRunDec::mOS2mMS(double mOS, std::pair<double,double>* mq, double asmu, double Mu,int nl, double fdelm){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN 
     }
     if(nl == 4 && (Nf<4||Nf>6)){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR 4 LOOPS and "<< nl << " FLAVORS" <<endl;
       RETURN  
     }
     double sum[5];
     double deltalight = deltamOS2mMS(mOS, mq, asmu/Pi, Mu, Nf-1, nl);
     sum[0]=(double) 1.;
     sum[1]=asmu*(this ->fMsFromOs1(Mu, mOS))/Pi;
     sum[2]=asmu*asmu*(this-> fMsFromOs2(Mu, mOS, Nf-1))/(Pi*Pi); 
     sum[3]=asmu*asmu*asmu*(this-> fMsFromOs3(Mu, mOS, Nf-1)+this->fZmM(Nf-1))/(Pi*Pi*Pi);
     sum[4]=asmu*asmu*asmu*asmu*(this->fMsFromOs4(Mu, mOS, Nf-1, fdelm))/(Pi*Pi*Pi*Pi);
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     erg += deltalight;
     return mOS*erg;       
}

// Function: double CRunDec::mMS2mSI(double mMS, double asmu, double mu, int nl)
double CRunDec::mMS2mSI(double mMS, double asmu, double mu, int nl){   
  double epsilonX = 1e-8;
  AsmMS asmq;
  asmq.Asexact  = asmu;
  asmq.mMSexact = mMS;
  for (;;) {
    double mbold = asmq.mMSexact;
    asmq = AsmMSrunexact(mMS, asmu, mu, mbold, nl);
    if (abs(asmq.mMSexact - mbold) < epsilonX) break;
  }
  return asmq.mMSexact;
}

// Coefficients of eq.(18) of [RunDec]
double CRunDec::fMsFromRi1(void){
     return (double) -4./3.;
}
double CRunDec::fMsFromRi2(void){
     double erg= (double) -995./72. + 19.*Zeta3/6. + 89.*Nf/144.;
     return erg;       
}
double CRunDec::fMsFromRi3(void){
     double erg= (double) -6663911./41472. + 408007.*Zeta3/6912. -185.*Zeta5/36.
                + (118325./7776. + 5*Zeta4/12. - 617.*Zeta3/216.)*Nf +
                (-4459./23328. - Zeta3/54.)*Nf*Nf;
     return erg;                    
}

// Function: double CRunDec::mRI2mMS(double mRI, double asmu, int nl)
double CRunDec::mRI2mMS(double mRI, double asmu, int nl){
     if(nl<0||nl>3){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double sum[4];
     sum[0]=(double) 1.;
     sum[1]=asmu*(this ->fMsFromRi1())/Pi;
     sum[2]=asmu*asmu*(this-> fMsFromRi2())/(Pi*Pi);
     sum[3]=asmu*asmu*asmu*(this-> fMsFromRi3())/(Pi*Pi*Pi);
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{
       for(int i=0; i<=nl; i++){
        erg+=sum[i];
       }
     }
     return mRI*erg;         
}

// Function: double CRunDec::mMS2mRGI(double mMS, double asmu, int nl)
double CRunDec::mMS2mRGI(double mMS, double asmu, int nl){
     if(nl<0||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl==0){
       return (double) mMS;
     }
     else{    
       double cAsmu= this->fSetcx(asmu/Pi, nl); 
       return (double) mMS/cAsmu;
     }
}

// Function: double CRunDec::mRGI2mMS(double mRGI, double asmu, int nl)
double CRunDec::mRGI2mMS(double mRGI, double asmu, int nl){
     if(nl<0||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN   
     }  
     if(nl==0){
       return (double) mRGI;
     }
     double cAsmu= this->fSetcx(asmu/Pi, nl); 
     return (double) mRGI*cAsmu;             
}

// Coefficients of eq.(17) of [RunDec]
double CRunDec::fOsFromMs1(double mu, double M){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) (cf + (3.*cf*lmM)/4.);
     return erg;
}

double CRunDec::fOsFromMs2(double mu, double M, double nl){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) ((1111.*ca*cf)/384. - (71.*cf*cf)/128. -
     (143.*cf*tr)/96. - (71.*cf*nl*tr)/96. +
     lmM*((185.*ca*cf)/96. - (9.*cf*cf)/32. - (13.*cf*tr)/24. - 
     (13.*cf*nl*tr)/24.) +
     lmM*lmM*((11.*ca*cf)/32. + (9.*cf*cf)/32. - (cf*tr)/8. - (cf*nl*tr)/8.) -
     (ca*cf*Zeta2)/2. + (15.*cf*cf*Zeta2)/8. + (3.*ca*cf*log(2)*Zeta2)/2. 
     - 3.*cf*cf*log(2)*Zeta2 +
     cf*tr*Zeta2 - (cf*nl*tr*Zeta2)/2. - (3.*ca*cf*Zeta3)/8. + 
     (3.*cf*cf*Zeta3)/4.);
     return erg;
}

double CRunDec::fOsFromMs3(double mu, double M, double nl){
     double lmM=log((mu*mu)/(M*M));
     double erg;

     erg= (double) (lmM*lmM*lmM*((121.*ca*ca*cf)/576. + (33.*ca*cf*cf)/128. +
     (9.*cf*cf*cf)/128. - (11.*ca*cf*tr)/72. - (3.*cf*cf*tr)/32. - 
     (11.*ca*cf*nl*tr)/72. -
     (3.*cf*cf*nl*tr)/32. + (cf*tr*tr)/36. + (cf*nl*tr*tr)/18. +
     (cf*nl*nl*tr*tr)/36.) + lmM*lmM*((2341.*ca*ca*cf)/1152. 
     + (21.*ca*cf*cf)/64. -
     (63.*cf*cf*cf)/128. - (373.*ca*cf*tr)/288. - (3.*cf*cf*tr)/32. -
     (373.*ca*cf*nl*tr)/288. - (3.*cf*cf*nl*tr)/32. + (13.*cf*tr*tr)/72. +
     (13.*cf*nl*tr*tr)/36. + (13.*cf*nl*nl*tr*tr)/72.) +
     lmM*((13243.*ca*ca*cf)/1728. - (4219.*ca*cf*cf)/1536. + 
     (495.*cf*cf*cf)/512. -
     (583.*ca*cf*tr)/108. - (307.*cf*cf*tr)/384. - (869.*ca*cf*nl*tr)/216. -
     (91.*cf*cf*nl*tr)/384. + (197.*cf*tr*tr)/216. + (143.*cf*nl*tr*tr)/108. +
     (89.*cf*nl*nl*tr*tr)/216. - (11.*ca*ca*cf*Zeta2)/12. + 
     (49.*ca*cf*cf*Zeta2)/16. +
     (45.*cf*cf*cf*Zeta2)/32. + (11.*ca*ca*cf*log(2)*Zeta2)/4. - 
     (35.*ca*cf*cf*log(2)*Zeta2)/8. -
     (9.*cf*cf*cf*log(2)*Zeta2)/4. + (13.*ca*cf*tr*Zeta2)/6. - 
     (cf*cf*tr*Zeta2)/2. -
     ca*cf*log(2)*tr*Zeta2 + 2.*cf*cf*log(2)*tr*Zeta2 - 
     (7.*ca*cf*nl*tr*Zeta2)/12. -
     (13.*cf*cf*nl*tr*Zeta2)/8. - ca*cf*log(2)*nl*tr*Zeta2 + 
     2.*cf*cf*log(2)*nl*tr*Zeta2 -
     (2.*cf*tr*tr*Zeta2)/3. - (cf*nl*tr*tr*Zeta2)/3. + 
     (cf*nl*nl*tr*tr*Zeta2)/3. -
     (11.*ca*ca*cf*Zeta3)/16. + (35.*ca*cf*cf*Zeta3)/32. + 
     (9.*cf*cf*cf*Zeta3)/16. -
     (ca*cf*tr*Zeta3)/2. + (cf*cf*tr*Zeta3)/4. - (ca*cf*nl*tr*Zeta3)/2. +
     (cf*cf*nl*tr*Zeta3)/4.));
     return erg;
}

double CRunDec::fOsFromMs4(double mu, double M, double nl, double err){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg = 3567.602784989066*err + 1727.2260148986106*lmM + 409.2429990574718*lmM*lmM + 66.93663194444443*lmM*lmM*lmM + 
           8.056278935185185*lmM*lmM*lmM*lmM + 
           nl*nl*nl*(-0.678141025604516*err - 0.3205521521864134*lmM - 0.0800290327210927*lmM*lmM - 
           0.010030864197530864*lmM*lmM*lmM) + 
           nl*(-745.7207145811878*err - 358.29765085086774*lmM - 87.39262571554698*lmM*lmM - 
           11.883873456790122*lmM*lmM*lmM - 1.2705439814814814*lmM*lmM*lmM*lmM) + 
           nl*nl*(43.396250117985666*err + 20.528466368867228*lmM + 4.971905254812516*lmM*lmM + 
           0.6304012345679011*lmM*lmM*lmM + 0.06655092592592593*lmM*lmM*lmM*lmM);
     return erg;
}

// Compute 2 Loop contribution to mOS2mMS at 3 loops for MS light quarks
double CRunDec::deltamOS2mMS(double mOS, std::pair<double,double>* mq, double asmu, double mu, int nlq, int nloops){
     double erg=0.0;
     if(!mq)
       return 0.0;

     for(int i = 0; i < 4; i++) {
       if(mq[i].first == 0.0)
        continue;

       double x = mq[i].first/mOS;
       double muf = mq[i].second;
       int nl = nlq-i;
       if(nloops >= 2) {
         erg += asmu*asmu*2./3.*(48.*x*x*x*x*log(x)*log(x) + 48.*x*x*log(x)+72.*x*x + 8.*Pi*Pi*(x*x*x*x-3.*x*x*x-3.*x)
                -48.*(x+1.)*(x+1.)*(x*x - x + 1.)*(log(x)*log(x+1) + PolyLog(2,-x))
                -48.*(x-1.)*(x-1.)*(x*x + x + 1.)*(log(x)*log(1-x) + PolyLog(2,x)))/96.;
       }
       if(nloops >= 3) {
         erg += asmu*asmu*asmu*(-21.8714*x - 4.348*x*x - 1.02211*x*x*x - 0.0493333*x*x*x*x + nl*(0.982667*x + 0.300333*x*x)
              + log(mu*mu/(mOS*mOS))*(-6.61056*x + 2.46511*x*x - 0.724333*x*x*x + nl*(0.534667*x - 0.22*x*x + 0.067*x*x*x))
              + log(x)*(16.9477*x - 1.10133*nl*x + 2.78756*x*x - 0.0343333*x*x*x)
              + 8./9.*(3./2.*log(muf*muf/(mq[i].first*mq[i].first))+2.)*x/24.*(24.*x*x*x*log(x)*log(x) + 12.*x*log(x)+ 24.*x
              + Pi*Pi*(4.*x*x*x-9.*x*x-3.) - 6.*(4.*x*x*x + 3.*x*x + 1.)*(log(x)*log(x+1) + PolyLog(2,-x))
              -6.*(x-1.)*(4.*x*x + x + 1.)*(log(x)*log(1-x) + PolyLog(2,x))));
       }
     }

     return erg;
}

// Compute 2 Loop contribution to mMS2mOS at 3 loops for MS light quarks
double CRunDec::deltamMS2mOS(double mMS, std::pair<double,double>* mq, double asmu, double mu, int nlq, int nloops){
     double erg=0.0;
     double lmu = log(mu*mu/(mMS*mMS));
     if(!mq)
       return 0.0;

     for(int i = 0; i < 4; i++) {
       if(mq[i].first == 0.0)
        continue;

       double x = mq[i].first/mMS;
       double muf = mq[i].second;
       int nl = nlq-i;
       if(nloops >= 2) {
         erg += -asmu*asmu*2./3.*(48.*x*x*x*x*log(x)*log(x) + 48.*x*x*log(x)+72.*x*x + 8.*Pi*Pi*(x*x*x*x-3.*x*x*x-3.*x)
                -48.*(x+1.)*(x+1.)*(x*x - x + 1.)*(log(x)*log(x+1) + PolyLog(2,-x))
                -48.*(x-1.)*(x-1.)*(x*x + x + 1.)*(log(x)*log(1-x) + PolyLog(2,x)))/96.;
       }
       if(nloops >= 3) {
         erg += -asmu*asmu*asmu*(-21.8714*x - 4.348*x*x - 1.02211*x*x*x - 0.0493333*x*x*x*x + nl*(0.982667*x + 0.300333*x*x)
              + lmu*(-6.61056*x + 2.46511*x*x - 0.724333*x*x*x + nl*(0.534667*x - 0.22*x*x + 0.067*x*x*x))
              + log(x)*(16.9477*x - 1.10133*nl*x + 2.78756*x*x - 0.0343333*x*x*x)
              + 8./9.*(3./2.*log(muf*muf/(mq[i].first*mq[i].first))+2.)*x/24.*(24.*x*x*x*log(x)*log(x) + 12.*x*log(x)+ 24.*x
              + Pi*Pi*(4.*x*x*x-9.*x*x-3.) - 6.*(4.*x*x*x + 3.*x*x + 1.)*(log(x)*log(x+1) + PolyLog(2,-x))
              -6.*(x-1.)*(4.*x*x + x + 1.)*(log(x)*log(1-x) + PolyLog(2,x))))
              + asmu*asmu*asmu*8./9.*(48.*(3.*lmu + 7.)*x*x*x*x*log(x)*log(x) + 144.*x*x*log(x)+312.*x*x 
              + 8.*Pi*Pi*(7.*x*x*x*x-15.*x*x*x-3.*x) - lmu*(-72.*x*x +
              + 12.*Pi*Pi*(-2.*x*x*x*x+3.*x*x*x - 3*x))
              - 48.*((7.*x+5.)*x*x*x + x + 3.*lmu/2.*(2.*x*x*x*x+x*x*x-x-2.) -1.)*(log(x)*log(x+1) + PolyLog(2,-x))
              - 48.*(x-1.)*(3.*lmu/2.*(2.*x*x*x+x*x+x+2.) + 7.*x*x*x + 2.*x*x + 2.*x + 1.)*(log(x)*log(1-x) + PolyLog(2,x)))/192.;
       }
     }

     return erg;
}
       
// z_m^inv
double CRunDec::fZmInvM(double nl){
     double erg;
     erg=(8481925./93312. + 
       (137.*nl)/216. + (652841.*Pi*Pi)/38880. - (nl*Pi*Pi)/27. - 
       (695.*Pi*Pi*Pi*Pi)/7776. - (575.*Pi*Pi*log(2))/162. - 
       (22.*Pi*Pi*log(2)*log(2))/81. - 
       (55.*log(2)*log(2)*log(2)*log(2))/162. - (220.*A4)/27. - 
       nl*nl*(-2353./23328. - (13.*Pi*Pi)/324. - (7.*Zeta3)/54.) + 
       (58.*Zeta3)/27. - 
       (1439.*Pi*Pi*Zeta3)/432. - nl*(246643./23328. + (967.*Pi*Pi)/648. - 
       (61.*Pi*Pi*Pi*Pi)/1944. + (11.*Pi*Pi*log(2))/81. - 
       (2.*Pi*Pi*log(2)*log(2))/81. - 
       log(2)*log(2)*log(2)*log(2)/81. - (8.*A4)/27. + 
       (241.*Zeta3)/72.) + 
       (1975.*Zeta5)/216.);
     return erg;
}

// Function: double CRunDec::mMS2mOS(double mMS, std::pair<double,double>* mq, double asmu,
//                           double mu,int nl)
double CRunDec::mMS2mOS(double mMS, std::pair<double,double>* mq, double asmu, double mu, int nl, double fdelm){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl == 4 && (Nf<4||Nf>6)){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR 4 LOOPS and "<< nl << " FLAVORS"<<endl;
       RETURN  
     }
     double deltalight = deltamMS2mOS(mMS, mq, asmu/Pi, mu, Nf-1, nl);
     double sum[5];
     sum[0]=(double) 1.;
     sum[1]=asmu*(this ->fOsFromMs1(mu, mMS))/Pi;
     sum[2]=asmu*asmu*((this-> fOsFromMs2(mu, mMS, Nf-1)))/(Pi*Pi);
     sum[3]=asmu*asmu*asmu*(this-> fOsFromMs3(mu, mMS, Nf-1)+
           this->fZmInvM(Nf-1))/(Pi*Pi*Pi);   
     sum[4]=asmu*asmu*asmu*asmu*(this->fOsFromMs4(mu, mMS, Nf-1, fdelm))/(Pi*Pi*Pi*Pi);  
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{ 
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     erg += deltalight;
     return mMS*erg;
}

// Function: double CRunDec::mMS2mOSmod(double mMS, std::pair<double,double>* mq, double asmu,
//                           double mu, int nf, int nloop, double fdelm), needed for several other mass relations
double CRunDec::mMS2mOSmod(double mMS, std::pair<double,double>* mq, double asmu, double mu,int nf, int nloop, double fdelm){
     if(nloop<0||nloop>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloop <<" LOOPS"<<endl;
       RETURN  
     }
     if(nloop == 4 && (nf<4||nf>6)){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR 4 LOOPS and "<< nf << " FLAVORS"<<endl;
       RETURN  
     }
     SetConstants(nf);
     double sum[5];
     double lmm = log((mu*mu)/(mMS*mMS));
     sum[0]=(double) 1.;


     sum[1]=asmu*(this ->fOsFromMs1(mu, mMS))/(Pi);
     sum[2]=asmu*asmu*((this-> fOsFromMs2(mu, mMS, Nf-1)))/(Pi*Pi) +
            asmu*asmu*(lmm)/(6.*Pi)*(this ->fOsFromMs1(mu, mMS))/(Pi);
     sum[3]=asmu*asmu*asmu*(this-> fOsFromMs3(mu, mMS, Nf-1)+this->fZmInvM(Nf-1))/(Pi*Pi*Pi)+
            asmu*asmu*asmu*2.*(lmm)/(6.*Pi)*((this-> fOsFromMs2(mu, mMS, Nf-1)))/(Pi*Pi)+
            asmu*asmu*asmu*(-11./72. + (11.*lmm)/24. + lmm*lmm/36.)/(Pi*Pi)*(this ->fOsFromMs1(mu, mMS))/(Pi);
     sum[4]=asmu*asmu*asmu*asmu*(this->fOsFromMs4(mu, mMS, Nf-1, fdelm))/(Pi*Pi*Pi*Pi) + 
			asmu*asmu*asmu*asmu*3.*(lmm)/(6.*Pi)*(this-> fOsFromMs3(mu, mMS, Nf-1)+this->fZmInvM(Nf-1))/(Pi*Pi*Pi) + 
            asmu*asmu*asmu*asmu*2.*(-11./72. + (11.*lmm)/24. + lmm*lmm/36.)/(Pi*Pi)*
                                  ((this-> fOsFromMs2(mu, mMS, Nf-1)))/(Pi*Pi) +
            asmu*asmu*asmu*asmu*(-564731./124416. + (2645.*lmm)/1728. + (167.*lmm*lmm)/576. + lmm*lmm*lmm/216. + 
   					(2633./31104. - (67.*lmm)/576. + lmm*lmm/36.)*(Nf-1) + (82043.*Zeta3)/27648.)/(Pi*Pi*Pi)
                    *(this ->fOsFromMs1(mu, mMS))/(Pi) + 
            asmu*asmu*asmu*asmu*(lmm)/(6.*Pi)*(lmm)/(6.*Pi)*((this-> fOsFromMs2(mu, mMS, Nf-1)))/(Pi*Pi);  
  
     double erg=0.0;
     if(nloop==0){
       erg=1;
     }
     else{ 
       for(int i=0; i<=nloop; i++){
         erg+=sum[i];
       }
     }
     return mMS*erg;
}      

// Coefficients of eq.(16) of [RunDec]
double CRunDec::fMumFromOs1(void){
     return (double) -cf;
}
  
double CRunDec::fMumFromOs2(void){
     double erg;
     erg= (double) ((-1111.*ca*cf)/384. + (199.*cf*cf)/128. + (143.*cf*tr)/96. +
     (71.*cf*(Nf-1)*tr)/96. + (ca*cf*Zeta2)/2. - (15.*cf*cf*Zeta2)/8. - 
     (3.*ca*cf*log(2)*Zeta2)/2. +
     3.*cf*cf*log(2)*Zeta2 - cf*tr*Zeta2 + (cf*(Nf-1)*tr*Zeta2)/2. + 
     (3.*ca*cf*Zeta3)/8. -
     (3.*cf*cf*Zeta3)/4.);     
     return erg;      
}

// z_m^SI
double CRunDec::fMumFromOs3(void){
     double erg;
     erg= (double) -7172965./93312. - 
	 (293.*(Nf-1))/216. - (618281.*Pi*Pi)/38880. - ((Nf-1)*Pi*Pi)/9. + 
     (695.*Pi*Pi*Pi*Pi)/7776. + (623.*Pi*Pi*log(2))/162. + 
     (22.*Pi*Pi*log(2)*log(2))/81. + 
     (55.*log(2)*log(2)*log(2)*log(2))/162. + (220.*A4)/27. + 
     (Nf-1)*(Nf-1)*(-2353./23328. - (13.*Pi*Pi)/324. - (7.*Zeta3)/54.) - 
     (70.*Zeta3)/27. + 
     (1439.*Pi*Pi*Zeta3)/432. + (Nf-1)*(246643./23328. + (967.*Pi*Pi)/648. - 
     (61.*Pi*Pi*Pi*Pi)/1944. + (11.*Pi*Pi*log(2))/81. - 
     (2*Pi*Pi*log(2)*log(2))/81. - 
     log(2)*log(2)*log(2)*log(2)/81. - (8.*A4)/27. + (241.*Zeta3)/72.) - 
     (1975.*Zeta5)/216.;
     return erg;   
}

double CRunDec::fMumFromOs4(double err){
     double erg;
     erg = err*(-3214.227044839041 + 692.4809215366435*(-1. + Nf) - 41.95978562498058*(-1. + Nf)*(-1. + Nf)
                + 0.678141025604516*(-1. + Nf)*(-1. + Nf)*(-1. + Nf));
     return erg;
}

// Function: double CRunDec::mOS2mSI(double mOS, std::pair<double,double>* mq, double asM, 
//                           int nl)
double CRunDec::mOS2mSI(double mOS, std::pair<double,double>* mq, double asM, int nl, double fdelm){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl == 4 && (Nf<4||Nf>6)){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR 4 LOOPS and "<< Nf << " FLAVORS"<<endl;
       RETURN  
     }

     double sum[5];
     double deltalight = deltamOS2mMS(mOS, mq, asM/Pi, mOS, Nf-1, nl);
     sum[0]=(double) 1.;
     sum[1]=asM*(this ->fMumFromOs1())/Pi;
     sum[2]=asM*asM*(this-> fMumFromOs2())/(Pi*Pi);
     sum[3]=asM*asM*asM*(this-> fMumFromOs3())/(Pi*Pi*Pi);     
     sum[4]=asM*asM*asM*asM*(this->fMumFromOs4(fdelm))/(Pi*Pi*Pi*Pi);
     

     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     erg += deltalight;
     return mOS*erg;       
       
}

// Function: double CRunDec::mOS2mMSrun(double mOS, std::pair<double,double>* mq, double asmu, 
//                           double mu, int nl)
double CRunDec::mOS2mMSrun(double mOS, std::pair<double,double>* mq, double asmu, double mu,
			   int nl){
     double asM=0.0;
     asM= this-> AlphasExact(asmu, mu, mOS, nl);
     double mum= this-> mOS2mSI(mOS, mq, asM, nl);
     double asmum= this-> AlphasExact(asmu, mu, mum, nl);
     double newM= this->mMS2mMS(mum, asmum, asmu, nl);
     return newM;       
}

// Function: double CRunDec::mMS2mOSrun(double mMS, std::pair<double,double>* mq, double asmu, 
//                           double mu, int nl)
double CRunDec::mMS2mOSrun(double mMS, std::pair<double,double>* mq, double asmu, double mu,
			   int nl){
     double mNeu = mMS2mSI(mMS, asmu, mu, nl);
     double asmNeu = AlphasExact(asmu, mu, mNeu, nl);
     return mMS2mOS(mNeu, mq, asmNeu, mNeu, nl);   
}

// Coefficients of eq.(19) of [RunDec]
double CRunDec::fRiFromMs(double alpha, double nl){
     double sum[4];
     sum[0]= 1.;
     sum[1]= (4.*alpha)/3.;
     sum[2]= alpha*alpha*((1123./72. - (89.*Nf)/144. - (19.*Zeta3)/6.));
     sum[3]= alpha*alpha*alpha*(6663911./41472. - (118325.*Nf)/7776. +
            (4459.*Nf*Nf)/23328. + 
            (4.*(1123./72. - (89.*Nf)/144. - (19.*Zeta3)/6.))/3. - 
            (408007.*Zeta3)/6912. + 
            (617.*Nf*Zeta3)/216. + (Nf*Nf*Zeta3)/54. - (4.*(-995./72. + 
            (89.*Nf)/144. + (19.*Zeta3)/6.))/3. - 
            (5.*Nf*Zeta4)/12. + (185.*Zeta5)/36.);
     double erg=0.0;
     if(nl==0){
       erg=1.;
     }
     else{
       erg=0.;
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     return (erg);        
}

// Function: double CRunDec::mMS2mRI(double mMS, double asmu, int nl)
double CRunDec::mMS2mRI(double mMS, double asmu, int nl){
     if(nl<0||nl>3){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     return (double) mMS*(this->fRiFromMs((asmu/Pi), nl)); 
}

// Coefficients needed for the transformation of mOS to mMSit
double CRunDec::fHelpmOS2mMSit(double mMS, double mOS, std::pair<double,double>* mq, double asmu,
                                double mu, int nl){
     double sum[4];
     double deltalight = deltamMS2mOS(mMS,mq,asmu/Pi,mu,Nf-1,nl);
     sum[0]=(double) 1.;
     sum[1]=asmu*(this ->fOsFromMs1(mu, mMS))/Pi;
     sum[2]=asmu*asmu*((this-> fOsFromMs2(mu, mMS, Nf-1)))/(Pi*Pi);
     sum[3]=asmu*asmu*asmu*(this-> fOsFromMs3(mu, mMS,Nf-1)+
           +this->fZmInvM(Nf-1))/(Pi*Pi*Pi);
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{ 
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     erg += deltalight;
     return erg;
}

// Function: double CRunDec::mOS2mMSit(double mOS, std::pair<double,double>* mq, double asmu, 
//                          double mu, int nl)
double CRunDec::mOS2mMSit(double mOS, std::pair<double,double>* mq, double asmu, double mu,
                           int nl){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN
     }
     double epsilonX= 1e-8;
     double x0=mOS-0.1*mOS; 
     double x1=mOS+0.1*mOS;
     double f0= (x0*(this->fHelpmOS2mMSit(x0,mOS, mq, asmu, mu, nl))-mOS);
     double f1= (x1*(this->fHelpmOS2mMSit(x1,mOS, mq, asmu, mu, nl))-mOS);
     if(f0*f1>0){
       cout<<"WARNING: No root can be calculatet!"<<endl;
       RETURN
     }
     double xTest;
     double fTest;
     do{
       xTest= (x0+x1)/2.;
       fTest= (xTest*(this->fHelpmOS2mMSit(xTest,mOS, mq, asmu, mu, nl))-mOS);
       if(f0*fTest<=0){x1= xTest;}
       else {x0= xTest;}
     }
     while(abs(x1-x0)>= epsilonX);
       double mNeu=xTest;
     return mNeu;               
}

double CRunDec::PSdelta(double asmu, double muf, double mu, int nl, int nloops) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     double lf = log(mu*mu/(muf*muf));
     double api = asmu;
     double ret[5];
     ret[0] = 0.0;
     ret[1] = api*4./3.;
     ret[2] = api*api*(97./9. + nl*(-22./27. - (2.*lf)/9.) + (11.*lf)/3.);
     ret[3] = api*api*api*(33623./216. + 3.*Pi*Pi - (3.*Pi*Pi*Pi*Pi)/16. + (610.*lf)/9.
              + (121.*lf*lf)/12. + nl*nl*(157./243. + (22.*lf)/81. + lf*lf/27.) +
              nl*(-7145./324. - (493.*lf)/54. - (11.*lf*lf)/9. - (13.*Zeta3)/9.) + (11.*Zeta3)/2.);
     ret[4] = api*api*api*api*(3567.723056629293 + (7271.*lf*lf)/24. + (1331.*lf*lf*lf)/48.
              + nl*nl*nl*(-2951./4374. - (157.*lf)/486. - (11.*lf*lf)/162. - 
              lf*lf*lf/162.) + nl*(-701.2303148875468 - (8485.*lf*lf)/144. - (121.*lf*lf*lf)/24.
              + lf*(-253189./864. - (3.*Pi*Pi)/2. + (3.*Pi*Pi*Pi*Pi)/32.
              - (44.*Zeta3)/3.)) + nl*nl*(1751971./46656. + Pi*Pi*Pi*Pi/135. + (773.*lf*lf)/216. +
              (11.*lf*lf*lf)/36. + lf*(15355./864. + (13.*Zeta3)/18.) + 
              (259.*Zeta3)/108.) + lf*(26125./18. + (99.*Pi*Pi)/4. - (99.*Pi*Pi*Pi*Pi)/64. + (363.*Zeta3)/8.));

     double delta = 0.0;
     for(int i = 0; i <= nloops; i++) {
       delta += ret[i];
     }
	 return delta;
}

// Function: double CRunDec::mOS2mPS(double mOS, std::pair<double,double>* mq, double asmu, double mu, double muf, int nl, int nloops)
double CRunDec::mOS2mPS(double mOS, std::pair<double,double>* mq, double asmu, double mu, double muf, int nl, int nloops) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }

	 return mOS - muf*PSdelta(asmu/Pi, muf, mu, nl, nloops);
}

// Function: double CRunDec::mMS2mPS(double mMS, std::pair<double,double>* mq, double asmu,
//                                   double mu, double muf, int nl, int nloops, double fdelm)
double CRunDec::mMS2mPS(double mMS, std::pair<double,double>* mq, double asmu,
                        double mu, double muf, int nl, int nloops, double fdelm) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }

     double delmuf = PSdelta(asmu/Pi, muf, mu, nl, nloops); 
     double exmOS = mMS2mOSmod(mMS, this->mq, asmu, mu, nl+1, nloops, fdelm);
     return (exmOS - muf*delmuf);
}

// Function: double CRunDec::mPS2mMS(double mPS, std::pair<double,double>* mq, double asmu,
//                                   double mu, double muf, int nl, int nloops, double fdelm)
double CRunDec::mPS2mMS(double mPS, std::pair<double,double>* mq, double asmu,
                        double mu, double muf, int nl, int nloops, double fdelm) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }

     double lowbound = mPS - mPS/4.;
     double highbound = mPS + mPS/4.;
     bool found = false;
     double f1 = mMS2mPS(lowbound, mq, asmu, mu, muf, nl, nloops, fdelm) - mPS;
     double f2 = mMS2mPS(highbound, mq, asmu, mu, muf, nl, nloops, fdelm) - mPS;
     for(int j = 0; j < 1000; j++) {
       if(f1*f2 < 0.0) {
         found = true;
         break;
       }
       if(fabs(f1) < fabs(f2)) {
         lowbound  += 1.5*(lowbound - highbound);
         if(lowbound < 0.0)
            lowbound = 0.0;
         f1 = mMS2mPS(lowbound, mq, asmu, mu, muf, nl, nloops, fdelm) - mPS;
       } else {
         highbound  -= 1.5*(lowbound - highbound);
         f2 = mMS2mPS(highbound, mq, asmu, mu, muf, nl, nloops, fdelm) - mPS;
       }
     }
     if(found) {
     double acc = 10e-10;
     double dx;
     double xmid;
     double mMS = f1 < 0.0 ? (dx=highbound-lowbound,lowbound) : (dx=lowbound-highbound,highbound);
     for(int j = 0; j < 1000; j++) {
       f2 = mMS2mPS(xmid = mMS+(dx *= 0.5), mq, asmu, mu, muf, nl, nloops, fdelm) - mPS;
       if(f2 <= 0.0)
         mMS = xmid;
       if(fabs(dx) < acc || f2 == 0.0)
         return mMS;
     }
     }
     return 0.0;
}

// Function: double CRunDec::mPS2mSI(double mPS, std::pair<double,double>* mq,
//                                   double (*as)(double), double muf, int nl, int nloops, double fdelm)
// The function pointer passed should contain the adress of a function computing alpha_s in dependence of mu
double CRunDec::mPS2mSI(double mPS, std::pair<double,double>* mq,
                        double (*as)(double), double muf, int nl, int nloops, double fdelm) {
    if(as == NULL) {
      cout << "Pointer to as == NULL! Aborting..." << endl;
      RETURN
    }
	double mMS1 = 0;
	double mMS = mPS;
	double acc = 10e-6; 
	while(fabs(mMS1 - mMS) > acc) {
		mMS1 = mMS; 
        mMS = mPS2mMS(mPS, mq, as(mMS1), mMS1, muf, nl, nloops, fdelm);
	}
	return mMS;
}


double CRunDec::E1p(double mOS, double asmu, double mu, int nl, int nloops) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     double lmm = log(mu/mOS);
     double lmm34 = log(3.*mu/(4.*asmu*mOS));
     double api = asmu/(Pi);
     double ret[5];

     ret[0] = 0.0;
     ret[1] = 1.0;
     ret[2] = api*(97./6. + nl*(-11./9. - (2.*lmm34)/
                 3.) + 11.*lmm34);
     ret[3] = api*api*(1793./12. + (2917./216. - (11.*nl)/18. + nl*nl/54.)*Pi*Pi - 
              (9.*Pi*Pi*Pi*Pi)/32. + (927.*lmm34)/4. + 
              (363.*lmm34*lmm34)/4. + nl*(-1693./72. - 
                (193.*lmm34)/6. - 11.*lmm34*lmm34 - (19.*Zeta3)/2.) + nl*nl*(77./108. + 
                lmm34 + lmm34*lmm34/3. + 
                (2.*Zeta3)/9.) + (275.*Zeta3)/4.);
     ret[4] = api*api*api*(1267919./1728. + Pi*Pi*Pi*Pi*(-723119./51840. + (11.*nl*nl)/
                 1080. - nl*nl*nl/4860. + nl*(59677./77760. + (3.*lmm34)/8.) - (99.*lmm34)/16.) + 
              (4521.*lmm34*lmm34)/2. + (1331.*
                lmm34*lmm34*lmm34)/2. + (114917.*Zeta3)/48. + 
              lmm34*(247675./96. + (3025.*Zeta3)/2.) + 
              Pi*Pi*(265.389067842508 + (865.*lmm)/18. + 
                (26897.*lmm34)/108. + nl*nl*(905./432. + 
                  (11.*lmm34)/9. - (11.*Zeta3)/6.) + 
                nl*nl*nl*(-19./486. - (2.*lmm34)/81. + 
                  Zeta3/27.) + nl*(-397591./7776. - (5095.*lmm34)/162. + (121.*Zeta3)/4.)) + (13432.614375 - 
                3289.906669391583*nl - (1000.*nl*nl*nl)/729. + nl*nl*
                 ((14002./81. - (416.*Zeta3)/3.)/3. + (3.*(12541./243. + (64.*Pi*Pi*Pi*Pi)/
                      135. + (368.*Zeta3)/3.))/4.))/32. + nl*(-52033./288. - 
                (10955.*lmm34*lmm34)/24. - 
                121.*lmm34*lmm34*lmm34 + lmm34*
                 (-166309./288. - (902.*Zeta3)/3.) - (8797.*Zeta3)/18. - 
                363.*Zeta5) + nl*nl*nl*(-98./729. - (5.*lmm34*lmm34)/
                 9. - (4.*lmm34*lmm34*lmm34)/27. + 
                lmm34*(-50./81. - (8.*Zeta3)/27.) - 
                (44.*Zeta3)/81. - (4.*Zeta5)/9.) + (3993.*Zeta5)/2. + 
              nl*nl*(3073./288. + (1027.*lmm34*lmm34)/36. + 
                (22.*lmm34*lmm34*lmm34)/3. + (3239.*Zeta3)/108. + 
                lmm34*(10351./288. + (158.*Zeta3)/9.) + 
                22.*Zeta5));


     double E = 0.0;
     for(int i = 0; i <= nloops; i++) {
       E += ret[i];
     }
	 return -E*(4.*asmu*asmu*mOS)/9.;
}

// Function: double CRunDec::mOS2m1S(double mOS, std::pair<double,double>* mq,
//                                   double asmu, double mu, int nl, int nloops)
double CRunDec::mOS2m1S(double mOS, std::pair<double,double>* mq, double asmu, double mu, int nl, int nloops) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
	 return mOS + 0.5*E1p(mOS, asmu, mu, nl, nloops);
}

// Function: double CRunDec::mMS2m1S(double mMS, std::pair<double,double>* mq,
//                                   double asmu, double mu, int nl, int nloops, double fdelm)
double CRunDec::mMS2m1S(double mMS, std::pair<double,double>* mq,
                        double asmu, double mu, int nl, int nloops, double fdelm) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     if(nl<3||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR 4 LOOPS and " << nl << " FLAVORS"<<endl;
       RETURN  
     }
     double lmm = log(mu*mu/(mMS*mMS));
     double log34 = log((3.*mu)/(4.*asmu*mMS));
     double m1S = mMS;
     switch(nloops) {
       case 1: m1S += asmu*mMS*(12. - 2.*asmu*Pi + 9.*lmm)/(9.*Pi);
               break;
       case 2: m1S += asmu*mMS*(-(2.*asmu)/9. + asmu*asmu*(-291. + 22.*nl - 198.*log34 + 
                      12.*nl*log34 - 3.*(8. + 6.*lmm))/(81.*Pi) + (4. + 3.*lmm)/(3.*Pi) +
                      (asmu*(2763. - 142.*nl + 96.*Pi*Pi - 16.*nl*Pi*Pi + 32.*Pi*Pi*ln2 + 
                      2036.*lmm - 104.*nl*lmm + 564.*lmm*lmm - 24.*nl*lmm*lmm - 48.*Zeta3))/(288.*Pi*Pi));
               break;
       case 3: m1S += asmu*mMS*(-(2.*asmu)/9. + asmu*asmu*(-291. + 22.*nl -  198.*log34 + 
                      12.*nl*log34)/(81.*Pi) + (4. + 3.*lmm)/(3.*Pi) - (2.*asmu*asmu*(4. + 3.*lmm))/(27.*Pi) + 
                      asmu*asmu*asmu*(-372. + 40.*nl - 792.*log34 + 
                      48.*nl*log34 - 279.*lmm + 30.*nl*lmm - 
                      594.*log34*lmm + 36.*nl*log34*lmm)/(243.*Pi*Pi)
                      + (asmu*(2763. - 142.*nl + 96.*Pi*Pi - 16.*nl*Pi*Pi + 32.*Pi*Pi*ln2 + 2036.*lmm - 104.*nl*lmm + 
                      564.*lmm*lmm - 24.*nl*lmm*lmm - 48.*Zeta3))/(288.*Pi*Pi) + asmu*asmu*asmu*(-2763.+ 142.*nl - 
                      96.*Pi*Pi + 16.*nl*Pi*Pi - 32.*Pi*Pi*ln2 -
                      2036.*lmm + 104.*nl*lmm - 564.*lmm*lmm + 
                      24.*nl*lmm*lmm + 48.*Zeta3)/(1296.*Pi*Pi)
                      + asmu*asmu*asmu*(-129096. + 20316.*nl - 616.*nl*nl - 
                      11668.*Pi*Pi + 528.*nl*Pi*Pi - 16.*nl*nl*Pi*Pi + 243.*Pi*Pi*Pi*Pi - 
                      200232.*log34 + 27792.*nl*log34 - 864.*nl*nl*log34 - 
                      78408.*log34*log34 + 9504.*nl*log34*log34 - 
                      288.*nl*nl*log34*log34 - 59400.*Zeta3 + 8208.*nl*Zeta3 - 
                      192.*nl*nl*Zeta3)/(3888*Pi*Pi) + 
                      (asmu*asmu*(42314585. - 4636940.*nl + 47060.*nl*nl + 7834092.*Pi*Pi - 
                      713520.*nl*Pi*Pi + 18720.*nl*nl*Pi*Pi - 41700.*Pi*Pi*Pi*Pi + 14640.*nl*Pi*Pi*Pi*Pi - 
                      1656000.*Pi*Pi*ln2 - 63360.*nl*Pi*Pi*ln2 - 126720.*Pi*Pi*ln2*ln2 + 
                      11520.*nl*Pi*Pi*ln2*ln2 - 158400.*ln2*ln2*ln2*ln2 + 5760.*nl*ln2*ln2*ln2*ln2 + 
                      33620760.*lmm - 3723120.*nl*lmm + 64080.*nl*nl*lmm + 1010880.*Pi*Pi*lmm - 220320.*nl*Pi*Pi*lmm + 8640.*nl*nl*Pi*Pi*lmm + 
                      336960.*Pi*Pi*ln2*lmm - 17280.*nl*Pi*Pi*ln2*lmm + 11726100.*lmm*lmm - 1247400.*nl*lmm*lmm + 28080.*nl*nl*lmm*lmm + 
                      2009880.*lmm*lmm*lmm - 185760.*nl*lmm*lmm*lmm + 4320.*nl*nl*lmm*lmm*lmm - 3801600.*A4 + 138240.*nl*A4 + 1002240.*Zeta3 - 
                      1561680.*nl*Zeta3 + 60480.*nl*nl*Zeta3 - 1554120.*Pi*Pi*Zeta3 - 894240.*lmm*Zeta3 - 362880.*nl*lmm*Zeta3 + 4266000.*Zeta5))/
                      (466560.*Pi*Pi*Pi));
                      break;
       case 4: m1S += - (2*asmu*asmu*mMS)/9 + (-291*asmu*asmu*asmu*mMS + 22*asmu*asmu*asmu*mMS*nl - 
                      198*asmu*asmu*asmu*mMS*log34 + 
                      12*asmu*asmu*asmu*mMS*nl*log34)/(81*Pi) + 
                      (asmu*mMS*(4 + 3*lmm))/(3*Pi) - 
                      (2*(4*asmu*asmu*asmu*mMS + 3*asmu*asmu*asmu*mMS*lmm))/(27*Pi) + 
                      (-372*asmu*asmu*asmu*asmu*mMS + 40*asmu*asmu*asmu*asmu*mMS*nl - 
                      792*asmu*asmu*asmu*asmu*mMS*log34 + 
                      48*asmu*asmu*asmu*asmu*mMS*nl*log34 - 
                      279*asmu*asmu*asmu*asmu*mMS*lmm + 30*asmu*asmu*asmu*asmu*mMS*nl*lmm - 
                      594*asmu*asmu*asmu*asmu*mMS*log34*lmm + 
                      36*asmu*asmu*asmu*asmu*mMS*nl*log34*lmm)/(243*Pi*Pi)
                      + 1.5139171475883256E-8*(-3.6643425048120157E6*asmu*asmu*asmu*asmu*mMS - 
                      3.670703784417217E9*asmu*asmu*asmu*asmu*asmu*mMS + 292307.7010004199*asmu*asmu*asmu*asmu*mMS*nl + 
                      6.48095550442438E8*asmu*asmu*asmu*asmu*asmu*mMS*nl - 
                      3.6315792483681194E7*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl + 
                      595364.9557740748*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*nl - 
                      2.245324461595123E8*asmu*asmu*asmu*asmu*asmu*mMS*log(mu/mMS) - 
                      3.062101652047988E9*asmu*asmu*asmu*asmu*asmu*mMS*log34 + 
                      5.902121099662828E8*asmu*asmu*asmu*asmu*asmu*mMS*nl*log34 - 
                      3.325469404123189E7*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*log34 + 
                      576205.4055744545*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*nl*log34 - 
                      1.1274215864457371E9*asmu*asmu*asmu*asmu*asmu*mMS*log34*log34 + 
                      2.230344299450968E8*asmu*asmu*asmu*asmu*asmu*mMS*nl*log34*log34 - 
                      1.3715687063042711E7*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*log34*log34 + 
                      263004.5457918066*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*nl*log34*log34 - 
                      3.1505314540400505E8*asmu*asmu*asmu*asmu*asmu*mMS*log34*log34*log34 + 
                      5.728239007345548E7*asmu*asmu*asmu*asmu*asmu*mMS*nl*log34*log34*log34 - 
                      3.4716600044518467E6*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*log34*log34*log34 + 
                      70134.54554448176*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*nl*log34*log34*log34 + 
                      7.228435160639514E7*asmu*asmu*asmu*asmu*mMS*lmm - 
                      1.061086986213349E8*asmu*asmu*asmu*asmu*asmu*mMS*lmm - 
                      9.657775840868242E6*asmu*asmu*asmu*asmu*mMS*nl*lmm + 
                      1.3033590170779591E7*asmu*asmu*asmu*asmu*asmu*mMS*nl*lmm + 
                      221297.1713327142*asmu*asmu*asmu*asmu*mMS*nl*nl*lmm - 
                      323617.3768070282*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*lmm - 
                      6.0602822464077026E7*asmu*asmu*asmu*asmu*asmu*mMS*log34*
                      lmm + 8.924620920535302E6*asmu*asmu*asmu*asmu*asmu*mMS*nl*log34*lmm - 
                      271771.36398486677*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*log34*
                      lmm - 4.2961792555091605E7*asmu*asmu*asmu*asmu*asmu*mMS*log34*log34*lmm + 
                      5.20749000667777E6*asmu*asmu*asmu*asmu*asmu*mMS*nl*log34*log34*
                      lmm - 157802.72747508393*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*log34*log34*lmm + 
                      3.6846948160874695E7*asmu*asmu*asmu*asmu*mMS*lmm*lmm - 
                      1.4084441353288308E7*asmu*asmu*asmu*asmu*asmu*mMS*lmm*lmm - 
                      4.890850849173104E6*asmu*asmu*asmu*asmu*mMS*nl*lmm*lmm + 
                      1.826785740978923E6*asmu*asmu*asmu*asmu*asmu*mMS*nl*lmm*lmm + 
                      108536.53795824476*asmu*asmu*asmu*asmu*mMS*nl*nl*lmm*lmm - 
                      50409.20461009625*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*lmm*lmm - 
                      1.0198001263077298E7*asmu*asmu*asmu*asmu*asmu*mMS*log34*
                      lmm*lmm + 1.0520181831672261E6*asmu*asmu*asmu*asmu*asmu*mMS*nl*log34*lmm*lmm - 
                      26300.454579180656*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*log34*
                      lmm*lmm + 9.281070040868293E6*asmu*asmu*asmu*asmu*mMS*
                      lmm*lmm*lmm - 2.039381082160633E6*asmu*asmu*asmu*asmu*asmu*mMS*
                      lmm*lmm*lmm - 919055.1091884745*asmu*asmu*asmu*asmu*mMS*nl*
                      lmm*lmm*lmm + 188486.5911507947*asmu*asmu*asmu*asmu*asmu*mMS*nl*
                      lmm*lmm*lmm + 20406.00584022232*asmu*asmu*asmu*asmu*mMS*nl*nl*
                      lmm*lmm*lmm - 4383.409096530109*asmu*asmu*asmu*asmu*asmu*mMS*nl*nl*
                      lmm*lmm*lmm + 1.3530751564824337E6*asmu*asmu*asmu*asmu*mMS*
                      lmm*lmm*lmm*lmm - 130284.49882603479*asmu*asmu*asmu*asmu*mMS*nl*
                      lmm*lmm*lmm*lmm + 3139.3855138803565*asmu*asmu*asmu*asmu*mMS*nl*nl*
                      lmm*lmm*lmm*lmm + 678107.270998157*asmu*asmu*asmu*asmu*mMS*
                      this->fOsFromMs4(mu,mMS,nl,fdelm)) + 
                      (asmu*asmu*mMS*(2763 - 142*nl + 96*Pi*Pi - 16*nl*Pi*Pi + 32*Pi*Pi*ln2 + 
                      2036*lmm - 104*nl*lmm + 564*lmm*lmm - 24*nl*lmm*lmm - 48*Zeta3))/
                      (288*Pi*Pi) + (-2763*asmu*asmu*asmu*asmu*mMS + 142*asmu*asmu*asmu*asmu*mMS*nl - 
                      96*asmu*asmu*asmu*asmu*mMS*Pi*Pi + 16*asmu*asmu*asmu*asmu*mMS*nl*Pi*Pi - 
                      32*asmu*asmu*asmu*asmu*mMS*Pi*Pi*ln2 - 2036*asmu*asmu*asmu*asmu*mMS*lmm + 
                      104*asmu*asmu*asmu*asmu*mMS*nl*lmm - 564*asmu*asmu*asmu*asmu*mMS*lmm*lmm + 
                      24*asmu*asmu*asmu*asmu*mMS*nl*lmm*lmm + 48*asmu*asmu*asmu*asmu*mMS*Zeta3)/(1296*Pi*Pi)
                      + (-129096*asmu*asmu*asmu*asmu*mMS + 20316*asmu*asmu*asmu*asmu*mMS*nl - 616*asmu*asmu*asmu*asmu*mMS*nl*nl - 
                      11668*asmu*asmu*asmu*asmu*mMS*Pi*Pi + 528*asmu*asmu*asmu*asmu*mMS*nl*Pi*Pi - 
                      16*asmu*asmu*asmu*asmu*mMS*nl*nl*Pi*Pi + 243*asmu*asmu*asmu*asmu*mMS*Pi*Pi*Pi*Pi - 
                      200232*asmu*asmu*asmu*asmu*mMS*log34 + 
                      27792*asmu*asmu*asmu*asmu*mMS*nl*log34 - 
                      864*asmu*asmu*asmu*asmu*mMS*nl*nl*log34 - 
                      78408*asmu*asmu*asmu*asmu*mMS*log34*log34 + 
                      9504*asmu*asmu*asmu*asmu*mMS*nl*log34*log34 - 
                      288*asmu*asmu*asmu*asmu*mMS*nl*nl*log34*log34 - 
                      59400*asmu*asmu*asmu*asmu*mMS*Zeta3 + 8208*asmu*asmu*asmu*asmu*mMS*nl*Zeta3 - 
                      192*asmu*asmu*asmu*asmu*mMS*nl*nl*Zeta3)/(3888*Pi*Pi) + 
                      (asmu*asmu*asmu*mMS*(42314585 - 4636940*nl + 47060*nl*nl + 7834092*Pi*Pi - 
                      713520*nl*Pi*Pi + 18720*nl*nl*Pi*Pi - 41700*Pi*Pi*Pi*Pi + 14640*nl*Pi*Pi*Pi*Pi - 
                      1656000*Pi*Pi*ln2 - 63360*nl*Pi*Pi*ln2 - 126720*Pi*Pi*ln2*ln2 + 
                      11520*nl*Pi*Pi*ln2*ln2 - 158400*ln2*ln2*ln2*ln2 + 5760*nl*ln2*ln2*ln2*ln2 + 
                      33620760*lmm - 3723120*nl*lmm + 
                      64080*nl*nl*lmm + 1010880*Pi*Pi*lmm - 
                      220320*nl*Pi*Pi*lmm + 8640*nl*nl*Pi*Pi*lmm + 
                      336960*Pi*Pi*ln2*lmm - 
                      17280*nl*Pi*Pi*ln2*lmm + 11726100*lmm*lmm - 
                      1247400*nl*lmm*lmm + 28080*nl*nl*lmm*lmm + 
                      2009880*lmm*lmm*lmm - 185760*nl*lmm*lmm*lmm + 
                      4320*nl*nl*lmm*lmm*lmm - 3801600*A4 + 
                      138240*nl*A4 + 1002240*Zeta3 - 1561680*nl*Zeta3 + 
                      60480*nl*nl*Zeta3 - 1554120*Pi*Pi*Zeta3 - 
                      894240*lmm*Zeta3 - 362880*nl*lmm*Zeta3 + 
                      4266000*Zeta5))/(466560*Pi*Pi*Pi);
                      break;
     }



     return m1S;  
}

// Function: double CRunDec::m1S2mMS(double m1S, std::pair<double,double>* mq,
//                                   double asmu, double mu, int nl, int nloops, double fdelm)
double CRunDec::m1S2mMS(double m1S, std::pair<double,double>* mq,
                        double asmu, double mu, int nl, int nloops, double fdelm) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }

     double lowbound = m1S - m1S/5.;
     double highbound = m1S + m1S/5.;
     bool found = false;
     double f1 = mMS2m1S(lowbound, mq, asmu, mu, nl, nloops, fdelm) - m1S;
     double f2 = mMS2m1S(highbound,mq, asmu, mu, nl, nloops, fdelm) - m1S;
     for(int j = 0; j < 1000; j++) {
       if(f1*f2 < 0.0) {
         found = true;
         break;
       }
       if(fabs(f1) < fabs(f2))
         f1 = mMS2m1S(lowbound += 1.5*(lowbound - highbound), mq, asmu, mu, nl, nloops, fdelm) - m1S;
       else
         f2 = mMS2m1S(highbound -= 1.5*(lowbound - highbound), mq, asmu, mu, nl, nloops, fdelm) - m1S;
     }
     if(found) {
     double acc = 10e-10;
     double dx;
     double xmid;
     double mMS = f1 < 0.0 ? (dx=highbound-lowbound,lowbound) : (dx=lowbound-highbound,highbound);
     for(int j = 0; j < 1000; j++) {
       f2 = mMS2m1S(xmid = mMS+(dx *= 0.5), mq, asmu, mu, nl, nloops, fdelm) - m1S;
       if(f2 <= 0.0)
         mMS = xmid;
       if(fabs(dx) < acc || f2 == 0.0)
         return mMS;
     }
     }
     return 0.0;
}

// Function: double CRunDec::m1S2mSI(double m1S, std::pair<double,double>* mq,
//                                   double (*as)(double), int nl, int nloops, double fdelm)
// The function pointer passed should contain the adress of a function computing alpha_s in dependence of mu
double CRunDec::m1S2mSI(double m1S, std::pair<double,double>* mq,
                        double (*as)(double), int nl, int nloops, double fdelm) {
    if(as == NULL) {
      cout << "Pointer to as == NULL! Aborting..." << endl;
      RETURN
    }
	double mMS1 = 0;
	double mMS = m1S;
	double acc = 10e-6; 
	while(fabs(mMS1 - mMS) > acc) {
		mMS1 = mMS; 
        mMS = m1S2mMS(m1S, mq, as(mMS1), mMS1, nl, nloops, fdelm);
	}
	return mMS;
}

double CRunDec::exOS2RS(double api, double mmu, double nnuf, int nnl, int nloops) {
	if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     double lmm = log(mmu*mmu/(nnuf*nnuf));
     double ret[5];
     double Nm[] = {1.0, 1.0, 1.0, 0.563, 0.547, 0.527};
     double nu[] = {1.0, 1.0, 1.0, 0.3951, 0.3696, 0.3289};
     double ctil[4][6];
     ctil[1][3] = -0.1638;
     ctil[1][4] = -0.1054;
     ctil[1][5] = 0.0238;
     ctil[2][3] = 0.2372;
     ctil[2][4] = 0.2736;
     ctil[2][5] = 0.3265;
     ctil[3][3] = 0.0217337;
     ctil[3][4] = 0.0423924;
     ctil[3][5] = 0.0595309;
     ctil[0][0] = 1;
     ctil[0][1] = 1;
     ctil[0][2] = 1;
     ctil[0][3] = 1;
     ctil[0][4] = 1;
     ctil[0][5] = 1;


     ret[0] = 0.0;
     ret[1] = api*nnuf*Pi*(1. + ctil[1][nnl] + ctil[2][nnl] + ctil[3][nnl])*Nm[nnl];
     ret[2] = api*api*(((11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(-1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(2. + nu[nnl])/tgamma(1. + nu[nnl]))*Nm[nnl])/
        2. - (nnuf*Pi*(1. + ctil[1][nnl] + ctil[2][nnl] + ctil[3][nnl])*
         (-33.*lmm + 2.*nnl*lmm)*Nm[nnl])/12.);

     ret[3] = api*api*api*(((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(2. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(3. + nu[nnl])/tgamma(1. + nu[nnl]))*Nm[nnl])/
        4. + nnuf*Pi*(1. + ctil[1][nnl] + ctil[2][nnl] + ctil[3][nnl])*
        lmm*((102. - (38.*nnl)/3.)/16. + 
         ((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*lmm)/16.)*Nm[nnl] + 
       ((11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(-1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(2. + nu[nnl])/tgamma(1. + nu[nnl]))*
         ((11.*lmm)/2. - (nnl*lmm)/3.)*Nm[nnl])/2.);

     ret[4] = api*api*api*api*(((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(2. + nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(3. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(4. + nu[nnl])/tgamma(1. + nu[nnl]))*Nm[nnl])/
        8. + ((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(2. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(3. + nu[nnl])/tgamma(1. + nu[nnl]))*
         ((33.*lmm)/4. - (nnl*lmm)/2.)*Nm[nnl])/4. + 
       nnuf*Pi*(1. + ctil[1][nnl] + ctil[2][nnl] + ctil[3][nnl])*
        lmm*((2857./2. - (5033.*nnl)/18. + (325.*nnl*nnl)/54.)/64. + 
         (5.*(102. - (38.*nnl)/3.)*(11. - (2.*nnl)/3.)*lmm)/128. + 
         ((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*lmm*lmm)/64.)*Nm[nnl] + 
       ((11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(-1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(2. + nu[nnl])/tgamma(1. + nu[nnl]))*
         (612.*lmm - 76.*nnl*lmm + 
          1089.*lmm*lmm - 132.*nnl*lmm*lmm + 
          4.*nnl*nnl*lmm*lmm)*Nm[nnl])/96.);


     double res = 0.0;
     for(int i = 0; i <= nloops; i++) {
       res += ret[i];
     }
	 return res;
}

double CRunDec::exOS2RSp(double api, double mmu, double nnuf, int nnl, int nloops) {
	if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     double lmm = log(mmu*mmu/(nnuf*nnuf));
     double ret[5];
     double Nm[] = {1.0, 1.0, 1.0, 0.563, 0.547, 0.527};
     double nu[] = {1.0, 1.0, 1.0, 0.3951, 0.3696, 0.3289};
     double ctil[4][6];
     ctil[1][3] = -0.1638;
     ctil[1][4] = -0.1054;
     ctil[1][5] = 0.0238;
     ctil[2][3] = 0.2372;
     ctil[2][4] = 0.2736;
     ctil[2][5] = 0.3265;
     ctil[3][3] = 0.0217337;
     ctil[3][4] = 0.0423924;
     ctil[3][5] = 0.0595309;
     ctil[0][0] = 1;
     ctil[0][1] = 1;
     ctil[0][2] = 1;
     ctil[0][3] = 1;
     ctil[0][4] = 1;
     ctil[0][5] = 1;


     ret[0] = 0.0;
     ret[1] = 0.0;
     ret[2] = api*api*((11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(-1. + nu[nnl]))/
         tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(nu[nnl]))/
         tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(1. + nu[nnl]))/
         tgamma(nu[nnl]) + tgamma(2. + nu[nnl])/tgamma(1. + nu[nnl]))*Nm[nnl])/2.;

     ret[3] = api*api*api*(((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(2. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(3. + nu[nnl])/tgamma(1. + nu[nnl]))*Nm[nnl])/
        4. + ((11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(-1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(2. + nu[nnl])/tgamma(1. + nu[nnl]))*
         ((11.*lmm)/2. - (nnl*lmm)/3.)*Nm[nnl])/2.);

     ret[4] = api*api*api*api*(((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(2. + nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(3. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(4. + nu[nnl])/tgamma(1. + nu[nnl]))*Nm[nnl])/
        8. + ((11. - (2.*nnl)/3.)*(11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(2. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(3. + nu[nnl])/tgamma(1. + nu[nnl]))*
         ((33.*lmm)/4. - (nnl*lmm)/2.)*Nm[nnl])/4. + 
       ((11. - (2.*nnl)/3.)*nnuf*Pi*((ctil[3][nnl]*tgamma(-1. + nu[nnl]))/
           tgamma(-2. + nu[nnl]) + (ctil[2][nnl]*tgamma(nu[nnl]))/
           tgamma(-1. + nu[nnl]) + (ctil[1][nnl]*tgamma(1. + nu[nnl]))/
           tgamma(nu[nnl]) + tgamma(2. + nu[nnl])/tgamma(1. + nu[nnl]))*
         (612.*lmm - 76.*nnl*lmm + 
          1089.*lmm*lmm - 132.*nnl*lmm*lmm + 
          4.*nnl*nnl*lmm*lmm)*Nm[nnl])/96.);

     double res = 0.0;
     for(int i = 0; i <= nloops; i++) {
       res += ret[i];
     }
	 return res;
}

// Function: CRunDec::mOS2mRS(double mOS, std::pair<double,double>* mq, double asmu,
//                            double mu, double nuf, int nl, int nloops, bool prime)
double CRunDec::mOS2mRS(double mOS, std::pair<double,double>* mq, double asmu,
                        double mu, double nuf, int nl, int nloops, bool prime) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     if(prime) {
	   return mOS - exOS2RSp(asmu/Pi, mu, nuf, nl, nloops);
     } else {
	   return mOS - exOS2RS(asmu/Pi, mu, nuf, nl, nloops);
     }
}

// Function: double CRunDec::mMS2mRS(double mMS, std::pair<double,double>* mq, double asmu,
//                                   double mu, double nuf, int nl, int nloops, double fdelm, bool prime)
double CRunDec::mMS2mRS(double mMS, std::pair<double,double>* mq, double asmu,
                        double mu, double nuf, int nl, int nloops, double fdelm, bool prime) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }
     if(prime) {
	   return mMS2mOSmod(mMS, mq, asmu, mu, nl+1, nloops, fdelm) - exOS2RSp(asmu/Pi, mu, nuf, nl, nloops);
     } else {
	   return mMS2mOSmod(mMS, mq, asmu, mu, nl+1, nloops, fdelm) - exOS2RS(asmu/Pi, mu, nuf, nl, nloops);
     }
}

// Function: double CRunDec::mRS2mMS(double mRS, std::pair<double,double>* mq, double asmu,
//                                   double mu, double muf, int nl, int nloops, double fdelm, bool prime)
double CRunDec::mRS2mMS(double mRS, std::pair<double,double>* mq, double asmu,
                        double mu, double muf, int nl, int nloops, double fdelm, bool prime) {
     if(nloops<0||nloops>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
       RETURN
     }

     double lowbound = mRS - mRS/5.;
     double highbound = mRS + mRS/5.;
     bool found = false;
     double f1 = mMS2mRS(lowbound, mq, asmu, mu, muf, nl, nloops, fdelm, prime) - mRS;
     double f2 = mMS2mRS(highbound, mq, asmu, mu, muf, nl, nloops, fdelm, prime) - mRS;
     for(int j = 0; j < 1000; j++) {
       if(f1*f2 < 0.0) {
         found = true;
         break;
       }
       if(fabs(f1) < fabs(f2))
         f1 = mMS2mRS(lowbound += 1.5*(lowbound - highbound), mq, asmu, mu, muf, nl, nloops, fdelm, prime) - mRS;
       else
         f2 = mMS2mRS(highbound -= 1.5*(lowbound - highbound), mq, asmu, mu, muf, nl, nloops, fdelm, prime) - mRS;
     }
     if(found) {
     double acc = 10e-10;
     double dx;
     double xmid;
     double mMS = f1 < 0.0 ? (dx=highbound-lowbound,lowbound) : (dx=lowbound-highbound,highbound);
     for(int j = 0; j < 1000; j++) {
       f2 = mMS2mRS(xmid = mMS+(dx *= 0.5), mq, asmu, mu, muf, nl, nloops, fdelm, prime) - mRS;
       if(f2 <= 0.0)
         mMS = xmid;
       if(fabs(dx) < acc || f2 == 0.0)
         return mMS;
     }
     }
     return 0.0;
}

// Function: double CRunDec::mRS2mSI(double mRS, std::pair<double,double>* mq, double (*as)(double),
//                                   double muf, int nl, int nloops, double fdelm, bool prime)
// The function pointer passed should contain the adress of a function computing alpha_s in dependence of mu
double CRunDec::mRS2mSI(double mRS, std::pair<double,double>* mq, double (*as)(double),
                        double muf, int nl, int nloops, double fdelm, bool prime) {
    if(as == NULL) {
      cout << "Pointer to as == NULL! Aborting..." << endl;
      RETURN
    }
	double mMS1 = 0;
	double mMS = mRS;
	double acc = 10e-6; 
	while(fabs(mMS1 - mMS) > acc) {
		mMS1 = mMS; 
        mMS = mRS2mMS(mRS, mq, as(mMS1), mMS1, muf, nl, nloops, fdelm, prime);
	}
	return mMS;
}

// Function: double CRunDec::mMS2mRGImod(double mMS, double asmu, int nl)
// See 'mMS2mRGI' but 'Alphas/Pi -> AlphaS*2*Beta0/Pi'
double CRunDec::mMS2mRGImod(double mMS, double asmu, int nl){
     if(nl<0||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl==0){
       return (double) mMS;
     }
     double bet0= 11./4. - Nf/6.;     
     double cAsmu= this->fSetcx((2*bet0*asmu)/Pi, nl); 
     return (double) mMS/cAsmu;     
}

double CRunDec::mMS2mkinA(double mbMSmus, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops) {
    if(nloops<0||nloops>3){
      cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
      RETURN
    }
    double ret[4];
    double lnmusms = log(mus*mus/(mbMSmus*mbMSmus));
    double lnmufmus = log(2*mufac/mus);
    double lnmusmcms = 0.0;
    double lnmcmusmcms = 0.0;
    double lnmcmsms = 0.0;
    double mcMSmusmc = 0.0;
    double musmc = 0.0;
    if(mcMSmusmcin > 0.0) {
       lnmusmcms = log(mus*mus/(mcMSmusmcin*mcMSmusmcin));
       lnmcmusmcms = log(musmcin*musmcin/(mcMSmusmcin*mcMSmusmcin));
       lnmcmsms = log(mcMSmusmcin/mbMSmus);
       mcMSmusmc = mcMSmusmcin;
       musmc = musmcin;
    }
    
    ret[0] = mbMSmus;
    ret[1] = apinlmus*((1.3333333333333333 + lnmusms)*mbMSmus - (16*mufac)/9. - (2*pow(mufac,2))/(3.*mbMSmus));
    ret[2] = apinlmus*apinlmus*(((19*pow(mcMSmusmc,6))/225. - (4*lnmcmsms*pow(mcMSmusmc,6))/45.)/pow(mbMSmus,5) + ((463*pow(mcMSmusmc,8))/39200. - (3*lnmcmsms*pow(mcMSmusmc,8))/140.)/pow(mbMSmus,7) + 
   ((997*pow(mcMSmusmc,10))/297675. - (8*lnmcmsms*pow(mcMSmusmc,10))/945.)/pow(mbMSmus,9) + ((1229*pow(mcMSmusmc,12))/940896. - (5*lnmcmsms*pow(mcMSmusmc,12))/1188.)/pow(mbMSmus,11) + 
   ((15371*pow(mcMSmusmc,14))/2.5050025e7 - (12*lnmcmsms*pow(mcMSmusmc,14))/5005.)/pow(mbMSmus,13) + ((529*pow(mcMSmusmc,16))/1.6224e6 - (7*lnmcmsms*pow(mcMSmusmc,16))/4680.)/pow(mbMSmus,15) + 
   ((16277*pow(mcMSmusmc,18))/8.6028075e7 - (16*lnmcmsms*pow(mcMSmusmc,18))/16065.)/pow(mbMSmus,17) + ((39163*pow(mcMSmusmc,20))/3.338528e8 - (9*lnmcmsms*pow(mcMSmusmc,20))/12920.)/pow(mbMSmus,19) + 
   ((39833*pow(mcMSmusmc,22))/5.20109667e8 - (20*lnmcmsms*pow(mcMSmusmc,22))/39501.)/pow(mbMSmus,21) + ((9727*pow(mcMSmusmc,24))/1.866312e8 - (11*lnmcmsms*pow(mcMSmusmc,24))/28980.)/pow(mbMSmus,23) + 
   (mcMSmusmc*pow(Pi,2))/6. + (pow(mcMSmusmc,3)*pow(Pi,2))/(6.*pow(mbMSmus,2)) + mufac*(-31.85185185185185 + lnmufmus*(9.777777777777779 - (16*NLOSKIN)/27.) + (128*NLOSKIN)/81. + (8*pow(Pi,2))/9.) + 
   ((-151*pow(mcMSmusmc,4))/216. + (13*lnmcmsms*pow(mcMSmusmc,4))/18. - (pow(lnmcmsms,2)*pow(mcMSmusmc,4))/3. - (pow(mcMSmusmc,4)*pow(Pi,2))/18.)/pow(mbMSmus,3) + 
   (-pow(mcMSmusmc,2) + pow(mufac,2)*(-9.222222222222221 + (2*lnmusms)/3. + lnmufmus*(3.6666666666666665 - (2*NLOSKIN)/9.) + (13*NLOSKIN)/27. + pow(Pi,2)/3.))/mbMSmus + 
   mbMSmus*(9.100694444444445 + (2*lnmusmcms)/9. + (6.708333333333333 + lnmusmcms/6.)*lnmusms + (15*pow(lnmusms,2))/8. + (5*pow(Pi,2))/18. + (ln2*pow(Pi,2))/9. + 
      NLMSOS*(-0.4930555555555556 - (13*lnmusms)/36. - pow(lnmusms,2)/12. - pow(Pi,2)/18.) - Zeta3/6.));
    ret[3] = apinlmus*apinlmus*apinlmus*((((9727*pow(mcMSmusmc,24))/2.799468e8 - (11*lnmcmsms*pow(mcMSmusmc,24))/43470.)*pow(mufac,2))/pow(mbMSmus,25) + (29783*mcMSmusmc*pow(Pi,2))/2430. - (1199*ln2*mcMSmusmc*pow(Pi,2))/81. - 
   (31*lnmcmsms*mcMSmusmc*pow(Pi,2))/18. + (lnmcmusmcms*mcMSmusmc*pow(Pi,2))/6. + (lnmusmcms*mcMSmusmc*pow(Pi,2))/18. + (31*lnmusms*mcMSmusmc*pow(Pi,2))/36. + (13*mcMSmusmc*pow(Pi,3))/162. + 
   NLMSOS*((-7*mcMSmusmc*pow(Pi,2))/27. + (2*ln2*mcMSmusmc*pow(Pi,2))/9. + (lnmcmsms*mcMSmusmc*pow(Pi,2))/9. - (lnmusms*mcMSmusmc*pow(Pi,2))/18.) + 
   ((149*pow(mcMSmusmc,5)*pow(Pi,2))/1080. - (13*lnmcmsms*pow(mcMSmusmc,5)*pow(Pi,2))/45. + (pow(mcMSmusmc,3)*pow(mufac,2)*pow(Pi,2))/9. - (pow(mcMSmusmc,5)*NLMSOS*pow(Pi,2))/20. + 
      (pow(mcMSmusmc,5)*pow(Pi,3))/240.)/pow(mbMSmus,4) + ((262769*pow(mcMSmusmc,7)*pow(Pi,2))/1.1907e7 - (29*lnmcmsms*pow(mcMSmusmc,7)*pow(Pi,2))/525. - 
      (5*pow(mcMSmusmc,7)*NLMSOS*pow(Pi,2))/756. - (25*pow(mcMSmusmc,7)*pow(Pi,3))/18144.)/pow(mbMSmus,6) + 
   ((416029*pow(mcMSmusmc,9)*pow(Pi,2))/6.001128e7 - (53*lnmcmsms*pow(mcMSmusmc,9)*pow(Pi,2))/2205. - (7*pow(mcMSmusmc,9)*NLMSOS*pow(Pi,2))/3240. - (77*pow(mcMSmusmc,9)*pow(Pi,3))/103680.)/
    pow(mbMSmus,8) + ((1055689*pow(mcMSmusmc,11)*pow(Pi,2))/3.4577928e8 - (85*lnmcmsms*pow(mcMSmusmc,11)*pow(Pi,2))/6237. - (3*pow(mcMSmusmc,11)*NLMSOS*pow(Pi,2))/3080. - 
      (17*pow(mcMSmusmc,11)*pow(Pi,3))/39424.)/pow(mbMSmus,10) + ((31786481*pow(mcMSmusmc,13)*pow(Pi,2))/1.967766372e10 - (125*lnmcmsms*pow(mcMSmusmc,13)*pow(Pi,2))/14157. - 
      (11*pow(mcMSmusmc,13)*NLMSOS*pow(Pi,2))/21060. - (1771*pow(mcMSmusmc,13)*pow(Pi,3))/6.469632e6)/pow(mbMSmus,12) + 
   ((17346493*pow(mcMSmusmc,15)*pow(Pi,2))/1.808754948e10 - (173*lnmcmsms*pow(mcMSmusmc,15)*pow(Pi,2))/27885. - (13*pow(mcMSmusmc,15)*NLMSOS*pow(Pi,2))/41580. - 
      (377*pow(mcMSmusmc,15)*pow(Pi,3))/2.02752e6)/pow(mbMSmus,14) + ((22757641*pow(mcMSmusmc,17)*pow(Pi,2))/3.6923796e10 - (229*lnmcmsms*pow(mcMSmusmc,17)*pow(Pi,2))/49725. - 
      (5*pow(mcMSmusmc,17)*NLMSOS*pow(Pi,2))/24752. - (1925*pow(mcMSmusmc,17)*pow(Pi,3))/1.4483456e7)/pow(mbMSmus,16) + 
   ((86836957*pow(mcMSmusmc,19)*pow(Pi,2))/2.0687188752e11 - (293*lnmcmsms*pow(mcMSmusmc,19)*pow(Pi,2))/82365. - (17*pow(mcMSmusmc,19)*NLMSOS*pow(Pi,2))/123120. - 
      (99671*pow(mcMSmusmc,19)*pow(Pi,3))/1.00859904e9)/pow(mbMSmus,18) + ((846435761*pow(mcMSmusmc,21)*pow(Pi,2))/2.83231951884e12 - (365*lnmcmsms*pow(mcMSmusmc,21)*pow(Pi,2))/128877. - 
      (19*pow(mcMSmusmc,21)*NLMSOS*pow(Pi,2))/192780. - (127699*pow(mcMSmusmc,21)*pow(Pi,3))/1.684537344e9)/pow(mbMSmus,20) + 
   ((171475369*pow(mcMSmusmc,23)*pow(Pi,2))/7.7816811996e11 - (445*lnmcmsms*pow(mcMSmusmc,23)*pow(Pi,2))/192717. - (7*pow(mcMSmusmc,23)*NLMSOS*pow(Pi,2))/96140. - 
      (81991*pow(mcMSmusmc,23)*pow(Pi,3))/1.374683136e9)/pow(mbMSmus,22) + ((106559417*pow(mcMSmusmc,25)*pow(Pi,2))/6.374388636e11 - (533*lnmcmsms*pow(mcMSmusmc,25)*pow(Pi,2))/277725. - 
      (23*pow(mcMSmusmc,25)*NLMSOS*pow(Pi,2))/415800. - (5698043*pow(mcMSmusmc,25)*pow(Pi,3))/1.189085184e11)/pow(mbMSmus,24) + 
   ((18607*pow(mcMSmusmc,3)*pow(Pi,2))/1620. - (1199*ln2*pow(mcMSmusmc,3)*pow(Pi,2))/81. - (271*lnmcmsms*pow(mcMSmusmc,3)*pow(Pi,2))/162. + (lnmcmusmcms*pow(mcMSmusmc,3)*pow(Pi,2))/2. + 
      (lnmusmcms*pow(mcMSmusmc,3)*pow(Pi,2))/18. + (19*lnmusms*pow(mcMSmusmc,3)*pow(Pi,2))/36. + (mcMSmusmc*pow(mufac,2)*pow(Pi,2))/9. + (7*pow(mcMSmusmc,3)*pow(Pi,3))/108. + 
      NLMSOS*((-7*pow(mcMSmusmc,3)*pow(Pi,2))/54. + (2*ln2*pow(mcMSmusmc,3)*pow(Pi,2))/9. + (lnmcmsms*pow(mcMSmusmc,3)*pow(Pi,2))/9. - (lnmusms*pow(mcMSmusmc,3)*pow(Pi,2))/18.))/pow(mbMSmus,2) + 
   ((-19654121*pow(mcMSmusmc,6))/1.7496e7 - (5329*pow(lnmcmsms,2)*pow(mcMSmusmc,6))/6480. + (94*lnmcmusmcms*pow(mcMSmusmc,6))/225. + (19*lnmusmcms*pow(mcMSmusmc,6))/675. + 
      (139*lnmusms*pow(mcMSmusmc,6))/1350. + (311*pow(mcMSmusmc,6)*pow(Pi,2))/2592. - (11*ln2*pow(mcMSmusmc,6)*pow(Pi,2))/648. + 
      pow(mufac,2)*((-151*pow(mcMSmusmc,4))/324. + (13*lnmcmsms*pow(mcMSmusmc,4))/27. - (2*pow(lnmcmsms,2)*pow(mcMSmusmc,4))/9. - (pow(mcMSmusmc,4)*pow(Pi,2))/27.) + 
      lnmcmsms*((321193*pow(mcMSmusmc,6))/233280. - (8*lnmcmusmcms*pow(mcMSmusmc,6))/15. - (4*lnmusmcms*pow(mcMSmusmc,6))/135. - (2*lnmusms*pow(mcMSmusmc,6))/135. - 
         (113*pow(mcMSmusmc,6)*pow(Pi,2))/1296.) + NLMSOS*((2729*pow(mcMSmusmc,6))/30375. - (19*lnmusms*pow(mcMSmusmc,6))/675. + 
         lnmcmsms*((-2*pow(mcMSmusmc,6))/75. + (4*lnmusms*pow(mcMSmusmc,6))/135.) + (4*pow(mcMSmusmc,6)*pow(Pi,2))/405.) - (29*pow(mcMSmusmc,6)*Zeta3)/648.)/pow(mbMSmus,5) + 
   ((-74049121*pow(mcMSmusmc,8))/1.185408e9 - (9451*pow(lnmcmsms,2)*pow(mcMSmusmc,8))/30240. + (179*lnmcmusmcms*pow(mcMSmusmc,8))/2450. + (463*lnmusmcms*pow(mcMSmusmc,8))/117600. - 
      (53*lnmusms*pow(mcMSmusmc,8))/235200. + ((38*pow(mcMSmusmc,6))/675. - (8*lnmcmsms*pow(mcMSmusmc,6))/135.)*pow(mufac,2) + (1789*pow(mcMSmusmc,8)*pow(Pi,2))/90720. - 
      (7*ln2*pow(mcMSmusmc,8)*pow(Pi,2))/1536. + lnmcmsms*((2975311*pow(mcMSmusmc,8))/1.27008e7 - (6*lnmcmusmcms*pow(mcMSmusmc,8))/35. - (lnmusmcms*pow(mcMSmusmc,8))/140. + 
         (11*lnmusms*pow(mcMSmusmc,8))/280. - (139*pow(mcMSmusmc,8)*pow(Pi,2))/4608.) + 
      NLMSOS*((174787*pow(mcMSmusmc,8))/1.2348e7 - (463*lnmusms*pow(mcMSmusmc,8))/117600. + lnmcmsms*((-37*pow(mcMSmusmc,8))/14700. + (lnmusms*pow(mcMSmusmc,8))/140.) + 
         (pow(mcMSmusmc,8)*pow(Pi,2))/420.) - (173*pow(mcMSmusmc,8)*Zeta3)/9216.)/pow(mbMSmus,7) + 
   ((13434791647*pow(mcMSmusmc,10))/5.184974592e12 - (382589*pow(lnmcmsms,2)*pow(mcMSmusmc,10))/2.3328e6 + (298*lnmcmusmcms*pow(mcMSmusmc,10))/11907. + (997*lnmusmcms*pow(mcMSmusmc,10))/893025. - 
      (7811*lnmusms*pow(mcMSmusmc,10))/1.78605e6 + ((463*pow(mcMSmusmc,8))/58800. - (lnmcmsms*pow(mcMSmusmc,8))/70.)*pow(mufac,2) + (359801*pow(mcMSmusmc,10)*pow(Pi,2))/4.89888e7 - 
      (17*ln2*pow(mcMSmusmc,10)*pow(Pi,2))/8640. + lnmcmsms*((3017558449*pow(mcMSmusmc,10))/4.1150592e10 - (16*lnmcmusmcms*pow(mcMSmusmc,10))/189. - (8*lnmusmcms*pow(mcMSmusmc,10))/2835. + 
         (92*lnmusms*pow(mcMSmusmc,10))/2835. - (17*pow(mcMSmusmc,10)*pow(Pi,2))/1080.) + 
      NLMSOS*((2816347*pow(mcMSmusmc,10))/5.6260575e8 - (997*lnmusms*pow(mcMSmusmc,10))/893025. + lnmcmsms*((-398*pow(mcMSmusmc,10))/893025. + (8*lnmusms*pow(mcMSmusmc,10))/2835.) + 
         (8*pow(mcMSmusmc,10)*pow(Pi,2))/8505.) - (187*pow(mcMSmusmc,10)*Zeta3)/17280.)/pow(mbMSmus,9) + 
   ((185657253148457*pow(mcMSmusmc,12))/2.46471470784e16 - (15261023*pow(lnmcmsms,2)*pow(mcMSmusmc,12))/1.49688e8 + (899*lnmcmusmcms*pow(mcMSmusmc,12))/78408. + 
      (1229*lnmusmcms*pow(mcMSmusmc,12))/2.822688e6 - (19255*lnmusms*pow(mcMSmusmc,12))/5.645376e6 + ((1994*pow(mcMSmusmc,10))/893025. - (16*lnmcmsms*pow(mcMSmusmc,10))/2835.)*pow(mufac,2) + 
      (84041429*pow(mcMSmusmc,12)*pow(Pi,2))/2.29920768e10 - (175*ln2*pow(mcMSmusmc,12)*pow(Pi,2))/165888. + 
      lnmcmsms*((202523874613*pow(mcMSmusmc,12))/6.22402704e12 - (5*lnmcmusmcms*pow(mcMSmusmc,12))/99. - (5*lnmusmcms*pow(mcMSmusmc,12))/3564. + (175*lnmusms*pow(mcMSmusmc,12))/7128. - 
         (40553*pow(mcMSmusmc,12)*pow(Pi,2))/4.1472e6) + NLMSOS*((326802499*pow(mcMSmusmc,12))/1.3692859488e11 - (1229*lnmusms*pow(mcMSmusmc,12))/2.822688e6 + 
         lnmcmsms*((-391*pow(mcMSmusmc,12))/4.939704e6 + (5*lnmusms*pow(mcMSmusmc,12))/3564.) + (5*pow(mcMSmusmc,12)*pow(Pi,2))/10692.) - (59231*pow(mcMSmusmc,12)*Zeta3)/8.2944e6)/pow(mbMSmus,11) + 
   ((396218064296685469*pow(mcMSmusmc,14))/5.896309609846656e19 - (845153*pow(lnmcmsms,2)*pow(mcMSmusmc,14))/1.2108096e7 + (3166*lnmcmusmcms*pow(mcMSmusmc,14))/511225. + 
      (15371*lnmusmcms*pow(mcMSmusmc,14))/7.5150075e7 - (362077*lnmusms*pow(mcMSmusmc,14))/1.5030015e8 + ((1229*pow(mcMSmusmc,12))/1.411344e6 - (5*lnmcmsms*pow(mcMSmusmc,12))/1782.)*pow(mufac,2) + 
      (82285201*pow(mcMSmusmc,14)*pow(Pi,2))/3.87459072e10 - (23*ln2*pow(mcMSmusmc,14)*pow(Pi,2))/35840. + 
      lnmcmsms*((64922338969*pow(mcMSmusmc,14))/3.718698984e12 - (24*lnmcmusmcms*pow(mcMSmusmc,14))/715. - (4*lnmusmcms*pow(mcMSmusmc,14))/5005. + (94*lnmusms*pow(mcMSmusmc,14))/5005. - 
         (2161*pow(mcMSmusmc,14)*pow(Pi,2))/322560.) + NLMSOS*((108352091581*pow(mcMSmusmc,14))/8.1243243081e13 - (15371*lnmusms*pow(mcMSmusmc,14))/7.5150075e7 + 
         lnmcmsms*((1153*pow(mcMSmusmc,14))/2.25450225e8 + (4*lnmusms*pow(mcMSmusmc,14))/5005.) + (4*pow(mcMSmusmc,14)*pow(Pi,2))/15015.) - (3287*pow(mcMSmusmc,14)*Zeta3)/645120.)/pow(mbMSmus,13) + 
   ((516851278553981044967.*pow(mcMSmusmc,16))/9.509568138760686e22 - (232939907*pow(lnmcmsms,2)*pow(mcMSmusmc,16))/4.576860288e9 + (283*lnmcmusmcms*pow(mcMSmusmc,16))/76050. + 
      (529*lnmusmcms*pow(mcMSmusmc,16))/4.8672e6 - (16651*lnmusms*pow(mcMSmusmc,16))/9.7344e6 + ((30742*pow(mcMSmusmc,14))/7.5150075e7 - (8*lnmcmsms*pow(mcMSmusmc,14))/5005.)*pow(mufac,2) + 
      (58806560951*pow(mcMSmusmc,16)*pow(Pi,2))/4.326189170688e13 - (1001*ln2*pow(mcMSmusmc,16)*pow(Pi,2))/2.359296e6 + 
      lnmcmsms*((12430141121803*pow(mcMSmusmc,16))/1.1780838381312e15 - (14*lnmcmusmcms*pow(mcMSmusmc,16))/585. - (7*lnmusmcms*pow(mcMSmusmc,16))/14040. + (413*lnmusms*pow(mcMSmusmc,16))/28080. - 
         (565351*pow(mcMSmusmc,16)*pow(Pi,2))/1.15605504e8) + NLMSOS*((214558103603*pow(mcMSmusmc,16))/2.60460712512e14 - (529*lnmusms*pow(mcMSmusmc,16))/4.8672e6 + 
         lnmcmsms*((355*pow(mcMSmusmc,16))/1.4455584e7 + (7*lnmusms*pow(mcMSmusmc,16))/14040.) + (7*pow(mcMSmusmc,16)*pow(Pi,2))/42120.) - (885457*pow(mcMSmusmc,16)*Zeta3)/2.31211008e8)/pow(mbMSmus,15)\
    + ((218348329088480198615629.*pow(mcMSmusmc,18))/5.00576874275692e25 - (1014809351*pow(lnmcmsms,2)*pow(mcMSmusmc,18))/2.615348736e10 + (7678*lnmcmusmcms*pow(mcMSmusmc,18))/3.186225e6 + 
      (16277*lnmusmcms*pow(mcMSmusmc,18))/2.58084225e8 - (641587*lnmusms*pow(mcMSmusmc,18))/5.1616845e8 + ((529*pow(mcMSmusmc,16))/2.4336e6 - (7*lnmcmsms*pow(mcMSmusmc,16))/7020.)*pow(mufac,2) + 
      (26463251891*pow(mcMSmusmc,18)*pow(Pi,2))/2.845499424768e13 - (4147*ln2*pow(mcMSmusmc,18)*pow(Pi,2))/1.3934592e7 + 
      lnmcmsms*((22606266853306541*pow(mcMSmusmc,18))/3.2684758005112013e18 - (32*lnmcmusmcms*pow(mcMSmusmc,18))/1785. - (16*lnmusmcms*pow(mcMSmusmc,18))/48195. + (568*lnmusms*pow(mcMSmusmc,18))/48195. - 
         (52013*pow(mcMSmusmc,18)*pow(Pi,2))/1.3934592e7) + NLMSOS*((20555048260909*pow(mcMSmusmc,18))/3.76818092235585e16 - (16277*lnmusms*pow(mcMSmusmc,18))/2.58084225e8 + 
         lnmcmsms*((197062*pow(mcMSmusmc,18))/7.381208835e9 + (16*lnmusms*pow(mcMSmusmc,18))/48195.) + (16*pow(mcMSmusmc,18)*pow(Pi,2))/144585.) - (83291*pow(mcMSmusmc,18)*Zeta3)/2.7869184e7)/
    pow(mbMSmus,17) + ((192872454677623233846986449.*pow(mcMSmusmc,20))/5.449931397868209e28 - (46137065941*pow(lnmcmsms,2)*pow(mcMSmusmc,20))/1.5084957888e12 + 
      (5507*lnmcmusmcms*pow(mcMSmusmc,20))/3.338528e6 + (39163*lnmusmcms*pow(mcMSmusmc,20))/1.0015584e9 - (1855169*lnmusms*pow(mcMSmusmc,20))/2.0031168e9 + 
      ((32554*pow(mcMSmusmc,18))/2.58084225e8 - (32*lnmcmsms*pow(mcMSmusmc,18))/48195.)*pow(mufac,2) + (9007367733163*pow(mcMSmusmc,20)*pow(Pi,2))/1.34810154565632e16 - 
      (143*ln2*pow(mcMSmusmc,20)*pow(Pi,2))/655360. + lnmcmsms*((844305629180580989*pow(mcMSmusmc,20))/1.7558329821198565e20 - (9*lnmcmusmcms*pow(mcMSmusmc,20))/646. - 
         (3*lnmusmcms*pow(mcMSmusmc,20))/12920. + (249*lnmusms*pow(mcMSmusmc,20))/25840. - (156353*pow(mcMSmusmc,20)*pow(Pi,2))/5.308416e7) + 
      NLMSOS*((22189567531163017*pow(mcMSmusmc,20))/5.834712481735737e19 - (39163*lnmusms*pow(mcMSmusmc,20))/1.0015584e9 + 
         lnmcmsms*((1211963*pow(mcMSmusmc,20))/5.012799792e10 + (3*lnmusms*pow(mcMSmusmc,20))/12920.) + (pow(mcMSmusmc,20)*pow(Pi,2))/12920.) - (254791*pow(mcMSmusmc,20)*Zeta3)/1.0616832e8)/
    pow(mbMSmus,19) + ((31442447067404839835736513193.*pow(mcMSmusmc,22))/1.0790864167779053e31 - (3078965960711*pow(lnmcmsms,2)*pow(mcMSmusmc,22))/1.24450902576e14 + 
      (5066*lnmcmusmcms*pow(mcMSmusmc,22))/4.298427e6 + (39833*lnmusmcms*pow(mcMSmusmc,22))/1.560329001e9 - (116005*lnmusms*pow(mcMSmusmc,22))/1.64245158e8 + 
      ((39163*pow(mcMSmusmc,20))/5.007792e8 - (3*lnmcmsms*pow(mcMSmusmc,20))/6460.)*pow(mufac,2) + (174200135864459*pow(mcMSmusmc,22)*pow(Pi,2))/3.495434721951744e17 - 
      (38675*ln2*pow(mcMSmusmc,22)*pow(Pi,2))/2.33570304e8 + lnmcmsms*((25315115748447270877.*pow(mcMSmusmc,22))/7.242811051244408e21 - (40*lnmcmusmcms*pow(mcMSmusmc,22))/3591. - 
         (20*lnmusmcms*pow(mcMSmusmc,22))/118503. + (50*lnmusms*pow(mcMSmusmc,22))/6237. - (13926181*pow(mcMSmusmc,22)*pow(Pi,2))/5.8392576e9) + 
      NLMSOS*((1640519393726677*pow(mcMSmusmc,22))/5.946258455651274e18 - (39833*lnmusms*pow(mcMSmusmc,22))/1.560329001e9 + 
         lnmcmsms*((14290513*pow(mcMSmusmc,22))/6.89665418442e11 + (20*lnmusms*pow(mcMSmusmc,22))/118503.) + (20*pow(mcMSmusmc,22)*pow(Pi,2))/355509.) - (23017987*pow(mcMSmusmc,22)*Zeta3)/1.16785152e10)/
    pow(mbMSmus,21) + ((5769036308265946032465571119823289.*pow(mcMSmusmc,24))/2.369996943791664e36 - (953082228281*pow(lnmcmsms,2)*pow(mcMSmusmc,24))/4.664604200256e13 + 
      (10163*lnmcmusmcms*pow(mcMSmusmc,24))/1.166445e7 + (9727*lnmusmcms*pow(mcMSmusmc,24))/5.598936e8 - (615749*lnmusms*pow(mcMSmusmc,24))/1.1197872e9 + 
      ((79666*pow(mcMSmusmc,22))/1.560329001e9 - (40*lnmcmsms*pow(mcMSmusmc,22))/118503.)*pow(mufac,2) + (240817793781176357*pow(mcMSmusmc,24)*pow(Pi,2))/6.28867544642696e20 - 
      (877591*ln2*pow(mcMSmusmc,24)*pow(Pi,2))/6.79477248e9 + lnmcmsms*((73910608414092571571.*pow(mcMSmusmc,24))/2.8097278338127475e22 - (22*lnmcmusmcms*pow(mcMSmusmc,24))/2415. - 
         (11*lnmusmcms*pow(mcMSmusmc,24))/86940. + (1177*lnmusms*pow(mcMSmusmc,24))/173880. - (1620816161*pow(mcMSmusmc,24)*pow(Pi,2))/8.2216747008e11) + 
      NLMSOS*((7801530877413386647*pow(mcMSmusmc,24))/3.7763267488425775e22 - (9727*lnmusms*pow(mcMSmusmc,24))/5.598936e8 + 
         lnmcmsms*((73801799*pow(mcMSmusmc,24))/4.23178780752e12 + (11*lnmusms*pow(mcMSmusmc,24))/86940.) + (11*pow(mcMSmusmc,24)*pow(Pi,2))/260820.) - (2710689767*pow(mcMSmusmc,24)*Zeta3)/1.64433494016e12
      )/pow(mbMSmus,23) + mufac*(-807.820987654321 + (20047*NLOSKIN)/243. - (1292*pow(NLOSKIN,2))/729. + pow(lnmufmus,2)*(-53.77777777777778 + (176*NLOSKIN)/27. - (16*pow(NLOSKIN,2))/81.) + 
      (1022*pow(Pi,2))/27. - (208*NLOSKIN*pow(Pi,2))/81. + (8*pow(NLOSKIN,2)*pow(Pi,2))/243. - (2*pow(Pi,4))/3. + 
      lnmufmus*(373.037037037037 - (3356*NLOSKIN)/81. + (256*pow(NLOSKIN,2))/243. - (88*pow(Pi,2))/9. + (16*NLOSKIN*pow(Pi,2))/27.) + 114*Zeta3 - (140*NLOSKIN*Zeta3)/27.) + 
   ((7777*pow(mcMSmusmc,4))/7776. + (20*A4*pow(mcMSmusmc,4))/27. + (5*pow(ln2,4)*pow(mcMSmusmc,4))/162. + (7*pow(lnmcmsms,3)*pow(mcMSmusmc,4))/27. - (2*pow(mcMSmusmc,2)*pow(mufac,2))/3. + 
      (37*pow(mcMSmusmc,4)*pow(Pi,2))/108. + (2*pow(ln2,2)*pow(mcMSmusmc,4)*pow(Pi,2))/81. + (271*pow(mcMSmusmc,4)*pow(Pi,4))/19440. + 
      lnmcmusmcms*((-56*pow(mcMSmusmc,4))/27. - (2*pow(mcMSmusmc,4)*pow(Pi,2))/9.) + lnmusms*((-2899*pow(mcMSmusmc,4))/1296. - (13*pow(mcMSmusmc,4)*pow(Pi,2))/108.) + 
      lnmusmcms*((-151*pow(mcMSmusmc,4))/648. - (pow(mcMSmusmc,4)*pow(Pi,2))/54.) + 
      pow(lnmcmsms,2)*((-94*pow(mcMSmusmc,4))/27. - (4*lnmcmusmcms*pow(mcMSmusmc,4))/3. - (lnmusmcms*pow(mcMSmusmc,4))/9. - (13*lnmusms*pow(mcMSmusmc,4))/18. + (pow(mcMSmusmc,4)*pow(Pi,2))/4.) + 
      ln2*((5*pow(mcMSmusmc,4)*pow(Pi,2))/144. - (lnmcmsms*pow(mcMSmusmc,4)*pow(Pi,2))/9.) - (2309*pow(mcMSmusmc,4)*Zeta3)/864. + 
      NLMSOS*((-1423*pow(mcMSmusmc,4))/3888. - (2*pow(lnmcmsms,3)*pow(mcMSmusmc,4))/27. + pow(lnmcmsms,2)*((13*pow(mcMSmusmc,4))/54. + (lnmusms*pow(mcMSmusmc,4))/9.) - 
         (13*pow(mcMSmusmc,4)*pow(Pi,2))/324. + lnmusms*((151*pow(mcMSmusmc,4))/648. + (pow(mcMSmusmc,4)*pow(Pi,2))/54.) + 
         lnmcmsms*(-pow(mcMSmusmc,4)/12. - (13*lnmusms*pow(mcMSmusmc,4))/54. + (pow(mcMSmusmc,4)*pow(Pi,2))/27.) + (pow(mcMSmusmc,4)*Zeta3)/3.) + 
      lnmcmsms*((2443*pow(mcMSmusmc,4))/648. + (20*lnmcmusmcms*pow(mcMSmusmc,4))/9. + (13*lnmusmcms*pow(mcMSmusmc,4))/54. + (241*lnmusms*pow(mcMSmusmc,4))/108. - (67*pow(mcMSmusmc,4)*pow(Pi,2))/144. + 
         (14*pow(mcMSmusmc,4)*Zeta3)/9.))/pow(mbMSmus,3) + mbMSmus*(80.65343149862825 - (212*A4)/27. - (53*pow(ln2,4))/162. - (4*lnmcmusmcms)/9. + pow(lnmusmcms,2)/27. + 
      (22.519675925925927 + (5*lnmusmcms)/8.)*pow(lnmusms,2) + (1693*pow(lnmusms,3))/432. + (594941*pow(Pi,2))/38880. - (20*pow(ln2,2)*pow(Pi,2))/81. - (451*pow(Pi,4))/7776. + 
      ln2*((-199*pow(Pi,2))/54. + (lnmusmcms*pow(Pi,2))/27. + (37*lnmusms*pow(Pi,2))/54.) + 
      NLMSOS*(-9.736839849108367 + (8*A4)/27. + pow(ln2,4)/81. + (-2.553240740740741 - lnmusmcms/36.)*pow(lnmusms,2) - (41*pow(lnmusms,3))/108. - (313*pow(Pi,2))/216. + (2*pow(ln2,2)*pow(Pi,2))/81. + 
         (61*pow(Pi,4))/1944. + lnmusmcms*(-0.16435185185185186 - pow(Pi,2)/54.) + ln2*((-11*pow(Pi,2))/81. - (lnmusms*pow(Pi,2))/27.) + 
         lnmusms*(-7.705246913580247 - (13*lnmusmcms)/108. - (47*pow(Pi,2))/108. - (7*Zeta3)/9.) - (667*Zeta3)/216.) + 
      lnmusms*(64.06558641975309 - lnmcmusmcms/3. + (109*lnmusmcms)/36. + pow(lnmusmcms,2)/36. + (185*pow(Pi,2))/108. - (97*Zeta3)/36.) + lnmusmcms*(4.08912037037037 + (5*pow(Pi,2))/54. - Zeta3/18.) + 
      pow(NLMSOS,2)*(0.1008659122085048 + (13*pow(lnmusms,2))/216. + pow(lnmusms,3)/108. + (13*pow(Pi,2))/324. + lnmusms*(0.13734567901234568 + pow(Pi,2)/54.) + (7*Zeta3)/54.) - (77*Zeta3)/72. - 
      (1439*pow(Pi,2)*Zeta3)/432. + (1975*Zeta5)/216.) + ((-68*pow(mcMSmusmc,2))/9. - 2*lnmcmusmcms*pow(mcMSmusmc,2) - (lnmusmcms*pow(mcMSmusmc,2))/3. - (25*lnmusms*pow(mcMSmusmc,2))/6. + 
      ((2*pow(mcMSmusmc,2))/9. + (lnmusms*pow(mcMSmusmc,2))/3.)*NLMSOS - (13*pow(mcMSmusmc,2)*pow(Pi,2))/12. + 
      lnmcmsms*(-4*pow(mcMSmusmc,2) - (3*pow(mcMSmusmc,2)*pow(Pi,2))/2. + (pow(mcMSmusmc,2)*pow(Pi,4))/12.) - (11*pow(mcMSmusmc,2)*Zeta3)/2. + (3*pow(mcMSmusmc,2)*pow(Pi,2)*Zeta3)/4. + 
      pow(mufac,2)*(-204.54166666666666 + (4*lnmusmcms)/27. + (7*pow(lnmusms,2))/12. + (13805*NLOSKIN)/648. - (209*pow(NLOSKIN,2))/486. + 
         pow(lnmufmus,2)*(-20.166666666666668 + (22*NLOSKIN)/9. - (2*pow(NLOSKIN,2))/27.) + (1307*pow(Pi,2))/108. + (2*ln2*pow(Pi,2))/27. - (23*NLOSKIN*pow(Pi,2))/27. + 
         (pow(NLOSKIN,2)*pow(Pi,2))/81. - pow(Pi,4)/4. + lnmusms*(12.805555555555555 + lnmusmcms/9. - (13*NLOSKIN)/27. - pow(Pi,2)/3.) + 
         NLMSOS*(-0.3287037037037037 - (13*lnmusms)/54. - pow(lnmusms,2)/18. - pow(Pi,2)/27.) + 
         lnmufmus*(114.83333333333333 + lnmusms*(-3.6666666666666665 + (2*NLOSKIN)/9.) - (691*NLOSKIN)/54. + (26*pow(NLOSKIN,2))/81. - (11*pow(Pi,2))/3. + (2*NLOSKIN*pow(Pi,2))/9.) + (1535*Zeta3)/36. - 
         (35*NLOSKIN*Zeta3)/18.) + (5*pow(mcMSmusmc,2)*Zeta5)/2.)/mbMSmus);
    
    double mkin = 0.0;
    for(int i = 0; i <= nloops; i++) {
       mkin += ret[i];
    }
    return mkin;
}

double CRunDec::mMS2mkinB(double mbMSmus, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops) {
    if(nloops<0||nloops>3){
      cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
      RETURN
    }
    double ret[4];
    double lnmusms = log(mus*mus/(mbMSmus*mbMSmus));
    double lnmufmus = log(2*mufac/mus);
    double lnmusmcms = 0.0;
    double lnmcmusmcms = 0.0;
    double lnmcmsms = 0.0;
    double lnmusmusmc = 0.0;
    double mcMSmusmc = 0.0;
    double musmc = 0.0;
    if(mcMSmusmcin > 0.0) {
       lnmusmusmc = log(mus*mus/(musmcin*musmcin));
       lnmusmcms = log(mus*mus/(mcMSmusmcin*mcMSmusmcin));
       lnmcmusmcms = log(musmcin*musmcin/(mcMSmusmcin*mcMSmusmcin));
       lnmcmsms = log(mcMSmusmcin/mbMSmus);
       mcMSmusmc = mcMSmusmcin;
       musmc = musmcin;
    }
    
    ret[0] = mbMSmus;
    ret[1] = apinlmus*((1.3333333333333333 + lnmusms)*mbMSmus - (16*mufac)/9. - (2*pow(mufac,2))/(3.*mbMSmus));
    ret[2] = apinlmus*apinlmus*(((19*pow(mcMSmusmc,6))/225. - (4*lnmcmsms*pow(mcMSmusmc,6))/45.)/pow(mbMSmus,5) + ((463*pow(mcMSmusmc,8))/39200. - (3*lnmcmsms*pow(mcMSmusmc,8))/140.)/pow(mbMSmus,7) + 
   ((997*pow(mcMSmusmc,10))/297675. - (8*lnmcmsms*pow(mcMSmusmc,10))/945.)/pow(mbMSmus,9) + ((1229*pow(mcMSmusmc,12))/940896. - (5*lnmcmsms*pow(mcMSmusmc,12))/1188.)/pow(mbMSmus,11) + 
   ((15371*pow(mcMSmusmc,14))/2.5050025e7 - (12*lnmcmsms*pow(mcMSmusmc,14))/5005.)/pow(mbMSmus,13) + ((529*pow(mcMSmusmc,16))/1.6224e6 - (7*lnmcmsms*pow(mcMSmusmc,16))/4680.)/pow(mbMSmus,15) + 
   ((16277*pow(mcMSmusmc,18))/8.6028075e7 - (16*lnmcmsms*pow(mcMSmusmc,18))/16065.)/pow(mbMSmus,17) + ((39163*pow(mcMSmusmc,20))/3.338528e8 - (9*lnmcmsms*pow(mcMSmusmc,20))/12920.)/pow(mbMSmus,19) + 
   ((39833*pow(mcMSmusmc,22))/5.20109667e8 - (20*lnmcmsms*pow(mcMSmusmc,22))/39501.)/pow(mbMSmus,21) + ((9727*pow(mcMSmusmc,24))/1.866312e8 - (11*lnmcmsms*pow(mcMSmusmc,24))/28980.)/pow(mbMSmus,23) + 
   (mcMSmusmc*pow(Pi,2))/6. + (pow(mcMSmusmc,3)*pow(Pi,2))/(6.*pow(mbMSmus,2)) + 
   mufac*(-31.85185185185185 + (8*lnmusmcms)/27. + lnmufmus*(9.777777777777779 - (16*NLOSKIN)/27.) + (128*NLOSKIN)/81. + (8*pow(Pi,2))/9.) + 
   ((-151*pow(mcMSmusmc,4))/216. + (13*lnmcmsms*pow(mcMSmusmc,4))/18. - (pow(lnmcmsms,2)*pow(mcMSmusmc,4))/3. - (pow(mcMSmusmc,4)*pow(Pi,2))/18.)/pow(mbMSmus,3) + 
   (-pow(mcMSmusmc,2) + pow(mufac,2)*(-9.222222222222221 + lnmusmcms/9. + (2*lnmusms)/3. + lnmufmus*(3.6666666666666665 - (2*NLOSKIN)/9.) + (13*NLOSKIN)/27. + pow(Pi,2)/3.))/mbMSmus + 
   mbMSmus*(9.100694444444445 + (161*lnmusms)/24. + (15*pow(lnmusms,2))/8. + (5*pow(Pi,2))/18. + (ln2*pow(Pi,2))/9. + NLMSOS*(-0.4930555555555556 - (13*lnmusms)/36. - pow(lnmusms,2)/12. - pow(Pi,2)/18.) - 
      Zeta3/6.));
    ret[3] = apinlmus*apinlmus*apinlmus*((((9727*pow(mcMSmusmc,24))/2.799468e8 - (11*lnmcmsms*pow(mcMSmusmc,24))/43470.)*pow(mufac,2))/pow(mbMSmus,25) + (29783*mcMSmusmc*pow(Pi,2))/2430. - (1199*ln2*mcMSmusmc*pow(Pi,2))/81. - 
   (31*lnmcmsms*mcMSmusmc*pow(Pi,2))/18. + (lnmcmusmcms*mcMSmusmc*pow(Pi,2))/6. + (31*lnmusms*mcMSmusmc*pow(Pi,2))/36. + (13*mcMSmusmc*pow(Pi,3))/162. + 
   NLMSOS*((-7*mcMSmusmc*pow(Pi,2))/27. + (2*ln2*mcMSmusmc*pow(Pi,2))/9. + (lnmcmsms*mcMSmusmc*pow(Pi,2))/9. - (lnmusms*mcMSmusmc*pow(Pi,2))/18.) + 
   ((149*pow(mcMSmusmc,5)*pow(Pi,2))/1080. - (13*lnmcmsms*pow(mcMSmusmc,5)*pow(Pi,2))/45. + (pow(mcMSmusmc,3)*pow(mufac,2)*pow(Pi,2))/9. - (pow(mcMSmusmc,5)*NLMSOS*pow(Pi,2))/20. + 
      (pow(mcMSmusmc,5)*pow(Pi,3))/240.)/pow(mbMSmus,4) + ((262769*pow(mcMSmusmc,7)*pow(Pi,2))/1.1907e7 - (29*lnmcmsms*pow(mcMSmusmc,7)*pow(Pi,2))/525. - 
      (5*pow(mcMSmusmc,7)*NLMSOS*pow(Pi,2))/756. - (25*pow(mcMSmusmc,7)*pow(Pi,3))/18144.)/pow(mbMSmus,6) + 
   ((416029*pow(mcMSmusmc,9)*pow(Pi,2))/6.001128e7 - (53*lnmcmsms*pow(mcMSmusmc,9)*pow(Pi,2))/2205. - (7*pow(mcMSmusmc,9)*NLMSOS*pow(Pi,2))/3240. - (77*pow(mcMSmusmc,9)*pow(Pi,3))/103680.)/
    pow(mbMSmus,8) + ((1055689*pow(mcMSmusmc,11)*pow(Pi,2))/3.4577928e8 - (85*lnmcmsms*pow(mcMSmusmc,11)*pow(Pi,2))/6237. - (3*pow(mcMSmusmc,11)*NLMSOS*pow(Pi,2))/3080. - 
      (17*pow(mcMSmusmc,11)*pow(Pi,3))/39424.)/pow(mbMSmus,10) + ((31786481*pow(mcMSmusmc,13)*pow(Pi,2))/1.967766372e10 - (125*lnmcmsms*pow(mcMSmusmc,13)*pow(Pi,2))/14157. - 
      (11*pow(mcMSmusmc,13)*NLMSOS*pow(Pi,2))/21060. - (1771*pow(mcMSmusmc,13)*pow(Pi,3))/6.469632e6)/pow(mbMSmus,12) + 
   ((17346493*pow(mcMSmusmc,15)*pow(Pi,2))/1.808754948e10 - (173*lnmcmsms*pow(mcMSmusmc,15)*pow(Pi,2))/27885. - (13*pow(mcMSmusmc,15)*NLMSOS*pow(Pi,2))/41580. - 
      (377*pow(mcMSmusmc,15)*pow(Pi,3))/2.02752e6)/pow(mbMSmus,14) + ((22757641*pow(mcMSmusmc,17)*pow(Pi,2))/3.6923796e10 - (229*lnmcmsms*pow(mcMSmusmc,17)*pow(Pi,2))/49725. - 
      (5*pow(mcMSmusmc,17)*NLMSOS*pow(Pi,2))/24752. - (1925*pow(mcMSmusmc,17)*pow(Pi,3))/1.4483456e7)/pow(mbMSmus,16) + 
   ((86836957*pow(mcMSmusmc,19)*pow(Pi,2))/2.0687188752e11 - (293*lnmcmsms*pow(mcMSmusmc,19)*pow(Pi,2))/82365. - (17*pow(mcMSmusmc,19)*NLMSOS*pow(Pi,2))/123120. - 
      (99671*pow(mcMSmusmc,19)*pow(Pi,3))/1.00859904e9)/pow(mbMSmus,18) + ((846435761*pow(mcMSmusmc,21)*pow(Pi,2))/2.83231951884e12 - (365*lnmcmsms*pow(mcMSmusmc,21)*pow(Pi,2))/128877. - 
      (19*pow(mcMSmusmc,21)*NLMSOS*pow(Pi,2))/192780. - (127699*pow(mcMSmusmc,21)*pow(Pi,3))/1.684537344e9)/pow(mbMSmus,20) + 
   ((171475369*pow(mcMSmusmc,23)*pow(Pi,2))/7.7816811996e11 - (445*lnmcmsms*pow(mcMSmusmc,23)*pow(Pi,2))/192717. - (7*pow(mcMSmusmc,23)*NLMSOS*pow(Pi,2))/96140. - 
      (81991*pow(mcMSmusmc,23)*pow(Pi,3))/1.374683136e9)/pow(mbMSmus,22) + ((106559417*pow(mcMSmusmc,25)*pow(Pi,2))/6.374388636e11 - (533*lnmcmsms*pow(mcMSmusmc,25)*pow(Pi,2))/277725. - 
      (23*pow(mcMSmusmc,25)*NLMSOS*pow(Pi,2))/415800. - (5698043*pow(mcMSmusmc,25)*pow(Pi,3))/1.189085184e11)/pow(mbMSmus,24) + 
   ((18607*pow(mcMSmusmc,3)*pow(Pi,2))/1620. - (1199*ln2*pow(mcMSmusmc,3)*pow(Pi,2))/81. - (271*lnmcmsms*pow(mcMSmusmc,3)*pow(Pi,2))/162. + (lnmcmusmcms*pow(mcMSmusmc,3)*pow(Pi,2))/2. + 
      (19*lnmusms*pow(mcMSmusmc,3)*pow(Pi,2))/36. + (mcMSmusmc*pow(mufac,2)*pow(Pi,2))/9. + (7*pow(mcMSmusmc,3)*pow(Pi,3))/108. + 
      NLMSOS*((-7*pow(mcMSmusmc,3)*pow(Pi,2))/54. + (2*ln2*pow(mcMSmusmc,3)*pow(Pi,2))/9. + (lnmcmsms*pow(mcMSmusmc,3)*pow(Pi,2))/9. - (lnmusms*pow(mcMSmusmc,3)*pow(Pi,2))/18.))/pow(mbMSmus,2) + 
   ((-19654121*pow(mcMSmusmc,6))/1.7496e7 - (5329*pow(lnmcmsms,2)*pow(mcMSmusmc,6))/6480. + (94*lnmcmusmcms*pow(mcMSmusmc,6))/225. + (139*lnmusms*pow(mcMSmusmc,6))/1350. + 
      (311*pow(mcMSmusmc,6)*pow(Pi,2))/2592. - (11*ln2*pow(mcMSmusmc,6)*pow(Pi,2))/648. + 
      pow(mufac,2)*((-151*pow(mcMSmusmc,4))/324. + (13*lnmcmsms*pow(mcMSmusmc,4))/27. - (2*pow(lnmcmsms,2)*pow(mcMSmusmc,4))/9. - (pow(mcMSmusmc,4)*pow(Pi,2))/27.) + 
      lnmcmsms*((321193*pow(mcMSmusmc,6))/233280. - (8*lnmcmusmcms*pow(mcMSmusmc,6))/15. - (2*lnmusms*pow(mcMSmusmc,6))/135. - (113*pow(mcMSmusmc,6)*pow(Pi,2))/1296.) + 
      NLMSOS*((2729*pow(mcMSmusmc,6))/30375. - (19*lnmusms*pow(mcMSmusmc,6))/675. + lnmcmsms*((-2*pow(mcMSmusmc,6))/75. + (4*lnmusms*pow(mcMSmusmc,6))/135.) + (4*pow(mcMSmusmc,6)*pow(Pi,2))/405.) - 
      (29*pow(mcMSmusmc,6)*Zeta3)/648.)/pow(mbMSmus,5) + ((-74049121*pow(mcMSmusmc,8))/1.185408e9 - (9451*pow(lnmcmsms,2)*pow(mcMSmusmc,8))/30240. + (179*lnmcmusmcms*pow(mcMSmusmc,8))/2450. - 
      (53*lnmusms*pow(mcMSmusmc,8))/235200. + ((38*pow(mcMSmusmc,6))/675. - (8*lnmcmsms*pow(mcMSmusmc,6))/135.)*pow(mufac,2) + (1789*pow(mcMSmusmc,8)*pow(Pi,2))/90720. - 
      (7*ln2*pow(mcMSmusmc,8)*pow(Pi,2))/1536. + lnmcmsms*((2975311*pow(mcMSmusmc,8))/1.27008e7 - (6*lnmcmusmcms*pow(mcMSmusmc,8))/35. + (11*lnmusms*pow(mcMSmusmc,8))/280. - 
         (139*pow(mcMSmusmc,8)*pow(Pi,2))/4608.) + NLMSOS*((174787*pow(mcMSmusmc,8))/1.2348e7 - (463*lnmusms*pow(mcMSmusmc,8))/117600. + 
         lnmcmsms*((-37*pow(mcMSmusmc,8))/14700. + (lnmusms*pow(mcMSmusmc,8))/140.) + (pow(mcMSmusmc,8)*pow(Pi,2))/420.) - (173*pow(mcMSmusmc,8)*Zeta3)/9216.)/pow(mbMSmus,7) + 
   ((13434791647*pow(mcMSmusmc,10))/5.184974592e12 - (382589*pow(lnmcmsms,2)*pow(mcMSmusmc,10))/2.3328e6 + (298*lnmcmusmcms*pow(mcMSmusmc,10))/11907. - (7811*lnmusms*pow(mcMSmusmc,10))/1.78605e6 + 
      ((463*pow(mcMSmusmc,8))/58800. - (lnmcmsms*pow(mcMSmusmc,8))/70.)*pow(mufac,2) + (359801*pow(mcMSmusmc,10)*pow(Pi,2))/4.89888e7 - (17*ln2*pow(mcMSmusmc,10)*pow(Pi,2))/8640. + 
      lnmcmsms*((3017558449*pow(mcMSmusmc,10))/4.1150592e10 - (16*lnmcmusmcms*pow(mcMSmusmc,10))/189. + (92*lnmusms*pow(mcMSmusmc,10))/2835. - (17*pow(mcMSmusmc,10)*pow(Pi,2))/1080.) + 
      NLMSOS*((2816347*pow(mcMSmusmc,10))/5.6260575e8 - (997*lnmusms*pow(mcMSmusmc,10))/893025. + lnmcmsms*((-398*pow(mcMSmusmc,10))/893025. + (8*lnmusms*pow(mcMSmusmc,10))/2835.) + 
         (8*pow(mcMSmusmc,10)*pow(Pi,2))/8505.) - (187*pow(mcMSmusmc,10)*Zeta3)/17280.)/pow(mbMSmus,9) + 
   ((185657253148457*pow(mcMSmusmc,12))/2.46471470784e16 - (15261023*pow(lnmcmsms,2)*pow(mcMSmusmc,12))/1.49688e8 + (899*lnmcmusmcms*pow(mcMSmusmc,12))/78408. - 
      (19255*lnmusms*pow(mcMSmusmc,12))/5.645376e6 + ((1994*pow(mcMSmusmc,10))/893025. - (16*lnmcmsms*pow(mcMSmusmc,10))/2835.)*pow(mufac,2) + (84041429*pow(mcMSmusmc,12)*pow(Pi,2))/2.29920768e10 - 
      (175*ln2*pow(mcMSmusmc,12)*pow(Pi,2))/165888. + lnmcmsms*((202523874613*pow(mcMSmusmc,12))/6.22402704e12 - (5*lnmcmusmcms*pow(mcMSmusmc,12))/99. + (175*lnmusms*pow(mcMSmusmc,12))/7128. - 
         (40553*pow(mcMSmusmc,12)*pow(Pi,2))/4.1472e6) + NLMSOS*((326802499*pow(mcMSmusmc,12))/1.3692859488e11 - (1229*lnmusms*pow(mcMSmusmc,12))/2.822688e6 + 
         lnmcmsms*((-391*pow(mcMSmusmc,12))/4.939704e6 + (5*lnmusms*pow(mcMSmusmc,12))/3564.) + (5*pow(mcMSmusmc,12)*pow(Pi,2))/10692.) - (59231*pow(mcMSmusmc,12)*Zeta3)/8.2944e6)/pow(mbMSmus,11) + 
   ((396218064296685469*pow(mcMSmusmc,14))/5.896309609846656e19 - (845153*pow(lnmcmsms,2)*pow(mcMSmusmc,14))/1.2108096e7 + (3166*lnmcmusmcms*pow(mcMSmusmc,14))/511225. - 
      (362077*lnmusms*pow(mcMSmusmc,14))/1.5030015e8 + ((1229*pow(mcMSmusmc,12))/1.411344e6 - (5*lnmcmsms*pow(mcMSmusmc,12))/1782.)*pow(mufac,2) + (82285201*pow(mcMSmusmc,14)*pow(Pi,2))/3.87459072e10 - 
      (23*ln2*pow(mcMSmusmc,14)*pow(Pi,2))/35840. + lnmcmsms*((64922338969*pow(mcMSmusmc,14))/3.718698984e12 - (24*lnmcmusmcms*pow(mcMSmusmc,14))/715. + (94*lnmusms*pow(mcMSmusmc,14))/5005. - 
         (2161*pow(mcMSmusmc,14)*pow(Pi,2))/322560.) + NLMSOS*((108352091581*pow(mcMSmusmc,14))/8.1243243081e13 - (15371*lnmusms*pow(mcMSmusmc,14))/7.5150075e7 + 
         lnmcmsms*((1153*pow(mcMSmusmc,14))/2.25450225e8 + (4*lnmusms*pow(mcMSmusmc,14))/5005.) + (4*pow(mcMSmusmc,14)*pow(Pi,2))/15015.) - (3287*pow(mcMSmusmc,14)*Zeta3)/645120.)/pow(mbMSmus,13) + 
   ((516851278553981044967.*pow(mcMSmusmc,16))/9.509568138760686e22 - (232939907*pow(lnmcmsms,2)*pow(mcMSmusmc,16))/4.576860288e9 + (283*lnmcmusmcms*pow(mcMSmusmc,16))/76050. - 
      (16651*lnmusms*pow(mcMSmusmc,16))/9.7344e6 + ((30742*pow(mcMSmusmc,14))/7.5150075e7 - (8*lnmcmsms*pow(mcMSmusmc,14))/5005.)*pow(mufac,2) + 
      (58806560951*pow(mcMSmusmc,16)*pow(Pi,2))/4.326189170688e13 - (1001*ln2*pow(mcMSmusmc,16)*pow(Pi,2))/2.359296e6 + 
      lnmcmsms*((12430141121803*pow(mcMSmusmc,16))/1.1780838381312e15 - (14*lnmcmusmcms*pow(mcMSmusmc,16))/585. + (413*lnmusms*pow(mcMSmusmc,16))/28080. - 
         (565351*pow(mcMSmusmc,16)*pow(Pi,2))/1.15605504e8) + NLMSOS*((214558103603*pow(mcMSmusmc,16))/2.60460712512e14 - (529*lnmusms*pow(mcMSmusmc,16))/4.8672e6 + 
         lnmcmsms*((355*pow(mcMSmusmc,16))/1.4455584e7 + (7*lnmusms*pow(mcMSmusmc,16))/14040.) + (7*pow(mcMSmusmc,16)*pow(Pi,2))/42120.) - (885457*pow(mcMSmusmc,16)*Zeta3)/2.31211008e8)/pow(mbMSmus,15)\
    + ((218348329088480198615629.*pow(mcMSmusmc,18))/5.00576874275692e25 - (1014809351*pow(lnmcmsms,2)*pow(mcMSmusmc,18))/2.615348736e10 + (7678*lnmcmusmcms*pow(mcMSmusmc,18))/3.186225e6 - 
      (641587*lnmusms*pow(mcMSmusmc,18))/5.1616845e8 + ((529*pow(mcMSmusmc,16))/2.4336e6 - (7*lnmcmsms*pow(mcMSmusmc,16))/7020.)*pow(mufac,2) + 
      (26463251891*pow(mcMSmusmc,18)*pow(Pi,2))/2.845499424768e13 - (4147*ln2*pow(mcMSmusmc,18)*pow(Pi,2))/1.3934592e7 + 
      lnmcmsms*((22606266853306541*pow(mcMSmusmc,18))/3.2684758005112013e18 - (32*lnmcmusmcms*pow(mcMSmusmc,18))/1785. + (568*lnmusms*pow(mcMSmusmc,18))/48195. - 
         (52013*pow(mcMSmusmc,18)*pow(Pi,2))/1.3934592e7) + NLMSOS*((20555048260909*pow(mcMSmusmc,18))/3.76818092235585e16 - (16277*lnmusms*pow(mcMSmusmc,18))/2.58084225e8 + 
         lnmcmsms*((197062*pow(mcMSmusmc,18))/7.381208835e9 + (16*lnmusms*pow(mcMSmusmc,18))/48195.) + (16*pow(mcMSmusmc,18)*pow(Pi,2))/144585.) - (83291*pow(mcMSmusmc,18)*Zeta3)/2.7869184e7)/
    pow(mbMSmus,17) + ((192872454677623233846986449.*pow(mcMSmusmc,20))/5.449931397868209e28 - (46137065941*pow(lnmcmsms,2)*pow(mcMSmusmc,20))/1.5084957888e12 + 
      (5507*lnmcmusmcms*pow(mcMSmusmc,20))/3.338528e6 - (1855169*lnmusms*pow(mcMSmusmc,20))/2.0031168e9 + ((32554*pow(mcMSmusmc,18))/2.58084225e8 - (32*lnmcmsms*pow(mcMSmusmc,18))/48195.)*pow(mufac,2) + 
      (9007367733163*pow(mcMSmusmc,20)*pow(Pi,2))/1.34810154565632e16 - (143*ln2*pow(mcMSmusmc,20)*pow(Pi,2))/655360. + 
      lnmcmsms*((844305629180580989*pow(mcMSmusmc,20))/1.7558329821198565e20 - (9*lnmcmusmcms*pow(mcMSmusmc,20))/646. + (249*lnmusms*pow(mcMSmusmc,20))/25840. - 
         (156353*pow(mcMSmusmc,20)*pow(Pi,2))/5.308416e7) + NLMSOS*((22189567531163017*pow(mcMSmusmc,20))/5.834712481735737e19 - (39163*lnmusms*pow(mcMSmusmc,20))/1.0015584e9 + 
         lnmcmsms*((1211963*pow(mcMSmusmc,20))/5.012799792e10 + (3*lnmusms*pow(mcMSmusmc,20))/12920.) + (pow(mcMSmusmc,20)*pow(Pi,2))/12920.) - (254791*pow(mcMSmusmc,20)*Zeta3)/1.0616832e8)/
    pow(mbMSmus,19) + ((31442447067404839835736513193.*pow(mcMSmusmc,22))/1.0790864167779053e31 - (3078965960711*pow(lnmcmsms,2)*pow(mcMSmusmc,22))/1.24450902576e14 + 
      (5066*lnmcmusmcms*pow(mcMSmusmc,22))/4.298427e6 - (116005*lnmusms*pow(mcMSmusmc,22))/1.64245158e8 + ((39163*pow(mcMSmusmc,20))/5.007792e8 - (3*lnmcmsms*pow(mcMSmusmc,20))/6460.)*pow(mufac,2) + 
      (174200135864459*pow(mcMSmusmc,22)*pow(Pi,2))/3.495434721951744e17 - (38675*ln2*pow(mcMSmusmc,22)*pow(Pi,2))/2.33570304e8 + 
      lnmcmsms*((25315115748447270877.*pow(mcMSmusmc,22))/7.242811051244408e21 - (40*lnmcmusmcms*pow(mcMSmusmc,22))/3591. + (50*lnmusms*pow(mcMSmusmc,22))/6237. - 
         (13926181*pow(mcMSmusmc,22)*pow(Pi,2))/5.8392576e9) + NLMSOS*((1640519393726677*pow(mcMSmusmc,22))/5.946258455651274e18 - (39833*lnmusms*pow(mcMSmusmc,22))/1.560329001e9 + 
         lnmcmsms*((14290513*pow(mcMSmusmc,22))/6.89665418442e11 + (20*lnmusms*pow(mcMSmusmc,22))/118503.) + (20*pow(mcMSmusmc,22)*pow(Pi,2))/355509.) - (23017987*pow(mcMSmusmc,22)*Zeta3)/1.16785152e10)/
    pow(mbMSmus,21) + ((5769036308265946032465571119823289.*pow(mcMSmusmc,24))/2.369996943791664e36 - (953082228281*pow(lnmcmsms,2)*pow(mcMSmusmc,24))/4.664604200256e13 + 
      (10163*lnmcmusmcms*pow(mcMSmusmc,24))/1.166445e7 - (615749*lnmusms*pow(mcMSmusmc,24))/1.1197872e9 + ((79666*pow(mcMSmusmc,22))/1.560329001e9 - (40*lnmcmsms*pow(mcMSmusmc,22))/118503.)*pow(mufac,2) + 
      (240817793781176357*pow(mcMSmusmc,24)*pow(Pi,2))/6.28867544642696e20 - (877591*ln2*pow(mcMSmusmc,24)*pow(Pi,2))/6.79477248e9 + 
      lnmcmsms*((73910608414092571571.*pow(mcMSmusmc,24))/2.8097278338127475e22 - (22*lnmcmusmcms*pow(mcMSmusmc,24))/2415. + (1177*lnmusms*pow(mcMSmusmc,24))/173880. - 
         (1620816161*pow(mcMSmusmc,24)*pow(Pi,2))/8.2216747008e11) + NLMSOS*((7801530877413386647*pow(mcMSmusmc,24))/3.7763267488425775e22 - (9727*lnmusms*pow(mcMSmusmc,24))/5.598936e8 + 
         lnmcmsms*((73801799*pow(mcMSmusmc,24))/4.23178780752e12 + (11*lnmusms*pow(mcMSmusmc,24))/86940.) + (11*pow(mcMSmusmc,24)*pow(Pi,2))/260820.) - (2710689767*pow(mcMSmusmc,24)*Zeta3)/1.64433494016e12
      )/pow(mbMSmus,23) + mufac*(-808.0925925925926 - (4*pow(lnmusmcms,2))/81. + (16*lnmusmusmc)/27. + (20047*NLOSKIN)/243. - (1292*pow(NLOSKIN,2))/729. + 
      pow(lnmufmus,2)*(-53.77777777777778 + (176*NLOSKIN)/27. - (16*pow(NLOSKIN,2))/81.) + (1022*pow(Pi,2))/27. - (208*NLOSKIN*pow(Pi,2))/81. + (8*pow(NLOSKIN,2)*pow(Pi,2))/243. - (2*pow(Pi,4))/3. + 
      lnmusmcms*(11.432098765432098 - (128*NLOSKIN)/243. - (8*pow(Pi,2))/27.) + lnmufmus*
       (373.037037037037 + lnmusmcms*(-3.259259259259259 + (16*NLOSKIN)/81.) - (3356*NLOSKIN)/81. + (256*pow(NLOSKIN,2))/243. - (88*pow(Pi,2))/9. + (16*NLOSKIN*pow(Pi,2))/27.) + 114*Zeta3 - 
      (140*NLOSKIN*Zeta3)/27.) + ((7777*pow(mcMSmusmc,4))/7776. + (20*A4*pow(mcMSmusmc,4))/27. + (5*pow(ln2,4)*pow(mcMSmusmc,4))/162. + (7*pow(lnmcmsms,3)*pow(mcMSmusmc,4))/27. - 
      (2*pow(mcMSmusmc,2)*pow(mufac,2))/3. + (37*pow(mcMSmusmc,4)*pow(Pi,2))/108. + (2*pow(ln2,2)*pow(mcMSmusmc,4)*pow(Pi,2))/81. + (271*pow(mcMSmusmc,4)*pow(Pi,4))/19440. + 
      lnmcmusmcms*((-56*pow(mcMSmusmc,4))/27. - (2*pow(mcMSmusmc,4)*pow(Pi,2))/9.) + lnmusms*((-2899*pow(mcMSmusmc,4))/1296. - (13*pow(mcMSmusmc,4)*pow(Pi,2))/108.) + 
      pow(lnmcmsms,2)*((-94*pow(mcMSmusmc,4))/27. - (4*lnmcmusmcms*pow(mcMSmusmc,4))/3. - (13*lnmusms*pow(mcMSmusmc,4))/18. + (pow(mcMSmusmc,4)*pow(Pi,2))/4.) + 
      ln2*((5*pow(mcMSmusmc,4)*pow(Pi,2))/144. - (lnmcmsms*pow(mcMSmusmc,4)*pow(Pi,2))/9.) - (2309*pow(mcMSmusmc,4)*Zeta3)/864. + 
      NLMSOS*((-1423*pow(mcMSmusmc,4))/3888. - (2*pow(lnmcmsms,3)*pow(mcMSmusmc,4))/27. + pow(lnmcmsms,2)*((13*pow(mcMSmusmc,4))/54. + (lnmusms*pow(mcMSmusmc,4))/9.) - 
         (13*pow(mcMSmusmc,4)*pow(Pi,2))/324. + lnmusms*((151*pow(mcMSmusmc,4))/648. + (pow(mcMSmusmc,4)*pow(Pi,2))/54.) + 
         lnmcmsms*(-pow(mcMSmusmc,4)/12. - (13*lnmusms*pow(mcMSmusmc,4))/54. + (pow(mcMSmusmc,4)*pow(Pi,2))/27.) + (pow(mcMSmusmc,4)*Zeta3)/3.) + 
      lnmcmsms*((2443*pow(mcMSmusmc,4))/648. + (20*lnmcmusmcms*pow(mcMSmusmc,4))/9. + (241*lnmusms*pow(mcMSmusmc,4))/108. - (67*pow(mcMSmusmc,4)*pow(Pi,2))/144. + (14*pow(mcMSmusmc,4)*Zeta3)/9.))/
    pow(mbMSmus,3) + mbMSmus*(80.85713520233196 - (212*A4)/27. - (53*pow(ln2,4))/162. - (4*lnmcmusmcms)/9. + (4*lnmusmcms)/9. + (19457*pow(lnmusms,2))/864. + (1693*pow(lnmusms,3))/432. - (4*lnmusmusmc)/9. + 
      (594941*pow(Pi,2))/38880. - (20*pow(ln2,2)*pow(Pi,2))/81. - (451*pow(Pi,4))/7776. + ln2*((-199*pow(Pi,2))/54. + (37*lnmusms*pow(Pi,2))/54.) + 
      NLMSOS*(-9.736839849108367 + (8*A4)/27. + pow(ln2,4)/81. - (1103*pow(lnmusms,2))/432. - (41*pow(lnmusms,3))/108. - (313*pow(Pi,2))/216. + (2*pow(ln2,2)*pow(Pi,2))/81. + (61*pow(Pi,4))/1944. + 
         ln2*((-11*pow(Pi,2))/81. - (lnmusms*pow(Pi,2))/27.) + lnmusms*(-7.705246913580247 - (47*pow(Pi,2))/108. - (7*Zeta3)/9.) - (667*Zeta3)/216.) + 
      lnmusms*(64.21836419753086 - lnmcmusmcms/3. + lnmusmcms/3. - lnmusmusmc/3. + (185*pow(Pi,2))/108. - (97*Zeta3)/36.) + 
      pow(NLMSOS,2)*(0.1008659122085048 + (13*pow(lnmusms,2))/216. + pow(lnmusms,3)/108. + (13*pow(Pi,2))/324. + lnmusms*(0.13734567901234568 + pow(Pi,2)/54.) + (7*Zeta3)/54.) - (77*Zeta3)/72. - 
      (1439*pow(Pi,2)*Zeta3)/432. + (1975*Zeta5)/216.) + ((-68*pow(mcMSmusmc,2))/9. - 2*lnmcmusmcms*pow(mcMSmusmc,2) - (25*lnmusms*pow(mcMSmusmc,2))/6. + 
      ((2*pow(mcMSmusmc,2))/9. + (lnmusms*pow(mcMSmusmc,2))/3.)*NLMSOS - (13*pow(mcMSmusmc,2)*pow(Pi,2))/12. + 
      lnmcmsms*(-4*pow(mcMSmusmc,2) - (3*pow(mcMSmusmc,2)*pow(Pi,2))/2. + (pow(mcMSmusmc,2)*pow(Pi,4))/12.) - (11*pow(mcMSmusmc,2)*Zeta3)/2. + (3*pow(mcMSmusmc,2)*pow(Pi,2)*Zeta3)/4. + 
      pow(mufac,2)*(-204.6435185185185 - pow(lnmusmcms,2)/54. + (7*pow(lnmusms,2))/12. + (2*lnmusmusmc)/9. + (13805*NLOSKIN)/648. - (209*pow(NLOSKIN,2))/486. + 
         pow(lnmufmus,2)*(-20.166666666666668 + (22*NLOSKIN)/9. - (2*pow(NLOSKIN,2))/27.) + (1307*pow(Pi,2))/108. + (2*ln2*pow(Pi,2))/27. - (23*NLOSKIN*pow(Pi,2))/27. + 
         (pow(NLOSKIN,2)*pow(Pi,2))/81. - pow(Pi,4)/4. + lnmusms*(12.805555555555555 - lnmusmcms/9. - (13*NLOSKIN)/27. - pow(Pi,2)/3.) + lnmusmcms*(3.5277777777777777 - (13*NLOSKIN)/81. - pow(Pi,2)/9.) + 
         NLMSOS*(-0.3287037037037037 - (13*lnmusms)/54. - pow(lnmusms,2)/18. - pow(Pi,2)/27.) + 
         lnmufmus*(114.83333333333333 + lnmusmcms*(-1.2222222222222223 + (2*NLOSKIN)/27.) + lnmusms*(-3.6666666666666665 + (2*NLOSKIN)/9.) - (691*NLOSKIN)/54. + (26*pow(NLOSKIN,2))/81. - (11*pow(Pi,2))/3. + 
            (2*NLOSKIN*pow(Pi,2))/9.) + (1535*Zeta3)/36. - (35*NLOSKIN*Zeta3)/18.) + (5*pow(mcMSmusmc,2)*Zeta5)/2.)/mbMSmus);
    
    double mkin = 0.0;
    for(int i = 0; i <= nloops; i++) {
       mkin += ret[i];
    }
    return mkin;
}

double CRunDec::mMS2mkinC(double mbMSmus, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops) {
    if(nloops<0||nloops>3){
      cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
      RETURN
    }
    double ret[4];
    double lnmusms = log(mus*mus/(mbMSmus*mbMSmus));
    double lnmufmus = log(2*mufac/mus);
    double lnmusmcms = 0.0;
    double lnmcmusmcms = 0.0;
    double lnmcmsms = 0.0;
    double mcMSmusmc = 0.0;
    double musmc = 0.0;
    if(mcMSmusmcin > 0.0) {
       lnmusmcms = log(mus*mus/(mcMSmusmcin*mcMSmusmcin));
       lnmcmusmcms = log(musmcin*musmcin/(mcMSmusmcin*mcMSmusmcin));
       lnmcmsms = log(mcMSmusmcin/mbMSmus);
       mcMSmusmc = mcMSmusmcin;
       musmc = musmcin;
    }
    
    ret[0] = mbMSmus;
    ret[1] = apinlmus*((1.3333333333333333 + lnmusms)*mbMSmus - (16*mufac)/9. - (2*pow(mufac,2))/(3.*mbMSmus));
    ret[2] = apinlmus*apinlmus*((2621*mbMSmus)/288. + (161*lnmusms*mbMSmus)/24. + (15*pow(lnmusms,2)*mbMSmus)/8. - pow(mcMSmusmc,2)/mbMSmus - (151*pow(mcMSmusmc,4))/(216.*pow(mbMSmus,3)) + 
   (13*lnmcmsms*pow(mcMSmusmc,4))/(18.*pow(mbMSmus,3)) - (pow(lnmcmsms,2)*pow(mcMSmusmc,4))/(3.*pow(mbMSmus,3)) + (19*pow(mcMSmusmc,6))/(225.*pow(mbMSmus,5)) - 
   (4*lnmcmsms*pow(mcMSmusmc,6))/(45.*pow(mbMSmus,5)) + (463*pow(mcMSmusmc,8))/(39200.*pow(mbMSmus,7)) - (3*lnmcmsms*pow(mcMSmusmc,8))/(140.*pow(mbMSmus,7)) + 
   (997*pow(mcMSmusmc,10))/(297675.*pow(mbMSmus,9)) - (8*lnmcmsms*pow(mcMSmusmc,10))/(945.*pow(mbMSmus,9)) + (1229*pow(mcMSmusmc,12))/(940896.*pow(mbMSmus,11)) - 
   (5*lnmcmsms*pow(mcMSmusmc,12))/(1188.*pow(mbMSmus,11)) + (15371*pow(mcMSmusmc,14))/(2.5050025e7*pow(mbMSmus,13)) - (12*lnmcmsms*pow(mcMSmusmc,14))/(5005.*pow(mbMSmus,13)) + 
   (529*pow(mcMSmusmc,16))/(1.6224e6*pow(mbMSmus,15)) - (7*lnmcmsms*pow(mcMSmusmc,16))/(4680.*pow(mbMSmus,15)) + (16277*pow(mcMSmusmc,18))/(8.6028075e7*pow(mbMSmus,17)) - 
   (16*lnmcmsms*pow(mcMSmusmc,18))/(16065.*pow(mbMSmus,17)) + (39163*pow(mcMSmusmc,20))/(3.338528e8*pow(mbMSmus,19)) - (9*lnmcmsms*pow(mcMSmusmc,20))/(12920.*pow(mbMSmus,19)) + 
   (39833*pow(mcMSmusmc,22))/(5.20109667e8*pow(mbMSmus,21)) - (20*lnmcmsms*pow(mcMSmusmc,22))/(39501.*pow(mbMSmus,21)) + (9727*pow(mcMSmusmc,24))/(1.866312e8*pow(mbMSmus,23)) - 
   (11*lnmcmsms*pow(mcMSmusmc,24))/(28980.*pow(mbMSmus,23)) - (860*mufac)/27. + (88*lnmufmus*mufac)/9. - (83*pow(mufac,2))/(9.*mbMSmus) + (11*lnmufmus*pow(mufac,2))/(3.*mbMSmus) + 
   (2*lnmusms*pow(mufac,2))/(3.*mbMSmus) - (71*mbMSmus*NLMSOS)/144. - (13*lnmusms*mbMSmus*NLMSOS)/36. - (pow(lnmusms,2)*mbMSmus*NLMSOS)/12. + (128*mufac*NLOSKIN)/81. - (16*lnmufmus*mufac*NLOSKIN)/27. + 
   (13*pow(mufac,2)*NLOSKIN)/(27.*mbMSmus) - (2*lnmufmus*pow(mufac,2)*NLOSKIN)/(9.*mbMSmus) + (5*mbMSmus*pow(Pi,2))/18. + (ln2*mbMSmus*pow(Pi,2))/9. + (mcMSmusmc*pow(Pi,2))/6. + 
   (pow(mcMSmusmc,3)*pow(Pi,2))/(6.*pow(mbMSmus,2)) - (pow(mcMSmusmc,4)*pow(Pi,2))/(18.*pow(mbMSmus,3)) + (8*mufac*pow(Pi,2))/9. + (pow(mufac,2)*pow(Pi,2))/(3.*mbMSmus) - 
   (mbMSmus*NLMSOS*pow(Pi,2))/18. - (mbMSmus*Zeta3)/6.);
    ret[3] = apinlmus*apinlmus*apinlmus*((((9727*pow(mcMSmusmc,24))/2.799468e8 - (11*lnmcmsms*pow(mcMSmusmc,24))/43470.)*pow(mufac,2))/pow(mbMSmus,25) + (29783*mcMSmusmc*pow(Pi,2))/2430. - (1199*ln2*mcMSmusmc*pow(Pi,2))/81. - 
   (31*lnmcmsms*mcMSmusmc*pow(Pi,2))/18. + (lnmcmusmcms*mcMSmusmc*pow(Pi,2))/6. + (31*lnmusms*mcMSmusmc*pow(Pi,2))/36. + (13*mcMSmusmc*pow(Pi,3))/162. + 
   NLMSOS*((-7*mcMSmusmc*pow(Pi,2))/27. + (2*ln2*mcMSmusmc*pow(Pi,2))/9. + (lnmcmsms*mcMSmusmc*pow(Pi,2))/9. - (lnmusms*mcMSmusmc*pow(Pi,2))/18.) + 
   ((149*pow(mcMSmusmc,5)*pow(Pi,2))/1080. - (13*lnmcmsms*pow(mcMSmusmc,5)*pow(Pi,2))/45. + (pow(mcMSmusmc,3)*pow(mufac,2)*pow(Pi,2))/9. - (pow(mcMSmusmc,5)*NLMSOS*pow(Pi,2))/20. + 
      (pow(mcMSmusmc,5)*pow(Pi,3))/240.)/pow(mbMSmus,4) + ((262769*pow(mcMSmusmc,7)*pow(Pi,2))/1.1907e7 - (29*lnmcmsms*pow(mcMSmusmc,7)*pow(Pi,2))/525. - 
      (5*pow(mcMSmusmc,7)*NLMSOS*pow(Pi,2))/756. - (25*pow(mcMSmusmc,7)*pow(Pi,3))/18144.)/pow(mbMSmus,6) + 
   ((416029*pow(mcMSmusmc,9)*pow(Pi,2))/6.001128e7 - (53*lnmcmsms*pow(mcMSmusmc,9)*pow(Pi,2))/2205. - (7*pow(mcMSmusmc,9)*NLMSOS*pow(Pi,2))/3240. - (77*pow(mcMSmusmc,9)*pow(Pi,3))/103680.)/
    pow(mbMSmus,8) + ((1055689*pow(mcMSmusmc,11)*pow(Pi,2))/3.4577928e8 - (85*lnmcmsms*pow(mcMSmusmc,11)*pow(Pi,2))/6237. - (3*pow(mcMSmusmc,11)*NLMSOS*pow(Pi,2))/3080. - 
      (17*pow(mcMSmusmc,11)*pow(Pi,3))/39424.)/pow(mbMSmus,10) + ((31786481*pow(mcMSmusmc,13)*pow(Pi,2))/1.967766372e10 - (125*lnmcmsms*pow(mcMSmusmc,13)*pow(Pi,2))/14157. - 
      (11*pow(mcMSmusmc,13)*NLMSOS*pow(Pi,2))/21060. - (1771*pow(mcMSmusmc,13)*pow(Pi,3))/6.469632e6)/pow(mbMSmus,12) + 
   ((17346493*pow(mcMSmusmc,15)*pow(Pi,2))/1.808754948e10 - (173*lnmcmsms*pow(mcMSmusmc,15)*pow(Pi,2))/27885. - (13*pow(mcMSmusmc,15)*NLMSOS*pow(Pi,2))/41580. - 
      (377*pow(mcMSmusmc,15)*pow(Pi,3))/2.02752e6)/pow(mbMSmus,14) + ((22757641*pow(mcMSmusmc,17)*pow(Pi,2))/3.6923796e10 - (229*lnmcmsms*pow(mcMSmusmc,17)*pow(Pi,2))/49725. - 
      (5*pow(mcMSmusmc,17)*NLMSOS*pow(Pi,2))/24752. - (1925*pow(mcMSmusmc,17)*pow(Pi,3))/1.4483456e7)/pow(mbMSmus,16) + 
   ((86836957*pow(mcMSmusmc,19)*pow(Pi,2))/2.0687188752e11 - (293*lnmcmsms*pow(mcMSmusmc,19)*pow(Pi,2))/82365. - (17*pow(mcMSmusmc,19)*NLMSOS*pow(Pi,2))/123120. - 
      (99671*pow(mcMSmusmc,19)*pow(Pi,3))/1.00859904e9)/pow(mbMSmus,18) + ((846435761*pow(mcMSmusmc,21)*pow(Pi,2))/2.83231951884e12 - (365*lnmcmsms*pow(mcMSmusmc,21)*pow(Pi,2))/128877. - 
      (19*pow(mcMSmusmc,21)*NLMSOS*pow(Pi,2))/192780. - (127699*pow(mcMSmusmc,21)*pow(Pi,3))/1.684537344e9)/pow(mbMSmus,20) + 
   ((171475369*pow(mcMSmusmc,23)*pow(Pi,2))/7.7816811996e11 - (445*lnmcmsms*pow(mcMSmusmc,23)*pow(Pi,2))/192717. - (7*pow(mcMSmusmc,23)*NLMSOS*pow(Pi,2))/96140. - 
      (81991*pow(mcMSmusmc,23)*pow(Pi,3))/1.374683136e9)/pow(mbMSmus,22) + ((106559417*pow(mcMSmusmc,25)*pow(Pi,2))/6.374388636e11 - (533*lnmcmsms*pow(mcMSmusmc,25)*pow(Pi,2))/277725. - 
      (23*pow(mcMSmusmc,25)*NLMSOS*pow(Pi,2))/415800. - (5698043*pow(mcMSmusmc,25)*pow(Pi,3))/1.189085184e11)/pow(mbMSmus,24) + 
   ((18607*pow(mcMSmusmc,3)*pow(Pi,2))/1620. - (1199*ln2*pow(mcMSmusmc,3)*pow(Pi,2))/81. - (271*lnmcmsms*pow(mcMSmusmc,3)*pow(Pi,2))/162. + (lnmcmusmcms*pow(mcMSmusmc,3)*pow(Pi,2))/2. + 
      (19*lnmusms*pow(mcMSmusmc,3)*pow(Pi,2))/36. + (mcMSmusmc*pow(mufac,2)*pow(Pi,2))/9. + (7*pow(mcMSmusmc,3)*pow(Pi,3))/108. + 
      NLMSOS*((-7*pow(mcMSmusmc,3)*pow(Pi,2))/54. + (2*ln2*pow(mcMSmusmc,3)*pow(Pi,2))/9. + (lnmcmsms*pow(mcMSmusmc,3)*pow(Pi,2))/9. - (lnmusms*pow(mcMSmusmc,3)*pow(Pi,2))/18.))/pow(mbMSmus,2) + 
   ((-19654121*pow(mcMSmusmc,6))/1.7496e7 - (5329*pow(lnmcmsms,2)*pow(mcMSmusmc,6))/6480. + (94*lnmcmusmcms*pow(mcMSmusmc,6))/225. + (139*lnmusms*pow(mcMSmusmc,6))/1350. + 
      (311*pow(mcMSmusmc,6)*pow(Pi,2))/2592. - (11*ln2*pow(mcMSmusmc,6)*pow(Pi,2))/648. + 
      pow(mufac,2)*((-151*pow(mcMSmusmc,4))/324. + (13*lnmcmsms*pow(mcMSmusmc,4))/27. - (2*pow(lnmcmsms,2)*pow(mcMSmusmc,4))/9. - (pow(mcMSmusmc,4)*pow(Pi,2))/27.) + 
      lnmcmsms*((321193*pow(mcMSmusmc,6))/233280. - (8*lnmcmusmcms*pow(mcMSmusmc,6))/15. - (2*lnmusms*pow(mcMSmusmc,6))/135. - (113*pow(mcMSmusmc,6)*pow(Pi,2))/1296.) + 
      NLMSOS*((2729*pow(mcMSmusmc,6))/30375. - (19*lnmusms*pow(mcMSmusmc,6))/675. + lnmcmsms*((-2*pow(mcMSmusmc,6))/75. + (4*lnmusms*pow(mcMSmusmc,6))/135.) + (4*pow(mcMSmusmc,6)*pow(Pi,2))/405.) - 
      (29*pow(mcMSmusmc,6)*Zeta3)/648.)/pow(mbMSmus,5) + ((-74049121*pow(mcMSmusmc,8))/1.185408e9 - (9451*pow(lnmcmsms,2)*pow(mcMSmusmc,8))/30240. + (179*lnmcmusmcms*pow(mcMSmusmc,8))/2450. - 
      (53*lnmusms*pow(mcMSmusmc,8))/235200. + ((38*pow(mcMSmusmc,6))/675. - (8*lnmcmsms*pow(mcMSmusmc,6))/135.)*pow(mufac,2) + (1789*pow(mcMSmusmc,8)*pow(Pi,2))/90720. - 
      (7*ln2*pow(mcMSmusmc,8)*pow(Pi,2))/1536. + lnmcmsms*((2975311*pow(mcMSmusmc,8))/1.27008e7 - (6*lnmcmusmcms*pow(mcMSmusmc,8))/35. + (11*lnmusms*pow(mcMSmusmc,8))/280. - 
         (139*pow(mcMSmusmc,8)*pow(Pi,2))/4608.) + NLMSOS*((174787*pow(mcMSmusmc,8))/1.2348e7 - (463*lnmusms*pow(mcMSmusmc,8))/117600. + 
         lnmcmsms*((-37*pow(mcMSmusmc,8))/14700. + (lnmusms*pow(mcMSmusmc,8))/140.) + (pow(mcMSmusmc,8)*pow(Pi,2))/420.) - (173*pow(mcMSmusmc,8)*Zeta3)/9216.)/pow(mbMSmus,7) + 
   ((13434791647*pow(mcMSmusmc,10))/5.184974592e12 - (382589*pow(lnmcmsms,2)*pow(mcMSmusmc,10))/2.3328e6 + (298*lnmcmusmcms*pow(mcMSmusmc,10))/11907. - (7811*lnmusms*pow(mcMSmusmc,10))/1.78605e6 + 
      ((463*pow(mcMSmusmc,8))/58800. - (lnmcmsms*pow(mcMSmusmc,8))/70.)*pow(mufac,2) + (359801*pow(mcMSmusmc,10)*pow(Pi,2))/4.89888e7 - (17*ln2*pow(mcMSmusmc,10)*pow(Pi,2))/8640. + 
      lnmcmsms*((3017558449*pow(mcMSmusmc,10))/4.1150592e10 - (16*lnmcmusmcms*pow(mcMSmusmc,10))/189. + (92*lnmusms*pow(mcMSmusmc,10))/2835. - (17*pow(mcMSmusmc,10)*pow(Pi,2))/1080.) + 
      NLMSOS*((2816347*pow(mcMSmusmc,10))/5.6260575e8 - (997*lnmusms*pow(mcMSmusmc,10))/893025. + lnmcmsms*((-398*pow(mcMSmusmc,10))/893025. + (8*lnmusms*pow(mcMSmusmc,10))/2835.) + 
         (8*pow(mcMSmusmc,10)*pow(Pi,2))/8505.) - (187*pow(mcMSmusmc,10)*Zeta3)/17280.)/pow(mbMSmus,9) + 
   ((185657253148457*pow(mcMSmusmc,12))/2.46471470784e16 - (15261023*pow(lnmcmsms,2)*pow(mcMSmusmc,12))/1.49688e8 + (899*lnmcmusmcms*pow(mcMSmusmc,12))/78408. - 
      (19255*lnmusms*pow(mcMSmusmc,12))/5.645376e6 + ((1994*pow(mcMSmusmc,10))/893025. - (16*lnmcmsms*pow(mcMSmusmc,10))/2835.)*pow(mufac,2) + (84041429*pow(mcMSmusmc,12)*pow(Pi,2))/2.29920768e10 - 
      (175*ln2*pow(mcMSmusmc,12)*pow(Pi,2))/165888. + lnmcmsms*((202523874613*pow(mcMSmusmc,12))/6.22402704e12 - (5*lnmcmusmcms*pow(mcMSmusmc,12))/99. + (175*lnmusms*pow(mcMSmusmc,12))/7128. - 
         (40553*pow(mcMSmusmc,12)*pow(Pi,2))/4.1472e6) + NLMSOS*((326802499*pow(mcMSmusmc,12))/1.3692859488e11 - (1229*lnmusms*pow(mcMSmusmc,12))/2.822688e6 + 
         lnmcmsms*((-391*pow(mcMSmusmc,12))/4.939704e6 + (5*lnmusms*pow(mcMSmusmc,12))/3564.) + (5*pow(mcMSmusmc,12)*pow(Pi,2))/10692.) - (59231*pow(mcMSmusmc,12)*Zeta3)/8.2944e6)/pow(mbMSmus,11) + 
   ((396218064296685469*pow(mcMSmusmc,14))/5.896309609846656e19 - (845153*pow(lnmcmsms,2)*pow(mcMSmusmc,14))/1.2108096e7 + (3166*lnmcmusmcms*pow(mcMSmusmc,14))/511225. - 
      (362077*lnmusms*pow(mcMSmusmc,14))/1.5030015e8 + ((1229*pow(mcMSmusmc,12))/1.411344e6 - (5*lnmcmsms*pow(mcMSmusmc,12))/1782.)*pow(mufac,2) + (82285201*pow(mcMSmusmc,14)*pow(Pi,2))/3.87459072e10 - 
      (23*ln2*pow(mcMSmusmc,14)*pow(Pi,2))/35840. + lnmcmsms*((64922338969*pow(mcMSmusmc,14))/3.718698984e12 - (24*lnmcmusmcms*pow(mcMSmusmc,14))/715. + (94*lnmusms*pow(mcMSmusmc,14))/5005. - 
         (2161*pow(mcMSmusmc,14)*pow(Pi,2))/322560.) + NLMSOS*((108352091581*pow(mcMSmusmc,14))/8.1243243081e13 - (15371*lnmusms*pow(mcMSmusmc,14))/7.5150075e7 + 
         lnmcmsms*((1153*pow(mcMSmusmc,14))/2.25450225e8 + (4*lnmusms*pow(mcMSmusmc,14))/5005.) + (4*pow(mcMSmusmc,14)*pow(Pi,2))/15015.) - (3287*pow(mcMSmusmc,14)*Zeta3)/645120.)/pow(mbMSmus,13) + 
   ((516851278553981044967.*pow(mcMSmusmc,16))/9.509568138760686e22 - (232939907*pow(lnmcmsms,2)*pow(mcMSmusmc,16))/4.576860288e9 + (283*lnmcmusmcms*pow(mcMSmusmc,16))/76050. - 
      (16651*lnmusms*pow(mcMSmusmc,16))/9.7344e6 + ((30742*pow(mcMSmusmc,14))/7.5150075e7 - (8*lnmcmsms*pow(mcMSmusmc,14))/5005.)*pow(mufac,2) + 
      (58806560951*pow(mcMSmusmc,16)*pow(Pi,2))/4.326189170688e13 - (1001*ln2*pow(mcMSmusmc,16)*pow(Pi,2))/2.359296e6 + 
      lnmcmsms*((12430141121803*pow(mcMSmusmc,16))/1.1780838381312e15 - (14*lnmcmusmcms*pow(mcMSmusmc,16))/585. + (413*lnmusms*pow(mcMSmusmc,16))/28080. - 
         (565351*pow(mcMSmusmc,16)*pow(Pi,2))/1.15605504e8) + NLMSOS*((214558103603*pow(mcMSmusmc,16))/2.60460712512e14 - (529*lnmusms*pow(mcMSmusmc,16))/4.8672e6 + 
         lnmcmsms*((355*pow(mcMSmusmc,16))/1.4455584e7 + (7*lnmusms*pow(mcMSmusmc,16))/14040.) + (7*pow(mcMSmusmc,16)*pow(Pi,2))/42120.) - (885457*pow(mcMSmusmc,16)*Zeta3)/2.31211008e8)/pow(mbMSmus,15)\
    + ((218348329088480198615629.*pow(mcMSmusmc,18))/5.00576874275692e25 - (1014809351*pow(lnmcmsms,2)*pow(mcMSmusmc,18))/2.615348736e10 + (7678*lnmcmusmcms*pow(mcMSmusmc,18))/3.186225e6 - 
      (641587*lnmusms*pow(mcMSmusmc,18))/5.1616845e8 + ((529*pow(mcMSmusmc,16))/2.4336e6 - (7*lnmcmsms*pow(mcMSmusmc,16))/7020.)*pow(mufac,2) + 
      (26463251891*pow(mcMSmusmc,18)*pow(Pi,2))/2.845499424768e13 - (4147*ln2*pow(mcMSmusmc,18)*pow(Pi,2))/1.3934592e7 + 
      lnmcmsms*((22606266853306541*pow(mcMSmusmc,18))/3.2684758005112013e18 - (32*lnmcmusmcms*pow(mcMSmusmc,18))/1785. + (568*lnmusms*pow(mcMSmusmc,18))/48195. - 
         (52013*pow(mcMSmusmc,18)*pow(Pi,2))/1.3934592e7) + NLMSOS*((20555048260909*pow(mcMSmusmc,18))/3.76818092235585e16 - (16277*lnmusms*pow(mcMSmusmc,18))/2.58084225e8 + 
         lnmcmsms*((197062*pow(mcMSmusmc,18))/7.381208835e9 + (16*lnmusms*pow(mcMSmusmc,18))/48195.) + (16*pow(mcMSmusmc,18)*pow(Pi,2))/144585.) - (83291*pow(mcMSmusmc,18)*Zeta3)/2.7869184e7)/
    pow(mbMSmus,17) + ((192872454677623233846986449.*pow(mcMSmusmc,20))/5.449931397868209e28 - (46137065941*pow(lnmcmsms,2)*pow(mcMSmusmc,20))/1.5084957888e12 + 
      (5507*lnmcmusmcms*pow(mcMSmusmc,20))/3.338528e6 - (1855169*lnmusms*pow(mcMSmusmc,20))/2.0031168e9 + ((32554*pow(mcMSmusmc,18))/2.58084225e8 - (32*lnmcmsms*pow(mcMSmusmc,18))/48195.)*pow(mufac,2) + 
      (9007367733163*pow(mcMSmusmc,20)*pow(Pi,2))/1.34810154565632e16 - (143*ln2*pow(mcMSmusmc,20)*pow(Pi,2))/655360. + 
      lnmcmsms*((844305629180580989*pow(mcMSmusmc,20))/1.7558329821198565e20 - (9*lnmcmusmcms*pow(mcMSmusmc,20))/646. + (249*lnmusms*pow(mcMSmusmc,20))/25840. - 
         (156353*pow(mcMSmusmc,20)*pow(Pi,2))/5.308416e7) + NLMSOS*((22189567531163017*pow(mcMSmusmc,20))/5.834712481735737e19 - (39163*lnmusms*pow(mcMSmusmc,20))/1.0015584e9 + 
         lnmcmsms*((1211963*pow(mcMSmusmc,20))/5.012799792e10 + (3*lnmusms*pow(mcMSmusmc,20))/12920.) + (pow(mcMSmusmc,20)*pow(Pi,2))/12920.) - (254791*pow(mcMSmusmc,20)*Zeta3)/1.0616832e8)/
    pow(mbMSmus,19) + ((31442447067404839835736513193.*pow(mcMSmusmc,22))/1.0790864167779053e31 - (3078965960711*pow(lnmcmsms,2)*pow(mcMSmusmc,22))/1.24450902576e14 + 
      (5066*lnmcmusmcms*pow(mcMSmusmc,22))/4.298427e6 - (116005*lnmusms*pow(mcMSmusmc,22))/1.64245158e8 + ((39163*pow(mcMSmusmc,20))/5.007792e8 - (3*lnmcmsms*pow(mcMSmusmc,20))/6460.)*pow(mufac,2) + 
      (174200135864459*pow(mcMSmusmc,22)*pow(Pi,2))/3.495434721951744e17 - (38675*ln2*pow(mcMSmusmc,22)*pow(Pi,2))/2.33570304e8 + 
      lnmcmsms*((25315115748447270877.*pow(mcMSmusmc,22))/7.242811051244408e21 - (40*lnmcmusmcms*pow(mcMSmusmc,22))/3591. + (50*lnmusms*pow(mcMSmusmc,22))/6237. - 
         (13926181*pow(mcMSmusmc,22)*pow(Pi,2))/5.8392576e9) + NLMSOS*((1640519393726677*pow(mcMSmusmc,22))/5.946258455651274e18 - (39833*lnmusms*pow(mcMSmusmc,22))/1.560329001e9 + 
         lnmcmsms*((14290513*pow(mcMSmusmc,22))/6.89665418442e11 + (20*lnmusms*pow(mcMSmusmc,22))/118503.) + (20*pow(mcMSmusmc,22)*pow(Pi,2))/355509.) - (23017987*pow(mcMSmusmc,22)*Zeta3)/1.16785152e10)/
    pow(mbMSmus,21) + ((5769036308265946032465571119823289.*pow(mcMSmusmc,24))/2.369996943791664e36 - (953082228281*pow(lnmcmsms,2)*pow(mcMSmusmc,24))/4.664604200256e13 + 
      (10163*lnmcmusmcms*pow(mcMSmusmc,24))/1.166445e7 - (615749*lnmusms*pow(mcMSmusmc,24))/1.1197872e9 + ((79666*pow(mcMSmusmc,22))/1.560329001e9 - (40*lnmcmsms*pow(mcMSmusmc,22))/118503.)*pow(mufac,2) + 
      (240817793781176357*pow(mcMSmusmc,24)*pow(Pi,2))/6.28867544642696e20 - (877591*ln2*pow(mcMSmusmc,24)*pow(Pi,2))/6.79477248e9 + 
      lnmcmsms*((73910608414092571571.*pow(mcMSmusmc,24))/2.8097278338127475e22 - (22*lnmcmusmcms*pow(mcMSmusmc,24))/2415. + (1177*lnmusms*pow(mcMSmusmc,24))/173880. - 
         (1620816161*pow(mcMSmusmc,24)*pow(Pi,2))/8.2216747008e11) + NLMSOS*((7801530877413386647*pow(mcMSmusmc,24))/3.7763267488425775e22 - (9727*lnmusms*pow(mcMSmusmc,24))/5.598936e8 + 
         lnmcmsms*((73801799*pow(mcMSmusmc,24))/4.23178780752e12 + (11*lnmusms*pow(mcMSmusmc,24))/86940.) + (11*pow(mcMSmusmc,24)*pow(Pi,2))/260820.) - (2710689767*pow(mcMSmusmc,24)*Zeta3)/1.64433494016e12
      )/pow(mbMSmus,23) + mufac*(-807.820987654321 + (20047*NLOSKIN)/243. - (1292*pow(NLOSKIN,2))/729. + pow(lnmufmus,2)*(-53.77777777777778 + (176*NLOSKIN)/27. - (16*pow(NLOSKIN,2))/81.) + 
      (1022*pow(Pi,2))/27. - (208*NLOSKIN*pow(Pi,2))/81. + (8*pow(NLOSKIN,2)*pow(Pi,2))/243. - (2*pow(Pi,4))/3. + 
      lnmufmus*(373.037037037037 - (3356*NLOSKIN)/81. + (256*pow(NLOSKIN,2))/243. - (88*pow(Pi,2))/9. + (16*NLOSKIN*pow(Pi,2))/27.) + 114*Zeta3 - (140*NLOSKIN*Zeta3)/27.) + 
   ((7777*pow(mcMSmusmc,4))/7776. + (20*A4*pow(mcMSmusmc,4))/27. + (5*pow(ln2,4)*pow(mcMSmusmc,4))/162. + (7*pow(lnmcmsms,3)*pow(mcMSmusmc,4))/27. - (2*pow(mcMSmusmc,2)*pow(mufac,2))/3. + 
      (37*pow(mcMSmusmc,4)*pow(Pi,2))/108. + (2*pow(ln2,2)*pow(mcMSmusmc,4)*pow(Pi,2))/81. + (271*pow(mcMSmusmc,4)*pow(Pi,4))/19440. + 
      lnmcmusmcms*((-56*pow(mcMSmusmc,4))/27. - (2*pow(mcMSmusmc,4)*pow(Pi,2))/9.) + lnmusms*((-2899*pow(mcMSmusmc,4))/1296. - (13*pow(mcMSmusmc,4)*pow(Pi,2))/108.) + 
      pow(lnmcmsms,2)*((-94*pow(mcMSmusmc,4))/27. - (4*lnmcmusmcms*pow(mcMSmusmc,4))/3. - (13*lnmusms*pow(mcMSmusmc,4))/18. + (pow(mcMSmusmc,4)*pow(Pi,2))/4.) + 
      ln2*((5*pow(mcMSmusmc,4)*pow(Pi,2))/144. - (lnmcmsms*pow(mcMSmusmc,4)*pow(Pi,2))/9.) - (2309*pow(mcMSmusmc,4)*Zeta3)/864. + 
      NLMSOS*((-1423*pow(mcMSmusmc,4))/3888. - (2*pow(lnmcmsms,3)*pow(mcMSmusmc,4))/27. + pow(lnmcmsms,2)*((13*pow(mcMSmusmc,4))/54. + (lnmusms*pow(mcMSmusmc,4))/9.) - 
         (13*pow(mcMSmusmc,4)*pow(Pi,2))/324. + lnmusms*((151*pow(mcMSmusmc,4))/648. + (pow(mcMSmusmc,4)*pow(Pi,2))/54.) + 
         lnmcmsms*(-pow(mcMSmusmc,4)/12. - (13*lnmusms*pow(mcMSmusmc,4))/54. + (pow(mcMSmusmc,4)*pow(Pi,2))/27.) + (pow(mcMSmusmc,4)*Zeta3)/3.) + 
      lnmcmsms*((2443*pow(mcMSmusmc,4))/648. + (20*lnmcmusmcms*pow(mcMSmusmc,4))/9. + (241*lnmusms*pow(mcMSmusmc,4))/108. - (67*pow(mcMSmusmc,4)*pow(Pi,2))/144. + (14*pow(mcMSmusmc,4)*Zeta3)/9.))/
    pow(mbMSmus,3) + mbMSmus*(80.85713520233196 - (212*A4)/27. - (53*pow(ln2,4))/162. + (19457*pow(lnmusms,2))/864. + (1693*pow(lnmusms,3))/432. + (594941*pow(Pi,2))/38880. - 
      (20*pow(ln2,2)*pow(Pi,2))/81. - (451*pow(Pi,4))/7776. + ln2*((-199*pow(Pi,2))/54. + (37*lnmusms*pow(Pi,2))/54.) + 
      NLMSOS*(-9.736839849108367 + (8*A4)/27. + pow(ln2,4)/81. - (1103*pow(lnmusms,2))/432. - (41*pow(lnmusms,3))/108. - (313*pow(Pi,2))/216. + (2*pow(ln2,2)*pow(Pi,2))/81. + (61*pow(Pi,4))/1944. + 
         ln2*((-11*pow(Pi,2))/81. - (lnmusms*pow(Pi,2))/27.) + lnmusms*(-7.705246913580247 - (47*pow(Pi,2))/108. - (7*Zeta3)/9.) - (667*Zeta3)/216.) + 
      lnmusms*(64.21836419753086 + (185*pow(Pi,2))/108. - (97*Zeta3)/36.) + pow(NLMSOS,2)*
       (0.1008659122085048 + (13*pow(lnmusms,2))/216. + pow(lnmusms,3)/108. + (13*pow(Pi,2))/324. + lnmusms*(0.13734567901234568 + pow(Pi,2)/54.) + (7*Zeta3)/54.) - (77*Zeta3)/72. - 
      (1439*pow(Pi,2)*Zeta3)/432. + (1975*Zeta5)/216.) + ((-68*pow(mcMSmusmc,2))/9. - 2*lnmcmusmcms*pow(mcMSmusmc,2) - (25*lnmusms*pow(mcMSmusmc,2))/6. + 
      ((2*pow(mcMSmusmc,2))/9. + (lnmusms*pow(mcMSmusmc,2))/3.)*NLMSOS - (13*pow(mcMSmusmc,2)*pow(Pi,2))/12. + 
      lnmcmsms*(-4*pow(mcMSmusmc,2) - (3*pow(mcMSmusmc,2)*pow(Pi,2))/2. + (pow(mcMSmusmc,2)*pow(Pi,4))/12.) - (11*pow(mcMSmusmc,2)*Zeta3)/2. + (3*pow(mcMSmusmc,2)*pow(Pi,2)*Zeta3)/4. + 
      pow(mufac,2)*(-204.54166666666666 + (7*pow(lnmusms,2))/12. + (13805*NLOSKIN)/648. - (209*pow(NLOSKIN,2))/486. + pow(lnmufmus,2)*(-20.166666666666668 + (22*NLOSKIN)/9. - (2*pow(NLOSKIN,2))/27.) + 
         (1307*pow(Pi,2))/108. + (2*ln2*pow(Pi,2))/27. - (23*NLOSKIN*pow(Pi,2))/27. + (pow(NLOSKIN,2)*pow(Pi,2))/81. - pow(Pi,4)/4. + lnmusms*(12.805555555555555 - (13*NLOSKIN)/27. - pow(Pi,2)/3.) + 
         NLMSOS*(-0.3287037037037037 - (13*lnmusms)/54. - pow(lnmusms,2)/18. - pow(Pi,2)/27.) + 
         lnmufmus*(114.83333333333333 + lnmusms*(-3.6666666666666665 + (2*NLOSKIN)/9.) - (691*NLOSKIN)/54. + (26*pow(NLOSKIN,2))/81. - (11*pow(Pi,2))/3. + (2*NLOSKIN*pow(Pi,2))/9.) + (1535*Zeta3)/36. - 
         (35*NLOSKIN*Zeta3)/18.) + (5*pow(mcMSmusmc,2)*Zeta5)/2.)/mbMSmus);
    
    double mkin = 0.0;
    for(int i = 0; i <= nloops; i++) {
       mkin += ret[i];
    }
    return mkin;
}

double CRunDec::mMS2mkinD(double mbMSmus, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops) {
    if(nloops<0||nloops>3){
      cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
      RETURN
    }
    double ret[4];
    double lnmusms = log(mus*mus/(mbMSmus*mbMSmus));
    double lnmufmus = log(2*mufac/mus);
    double lnmusmcms = 0.0;
    double lnmcmusmcms = 0.0;
    double lnmcmsms = 0.0;
    double mcMSmusmc = 0.0;
    double musmc = 0.0;
    if(mcMSmusmcin > 0.0) {
       lnmusmcms = log(mus*mus/(mcMSmusmcin*mcMSmusmcin));
       lnmcmusmcms = log(musmcin*musmcin/(mcMSmusmcin*mcMSmusmcin));
       lnmcmsms = log(mcMSmusmcin/mbMSmus);
       mcMSmusmc = mcMSmusmcin;
       musmc = musmcin;
    }
    
    ret[0] = mbMSmus;
    ret[1] = apinlmus*((1.3333333333333333 + lnmusms)*mbMSmus - (16*mufac)/9. - (2*pow(mufac,2))/(3.*mbMSmus));
    ret[2] = apinlmus*apinlmus*((pow(mufac,2)*(-9.222222222222221 + (2*lnmusms)/3. + lnmufmus*(3.6666666666666665 - (2*NLOSKIN)/9.) + (13*NLOSKIN)/27. + pow(Pi,2)/3.))/mbMSmus + 
   mufac*(-31.85185185185185 + lnmufmus*(9.777777777777779 - (16*NLOSKIN)/27.) + (128*NLOSKIN)/81. + (8*pow(Pi,2))/9.) + 
   mbMSmus*(9.59375 + (509*lnmusms)/72. + (47*pow(lnmusms,2))/24. + pow(Pi,2)/3. + (ln2*pow(Pi,2))/9. + NLMSOS*(-0.4930555555555556 - (13*lnmusms)/36. - pow(lnmusms,2)/12. - pow(Pi,2)/18.) - Zeta3/6.));
    ret[3] = apinlmus*apinlmus*apinlmus*(mufac*(-807.820987654321 + (20047*NLOSKIN)/243. - (1292*pow(NLOSKIN,2))/729. + pow(lnmufmus,2)*(-53.77777777777778 + (176*NLOSKIN)/27. - (16*pow(NLOSKIN,2))/81.) + (1022*pow(Pi,2))/27. - 
      (208*NLOSKIN*pow(Pi,2))/81. + (8*pow(NLOSKIN,2)*pow(Pi,2))/243. - (2*pow(Pi,4))/3. + 
      lnmufmus*(373.037037037037 - (3356*NLOSKIN)/81. + (256*pow(NLOSKIN,2))/243. - (88*pow(Pi,2))/9. + (16*NLOSKIN*pow(Pi,2))/27.) + 114*Zeta3 - (140*NLOSKIN*Zeta3)/27.) + 
   (pow(mufac,2)*(-204.21296296296296 + (23*pow(lnmusms,2))/36. + (13805*NLOSKIN)/648. - (209*pow(NLOSKIN,2))/486. + pow(lnmufmus,2)*(-20.166666666666668 + (22*NLOSKIN)/9. - (2*pow(NLOSKIN,2))/27.) + 
        (437*pow(Pi,2))/36. + (2*ln2*pow(Pi,2))/27. - (23*NLOSKIN*pow(Pi,2))/27. + (pow(NLOSKIN,2)*pow(Pi,2))/81. - pow(Pi,4)/4. + lnmusms*(13.046296296296296 - (13*NLOSKIN)/27. - pow(Pi,2)/3.) + 
        NLMSOS*(-0.3287037037037037 - (13*lnmusms)/54. - pow(lnmusms,2)/18. - pow(Pi,2)/27.) + 
        lnmufmus*(114.83333333333333 + lnmusms*(-3.6666666666666665 + (2*NLOSKIN)/9.) - (691*NLOSKIN)/54. + (26*pow(NLOSKIN,2))/81. - (11*pow(Pi,2))/3. + (2*NLOSKIN*pow(Pi,2))/9.) + (1535*Zeta3)/36. - 
        (35*NLOSKIN*Zeta3)/18.))/mbMSmus + mbMSmus*(90.69484096364883 - (220*A4)/27. - (55*pow(ln2,4))/162. + (21715*pow(lnmusms,2))/864. + (1861*pow(lnmusms,3))/432. + (652841*pow(Pi,2))/38880. - 
      (22*pow(ln2,2)*pow(Pi,2))/81. - (695*pow(Pi,4))/7776. + ln2*((-575*pow(Pi,2))/162. + (13*lnmusms*pow(Pi,2))/18.) + 
      NLMSOS*(-9.938571673525377 + (8*A4)/27. + pow(ln2,4)/81. - (385*pow(lnmusms,2))/144. - (43*pow(lnmusms,3))/108. - (991*pow(Pi,2))/648. + (2*pow(ln2,2)*pow(Pi,2))/81. + (61*pow(Pi,4))/1944. + 
         ln2*((-11*pow(Pi,2))/81. - (lnmusms*pow(Pi,2))/27.) + lnmusms*(-7.979938271604938 - (17*pow(Pi,2))/36. - (7*Zeta3)/9.) - (241*Zeta3)/72.) + 
      lnmusms*(72.06095679012346 + (13*pow(Pi,2))/6. - (23*Zeta3)/12.) + pow(NLMSOS,2)*
       (0.1008659122085048 + (13*pow(lnmusms,2))/216. + pow(lnmusms,3)/108. + (13*pow(Pi,2))/324. + lnmusms*(0.13734567901234568 + pow(Pi,2)/54.) + (7*Zeta3)/54.) + (58*Zeta3)/27. - 
      (1439*pow(Pi,2)*Zeta3)/432. + (1975*Zeta5)/216.));
    
    double mkin = 0.0;
    for(int i = 0; i <= nloops; i++) {
       mkin += ret[i];
    }
    return mkin;
}

double CRunDec::mkin2mMSA(double mkin, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops) {
    if(nloops<0||nloops>3){
      cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
      RETURN
    }
    double ret[5];
    double lnmusmkin = log(mus*mus/(mkin*mkin));
    double lnmufmus = log(2*mufac/mus);
    double lnmusmcms = 0.0;
    double lnmcmsmkin = 0.0;
    double lnmcmusmcms = 0.0;  
    double mcMSmusmc = 0.0;
    double musmc = 0.0;
    if(mcMSmusmcin > 0.0) {
       lnmusmcms = log(mus*mus/(mcMSmusmcin*mcMSmusmcin));
       lnmcmsmkin = log(mcMSmusmcin/mkin);
       lnmcmusmcms = log(musmcin*musmcin/(mcMSmusmcin*mcMSmusmcin));
       mcMSmusmc = mcMSmusmcin;
       musmc = musmcin;
    }
    
    ret[0] = mkin;
    ret[1] = apinlmus*((-4*mkin)/3. - lnmusmkin*mkin + (16*mufac)/9. + (2*pow(mufac,2))/(3.*mkin));
    ret[2] = apinlmus*apinlmus*((-9727*pow(mcMSmusmc,24))/(1.866312e8*pow(mkin,23)) + (11*lnmcmsmkin*pow(mcMSmusmc,24))/(28980.*pow(mkin,23)) - (39833*pow(mcMSmusmc,22))/(5.20109667e8*pow(mkin,21)) + 
   (20*lnmcmsmkin*pow(mcMSmusmc,22))/(39501.*pow(mkin,21)) - (39163*pow(mcMSmusmc,20))/(3.338528e8*pow(mkin,19)) + (9*lnmcmsmkin*pow(mcMSmusmc,20))/(12920.*pow(mkin,19)) - 
   (16277*pow(mcMSmusmc,18))/(8.6028075e7*pow(mkin,17)) + (16*lnmcmsmkin*pow(mcMSmusmc,18))/(16065.*pow(mkin,17)) - (529*pow(mcMSmusmc,16))/(1.6224e6*pow(mkin,15)) + 
   (7*lnmcmsmkin*pow(mcMSmusmc,16))/(4680.*pow(mkin,15)) - (15371*pow(mcMSmusmc,14))/(2.5050025e7*pow(mkin,13)) + (12*lnmcmsmkin*pow(mcMSmusmc,14))/(5005.*pow(mkin,13)) - 
   (1229*pow(mcMSmusmc,12))/(940896.*pow(mkin,11)) + (5*lnmcmsmkin*pow(mcMSmusmc,12))/(1188.*pow(mkin,11)) - (997*pow(mcMSmusmc,10))/(297675.*pow(mkin,9)) + 
   (8*lnmcmsmkin*pow(mcMSmusmc,10))/(945.*pow(mkin,9)) - (463*pow(mcMSmusmc,8))/(39200.*pow(mkin,7)) + (3*lnmcmsmkin*pow(mcMSmusmc,8))/(140.*pow(mkin,7)) - (19*pow(mcMSmusmc,6))/(225.*pow(mkin,5)) + 
   (4*lnmcmsmkin*pow(mcMSmusmc,6))/(45.*pow(mkin,5)) + (151*pow(mcMSmusmc,4))/(216.*pow(mkin,3)) - (13*lnmcmsmkin*pow(mcMSmusmc,4))/(18.*pow(mkin,3)) + 
   (pow(lnmcmsmkin,2)*pow(mcMSmusmc,4))/(3.*pow(mkin,3)) + pow(mcMSmusmc,2)/mkin - (959*mkin)/96. - (2*lnmusmcms*mkin)/9. - (145*lnmusmkin*mkin)/24. - (lnmusmcms*lnmusmkin*mkin)/6. - 
   (7*pow(lnmusmkin,2)*mkin)/8. + (892*mufac)/27. - (88*lnmufmus*mufac)/9. - (16*lnmusmkin*mufac)/9. + (95*pow(mufac,2))/(9.*mkin) - (11*lnmufmus*pow(mufac,2))/(3.*mkin) - 
   (2*lnmusmkin*pow(mufac,2))/(3.*mkin) + (71*mkin*NLMSOS)/144. + (13*lnmusmkin*mkin*NLMSOS)/36. + (pow(lnmusmkin,2)*mkin*NLMSOS)/12. - (128*mufac*NLOSKIN)/81. + (16*lnmufmus*mufac*NLOSKIN)/27. - 
   (13*pow(mufac,2)*NLOSKIN)/(27.*mkin) + (2*lnmufmus*pow(mufac,2)*NLOSKIN)/(9.*mkin) - (mcMSmusmc*pow(Pi,2))/6. + (pow(mcMSmusmc,4)*pow(Pi,2))/(18.*pow(mkin,3)) - 
   (pow(mcMSmusmc,3)*pow(Pi,2))/(6.*pow(mkin,2)) - (5*mkin*pow(Pi,2))/18. - (ln2*mkin*pow(Pi,2))/9. - (8*mufac*pow(Pi,2))/9. - (pow(mufac,2)*pow(Pi,2))/(3.*mkin) + (mkin*NLMSOS*pow(Pi,2))/18. + 
   (mkin*Zeta3)/6.);
    ret[3] = apinlmus*apinlmus*apinlmus*((-8439929370090150606400742890600889.*pow(mcMSmusmc,24))/(2.369996943791664e36*pow(mkin,23)) + (260257620268030848013.*lnmcmsmkin*pow(mcMSmusmc,24))/(2.8097278338127475e22*pow(mkin,23)) + 
   (953082228281*pow(lnmcmsmkin,2)*pow(mcMSmusmc,24))/(4.664604200256e13*pow(mkin,23)) - (10163*lnmcmusmcms*pow(mcMSmusmc,24))/(1.166445e7*pow(mkin,23)) + 
   (22*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,24))/(2415.*pow(mkin,23)) - (9727*lnmusmcms*pow(mcMSmusmc,24))/(5.598936e8*pow(mkin,23)) + (11*lnmcmsmkin*lnmusmcms*pow(mcMSmusmc,24))/(86940.*pow(mkin,23)) - 
   (9727*lnmusmkin*pow(mcMSmusmc,24))/(4.4791488e7*pow(mkin,23)) + (55*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,24))/(34776.*pow(mkin,23)) - 
   (47848557349149553062713793193.*pow(mcMSmusmc,22))/(1.0790864167779053e31*pow(mkin,21)) + (79809944445306329123.*lnmcmsmkin*pow(mcMSmusmc,22))/(7.242811051244408e21*pow(mkin,21)) + 
   (3078965960711*pow(lnmcmsmkin,2)*pow(mcMSmusmc,22))/(1.24450902576e14*pow(mkin,21)) - (5066*lnmcmusmcms*pow(mcMSmusmc,22))/(4.298427e6*pow(mkin,21)) + 
   (40*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,22))/(3591.*pow(mkin,21)) - (39833*lnmusmcms*pow(mcMSmusmc,22))/(1.560329001e9*pow(mkin,21)) + 
   (20*lnmcmsmkin*lnmusmcms*pow(mcMSmusmc,22))/(118503.*pow(mkin,21)) - (995825*lnmusmkin*pow(mcMSmusmc,22))/(3.120658002e9*pow(mkin,21)) + 
   (250*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,22))/(118503.*pow(mkin,21)) - (308474675741863370556669649.*pow(mcMSmusmc,20))/(5.449931397868209e28*pow(mkin,19)) + 
   (2335763394163431811*lnmcmsmkin*pow(mcMSmusmc,20))/(1.7558329821198565e20*pow(mkin,19)) + (46137065941*pow(lnmcmsmkin,2)*pow(mcMSmusmc,20))/(1.5084957888e12*pow(mkin,19)) - 
   (5507*lnmcmusmcms*pow(mcMSmusmc,20))/(3.338528e6*pow(mkin,19)) + (9*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,20))/(646.*pow(mkin,19)) - (39163*lnmusmcms*pow(mcMSmusmc,20))/(1.0015584e9*pow(mkin,19)) + 
   (3*lnmcmsmkin*lnmusmcms*pow(mcMSmusmc,20))/(12920.*pow(mkin,19)) - (39163*lnmusmkin*pow(mcMSmusmc,20))/(8.0124672e7*pow(mkin,19)) + (15*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,20))/(5168.*pow(mkin,19)) - 
   (372869427163814153380429.*pow(mcMSmusmc,18))/(5.00576874275692e25*pow(mkin,17)) + (53349597791833939*lnmcmsmkin*pow(mcMSmusmc,18))/(3.2684758005112013e18*pow(mkin,17)) + 
   (1014809351*pow(lnmcmsmkin,2)*pow(mcMSmusmc,18))/(2.615348736e10*pow(mkin,17)) - (7678*lnmcmusmcms*pow(mcMSmusmc,18))/(3.186225e6*pow(mkin,17)) + 
   (32*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,18))/(1785.*pow(mkin,17)) - (16277*lnmusmcms*pow(mcMSmusmc,18))/(2.58084225e8*pow(mkin,17)) + 
   (16*lnmcmsmkin*lnmusmcms*pow(mcMSmusmc,18))/(48195.*pow(mkin,17)) - (16277*lnmusmkin*pow(mcMSmusmc,18))/(2.0646738e7*pow(mkin,17)) + (40*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,18))/(9639.*pow(mkin,17)) - 
   (968011304335065184487.*pow(mcMSmusmc,16))/(9.509568138760686e22*pow(mkin,15)) + (23986410569717*lnmcmsmkin*pow(mcMSmusmc,16))/(1.1780838381312e15*pow(mkin,15)) + 
   (232939907*pow(lnmcmsmkin,2)*pow(mcMSmusmc,16))/(4.576860288e9*pow(mkin,15)) - (283*lnmcmusmcms*pow(mcMSmusmc,16))/(76050.*pow(mkin,15)) + 
   (14*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,16))/(585.*pow(mkin,15)) - (529*lnmusmcms*pow(mcMSmusmc,16))/(4.8672e6*pow(mkin,15)) + (7*lnmcmsmkin*lnmusmcms*pow(mcMSmusmc,16))/(14040.*pow(mkin,15)) - 
   (529*lnmusmkin*pow(mcMSmusmc,16))/(389376.*pow(mkin,15)) + (35*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,16))/(5616.*pow(mkin,15)) - 
   (858973162016800669*pow(mcMSmusmc,14))/(5.896309609846656e19*pow(mkin,13)) + (95564969831*lnmcmsmkin*pow(mcMSmusmc,14))/(3.718698984e12*pow(mkin,13)) + 
   (845153*pow(lnmcmsmkin,2)*pow(mcMSmusmc,14))/(1.2108096e7*pow(mkin,13)) - (3166*lnmcmusmcms*pow(mcMSmusmc,14))/(511225.*pow(mkin,13)) + 
   (24*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,14))/(715.*pow(mkin,13)) - (15371*lnmusmcms*pow(mcMSmusmc,14))/(7.5150075e7*pow(mkin,13)) + (4*lnmcmsmkin*lnmusmcms*pow(mcMSmusmc,14))/(5005.*pow(mkin,13)) - 
   (15371*lnmusmkin*pow(mcMSmusmc,14))/(6.012006e6*pow(mkin,13)) + (10*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,14))/(1001.*pow(mkin,13)) - 
   (540989122348457*pow(mcMSmusmc,12))/(2.46471470784e16*pow(mkin,11)) + (199138925387*lnmcmsmkin*pow(mcMSmusmc,12))/(6.22402704e12*pow(mkin,11)) + 
   (15261023*pow(lnmcmsmkin,2)*pow(mcMSmusmc,12))/(1.49688e8*pow(mkin,11)) - (899*lnmcmusmcms*pow(mcMSmusmc,12))/(78408.*pow(mkin,11)) + 
   (5*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,12))/(99.*pow(mkin,11)) - (1229*lnmusmcms*pow(mcMSmusmc,12))/(2.822688e6*pow(mkin,11)) + (5*lnmcmsmkin*lnmusmcms*pow(mcMSmusmc,12))/(3564.*pow(mkin,11)) - 
   (30725*lnmusmkin*pow(mcMSmusmc,12))/(5.645376e6*pow(mkin,11)) + (125*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,12))/(7128.*pow(mkin,11)) - (174878652127*pow(mcMSmusmc,10))/(5.184974592e12*pow(mkin,9)) + 
   (1395062351*lnmcmsmkin*pow(mcMSmusmc,10))/(4.1150592e10*pow(mkin,9)) + (382589*pow(lnmcmsmkin,2)*pow(mcMSmusmc,10))/(2.3328e6*pow(mkin,9)) - 
   (298*lnmcmusmcms*pow(mcMSmusmc,10))/(11907.*pow(mkin,9)) + (16*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,10))/(189.*pow(mkin,9)) - (997*lnmusmcms*pow(mcMSmusmc,10))/(893025.*pow(mkin,9)) + 
   (8*lnmcmsmkin*lnmusmcms*pow(mcMSmusmc,10))/(2835.*pow(mkin,9)) - (997*lnmusmkin*pow(mcMSmusmc,10))/(71442.*pow(mkin,9)) + (20*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,10))/(567.*pow(mkin,9)) - 
   (32093279*pow(mcMSmusmc,8))/(1.185408e9*pow(mkin,7)) - (253711*lnmcmsmkin*pow(mcMSmusmc,8))/(1.27008e7*pow(mkin,7)) + (9451*pow(lnmcmsmkin,2)*pow(mcMSmusmc,8))/(30240.*pow(mkin,7)) - 
   (179*lnmcmusmcms*pow(mcMSmusmc,8))/(2450.*pow(mkin,7)) + (6*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,8))/(35.*pow(mkin,7)) - (463*lnmusmcms*pow(mcMSmusmc,8))/(117600.*pow(mkin,7)) + 
   (lnmcmsmkin*lnmusmcms*pow(mcMSmusmc,8))/(140.*pow(mkin,7)) - (463*lnmusmkin*pow(mcMSmusmc,8))/(9408.*pow(mkin,7)) + (5*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,8))/(56.*pow(mkin,7)) + 
   (10893161*pow(mcMSmusmc,6))/(1.7496e7*pow(mkin,5)) - (169129*lnmcmsmkin*pow(mcMSmusmc,6))/(233280.*pow(mkin,5)) + (5329*pow(lnmcmsmkin,2)*pow(mcMSmusmc,6))/(6480.*pow(mkin,5)) - 
   (94*lnmcmusmcms*pow(mcMSmusmc,6))/(225.*pow(mkin,5)) + (8*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,6))/(15.*pow(mkin,5)) - (19*lnmusmcms*pow(mcMSmusmc,6))/(675.*pow(mkin,5)) + 
   (4*lnmcmsmkin*lnmusmcms*pow(mcMSmusmc,6))/(135.*pow(mkin,5)) - (19*lnmusmkin*pow(mcMSmusmc,6))/(54.*pow(mkin,5)) + (10*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,6))/(27.*pow(mkin,5)) + 
   (10103*pow(mcMSmusmc,4))/(7776.*pow(mkin,3)) - (20*A4*pow(mcMSmusmc,4))/(27.*pow(mkin,3)) - (5*pow(ln2,4)*pow(mcMSmusmc,4))/(162.*pow(mkin,3)) - 
   (4051*lnmcmsmkin*pow(mcMSmusmc,4))/(648.*pow(mkin,3)) + (136*pow(lnmcmsmkin,2)*pow(mcMSmusmc,4))/(27.*pow(mkin,3)) - (7*pow(lnmcmsmkin,3)*pow(mcMSmusmc,4))/(27.*pow(mkin,3)) + 
   (56*lnmcmusmcms*pow(mcMSmusmc,4))/(27.*pow(mkin,3)) - (20*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,4))/(9.*pow(mkin,3)) + (4*pow(lnmcmsmkin,2)*lnmcmusmcms*pow(mcMSmusmc,4))/(3.*pow(mkin,3)) + 
   (151*lnmusmcms*pow(mcMSmusmc,4))/(648.*pow(mkin,3)) - (13*lnmcmsmkin*lnmusmcms*pow(mcMSmusmc,4))/(54.*pow(mkin,3)) + (pow(lnmcmsmkin,2)*lnmusmcms*pow(mcMSmusmc,4))/(9.*pow(mkin,3)) + 
   (3775*lnmusmkin*pow(mcMSmusmc,4))/(1296.*pow(mkin,3)) - (325*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,4))/(108.*pow(mkin,3)) + (25*pow(lnmcmsmkin,2)*lnmusmkin*pow(mcMSmusmc,4))/(18.*pow(mkin,3)) + 
   (86*pow(mcMSmusmc,2))/(9.*mkin) + (4*lnmcmsmkin*pow(mcMSmusmc,2))/mkin + (2*lnmcmusmcms*pow(mcMSmusmc,2))/mkin + (lnmusmcms*pow(mcMSmusmc,2))/(3.*mkin) + (25*lnmusmkin*pow(mcMSmusmc,2))/(6.*mkin) - 
   (8518453*mkin)/93312. + (212*A4*mkin)/27. + (53*pow(ln2,4)*mkin)/162. + (4*lnmcmusmcms*mkin)/9. - (421*lnmusmcms*mkin)/96. - (pow(lnmusmcms,2)*mkin)/27. - (9851*lnmusmkin*mkin)/162. + 
   (lnmcmusmcms*lnmusmkin*mkin)/3. - (101*lnmusmcms*lnmusmkin*mkin)/36. - (pow(lnmusmcms,2)*lnmusmkin*mkin)/36. - (12401*pow(lnmusmkin,2)*mkin)/864. - (7*lnmusmcms*pow(lnmusmkin,2)*mkin)/24. - 
   (505*pow(lnmusmkin,3)*mkin)/432. + (7495*mufac)/9. - (3416*lnmufmus*mufac)/9. + (484*pow(lnmufmus,2)*mufac)/9. + (16*lnmusmcms*mufac)/81. - (982*lnmusmkin*mufac)/27. + (88*lnmufmus*lnmusmkin*mufac)/9. - 
   (8*lnmusmcms*lnmusmkin*mufac)/27. - (14*pow(lnmusmkin,2)*mufac)/9. + (289*pow(mcMSmusmc,24)*mufac)/(198450.*pow(mkin,24)) - (44*lnmcmsmkin*pow(mcMSmusmc,24)*mufac)/(2835.*pow(mkin,24)) + 
   (62384*pow(mcMSmusmc,22)*mufac)/(3.1843449e7*pow(mkin,22)) - (320*lnmcmsmkin*pow(mcMSmusmc,22)*mufac)/(16929.*pow(mkin,22)) + (1417*pow(mcMSmusmc,20)*mufac)/(520200.*pow(mkin,20)) - 
   (2*lnmcmsmkin*pow(mcMSmusmc,20)*mufac)/(85.*pow(mkin,20)) + (10576*pow(mcMSmusmc,18)*mufac)/(2.679075e6*pow(mkin,18)) - (256*lnmcmsmkin*pow(mcMSmusmc,18)*mufac)/(8505.*pow(mkin,18)) + 
   (661*pow(mcMSmusmc,16)*mufac)/(109512.*pow(mkin,16)) - (14*lnmcmsmkin*pow(mcMSmusmc,16)*mufac)/(351.*pow(mkin,16)) + (13232*pow(mcMSmusmc,14)*mufac)/(1.334025e6*pow(mkin,14)) - 
   (64*lnmcmsmkin*pow(mcMSmusmc,14)*mufac)/(1155.*pow(mkin,14)) + (79*pow(mcMSmusmc,12)*mufac)/(4374.*pow(mkin,12)) - (20*lnmcmsmkin*pow(mcMSmusmc,12)*mufac)/(243.*pow(mkin,12)) + 
   (3824*pow(mcMSmusmc,10)*mufac)/(99225.*pow(mkin,10)) - (128*lnmcmsmkin*pow(mcMSmusmc,10)*mufac)/(945.*pow(mkin,10)) + (49*pow(mcMSmusmc,8)*mufac)/(450.*pow(mkin,8)) - 
   (4*lnmcmsmkin*pow(mcMSmusmc,8)*mufac)/(15.*pow(mkin,8)) + (16*pow(mcMSmusmc,6)*mufac)/(27.*pow(mkin,6)) - (64*lnmcmsmkin*pow(mcMSmusmc,6)*mufac)/(81.*pow(mkin,6)) - 
   (22*pow(mcMSmusmc,4)*mufac)/(9.*pow(mkin,4)) + (8*lnmcmsmkin*pow(mcMSmusmc,4)*mufac)/(3.*pow(mkin,4)) - (16*pow(lnmcmsmkin,2)*pow(mcMSmusmc,4)*mufac)/(9.*pow(mkin,4)) - 
   (16*pow(mcMSmusmc,2)*mufac)/(9.*pow(mkin,2)) + (289*pow(mcMSmusmc,24)*pow(mufac,2))/(529200.*pow(mkin,25)) - (11*lnmcmsmkin*pow(mcMSmusmc,24)*pow(mufac,2))/(1890.*pow(mkin,25)) + 
   (7798*pow(mcMSmusmc,22)*pow(mufac,2))/(1.0614483e7*pow(mkin,23)) - (40*lnmcmsmkin*pow(mcMSmusmc,22)*pow(mufac,2))/(5643.*pow(mkin,23)) + 
   (1417*pow(mcMSmusmc,20)*pow(mufac,2))/(1.3872e6*pow(mkin,21)) - (3*lnmcmsmkin*pow(mcMSmusmc,20)*pow(mufac,2))/(340.*pow(mkin,21)) + 
   (1322*pow(mcMSmusmc,18)*pow(mufac,2))/(893025.*pow(mkin,19)) - (32*lnmcmsmkin*pow(mcMSmusmc,18)*pow(mufac,2))/(2835.*pow(mkin,19)) + 
   (661*pow(mcMSmusmc,16)*pow(mufac,2))/(292032.*pow(mkin,17)) - (7*lnmcmsmkin*pow(mcMSmusmc,16)*pow(mufac,2))/(468.*pow(mkin,17)) + (1654*pow(mcMSmusmc,14)*pow(mufac,2))/(444675.*pow(mkin,15)) - 
   (8*lnmcmsmkin*pow(mcMSmusmc,14)*pow(mufac,2))/(385.*pow(mkin,15)) + (79*pow(mcMSmusmc,12)*pow(mufac,2))/(11664.*pow(mkin,13)) - 
   (5*lnmcmsmkin*pow(mcMSmusmc,12)*pow(mufac,2))/(162.*pow(mkin,13)) + (478*pow(mcMSmusmc,10)*pow(mufac,2))/(33075.*pow(mkin,11)) - 
   (16*lnmcmsmkin*pow(mcMSmusmc,10)*pow(mufac,2))/(315.*pow(mkin,11)) + (49*pow(mcMSmusmc,8)*pow(mufac,2))/(1200.*pow(mkin,9)) - (lnmcmsmkin*pow(mcMSmusmc,8)*pow(mufac,2))/(10.*pow(mkin,9)) + 
   (2*pow(mcMSmusmc,6)*pow(mufac,2))/(9.*pow(mkin,7)) - (8*lnmcmsmkin*pow(mcMSmusmc,6)*pow(mufac,2))/(27.*pow(mkin,7)) - (11*pow(mcMSmusmc,4)*pow(mufac,2))/(12.*pow(mkin,5)) + 
   (lnmcmsmkin*pow(mcMSmusmc,4)*pow(mufac,2))/pow(mkin,5) - (2*pow(lnmcmsmkin,2)*pow(mcMSmusmc,4)*pow(mufac,2))/(3.*pow(mkin,5)) - (2*pow(mcMSmusmc,2)*pow(mufac,2))/(3.*pow(mkin,3)) + 
   (151763*pow(mufac,2))/(648.*mkin) - (733*lnmufmus*pow(mufac,2))/(6.*mkin) + (121*pow(lnmufmus,2)*pow(mufac,2))/(6.*mkin) + (2*lnmusmcms*pow(mufac,2))/(27.*mkin) - 
   (425*lnmusmkin*pow(mufac,2))/(36.*mkin) + (11*lnmufmus*lnmusmkin*pow(mufac,2))/(3.*mkin) - (lnmusmcms*lnmusmkin*pow(mufac,2))/(9.*mkin) - (7*pow(lnmusmkin,2)*pow(mufac,2))/(12.*mkin) - 
   (7801530877413386647*pow(mcMSmusmc,24)*NLMSOS)/(3.7763267488425775e22*pow(mkin,23)) - (73801799*lnmcmsmkin*pow(mcMSmusmc,24)*NLMSOS)/(4.23178780752e12*pow(mkin,23)) + 
   (9727*lnmusmkin*pow(mcMSmusmc,24)*NLMSOS)/(5.598936e8*pow(mkin,23)) - (11*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,24)*NLMSOS)/(86940.*pow(mkin,23)) - 
   (1640519393726677*pow(mcMSmusmc,22)*NLMSOS)/(5.946258455651274e18*pow(mkin,21)) - (14290513*lnmcmsmkin*pow(mcMSmusmc,22)*NLMSOS)/(6.89665418442e11*pow(mkin,21)) + 
   (39833*lnmusmkin*pow(mcMSmusmc,22)*NLMSOS)/(1.560329001e9*pow(mkin,21)) - (20*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,22)*NLMSOS)/(118503.*pow(mkin,21)) - 
   (22189567531163017*pow(mcMSmusmc,20)*NLMSOS)/(5.834712481735737e19*pow(mkin,19)) - (1211963*lnmcmsmkin*pow(mcMSmusmc,20)*NLMSOS)/(5.012799792e10*pow(mkin,19)) + 
   (39163*lnmusmkin*pow(mcMSmusmc,20)*NLMSOS)/(1.0015584e9*pow(mkin,19)) - (3*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,20)*NLMSOS)/(12920.*pow(mkin,19)) - 
   (20555048260909*pow(mcMSmusmc,18)*NLMSOS)/(3.76818092235585e16*pow(mkin,17)) - (197062*lnmcmsmkin*pow(mcMSmusmc,18)*NLMSOS)/(7.381208835e9*pow(mkin,17)) + 
   (16277*lnmusmkin*pow(mcMSmusmc,18)*NLMSOS)/(2.58084225e8*pow(mkin,17)) - (16*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,18)*NLMSOS)/(48195.*pow(mkin,17)) - 
   (214558103603*pow(mcMSmusmc,16)*NLMSOS)/(2.60460712512e14*pow(mkin,15)) - (355*lnmcmsmkin*pow(mcMSmusmc,16)*NLMSOS)/(1.4455584e7*pow(mkin,15)) + 
   (529*lnmusmkin*pow(mcMSmusmc,16)*NLMSOS)/(4.8672e6*pow(mkin,15)) - (7*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,16)*NLMSOS)/(14040.*pow(mkin,15)) - 
   (108352091581*pow(mcMSmusmc,14)*NLMSOS)/(8.1243243081e13*pow(mkin,13)) - (1153*lnmcmsmkin*pow(mcMSmusmc,14)*NLMSOS)/(2.25450225e8*pow(mkin,13)) + 
   (15371*lnmusmkin*pow(mcMSmusmc,14)*NLMSOS)/(7.5150075e7*pow(mkin,13)) - (4*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,14)*NLMSOS)/(5005.*pow(mkin,13)) - 
   (326802499*pow(mcMSmusmc,12)*NLMSOS)/(1.3692859488e11*pow(mkin,11)) + (391*lnmcmsmkin*pow(mcMSmusmc,12)*NLMSOS)/(4.939704e6*pow(mkin,11)) + 
   (1229*lnmusmkin*pow(mcMSmusmc,12)*NLMSOS)/(2.822688e6*pow(mkin,11)) - (5*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,12)*NLMSOS)/(3564.*pow(mkin,11)) - 
   (2816347*pow(mcMSmusmc,10)*NLMSOS)/(5.6260575e8*pow(mkin,9)) + (398*lnmcmsmkin*pow(mcMSmusmc,10)*NLMSOS)/(893025.*pow(mkin,9)) + (997*lnmusmkin*pow(mcMSmusmc,10)*NLMSOS)/(893025.*pow(mkin,9)) - 
   (8*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,10)*NLMSOS)/(2835.*pow(mkin,9)) - (174787*pow(mcMSmusmc,8)*NLMSOS)/(1.2348e7*pow(mkin,7)) + (37*lnmcmsmkin*pow(mcMSmusmc,8)*NLMSOS)/(14700.*pow(mkin,7)) + 
   (463*lnmusmkin*pow(mcMSmusmc,8)*NLMSOS)/(117600.*pow(mkin,7)) - (lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,8)*NLMSOS)/(140.*pow(mkin,7)) - (2729*pow(mcMSmusmc,6)*NLMSOS)/(30375.*pow(mkin,5)) + 
   (2*lnmcmsmkin*pow(mcMSmusmc,6)*NLMSOS)/(75.*pow(mkin,5)) + (19*lnmusmkin*pow(mcMSmusmc,6)*NLMSOS)/(675.*pow(mkin,5)) - (4*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,6)*NLMSOS)/(135.*pow(mkin,5)) + 
   (1423*pow(mcMSmusmc,4)*NLMSOS)/(3888.*pow(mkin,3)) + (lnmcmsmkin*pow(mcMSmusmc,4)*NLMSOS)/(12.*pow(mkin,3)) - (13*pow(lnmcmsmkin,2)*pow(mcMSmusmc,4)*NLMSOS)/(54.*pow(mkin,3)) + 
   (2*pow(lnmcmsmkin,3)*pow(mcMSmusmc,4)*NLMSOS)/(27.*pow(mkin,3)) - (151*lnmusmkin*pow(mcMSmusmc,4)*NLMSOS)/(648.*pow(mkin,3)) + (13*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,4)*NLMSOS)/(54.*pow(mkin,3)) - 
   (pow(lnmcmsmkin,2)*lnmusmkin*pow(mcMSmusmc,4)*NLMSOS)/(9.*pow(mkin,3)) - (2*pow(mcMSmusmc,2)*NLMSOS)/(9.*mkin) - (lnmusmkin*pow(mcMSmusmc,2)*NLMSOS)/(3.*mkin) + (241937*mkin*NLMSOS)/23328. - 
   (8*A4*mkin*NLMSOS)/27. - (pow(ln2,4)*mkin*NLMSOS)/81. + (71*lnmusmcms*mkin*NLMSOS)/432. + (2477*lnmusmkin*mkin*NLMSOS)/324. + (13*lnmusmcms*lnmusmkin*mkin*NLMSOS)/108. + 
   (911*pow(lnmusmkin,2)*mkin*NLMSOS)/432. + (lnmusmcms*pow(lnmusmkin,2)*mkin*NLMSOS)/36. + (23*pow(lnmusmkin,3)*mkin*NLMSOS)/108. - (11*mufac*NLMSOS)/27. + (4*lnmusmkin*mufac*NLMSOS)/81. + 
   (4*pow(lnmusmkin,2)*mufac*NLMSOS)/27. - (11*pow(mufac,2)*NLMSOS)/(72.*mkin) + (lnmusmkin*pow(mufac,2)*NLMSOS)/(54.*mkin) + (pow(lnmusmkin,2)*pow(mufac,2)*NLMSOS)/(18.*mkin) - 
   (2353*mkin*pow(NLMSOS,2))/23328. - (89*lnmusmkin*mkin*pow(NLMSOS,2))/648. - (13*pow(lnmusmkin,2)*mkin*pow(NLMSOS,2))/216. - (pow(lnmusmkin,3)*mkin*pow(NLMSOS,2))/108. - (20303*mufac*NLOSKIN)/243. + 
   (3388*lnmufmus*mufac*NLOSKIN)/81. - (176*pow(lnmufmus,2)*mufac*NLOSKIN)/27. + (128*lnmusmkin*mufac*NLOSKIN)/81. - (16*lnmufmus*lnmusmkin*mufac*NLOSKIN)/27. - (14429*pow(mufac,2)*NLOSKIN)/(648.*mkin) + 
   (715*lnmufmus*pow(mufac,2)*NLOSKIN)/(54.*mkin) - (22*pow(lnmufmus,2)*pow(mufac,2)*NLOSKIN)/(9.*mkin) + (13*lnmusmkin*pow(mufac,2)*NLOSKIN)/(27.*mkin) - 
   (2*lnmufmus*lnmusmkin*pow(mufac,2)*NLOSKIN)/(9.*mkin) + (1292*mufac*pow(NLOSKIN,2))/729. - (256*lnmufmus*mufac*pow(NLOSKIN,2))/243. + (16*pow(lnmufmus,2)*mufac*pow(NLOSKIN,2))/81. + 
   (209*pow(mufac,2)*pow(NLOSKIN,2))/(486.*mkin) - (26*lnmufmus*pow(mufac,2)*pow(NLOSKIN,2))/(81.*mkin) + (2*pow(lnmufmus,2)*pow(mufac,2)*pow(NLOSKIN,2))/(27.*mkin) - 
   (30053*mcMSmusmc*pow(Pi,2))/2430. + (1199*ln2*mcMSmusmc*pow(Pi,2))/81. + (31*lnmcmsmkin*mcMSmusmc*pow(Pi,2))/18. - (lnmcmusmcms*mcMSmusmc*pow(Pi,2))/6. - (lnmusmcms*mcMSmusmc*pow(Pi,2))/18. - 
   (25*lnmusmkin*mcMSmusmc*pow(Pi,2))/36. - (106559417*pow(mcMSmusmc,25)*pow(Pi,2))/(6.374388636e11*pow(mkin,24)) + (533*lnmcmsmkin*pow(mcMSmusmc,25)*pow(Pi,2))/(277725.*pow(mkin,24)) - 
   (240817793781176357*pow(mcMSmusmc,24)*pow(Pi,2))/(6.28867544642696e20*pow(mkin,23)) + (877591*ln2*pow(mcMSmusmc,24)*pow(Pi,2))/(6.79477248e9*pow(mkin,23)) + 
   (1620816161*lnmcmsmkin*pow(mcMSmusmc,24)*pow(Pi,2))/(8.2216747008e11*pow(mkin,23)) - (171475369*pow(mcMSmusmc,23)*pow(Pi,2))/(7.7816811996e11*pow(mkin,22)) + 
   (445*lnmcmsmkin*pow(mcMSmusmc,23)*pow(Pi,2))/(192717.*pow(mkin,22)) - (174200135864459*pow(mcMSmusmc,22)*pow(Pi,2))/(3.495434721951744e17*pow(mkin,21)) + 
   (38675*ln2*pow(mcMSmusmc,22)*pow(Pi,2))/(2.33570304e8*pow(mkin,21)) + (13926181*lnmcmsmkin*pow(mcMSmusmc,22)*pow(Pi,2))/(5.8392576e9*pow(mkin,21)) - 
   (846435761*pow(mcMSmusmc,21)*pow(Pi,2))/(2.83231951884e12*pow(mkin,20)) + (365*lnmcmsmkin*pow(mcMSmusmc,21)*pow(Pi,2))/(128877.*pow(mkin,20)) - 
   (9007367733163*pow(mcMSmusmc,20)*pow(Pi,2))/(1.34810154565632e16*pow(mkin,19)) + (143*ln2*pow(mcMSmusmc,20)*pow(Pi,2))/(655360.*pow(mkin,19)) + 
   (156353*lnmcmsmkin*pow(mcMSmusmc,20)*pow(Pi,2))/(5.308416e7*pow(mkin,19)) - (86836957*pow(mcMSmusmc,19)*pow(Pi,2))/(2.0687188752e11*pow(mkin,18)) + 
   (293*lnmcmsmkin*pow(mcMSmusmc,19)*pow(Pi,2))/(82365.*pow(mkin,18)) - (26463251891*pow(mcMSmusmc,18)*pow(Pi,2))/(2.845499424768e13*pow(mkin,17)) + 
   (4147*ln2*pow(mcMSmusmc,18)*pow(Pi,2))/(1.3934592e7*pow(mkin,17)) + (52013*lnmcmsmkin*pow(mcMSmusmc,18)*pow(Pi,2))/(1.3934592e7*pow(mkin,17)) - 
   (22757641*pow(mcMSmusmc,17)*pow(Pi,2))/(3.6923796e10*pow(mkin,16)) + (229*lnmcmsmkin*pow(mcMSmusmc,17)*pow(Pi,2))/(49725.*pow(mkin,16)) - 
   (58806560951*pow(mcMSmusmc,16)*pow(Pi,2))/(4.326189170688e13*pow(mkin,15)) + (1001*ln2*pow(mcMSmusmc,16)*pow(Pi,2))/(2.359296e6*pow(mkin,15)) + 
   (565351*lnmcmsmkin*pow(mcMSmusmc,16)*pow(Pi,2))/(1.15605504e8*pow(mkin,15)) - (17346493*pow(mcMSmusmc,15)*pow(Pi,2))/(1.808754948e10*pow(mkin,14)) + 
   (173*lnmcmsmkin*pow(mcMSmusmc,15)*pow(Pi,2))/(27885.*pow(mkin,14)) - (82285201*pow(mcMSmusmc,14)*pow(Pi,2))/(3.87459072e10*pow(mkin,13)) + 
   (23*ln2*pow(mcMSmusmc,14)*pow(Pi,2))/(35840.*pow(mkin,13)) + (2161*lnmcmsmkin*pow(mcMSmusmc,14)*pow(Pi,2))/(322560.*pow(mkin,13)) - 
   (31786481*pow(mcMSmusmc,13)*pow(Pi,2))/(1.967766372e10*pow(mkin,12)) + (125*lnmcmsmkin*pow(mcMSmusmc,13)*pow(Pi,2))/(14157.*pow(mkin,12)) - 
   (84041429*pow(mcMSmusmc,12)*pow(Pi,2))/(2.29920768e10*pow(mkin,11)) + (175*ln2*pow(mcMSmusmc,12)*pow(Pi,2))/(165888.*pow(mkin,11)) + 
   (40553*lnmcmsmkin*pow(mcMSmusmc,12)*pow(Pi,2))/(4.1472e6*pow(mkin,11)) - (1055689*pow(mcMSmusmc,11)*pow(Pi,2))/(3.4577928e8*pow(mkin,10)) + 
   (85*lnmcmsmkin*pow(mcMSmusmc,11)*pow(Pi,2))/(6237.*pow(mkin,10)) - (359801*pow(mcMSmusmc,10)*pow(Pi,2))/(4.89888e7*pow(mkin,9)) + (17*ln2*pow(mcMSmusmc,10)*pow(Pi,2))/(8640.*pow(mkin,9)) + 
   (17*lnmcmsmkin*pow(mcMSmusmc,10)*pow(Pi,2))/(1080.*pow(mkin,9)) - (416029*pow(mcMSmusmc,9)*pow(Pi,2))/(6.001128e7*pow(mkin,8)) + (53*lnmcmsmkin*pow(mcMSmusmc,9)*pow(Pi,2))/(2205.*pow(mkin,8)) - 
   (1789*pow(mcMSmusmc,8)*pow(Pi,2))/(90720.*pow(mkin,7)) + (7*ln2*pow(mcMSmusmc,8)*pow(Pi,2))/(1536.*pow(mkin,7)) + (139*lnmcmsmkin*pow(mcMSmusmc,8)*pow(Pi,2))/(4608.*pow(mkin,7)) - 
   (262769*pow(mcMSmusmc,7)*pow(Pi,2))/(1.1907e7*pow(mkin,6)) + (29*lnmcmsmkin*pow(mcMSmusmc,7)*pow(Pi,2))/(525.*pow(mkin,6)) - (311*pow(mcMSmusmc,6)*pow(Pi,2))/(2592.*pow(mkin,5)) + 
   (11*ln2*pow(mcMSmusmc,6)*pow(Pi,2))/(648.*pow(mkin,5)) + (113*lnmcmsmkin*pow(mcMSmusmc,6)*pow(Pi,2))/(1296.*pow(mkin,5)) - (149*pow(mcMSmusmc,5)*pow(Pi,2))/(1080.*pow(mkin,4)) + 
   (13*lnmcmsmkin*pow(mcMSmusmc,5)*pow(Pi,2))/(45.*pow(mkin,4)) - (pow(mcMSmusmc,4)*pow(Pi,2))/(12.*pow(mkin,3)) - (5*ln2*pow(mcMSmusmc,4)*pow(Pi,2))/(144.*pow(mkin,3)) - 
   (2*pow(ln2,2)*pow(mcMSmusmc,4)*pow(Pi,2))/(81.*pow(mkin,3)) + (67*lnmcmsmkin*pow(mcMSmusmc,4)*pow(Pi,2))/(144.*pow(mkin,3)) + (ln2*lnmcmsmkin*pow(mcMSmusmc,4)*pow(Pi,2))/(9.*pow(mkin,3)) - 
   (pow(lnmcmsmkin,2)*pow(mcMSmusmc,4)*pow(Pi,2))/(4.*pow(mkin,3)) + (2*lnmcmusmcms*pow(mcMSmusmc,4)*pow(Pi,2))/(9.*pow(mkin,3)) + (lnmusmcms*pow(mcMSmusmc,4)*pow(Pi,2))/(54.*pow(mkin,3)) + 
   (25*lnmusmkin*pow(mcMSmusmc,4)*pow(Pi,2))/(108.*pow(mkin,3)) - (19507*pow(mcMSmusmc,3)*pow(Pi,2))/(1620.*pow(mkin,2)) + (1199*ln2*pow(mcMSmusmc,3)*pow(Pi,2))/(81.*pow(mkin,2)) + 
   (271*lnmcmsmkin*pow(mcMSmusmc,3)*pow(Pi,2))/(162.*pow(mkin,2)) - (lnmcmusmcms*pow(mcMSmusmc,3)*pow(Pi,2))/(2.*pow(mkin,2)) - (lnmusmcms*pow(mcMSmusmc,3)*pow(Pi,2))/(18.*pow(mkin,2)) - 
   (25*lnmusmkin*pow(mcMSmusmc,3)*pow(Pi,2))/(36.*pow(mkin,2)) + (13*pow(mcMSmusmc,2)*pow(Pi,2))/(12.*mkin) + (3*lnmcmsmkin*pow(mcMSmusmc,2)*pow(Pi,2))/(2.*mkin) - (587741*mkin*pow(Pi,2))/38880. + 
   (203*ln2*mkin*pow(Pi,2))/54. + (20*pow(ln2,2)*mkin*pow(Pi,2))/81. - (5*lnmusmcms*mkin*pow(Pi,2))/54. - (ln2*lnmusmcms*mkin*pow(Pi,2))/27. - (125*lnmusmkin*mkin*pow(Pi,2))/108. - 
   (25*ln2*lnmusmkin*mkin*pow(Pi,2))/54. - (3154*mufac*pow(Pi,2))/81. - (16*ln2*mufac*pow(Pi,2))/81. + (88*lnmufmus*mufac*pow(Pi,2))/9. + (8*lnmusmkin*mufac*pow(Pi,2))/9. - 
   (8*pow(mcMSmusmc,4)*mufac*pow(Pi,2))/(27.*pow(mkin,4)) + (16*pow(mcMSmusmc,3)*mufac*pow(Pi,2))/(27.*pow(mkin,3)) - (pow(mcMSmusmc,4)*pow(mufac,2)*pow(Pi,2))/(9.*pow(mkin,5)) + 
   (2*pow(mcMSmusmc,3)*pow(mufac,2)*pow(Pi,2))/(9.*pow(mkin,4)) - (1379*pow(mufac,2)*pow(Pi,2))/(108.*mkin) - (2*ln2*pow(mufac,2)*pow(Pi,2))/(27.*mkin) + 
   (11*lnmufmus*pow(mufac,2)*pow(Pi,2))/(3.*mkin) + (lnmusmkin*pow(mufac,2)*pow(Pi,2))/(3.*mkin) + (7*mcMSmusmc*NLMSOS*pow(Pi,2))/27. - (2*ln2*mcMSmusmc*NLMSOS*pow(Pi,2))/9. - 
   (lnmcmsmkin*mcMSmusmc*NLMSOS*pow(Pi,2))/9. + (lnmusmkin*mcMSmusmc*NLMSOS*pow(Pi,2))/18. + (23*pow(mcMSmusmc,25)*NLMSOS*pow(Pi,2))/(415800.*pow(mkin,24)) - 
   (11*pow(mcMSmusmc,24)*NLMSOS*pow(Pi,2))/(260820.*pow(mkin,23)) + (7*pow(mcMSmusmc,23)*NLMSOS*pow(Pi,2))/(96140.*pow(mkin,22)) - (20*pow(mcMSmusmc,22)*NLMSOS*pow(Pi,2))/(355509.*pow(mkin,21)) + 
   (19*pow(mcMSmusmc,21)*NLMSOS*pow(Pi,2))/(192780.*pow(mkin,20)) - (pow(mcMSmusmc,20)*NLMSOS*pow(Pi,2))/(12920.*pow(mkin,19)) + (17*pow(mcMSmusmc,19)*NLMSOS*pow(Pi,2))/(123120.*pow(mkin,18)) - 
   (16*pow(mcMSmusmc,18)*NLMSOS*pow(Pi,2))/(144585.*pow(mkin,17)) + (5*pow(mcMSmusmc,17)*NLMSOS*pow(Pi,2))/(24752.*pow(mkin,16)) - (7*pow(mcMSmusmc,16)*NLMSOS*pow(Pi,2))/(42120.*pow(mkin,15)) + 
   (13*pow(mcMSmusmc,15)*NLMSOS*pow(Pi,2))/(41580.*pow(mkin,14)) - (4*pow(mcMSmusmc,14)*NLMSOS*pow(Pi,2))/(15015.*pow(mkin,13)) + (11*pow(mcMSmusmc,13)*NLMSOS*pow(Pi,2))/(21060.*pow(mkin,12)) - 
   (5*pow(mcMSmusmc,12)*NLMSOS*pow(Pi,2))/(10692.*pow(mkin,11)) + (3*pow(mcMSmusmc,11)*NLMSOS*pow(Pi,2))/(3080.*pow(mkin,10)) - (8*pow(mcMSmusmc,10)*NLMSOS*pow(Pi,2))/(8505.*pow(mkin,9)) + 
   (7*pow(mcMSmusmc,9)*NLMSOS*pow(Pi,2))/(3240.*pow(mkin,8)) - (pow(mcMSmusmc,8)*NLMSOS*pow(Pi,2))/(420.*pow(mkin,7)) + (5*pow(mcMSmusmc,7)*NLMSOS*pow(Pi,2))/(756.*pow(mkin,6)) - 
   (4*pow(mcMSmusmc,6)*NLMSOS*pow(Pi,2))/(405.*pow(mkin,5)) + (pow(mcMSmusmc,5)*NLMSOS*pow(Pi,2))/(20.*pow(mkin,4)) + (13*pow(mcMSmusmc,4)*NLMSOS*pow(Pi,2))/(324.*pow(mkin,3)) - 
   (lnmcmsmkin*pow(mcMSmusmc,4)*NLMSOS*pow(Pi,2))/(27.*pow(mkin,3)) - (lnmusmkin*pow(mcMSmusmc,4)*NLMSOS*pow(Pi,2))/(54.*pow(mkin,3)) + (7*pow(mcMSmusmc,3)*NLMSOS*pow(Pi,2))/(54.*pow(mkin,2)) - 
   (2*ln2*pow(mcMSmusmc,3)*NLMSOS*pow(Pi,2))/(9.*pow(mkin,2)) - (lnmcmsmkin*pow(mcMSmusmc,3)*NLMSOS*pow(Pi,2))/(9.*pow(mkin,2)) + (lnmusmkin*pow(mcMSmusmc,3)*NLMSOS*pow(Pi,2))/(18.*pow(mkin,2)) + 
   (305*mkin*NLMSOS*pow(Pi,2))/216. + (11*ln2*mkin*NLMSOS*pow(Pi,2))/81. - (2*pow(ln2,2)*mkin*NLMSOS*pow(Pi,2))/81. + (lnmusmcms*mkin*NLMSOS*pow(Pi,2))/54. + (35*lnmusmkin*mkin*NLMSOS*pow(Pi,2))/108. + 
   (ln2*lnmusmkin*mkin*NLMSOS*pow(Pi,2))/27. + (8*mufac*NLMSOS*pow(Pi,2))/81. + (pow(mufac,2)*NLMSOS*pow(Pi,2))/(27.*mkin) - (13*mkin*pow(NLMSOS,2)*pow(Pi,2))/324. - 
   (lnmusmkin*mkin*pow(NLMSOS,2)*pow(Pi,2))/54. + (208*mufac*NLOSKIN*pow(Pi,2))/81. - (16*lnmufmus*mufac*NLOSKIN*pow(Pi,2))/27. + (23*pow(mufac,2)*NLOSKIN*pow(Pi,2))/(27.*mkin) - 
   (2*lnmufmus*pow(mufac,2)*NLOSKIN*pow(Pi,2))/(9.*mkin) - (8*mufac*pow(NLOSKIN,2)*pow(Pi,2))/243. - (pow(mufac,2)*pow(NLOSKIN,2)*pow(Pi,2))/(81.*mkin) - (13*mcMSmusmc*pow(Pi,3))/162. + 
   (5698043*pow(mcMSmusmc,25)*pow(Pi,3))/(1.189085184e11*pow(mkin,24)) + (81991*pow(mcMSmusmc,23)*pow(Pi,3))/(1.374683136e9*pow(mkin,22)) + 
   (127699*pow(mcMSmusmc,21)*pow(Pi,3))/(1.684537344e9*pow(mkin,20)) + (99671*pow(mcMSmusmc,19)*pow(Pi,3))/(1.00859904e9*pow(mkin,18)) + 
   (1925*pow(mcMSmusmc,17)*pow(Pi,3))/(1.4483456e7*pow(mkin,16)) + (377*pow(mcMSmusmc,15)*pow(Pi,3))/(2.02752e6*pow(mkin,14)) + (1771*pow(mcMSmusmc,13)*pow(Pi,3))/(6.469632e6*pow(mkin,12)) + 
   (17*pow(mcMSmusmc,11)*pow(Pi,3))/(39424.*pow(mkin,10)) + (77*pow(mcMSmusmc,9)*pow(Pi,3))/(103680.*pow(mkin,8)) + (25*pow(mcMSmusmc,7)*pow(Pi,3))/(18144.*pow(mkin,6)) - 
   (pow(mcMSmusmc,5)*pow(Pi,3))/(240.*pow(mkin,4)) - (7*pow(mcMSmusmc,3)*pow(Pi,3))/(108.*pow(mkin,2)) - (271*pow(mcMSmusmc,4)*pow(Pi,4))/(19440.*pow(mkin,3)) - 
   (lnmcmsmkin*pow(mcMSmusmc,2)*pow(Pi,4))/(12.*mkin) + (451*mkin*pow(Pi,4))/7776. + (2*mufac*pow(Pi,4))/3. + (pow(mufac,2)*pow(Pi,4))/(4.*mkin) - (61*mkin*NLMSOS*pow(Pi,4))/1944. + 
   (2710689767*pow(mcMSmusmc,24)*Zeta3)/(1.64433494016e12*pow(mkin,23)) + (23017987*pow(mcMSmusmc,22)*Zeta3)/(1.16785152e10*pow(mkin,21)) + (254791*pow(mcMSmusmc,20)*Zeta3)/(1.0616832e8*pow(mkin,19)) + 
   (83291*pow(mcMSmusmc,18)*Zeta3)/(2.7869184e7*pow(mkin,17)) + (885457*pow(mcMSmusmc,16)*Zeta3)/(2.31211008e8*pow(mkin,15)) + (3287*pow(mcMSmusmc,14)*Zeta3)/(645120.*pow(mkin,13)) + 
   (59231*pow(mcMSmusmc,12)*Zeta3)/(8.2944e6*pow(mkin,11)) + (187*pow(mcMSmusmc,10)*Zeta3)/(17280.*pow(mkin,9)) + (173*pow(mcMSmusmc,8)*Zeta3)/(9216.*pow(mkin,7)) + 
   (29*pow(mcMSmusmc,6)*Zeta3)/(648.*pow(mkin,5)) + (2309*pow(mcMSmusmc,4)*Zeta3)/(864.*pow(mkin,3)) - (14*lnmcmsmkin*pow(mcMSmusmc,4)*Zeta3)/(9.*pow(mkin,3)) + (11*pow(mcMSmusmc,2)*Zeta3)/(2.*mkin) + 
   (23*mkin*Zeta3)/24. + (lnmusmcms*mkin*Zeta3)/18. + (85*lnmusmkin*mkin*Zeta3)/36. - (3070*mufac*Zeta3)/27. - (1535*pow(mufac,2)*Zeta3)/(36.*mkin) - (pow(mcMSmusmc,4)*NLMSOS*Zeta3)/(3.*pow(mkin,3)) + 
   (667*mkin*NLMSOS*Zeta3)/216. + (7*lnmusmkin*mkin*NLMSOS*Zeta3)/9. - (7*mkin*pow(NLMSOS,2)*Zeta3)/54. + (140*mufac*NLOSKIN*Zeta3)/27. + (35*pow(mufac,2)*NLOSKIN*Zeta3)/(18.*mkin) - 
   (3*pow(mcMSmusmc,2)*pow(Pi,2)*Zeta3)/(4.*mkin) + (1439*mkin*pow(Pi,2)*Zeta3)/432. - (5*pow(mcMSmusmc,2)*Zeta5)/(2.*mkin) - (1975*mkin*Zeta5)/216.);
    
    
    double mMS = 0.0;
    for(int i = 0; i <= nloops; i++) {
       mMS += ret[i];
    }    
    return mMS;
}

double CRunDec::mkin2mMSB(double mkin, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops) {
    if(nloops<0||nloops>3){
      cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
      RETURN
    }
    double ret[5];
    double lnmusmkin = log(mus*mus/(mkin*mkin));
    double lnmufmus = log(2*mufac/mus);
    double lnmusmcms = 0.0;
    double lnmcmsmkin = 0.0;
    double lnmcmusmcms = 0.0;
    double lnmusmusmc = 0.0;  
    double mcMSmusmc = 0.0;
    double musmc = 0.0;
    if(mcMSmusmcin > 0.0) {
       lnmusmcms = log(mus*mus/(mcMSmusmcin*mcMSmusmcin));
       lnmusmusmc = log(mus*mus/(musmcin*musmcin));
       lnmcmsmkin = log(mcMSmusmcin/mkin);
       lnmcmusmcms = log(musmcin*musmcin/(mcMSmusmcin*mcMSmusmcin));
       mcMSmusmc = mcMSmusmcin;
       musmc = musmcin;
    }
    
    ret[0] = mkin;
    ret[1] = apinlmus*((-4*mkin)/3. - lnmusmkin*mkin + (16*mufac)/9. + (2*pow(mufac,2))/(3.*mkin));
    ret[2] = apinlmus*apinlmus*((-9727*pow(mcMSmusmc,24))/(1.866312e8*pow(mkin,23)) + (11*lnmcmsmkin*pow(mcMSmusmc,24))/(28980.*pow(mkin,23)) - (39833*pow(mcMSmusmc,22))/(5.20109667e8*pow(mkin,21)) + 
   (20*lnmcmsmkin*pow(mcMSmusmc,22))/(39501.*pow(mkin,21)) - (39163*pow(mcMSmusmc,20))/(3.338528e8*pow(mkin,19)) + (9*lnmcmsmkin*pow(mcMSmusmc,20))/(12920.*pow(mkin,19)) - 
   (16277*pow(mcMSmusmc,18))/(8.6028075e7*pow(mkin,17)) + (16*lnmcmsmkin*pow(mcMSmusmc,18))/(16065.*pow(mkin,17)) - (529*pow(mcMSmusmc,16))/(1.6224e6*pow(mkin,15)) + 
   (7*lnmcmsmkin*pow(mcMSmusmc,16))/(4680.*pow(mkin,15)) - (15371*pow(mcMSmusmc,14))/(2.5050025e7*pow(mkin,13)) + (12*lnmcmsmkin*pow(mcMSmusmc,14))/(5005.*pow(mkin,13)) - 
   (1229*pow(mcMSmusmc,12))/(940896.*pow(mkin,11)) + (5*lnmcmsmkin*pow(mcMSmusmc,12))/(1188.*pow(mkin,11)) - (997*pow(mcMSmusmc,10))/(297675.*pow(mkin,9)) + 
   (8*lnmcmsmkin*pow(mcMSmusmc,10))/(945.*pow(mkin,9)) - (463*pow(mcMSmusmc,8))/(39200.*pow(mkin,7)) + (3*lnmcmsmkin*pow(mcMSmusmc,8))/(140.*pow(mkin,7)) - (19*pow(mcMSmusmc,6))/(225.*pow(mkin,5)) + 
   (4*lnmcmsmkin*pow(mcMSmusmc,6))/(45.*pow(mkin,5)) + (151*pow(mcMSmusmc,4))/(216.*pow(mkin,3)) - (13*lnmcmsmkin*pow(mcMSmusmc,4))/(18.*pow(mkin,3)) + 
   (pow(lnmcmsmkin,2)*pow(mcMSmusmc,4))/(3.*pow(mkin,3)) + pow(mcMSmusmc,2)/mkin - (959*mkin)/96. - (145*lnmusmkin*mkin)/24. - (7*pow(lnmusmkin,2)*mkin)/8. + (892*mufac)/27. - (88*lnmufmus*mufac)/9. - 
   (8*lnmusmcms*mufac)/27. - (16*lnmusmkin*mufac)/9. + (95*pow(mufac,2))/(9.*mkin) - (11*lnmufmus*pow(mufac,2))/(3.*mkin) - (lnmusmcms*pow(mufac,2))/(9.*mkin) - (2*lnmusmkin*pow(mufac,2))/(3.*mkin) + 
   (71*mkin*NLMSOS)/144. + (13*lnmusmkin*mkin*NLMSOS)/36. + (pow(lnmusmkin,2)*mkin*NLMSOS)/12. - (128*mufac*NLOSKIN)/81. + (16*lnmufmus*mufac*NLOSKIN)/27. - (13*pow(mufac,2)*NLOSKIN)/(27.*mkin) + 
   (2*lnmufmus*pow(mufac,2)*NLOSKIN)/(9.*mkin) - (mcMSmusmc*pow(Pi,2))/6. + (pow(mcMSmusmc,4)*pow(Pi,2))/(18.*pow(mkin,3)) - (pow(mcMSmusmc,3)*pow(Pi,2))/(6.*pow(mkin,2)) - 
   (5*mkin*pow(Pi,2))/18. - (ln2*mkin*pow(Pi,2))/9. - (8*mufac*pow(Pi,2))/9. - (pow(mufac,2)*pow(Pi,2))/(3.*mkin) + (mkin*NLMSOS*pow(Pi,2))/18. + (mkin*Zeta3)/6.);
    ret[3] = apinlmus*apinlmus*apinlmus*((-8439929370090150606400742890600889.*pow(mcMSmusmc,24))/(2.369996943791664e36*pow(mkin,23)) + (260257620268030848013.*lnmcmsmkin*pow(mcMSmusmc,24))/(2.8097278338127475e22*pow(mkin,23)) + 
   (953082228281*pow(lnmcmsmkin,2)*pow(mcMSmusmc,24))/(4.664604200256e13*pow(mkin,23)) - (10163*lnmcmusmcms*pow(mcMSmusmc,24))/(1.166445e7*pow(mkin,23)) + 
   (22*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,24))/(2415.*pow(mkin,23)) - (9727*lnmusmkin*pow(mcMSmusmc,24))/(4.4791488e7*pow(mkin,23)) + 
   (55*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,24))/(34776.*pow(mkin,23)) - (47848557349149553062713793193.*pow(mcMSmusmc,22))/(1.0790864167779053e31*pow(mkin,21)) + 
   (79809944445306329123.*lnmcmsmkin*pow(mcMSmusmc,22))/(7.242811051244408e21*pow(mkin,21)) + (3078965960711*pow(lnmcmsmkin,2)*pow(mcMSmusmc,22))/(1.24450902576e14*pow(mkin,21)) - 
   (5066*lnmcmusmcms*pow(mcMSmusmc,22))/(4.298427e6*pow(mkin,21)) + (40*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,22))/(3591.*pow(mkin,21)) - 
   (995825*lnmusmkin*pow(mcMSmusmc,22))/(3.120658002e9*pow(mkin,21)) + (250*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,22))/(118503.*pow(mkin,21)) - 
   (308474675741863370556669649.*pow(mcMSmusmc,20))/(5.449931397868209e28*pow(mkin,19)) + (2335763394163431811.*lnmcmsmkin*pow(mcMSmusmc,20))/(1.7558329821198565e20*pow(mkin,19)) + 
   (46137065941*pow(lnmcmsmkin,2)*pow(mcMSmusmc,20))/(1.5084957888e12*pow(mkin,19)) - (5507*lnmcmusmcms*pow(mcMSmusmc,20))/(3.338528e6*pow(mkin,19)) + 
   (9*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,20))/(646.*pow(mkin,19)) - (39163*lnmusmkin*pow(mcMSmusmc,20))/(8.0124672e7*pow(mkin,19)) + (15*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,20))/(5168.*pow(mkin,19)) - 
   (372869427163814153380429.*pow(mcMSmusmc,18))/(5.00576874275692e25*pow(mkin,17)) + (53349597791833939*lnmcmsmkin*pow(mcMSmusmc,18))/(3.2684758005112013e18*pow(mkin,17)) + 
   (1014809351*pow(lnmcmsmkin,2)*pow(mcMSmusmc,18))/(2.615348736e10*pow(mkin,17)) - (7678*lnmcmusmcms*pow(mcMSmusmc,18))/(3.186225e6*pow(mkin,17)) + 
   (32*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,18))/(1785.*pow(mkin,17)) - (16277*lnmusmkin*pow(mcMSmusmc,18))/(2.0646738e7*pow(mkin,17)) + 
   (40*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,18))/(9639.*pow(mkin,17)) - (968011304335065184487.*pow(mcMSmusmc,16))/(9.509568138760686e22*pow(mkin,15)) + 
   (23986410569717*lnmcmsmkin*pow(mcMSmusmc,16))/(1.1780838381312e15*pow(mkin,15)) + (232939907*pow(lnmcmsmkin,2)*pow(mcMSmusmc,16))/(4.576860288e9*pow(mkin,15)) - 
   (283*lnmcmusmcms*pow(mcMSmusmc,16))/(76050.*pow(mkin,15)) + (14*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,16))/(585.*pow(mkin,15)) - (529*lnmusmkin*pow(mcMSmusmc,16))/(389376.*pow(mkin,15)) + 
   (35*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,16))/(5616.*pow(mkin,15)) - (858973162016800669*pow(mcMSmusmc,14))/(5.896309609846656e19*pow(mkin,13)) + 
   (95564969831*lnmcmsmkin*pow(mcMSmusmc,14))/(3.718698984e12*pow(mkin,13)) + (845153*pow(lnmcmsmkin,2)*pow(mcMSmusmc,14))/(1.2108096e7*pow(mkin,13)) - 
   (3166*lnmcmusmcms*pow(mcMSmusmc,14))/(511225.*pow(mkin,13)) + (24*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,14))/(715.*pow(mkin,13)) - (15371*lnmusmkin*pow(mcMSmusmc,14))/(6.012006e6*pow(mkin,13)) + 
   (10*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,14))/(1001.*pow(mkin,13)) - (540989122348457*pow(mcMSmusmc,12))/(2.46471470784e16*pow(mkin,11)) + 
   (199138925387*lnmcmsmkin*pow(mcMSmusmc,12))/(6.22402704e12*pow(mkin,11)) + (15261023*pow(lnmcmsmkin,2)*pow(mcMSmusmc,12))/(1.49688e8*pow(mkin,11)) - 
   (899*lnmcmusmcms*pow(mcMSmusmc,12))/(78408.*pow(mkin,11)) + (5*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,12))/(99.*pow(mkin,11)) - (30725*lnmusmkin*pow(mcMSmusmc,12))/(5.645376e6*pow(mkin,11)) + 
   (125*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,12))/(7128.*pow(mkin,11)) - (174878652127*pow(mcMSmusmc,10))/(5.184974592e12*pow(mkin,9)) + 
   (1395062351*lnmcmsmkin*pow(mcMSmusmc,10))/(4.1150592e10*pow(mkin,9)) + (382589*pow(lnmcmsmkin,2)*pow(mcMSmusmc,10))/(2.3328e6*pow(mkin,9)) - 
   (298*lnmcmusmcms*pow(mcMSmusmc,10))/(11907.*pow(mkin,9)) + (16*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,10))/(189.*pow(mkin,9)) - (997*lnmusmkin*pow(mcMSmusmc,10))/(71442.*pow(mkin,9)) + 
   (20*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,10))/(567.*pow(mkin,9)) - (32093279*pow(mcMSmusmc,8))/(1.185408e9*pow(mkin,7)) - (253711*lnmcmsmkin*pow(mcMSmusmc,8))/(1.27008e7*pow(mkin,7)) + 
   (9451*pow(lnmcmsmkin,2)*pow(mcMSmusmc,8))/(30240.*pow(mkin,7)) - (179*lnmcmusmcms*pow(mcMSmusmc,8))/(2450.*pow(mkin,7)) + (6*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,8))/(35.*pow(mkin,7)) - 
   (463*lnmusmkin*pow(mcMSmusmc,8))/(9408.*pow(mkin,7)) + (5*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,8))/(56.*pow(mkin,7)) + (10893161*pow(mcMSmusmc,6))/(1.7496e7*pow(mkin,5)) - 
   (169129*lnmcmsmkin*pow(mcMSmusmc,6))/(233280.*pow(mkin,5)) + (5329*pow(lnmcmsmkin,2)*pow(mcMSmusmc,6))/(6480.*pow(mkin,5)) - (94*lnmcmusmcms*pow(mcMSmusmc,6))/(225.*pow(mkin,5)) + 
   (8*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,6))/(15.*pow(mkin,5)) - (19*lnmusmkin*pow(mcMSmusmc,6))/(54.*pow(mkin,5)) + (10*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,6))/(27.*pow(mkin,5)) + 
   (10103*pow(mcMSmusmc,4))/(7776.*pow(mkin,3)) - (20*A4*pow(mcMSmusmc,4))/(27.*pow(mkin,3)) - (5*pow(ln2,4)*pow(mcMSmusmc,4))/(162.*pow(mkin,3)) - 
   (4051*lnmcmsmkin*pow(mcMSmusmc,4))/(648.*pow(mkin,3)) + (136*pow(lnmcmsmkin,2)*pow(mcMSmusmc,4))/(27.*pow(mkin,3)) - (7*pow(lnmcmsmkin,3)*pow(mcMSmusmc,4))/(27.*pow(mkin,3)) + 
   (56*lnmcmusmcms*pow(mcMSmusmc,4))/(27.*pow(mkin,3)) - (20*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,4))/(9.*pow(mkin,3)) + (4*pow(lnmcmsmkin,2)*lnmcmusmcms*pow(mcMSmusmc,4))/(3.*pow(mkin,3)) + 
   (3775*lnmusmkin*pow(mcMSmusmc,4))/(1296.*pow(mkin,3)) - (325*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,4))/(108.*pow(mkin,3)) + (25*pow(lnmcmsmkin,2)*lnmusmkin*pow(mcMSmusmc,4))/(18.*pow(mkin,3)) + 
   (86*pow(mcMSmusmc,2))/(9.*mkin) + (4*lnmcmsmkin*pow(mcMSmusmc,2))/mkin + (2*lnmcmusmcms*pow(mcMSmusmc,2))/mkin + (25*lnmusmkin*pow(mcMSmusmc,2))/(6.*mkin) - (8537461*mkin)/93312. + (212*A4*mkin)/27. + 
   (53*pow(ln2,4)*mkin)/162. + (4*lnmcmusmcms*mkin)/9. - (4*lnmusmcms*mkin)/9. - (39503*lnmusmkin*mkin)/648. + (lnmcmusmcms*lnmusmkin*mkin)/3. - (lnmusmcms*lnmusmkin*mkin)/3. - 
   (12401*pow(lnmusmkin,2)*mkin)/864. - (505*pow(lnmusmkin,3)*mkin)/432. + (4*lnmusmusmc*mkin)/9. + (lnmusmkin*lnmusmusmc*mkin)/3. + (67477*mufac)/81. - (3416*lnmufmus*mufac)/9. + 
   (484*pow(lnmufmus,2)*mufac)/9. - (314*lnmusmcms*mufac)/27. + (88*lnmufmus*lnmusmcms*mufac)/27. + (4*pow(lnmusmcms,2)*mufac)/81. - (982*lnmusmkin*mufac)/27. + (88*lnmufmus*lnmusmkin*mufac)/9. + 
   (8*lnmusmcms*lnmusmkin*mufac)/27. - (14*pow(lnmusmkin,2)*mufac)/9. - (16*lnmusmusmc*mufac)/27. + (289*pow(mcMSmusmc,24)*mufac)/(198450.*pow(mkin,24)) - 
   (44*lnmcmsmkin*pow(mcMSmusmc,24)*mufac)/(2835.*pow(mkin,24)) + (62384*pow(mcMSmusmc,22)*mufac)/(3.1843449e7*pow(mkin,22)) - (320*lnmcmsmkin*pow(mcMSmusmc,22)*mufac)/(16929.*pow(mkin,22)) + 
   (1417*pow(mcMSmusmc,20)*mufac)/(520200.*pow(mkin,20)) - (2*lnmcmsmkin*pow(mcMSmusmc,20)*mufac)/(85.*pow(mkin,20)) + (10576*pow(mcMSmusmc,18)*mufac)/(2.679075e6*pow(mkin,18)) - 
   (256*lnmcmsmkin*pow(mcMSmusmc,18)*mufac)/(8505.*pow(mkin,18)) + (661*pow(mcMSmusmc,16)*mufac)/(109512.*pow(mkin,16)) - (14*lnmcmsmkin*pow(mcMSmusmc,16)*mufac)/(351.*pow(mkin,16)) + 
   (13232*pow(mcMSmusmc,14)*mufac)/(1.334025e6*pow(mkin,14)) - (64*lnmcmsmkin*pow(mcMSmusmc,14)*mufac)/(1155.*pow(mkin,14)) + (79*pow(mcMSmusmc,12)*mufac)/(4374.*pow(mkin,12)) - 
   (20*lnmcmsmkin*pow(mcMSmusmc,12)*mufac)/(243.*pow(mkin,12)) + (3824*pow(mcMSmusmc,10)*mufac)/(99225.*pow(mkin,10)) - (128*lnmcmsmkin*pow(mcMSmusmc,10)*mufac)/(945.*pow(mkin,10)) + 
   (49*pow(mcMSmusmc,8)*mufac)/(450.*pow(mkin,8)) - (4*lnmcmsmkin*pow(mcMSmusmc,8)*mufac)/(15.*pow(mkin,8)) + (16*pow(mcMSmusmc,6)*mufac)/(27.*pow(mkin,6)) - 
   (64*lnmcmsmkin*pow(mcMSmusmc,6)*mufac)/(81.*pow(mkin,6)) - (22*pow(mcMSmusmc,4)*mufac)/(9.*pow(mkin,4)) + (8*lnmcmsmkin*pow(mcMSmusmc,4)*mufac)/(3.*pow(mkin,4)) - 
   (16*pow(lnmcmsmkin,2)*pow(mcMSmusmc,4)*mufac)/(9.*pow(mkin,4)) - (16*pow(mcMSmusmc,2)*mufac)/(9.*pow(mkin,2)) + (289*pow(mcMSmusmc,24)*pow(mufac,2))/(529200.*pow(mkin,25)) - 
   (11*lnmcmsmkin*pow(mcMSmusmc,24)*pow(mufac,2))/(1890.*pow(mkin,25)) + (7798*pow(mcMSmusmc,22)*pow(mufac,2))/(1.0614483e7*pow(mkin,23)) - 
   (40*lnmcmsmkin*pow(mcMSmusmc,22)*pow(mufac,2))/(5643.*pow(mkin,23)) + (1417*pow(mcMSmusmc,20)*pow(mufac,2))/(1.3872e6*pow(mkin,21)) - 
   (3*lnmcmsmkin*pow(mcMSmusmc,20)*pow(mufac,2))/(340.*pow(mkin,21)) + (1322*pow(mcMSmusmc,18)*pow(mufac,2))/(893025.*pow(mkin,19)) - 
   (32*lnmcmsmkin*pow(mcMSmusmc,18)*pow(mufac,2))/(2835.*pow(mkin,19)) + (661*pow(mcMSmusmc,16)*pow(mufac,2))/(292032.*pow(mkin,17)) - 
   (7*lnmcmsmkin*pow(mcMSmusmc,16)*pow(mufac,2))/(468.*pow(mkin,17)) + (1654*pow(mcMSmusmc,14)*pow(mufac,2))/(444675.*pow(mkin,15)) - 
   (8*lnmcmsmkin*pow(mcMSmusmc,14)*pow(mufac,2))/(385.*pow(mkin,15)) + (79*pow(mcMSmusmc,12)*pow(mufac,2))/(11664.*pow(mkin,13)) - 
   (5*lnmcmsmkin*pow(mcMSmusmc,12)*pow(mufac,2))/(162.*pow(mkin,13)) + (478*pow(mcMSmusmc,10)*pow(mufac,2))/(33075.*pow(mkin,11)) - 
   (16*lnmcmsmkin*pow(mcMSmusmc,10)*pow(mufac,2))/(315.*pow(mkin,11)) + (49*pow(mcMSmusmc,8)*pow(mufac,2))/(1200.*pow(mkin,9)) - (lnmcmsmkin*pow(mcMSmusmc,8)*pow(mufac,2))/(10.*pow(mkin,9)) + 
   (2*pow(mcMSmusmc,6)*pow(mufac,2))/(9.*pow(mkin,7)) - (8*lnmcmsmkin*pow(mcMSmusmc,6)*pow(mufac,2))/(27.*pow(mkin,7)) - (11*pow(mcMSmusmc,4)*pow(mufac,2))/(12.*pow(mkin,5)) + 
   (lnmcmsmkin*pow(mcMSmusmc,4)*pow(mufac,2))/pow(mkin,5) - (2*pow(lnmcmsmkin,2)*pow(mcMSmusmc,4)*pow(mufac,2))/(3.*pow(mkin,5)) - (2*pow(mcMSmusmc,2)*pow(mufac,2))/(3.*pow(mkin,3)) + 
   (151829*pow(mufac,2))/(648.*mkin) - (733*lnmufmus*pow(mufac,2))/(6.*mkin) + (121*pow(lnmufmus,2)*pow(mufac,2))/(6.*mkin) - (15*lnmusmcms*pow(mufac,2))/(4.*mkin) + 
   (11*lnmufmus*lnmusmcms*pow(mufac,2))/(9.*mkin) + (pow(lnmusmcms,2)*pow(mufac,2))/(54.*mkin) - (425*lnmusmkin*pow(mufac,2))/(36.*mkin) + (11*lnmufmus*lnmusmkin*pow(mufac,2))/(3.*mkin) + 
   (lnmusmcms*lnmusmkin*pow(mufac,2))/(9.*mkin) - (7*pow(lnmusmkin,2)*pow(mufac,2))/(12.*mkin) - (2*lnmusmusmc*pow(mufac,2))/(9.*mkin) - 
   (7801530877413386647*pow(mcMSmusmc,24)*NLMSOS)/(3.7763267488425775e22*pow(mkin,23)) - (73801799*lnmcmsmkin*pow(mcMSmusmc,24)*NLMSOS)/(4.23178780752e12*pow(mkin,23)) + 
   (9727*lnmusmkin*pow(mcMSmusmc,24)*NLMSOS)/(5.598936e8*pow(mkin,23)) - (11*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,24)*NLMSOS)/(86940.*pow(mkin,23)) - 
   (1640519393726677*pow(mcMSmusmc,22)*NLMSOS)/(5.946258455651274e18*pow(mkin,21)) - (14290513*lnmcmsmkin*pow(mcMSmusmc,22)*NLMSOS)/(6.89665418442e11*pow(mkin,21)) + 
   (39833*lnmusmkin*pow(mcMSmusmc,22)*NLMSOS)/(1.560329001e9*pow(mkin,21)) - (20*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,22)*NLMSOS)/(118503.*pow(mkin,21)) - 
   (22189567531163017*pow(mcMSmusmc,20)*NLMSOS)/(5.834712481735737e19*pow(mkin,19)) - (1211963*lnmcmsmkin*pow(mcMSmusmc,20)*NLMSOS)/(5.012799792e10*pow(mkin,19)) + 
   (39163*lnmusmkin*pow(mcMSmusmc,20)*NLMSOS)/(1.0015584e9*pow(mkin,19)) - (3*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,20)*NLMSOS)/(12920.*pow(mkin,19)) - 
   (20555048260909*pow(mcMSmusmc,18)*NLMSOS)/(3.76818092235585e16*pow(mkin,17)) - (197062*lnmcmsmkin*pow(mcMSmusmc,18)*NLMSOS)/(7.381208835e9*pow(mkin,17)) + 
   (16277*lnmusmkin*pow(mcMSmusmc,18)*NLMSOS)/(2.58084225e8*pow(mkin,17)) - (16*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,18)*NLMSOS)/(48195.*pow(mkin,17)) - 
   (214558103603*pow(mcMSmusmc,16)*NLMSOS)/(2.60460712512e14*pow(mkin,15)) - (355*lnmcmsmkin*pow(mcMSmusmc,16)*NLMSOS)/(1.4455584e7*pow(mkin,15)) + 
   (529*lnmusmkin*pow(mcMSmusmc,16)*NLMSOS)/(4.8672e6*pow(mkin,15)) - (7*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,16)*NLMSOS)/(14040.*pow(mkin,15)) - 
   (108352091581*pow(mcMSmusmc,14)*NLMSOS)/(8.1243243081e13*pow(mkin,13)) - (1153*lnmcmsmkin*pow(mcMSmusmc,14)*NLMSOS)/(2.25450225e8*pow(mkin,13)) + 
   (15371*lnmusmkin*pow(mcMSmusmc,14)*NLMSOS)/(7.5150075e7*pow(mkin,13)) - (4*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,14)*NLMSOS)/(5005.*pow(mkin,13)) - 
   (326802499*pow(mcMSmusmc,12)*NLMSOS)/(1.3692859488e11*pow(mkin,11)) + (391*lnmcmsmkin*pow(mcMSmusmc,12)*NLMSOS)/(4.939704e6*pow(mkin,11)) + 
   (1229*lnmusmkin*pow(mcMSmusmc,12)*NLMSOS)/(2.822688e6*pow(mkin,11)) - (5*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,12)*NLMSOS)/(3564.*pow(mkin,11)) - 
   (2816347*pow(mcMSmusmc,10)*NLMSOS)/(5.6260575e8*pow(mkin,9)) + (398*lnmcmsmkin*pow(mcMSmusmc,10)*NLMSOS)/(893025.*pow(mkin,9)) + (997*lnmusmkin*pow(mcMSmusmc,10)*NLMSOS)/(893025.*pow(mkin,9)) - 
   (8*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,10)*NLMSOS)/(2835.*pow(mkin,9)) - (174787*pow(mcMSmusmc,8)*NLMSOS)/(1.2348e7*pow(mkin,7)) + (37*lnmcmsmkin*pow(mcMSmusmc,8)*NLMSOS)/(14700.*pow(mkin,7)) + 
   (463*lnmusmkin*pow(mcMSmusmc,8)*NLMSOS)/(117600.*pow(mkin,7)) - (lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,8)*NLMSOS)/(140.*pow(mkin,7)) - (2729*pow(mcMSmusmc,6)*NLMSOS)/(30375.*pow(mkin,5)) + 
   (2*lnmcmsmkin*pow(mcMSmusmc,6)*NLMSOS)/(75.*pow(mkin,5)) + (19*lnmusmkin*pow(mcMSmusmc,6)*NLMSOS)/(675.*pow(mkin,5)) - (4*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,6)*NLMSOS)/(135.*pow(mkin,5)) + 
   (1423*pow(mcMSmusmc,4)*NLMSOS)/(3888.*pow(mkin,3)) + (lnmcmsmkin*pow(mcMSmusmc,4)*NLMSOS)/(12.*pow(mkin,3)) - (13*pow(lnmcmsmkin,2)*pow(mcMSmusmc,4)*NLMSOS)/(54.*pow(mkin,3)) + 
   (2*pow(lnmcmsmkin,3)*pow(mcMSmusmc,4)*NLMSOS)/(27.*pow(mkin,3)) - (151*lnmusmkin*pow(mcMSmusmc,4)*NLMSOS)/(648.*pow(mkin,3)) + (13*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,4)*NLMSOS)/(54.*pow(mkin,3)) - 
   (pow(lnmcmsmkin,2)*lnmusmkin*pow(mcMSmusmc,4)*NLMSOS)/(9.*pow(mkin,3)) - (2*pow(mcMSmusmc,2)*NLMSOS)/(9.*mkin) - (lnmusmkin*pow(mcMSmusmc,2)*NLMSOS)/(3.*mkin) + (241937*mkin*NLMSOS)/23328. - 
   (8*A4*mkin*NLMSOS)/27. - (pow(ln2,4)*mkin*NLMSOS)/81. + (2477*lnmusmkin*mkin*NLMSOS)/324. + (911*pow(lnmusmkin,2)*mkin*NLMSOS)/432. + (23*pow(lnmusmkin,3)*mkin*NLMSOS)/108. - (11*mufac*NLMSOS)/27. + 
   (4*lnmusmkin*mufac*NLMSOS)/81. + (4*pow(lnmusmkin,2)*mufac*NLMSOS)/27. - (11*pow(mufac,2)*NLMSOS)/(72.*mkin) + (lnmusmkin*pow(mufac,2)*NLMSOS)/(54.*mkin) + 
   (pow(lnmusmkin,2)*pow(mufac,2)*NLMSOS)/(18.*mkin) - (2353*mkin*pow(NLMSOS,2))/23328. - (89*lnmusmkin*mkin*pow(NLMSOS,2))/648. - (13*pow(lnmusmkin,2)*mkin*pow(NLMSOS,2))/216. - 
   (pow(lnmusmkin,3)*mkin*pow(NLMSOS,2))/108. - (20303*mufac*NLOSKIN)/243. + (3388*lnmufmus*mufac*NLOSKIN)/81. - (176*pow(lnmufmus,2)*mufac*NLOSKIN)/27. + (128*lnmusmcms*mufac*NLOSKIN)/243. - 
   (16*lnmufmus*lnmusmcms*mufac*NLOSKIN)/81. + (128*lnmusmkin*mufac*NLOSKIN)/81. - (16*lnmufmus*lnmusmkin*mufac*NLOSKIN)/27. - (14429*pow(mufac,2)*NLOSKIN)/(648.*mkin) + 
   (715*lnmufmus*pow(mufac,2)*NLOSKIN)/(54.*mkin) - (22*pow(lnmufmus,2)*pow(mufac,2)*NLOSKIN)/(9.*mkin) + (13*lnmusmcms*pow(mufac,2)*NLOSKIN)/(81.*mkin) - 
   (2*lnmufmus*lnmusmcms*pow(mufac,2)*NLOSKIN)/(27.*mkin) + (13*lnmusmkin*pow(mufac,2)*NLOSKIN)/(27.*mkin) - (2*lnmufmus*lnmusmkin*pow(mufac,2)*NLOSKIN)/(9.*mkin) + (1292*mufac*pow(NLOSKIN,2))/729. - 
   (256*lnmufmus*mufac*pow(NLOSKIN,2))/243. + (16*pow(lnmufmus,2)*mufac*pow(NLOSKIN,2))/81. + (209*pow(mufac,2)*pow(NLOSKIN,2))/(486.*mkin) - (26*lnmufmus*pow(mufac,2)*pow(NLOSKIN,2))/(81.*mkin) + 
   (2*pow(lnmufmus,2)*pow(mufac,2)*pow(NLOSKIN,2))/(27.*mkin) - (30053*mcMSmusmc*pow(Pi,2))/2430. + (1199*ln2*mcMSmusmc*pow(Pi,2))/81. + (31*lnmcmsmkin*mcMSmusmc*pow(Pi,2))/18. - 
   (lnmcmusmcms*mcMSmusmc*pow(Pi,2))/6. - (25*lnmusmkin*mcMSmusmc*pow(Pi,2))/36. - (106559417*pow(mcMSmusmc,25)*pow(Pi,2))/(6.374388636e11*pow(mkin,24)) + 
   (533*lnmcmsmkin*pow(mcMSmusmc,25)*pow(Pi,2))/(277725.*pow(mkin,24)) - (240817793781176357*pow(mcMSmusmc,24)*pow(Pi,2))/(6.28867544642696e20*pow(mkin,23)) + 
   (877591*ln2*pow(mcMSmusmc,24)*pow(Pi,2))/(6.79477248e9*pow(mkin,23)) + (1620816161*lnmcmsmkin*pow(mcMSmusmc,24)*pow(Pi,2))/(8.2216747008e11*pow(mkin,23)) - 
   (171475369*pow(mcMSmusmc,23)*pow(Pi,2))/(7.7816811996e11*pow(mkin,22)) + (445*lnmcmsmkin*pow(mcMSmusmc,23)*pow(Pi,2))/(192717.*pow(mkin,22)) - 
   (174200135864459*pow(mcMSmusmc,22)*pow(Pi,2))/(3.495434721951744e17*pow(mkin,21)) + (38675*ln2*pow(mcMSmusmc,22)*pow(Pi,2))/(2.33570304e8*pow(mkin,21)) + 
   (13926181*lnmcmsmkin*pow(mcMSmusmc,22)*pow(Pi,2))/(5.8392576e9*pow(mkin,21)) - (846435761*pow(mcMSmusmc,21)*pow(Pi,2))/(2.83231951884e12*pow(mkin,20)) + 
   (365*lnmcmsmkin*pow(mcMSmusmc,21)*pow(Pi,2))/(128877.*pow(mkin,20)) - (9007367733163*pow(mcMSmusmc,20)*pow(Pi,2))/(1.34810154565632e16*pow(mkin,19)) + 
   (143*ln2*pow(mcMSmusmc,20)*pow(Pi,2))/(655360.*pow(mkin,19)) + (156353*lnmcmsmkin*pow(mcMSmusmc,20)*pow(Pi,2))/(5.308416e7*pow(mkin,19)) - 
   (86836957*pow(mcMSmusmc,19)*pow(Pi,2))/(2.0687188752e11*pow(mkin,18)) + (293*lnmcmsmkin*pow(mcMSmusmc,19)*pow(Pi,2))/(82365.*pow(mkin,18)) - 
   (26463251891*pow(mcMSmusmc,18)*pow(Pi,2))/(2.845499424768e13*pow(mkin,17)) + (4147*ln2*pow(mcMSmusmc,18)*pow(Pi,2))/(1.3934592e7*pow(mkin,17)) + 
   (52013*lnmcmsmkin*pow(mcMSmusmc,18)*pow(Pi,2))/(1.3934592e7*pow(mkin,17)) - (22757641*pow(mcMSmusmc,17)*pow(Pi,2))/(3.6923796e10*pow(mkin,16)) + 
   (229*lnmcmsmkin*pow(mcMSmusmc,17)*pow(Pi,2))/(49725.*pow(mkin,16)) - (58806560951*pow(mcMSmusmc,16)*pow(Pi,2))/(4.326189170688e13*pow(mkin,15)) + 
   (1001*ln2*pow(mcMSmusmc,16)*pow(Pi,2))/(2.359296e6*pow(mkin,15)) + (565351*lnmcmsmkin*pow(mcMSmusmc,16)*pow(Pi,2))/(1.15605504e8*pow(mkin,15)) - 
   (17346493*pow(mcMSmusmc,15)*pow(Pi,2))/(1.808754948e10*pow(mkin,14)) + (173*lnmcmsmkin*pow(mcMSmusmc,15)*pow(Pi,2))/(27885.*pow(mkin,14)) - 
   (82285201*pow(mcMSmusmc,14)*pow(Pi,2))/(3.87459072e10*pow(mkin,13)) + (23*ln2*pow(mcMSmusmc,14)*pow(Pi,2))/(35840.*pow(mkin,13)) + 
   (2161*lnmcmsmkin*pow(mcMSmusmc,14)*pow(Pi,2))/(322560.*pow(mkin,13)) - (31786481*pow(mcMSmusmc,13)*pow(Pi,2))/(1.967766372e10*pow(mkin,12)) + 
   (125*lnmcmsmkin*pow(mcMSmusmc,13)*pow(Pi,2))/(14157.*pow(mkin,12)) - (84041429*pow(mcMSmusmc,12)*pow(Pi,2))/(2.29920768e10*pow(mkin,11)) + 
   (175*ln2*pow(mcMSmusmc,12)*pow(Pi,2))/(165888.*pow(mkin,11)) + (40553*lnmcmsmkin*pow(mcMSmusmc,12)*pow(Pi,2))/(4.1472e6*pow(mkin,11)) - 
   (1055689*pow(mcMSmusmc,11)*pow(Pi,2))/(3.4577928e8*pow(mkin,10)) + (85*lnmcmsmkin*pow(mcMSmusmc,11)*pow(Pi,2))/(6237.*pow(mkin,10)) - 
   (359801*pow(mcMSmusmc,10)*pow(Pi,2))/(4.89888e7*pow(mkin,9)) + (17*ln2*pow(mcMSmusmc,10)*pow(Pi,2))/(8640.*pow(mkin,9)) + (17*lnmcmsmkin*pow(mcMSmusmc,10)*pow(Pi,2))/(1080.*pow(mkin,9)) - 
   (416029*pow(mcMSmusmc,9)*pow(Pi,2))/(6.001128e7*pow(mkin,8)) + (53*lnmcmsmkin*pow(mcMSmusmc,9)*pow(Pi,2))/(2205.*pow(mkin,8)) - (1789*pow(mcMSmusmc,8)*pow(Pi,2))/(90720.*pow(mkin,7)) + 
   (7*ln2*pow(mcMSmusmc,8)*pow(Pi,2))/(1536.*pow(mkin,7)) + (139*lnmcmsmkin*pow(mcMSmusmc,8)*pow(Pi,2))/(4608.*pow(mkin,7)) - (262769*pow(mcMSmusmc,7)*pow(Pi,2))/(1.1907e7*pow(mkin,6)) + 
   (29*lnmcmsmkin*pow(mcMSmusmc,7)*pow(Pi,2))/(525.*pow(mkin,6)) - (311*pow(mcMSmusmc,6)*pow(Pi,2))/(2592.*pow(mkin,5)) + (11*ln2*pow(mcMSmusmc,6)*pow(Pi,2))/(648.*pow(mkin,5)) + 
   (113*lnmcmsmkin*pow(mcMSmusmc,6)*pow(Pi,2))/(1296.*pow(mkin,5)) - (149*pow(mcMSmusmc,5)*pow(Pi,2))/(1080.*pow(mkin,4)) + (13*lnmcmsmkin*pow(mcMSmusmc,5)*pow(Pi,2))/(45.*pow(mkin,4)) - 
   (pow(mcMSmusmc,4)*pow(Pi,2))/(12.*pow(mkin,3)) - (5*ln2*pow(mcMSmusmc,4)*pow(Pi,2))/(144.*pow(mkin,3)) - (2*pow(ln2,2)*pow(mcMSmusmc,4)*pow(Pi,2))/(81.*pow(mkin,3)) + 
   (67*lnmcmsmkin*pow(mcMSmusmc,4)*pow(Pi,2))/(144.*pow(mkin,3)) + (ln2*lnmcmsmkin*pow(mcMSmusmc,4)*pow(Pi,2))/(9.*pow(mkin,3)) - 
   (pow(lnmcmsmkin,2)*pow(mcMSmusmc,4)*pow(Pi,2))/(4.*pow(mkin,3)) + (2*lnmcmusmcms*pow(mcMSmusmc,4)*pow(Pi,2))/(9.*pow(mkin,3)) + (25*lnmusmkin*pow(mcMSmusmc,4)*pow(Pi,2))/(108.*pow(mkin,3)) - 
   (19507*pow(mcMSmusmc,3)*pow(Pi,2))/(1620.*pow(mkin,2)) + (1199*ln2*pow(mcMSmusmc,3)*pow(Pi,2))/(81.*pow(mkin,2)) + (271*lnmcmsmkin*pow(mcMSmusmc,3)*pow(Pi,2))/(162.*pow(mkin,2)) - 
   (lnmcmusmcms*pow(mcMSmusmc,3)*pow(Pi,2))/(2.*pow(mkin,2)) - (25*lnmusmkin*pow(mcMSmusmc,3)*pow(Pi,2))/(36.*pow(mkin,2)) + (13*pow(mcMSmusmc,2)*pow(Pi,2))/(12.*mkin) + 
   (3*lnmcmsmkin*pow(mcMSmusmc,2)*pow(Pi,2))/(2.*mkin) - (587741*mkin*pow(Pi,2))/38880. + (203*ln2*mkin*pow(Pi,2))/54. + (20*pow(ln2,2)*mkin*pow(Pi,2))/81. - (125*lnmusmkin*mkin*pow(Pi,2))/108. - 
   (25*ln2*lnmusmkin*mkin*pow(Pi,2))/54. - (3154*mufac*pow(Pi,2))/81. - (16*ln2*mufac*pow(Pi,2))/81. + (88*lnmufmus*mufac*pow(Pi,2))/9. + (8*lnmusmcms*mufac*pow(Pi,2))/27. + 
   (8*lnmusmkin*mufac*pow(Pi,2))/9. - (8*pow(mcMSmusmc,4)*mufac*pow(Pi,2))/(27.*pow(mkin,4)) + (16*pow(mcMSmusmc,3)*mufac*pow(Pi,2))/(27.*pow(mkin,3)) - 
   (pow(mcMSmusmc,4)*pow(mufac,2)*pow(Pi,2))/(9.*pow(mkin,5)) + (2*pow(mcMSmusmc,3)*pow(mufac,2)*pow(Pi,2))/(9.*pow(mkin,4)) - (1379*pow(mufac,2)*pow(Pi,2))/(108.*mkin) - 
   (2*ln2*pow(mufac,2)*pow(Pi,2))/(27.*mkin) + (11*lnmufmus*pow(mufac,2)*pow(Pi,2))/(3.*mkin) + (lnmusmcms*pow(mufac,2)*pow(Pi,2))/(9.*mkin) + (lnmusmkin*pow(mufac,2)*pow(Pi,2))/(3.*mkin) + 
   (7*mcMSmusmc*NLMSOS*pow(Pi,2))/27. - (2*ln2*mcMSmusmc*NLMSOS*pow(Pi,2))/9. - (lnmcmsmkin*mcMSmusmc*NLMSOS*pow(Pi,2))/9. + (lnmusmkin*mcMSmusmc*NLMSOS*pow(Pi,2))/18. + 
   (23*pow(mcMSmusmc,25)*NLMSOS*pow(Pi,2))/(415800.*pow(mkin,24)) - (11*pow(mcMSmusmc,24)*NLMSOS*pow(Pi,2))/(260820.*pow(mkin,23)) + (7*pow(mcMSmusmc,23)*NLMSOS*pow(Pi,2))/(96140.*pow(mkin,22)) - 
   (20*pow(mcMSmusmc,22)*NLMSOS*pow(Pi,2))/(355509.*pow(mkin,21)) + (19*pow(mcMSmusmc,21)*NLMSOS*pow(Pi,2))/(192780.*pow(mkin,20)) - (pow(mcMSmusmc,20)*NLMSOS*pow(Pi,2))/(12920.*pow(mkin,19)) + 
   (17*pow(mcMSmusmc,19)*NLMSOS*pow(Pi,2))/(123120.*pow(mkin,18)) - (16*pow(mcMSmusmc,18)*NLMSOS*pow(Pi,2))/(144585.*pow(mkin,17)) + (5*pow(mcMSmusmc,17)*NLMSOS*pow(Pi,2))/(24752.*pow(mkin,16)) - 
   (7*pow(mcMSmusmc,16)*NLMSOS*pow(Pi,2))/(42120.*pow(mkin,15)) + (13*pow(mcMSmusmc,15)*NLMSOS*pow(Pi,2))/(41580.*pow(mkin,14)) - (4*pow(mcMSmusmc,14)*NLMSOS*pow(Pi,2))/(15015.*pow(mkin,13)) + 
   (11*pow(mcMSmusmc,13)*NLMSOS*pow(Pi,2))/(21060.*pow(mkin,12)) - (5*pow(mcMSmusmc,12)*NLMSOS*pow(Pi,2))/(10692.*pow(mkin,11)) + (3*pow(mcMSmusmc,11)*NLMSOS*pow(Pi,2))/(3080.*pow(mkin,10)) - 
   (8*pow(mcMSmusmc,10)*NLMSOS*pow(Pi,2))/(8505.*pow(mkin,9)) + (7*pow(mcMSmusmc,9)*NLMSOS*pow(Pi,2))/(3240.*pow(mkin,8)) - (pow(mcMSmusmc,8)*NLMSOS*pow(Pi,2))/(420.*pow(mkin,7)) + 
   (5*pow(mcMSmusmc,7)*NLMSOS*pow(Pi,2))/(756.*pow(mkin,6)) - (4*pow(mcMSmusmc,6)*NLMSOS*pow(Pi,2))/(405.*pow(mkin,5)) + (pow(mcMSmusmc,5)*NLMSOS*pow(Pi,2))/(20.*pow(mkin,4)) + 
   (13*pow(mcMSmusmc,4)*NLMSOS*pow(Pi,2))/(324.*pow(mkin,3)) - (lnmcmsmkin*pow(mcMSmusmc,4)*NLMSOS*pow(Pi,2))/(27.*pow(mkin,3)) - (lnmusmkin*pow(mcMSmusmc,4)*NLMSOS*pow(Pi,2))/(54.*pow(mkin,3)) + 
   (7*pow(mcMSmusmc,3)*NLMSOS*pow(Pi,2))/(54.*pow(mkin,2)) - (2*ln2*pow(mcMSmusmc,3)*NLMSOS*pow(Pi,2))/(9.*pow(mkin,2)) - (lnmcmsmkin*pow(mcMSmusmc,3)*NLMSOS*pow(Pi,2))/(9.*pow(mkin,2)) + 
   (lnmusmkin*pow(mcMSmusmc,3)*NLMSOS*pow(Pi,2))/(18.*pow(mkin,2)) + (305*mkin*NLMSOS*pow(Pi,2))/216. + (11*ln2*mkin*NLMSOS*pow(Pi,2))/81. - (2*pow(ln2,2)*mkin*NLMSOS*pow(Pi,2))/81. + 
   (35*lnmusmkin*mkin*NLMSOS*pow(Pi,2))/108. + (ln2*lnmusmkin*mkin*NLMSOS*pow(Pi,2))/27. + (8*mufac*NLMSOS*pow(Pi,2))/81. + (pow(mufac,2)*NLMSOS*pow(Pi,2))/(27.*mkin) - 
   (13*mkin*pow(NLMSOS,2)*pow(Pi,2))/324. - (lnmusmkin*mkin*pow(NLMSOS,2)*pow(Pi,2))/54. + (208*mufac*NLOSKIN*pow(Pi,2))/81. - (16*lnmufmus*mufac*NLOSKIN*pow(Pi,2))/27. + 
   (23*pow(mufac,2)*NLOSKIN*pow(Pi,2))/(27.*mkin) - (2*lnmufmus*pow(mufac,2)*NLOSKIN*pow(Pi,2))/(9.*mkin) - (8*mufac*pow(NLOSKIN,2)*pow(Pi,2))/243. - 
   (pow(mufac,2)*pow(NLOSKIN,2)*pow(Pi,2))/(81.*mkin) - (13*mcMSmusmc*pow(Pi,3))/162. + (5698043*pow(mcMSmusmc,25)*pow(Pi,3))/(1.189085184e11*pow(mkin,24)) + 
   (81991*pow(mcMSmusmc,23)*pow(Pi,3))/(1.374683136e9*pow(mkin,22)) + (127699*pow(mcMSmusmc,21)*pow(Pi,3))/(1.684537344e9*pow(mkin,20)) + 
   (99671*pow(mcMSmusmc,19)*pow(Pi,3))/(1.00859904e9*pow(mkin,18)) + (1925*pow(mcMSmusmc,17)*pow(Pi,3))/(1.4483456e7*pow(mkin,16)) + (377*pow(mcMSmusmc,15)*pow(Pi,3))/(2.02752e6*pow(mkin,14)) + 
   (1771*pow(mcMSmusmc,13)*pow(Pi,3))/(6.469632e6*pow(mkin,12)) + (17*pow(mcMSmusmc,11)*pow(Pi,3))/(39424.*pow(mkin,10)) + (77*pow(mcMSmusmc,9)*pow(Pi,3))/(103680.*pow(mkin,8)) + 
   (25*pow(mcMSmusmc,7)*pow(Pi,3))/(18144.*pow(mkin,6)) - (pow(mcMSmusmc,5)*pow(Pi,3))/(240.*pow(mkin,4)) - (7*pow(mcMSmusmc,3)*pow(Pi,3))/(108.*pow(mkin,2)) - 
   (271*pow(mcMSmusmc,4)*pow(Pi,4))/(19440.*pow(mkin,3)) - (lnmcmsmkin*pow(mcMSmusmc,2)*pow(Pi,4))/(12.*mkin) + (451*mkin*pow(Pi,4))/7776. + (2*mufac*pow(Pi,4))/3. + 
   (pow(mufac,2)*pow(Pi,4))/(4.*mkin) - (61*mkin*NLMSOS*pow(Pi,4))/1944. + (2710689767*pow(mcMSmusmc,24)*Zeta3)/(1.64433494016e12*pow(mkin,23)) + 
   (23017987*pow(mcMSmusmc,22)*Zeta3)/(1.16785152e10*pow(mkin,21)) + (254791*pow(mcMSmusmc,20)*Zeta3)/(1.0616832e8*pow(mkin,19)) + (83291*pow(mcMSmusmc,18)*Zeta3)/(2.7869184e7*pow(mkin,17)) + 
   (885457*pow(mcMSmusmc,16)*Zeta3)/(2.31211008e8*pow(mkin,15)) + (3287*pow(mcMSmusmc,14)*Zeta3)/(645120.*pow(mkin,13)) + (59231*pow(mcMSmusmc,12)*Zeta3)/(8.2944e6*pow(mkin,11)) + 
   (187*pow(mcMSmusmc,10)*Zeta3)/(17280.*pow(mkin,9)) + (173*pow(mcMSmusmc,8)*Zeta3)/(9216.*pow(mkin,7)) + (29*pow(mcMSmusmc,6)*Zeta3)/(648.*pow(mkin,5)) + 
   (2309*pow(mcMSmusmc,4)*Zeta3)/(864.*pow(mkin,3)) - (14*lnmcmsmkin*pow(mcMSmusmc,4)*Zeta3)/(9.*pow(mkin,3)) + (11*pow(mcMSmusmc,2)*Zeta3)/(2.*mkin) + (23*mkin*Zeta3)/24. + 
   (85*lnmusmkin*mkin*Zeta3)/36. - (3070*mufac*Zeta3)/27. - (1535*pow(mufac,2)*Zeta3)/(36.*mkin) - (pow(mcMSmusmc,4)*NLMSOS*Zeta3)/(3.*pow(mkin,3)) + (667*mkin*NLMSOS*Zeta3)/216. + 
   (7*lnmusmkin*mkin*NLMSOS*Zeta3)/9. - (7*mkin*pow(NLMSOS,2)*Zeta3)/54. + (140*mufac*NLOSKIN*Zeta3)/27. + (35*pow(mufac,2)*NLOSKIN*Zeta3)/(18.*mkin) - (3*pow(mcMSmusmc,2)*pow(Pi,2)*Zeta3)/(4.*mkin) + 
   (1439*mkin*pow(Pi,2)*Zeta3)/432. - (5*pow(mcMSmusmc,2)*Zeta5)/(2.*mkin) - (1975*mkin*Zeta5)/216.);
    
    
    double mMS = 0.0;
    for(int i = 0; i <= nloops; i++) {
       mMS += ret[i];
    }    
    return mMS;
}

double CRunDec::mkin2mMSC(double mkin, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops) {
    if(nloops<0||nloops>3){
      cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
      RETURN
    }
    double ret[5];
    double lnmusmkin = log(mus*mus/(mkin*mkin));
    double lnmufmus = log(2*mufac/mus);
    double lnmusmcms = 0.0;
    double lnmcmsmkin = 0.0;
    double lnmcmusmcms = 0.0;  
    double mcMSmusmc = 0.0;
    double musmc = 0.0;
    if(mcMSmusmcin > 0.0) {
       lnmusmcms = log(mus*mus/(mcMSmusmcin*mcMSmusmcin));
       lnmcmsmkin = log(mcMSmusmcin/mkin);
       lnmcmusmcms = log(musmcin*musmcin/(mcMSmusmcin*mcMSmusmcin));
       mcMSmusmc = mcMSmusmcin;
       musmc = musmcin;
    }
    
    ret[0] = mkin;
    ret[1] = apinlmus*((-4*mkin)/3. - lnmusmkin*mkin + (16*mufac)/9. + (2*pow(mufac,2))/(3.*mkin));
    ret[2] = apinlmus*apinlmus*((-9727*pow(mcMSmusmc,24))/(1.866312e8*pow(mkin,23)) + (11*lnmcmsmkin*pow(mcMSmusmc,24))/(28980.*pow(mkin,23)) - (39833*pow(mcMSmusmc,22))/(5.20109667e8*pow(mkin,21)) + 
   (20*lnmcmsmkin*pow(mcMSmusmc,22))/(39501.*pow(mkin,21)) - (39163*pow(mcMSmusmc,20))/(3.338528e8*pow(mkin,19)) + (9*lnmcmsmkin*pow(mcMSmusmc,20))/(12920.*pow(mkin,19)) - 
   (16277*pow(mcMSmusmc,18))/(8.6028075e7*pow(mkin,17)) + (16*lnmcmsmkin*pow(mcMSmusmc,18))/(16065.*pow(mkin,17)) - (529*pow(mcMSmusmc,16))/(1.6224e6*pow(mkin,15)) + 
   (7*lnmcmsmkin*pow(mcMSmusmc,16))/(4680.*pow(mkin,15)) - (15371*pow(mcMSmusmc,14))/(2.5050025e7*pow(mkin,13)) + (12*lnmcmsmkin*pow(mcMSmusmc,14))/(5005.*pow(mkin,13)) - 
   (1229*pow(mcMSmusmc,12))/(940896.*pow(mkin,11)) + (5*lnmcmsmkin*pow(mcMSmusmc,12))/(1188.*pow(mkin,11)) - (997*pow(mcMSmusmc,10))/(297675.*pow(mkin,9)) + 
   (8*lnmcmsmkin*pow(mcMSmusmc,10))/(945.*pow(mkin,9)) - (463*pow(mcMSmusmc,8))/(39200.*pow(mkin,7)) + (3*lnmcmsmkin*pow(mcMSmusmc,8))/(140.*pow(mkin,7)) - (19*pow(mcMSmusmc,6))/(225.*pow(mkin,5)) + 
   (4*lnmcmsmkin*pow(mcMSmusmc,6))/(45.*pow(mkin,5)) + (151*pow(mcMSmusmc,4))/(216.*pow(mkin,3)) - (13*lnmcmsmkin*pow(mcMSmusmc,4))/(18.*pow(mkin,3)) + 
   (pow(lnmcmsmkin,2)*pow(mcMSmusmc,4))/(3.*pow(mkin,3)) + pow(mcMSmusmc,2)/mkin - (959*mkin)/96. - (145*lnmusmkin*mkin)/24. - (7*pow(lnmusmkin,2)*mkin)/8. + (892*mufac)/27. - (88*lnmufmus*mufac)/9. - 
   (16*lnmusmkin*mufac)/9. + (95*pow(mufac,2))/(9.*mkin) - (11*lnmufmus*pow(mufac,2))/(3.*mkin) - (2*lnmusmkin*pow(mufac,2))/(3.*mkin) + (71*mkin*NLMSOS)/144. + (13*lnmusmkin*mkin*NLMSOS)/36. + 
   (pow(lnmusmkin,2)*mkin*NLMSOS)/12. - (128*mufac*NLOSKIN)/81. + (16*lnmufmus*mufac*NLOSKIN)/27. - (13*pow(mufac,2)*NLOSKIN)/(27.*mkin) + (2*lnmufmus*pow(mufac,2)*NLOSKIN)/(9.*mkin) - 
   (mcMSmusmc*pow(Pi,2))/6. + (pow(mcMSmusmc,4)*pow(Pi,2))/(18.*pow(mkin,3)) - (pow(mcMSmusmc,3)*pow(Pi,2))/(6.*pow(mkin,2)) - (5*mkin*pow(Pi,2))/18. - (ln2*mkin*pow(Pi,2))/9. - 
   (8*mufac*pow(Pi,2))/9. - (pow(mufac,2)*pow(Pi,2))/(3.*mkin) + (mkin*NLMSOS*pow(Pi,2))/18. + (mkin*Zeta3)/6.);
    ret[3] = apinlmus*apinlmus*apinlmus*((-8439929370090150606400742890600889.*pow(mcMSmusmc,24))/(2.369996943791664e36*pow(mkin,23)) + (260257620268030848013.*lnmcmsmkin*pow(mcMSmusmc,24))/(2.8097278338127475e22*pow(mkin,23)) + 
   (953082228281*pow(lnmcmsmkin,2)*pow(mcMSmusmc,24))/(4.664604200256e13*pow(mkin,23)) - (10163*lnmcmusmcms*pow(mcMSmusmc,24))/(1.166445e7*pow(mkin,23)) + 
   (22*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,24))/(2415.*pow(mkin,23)) - (9727*lnmusmkin*pow(mcMSmusmc,24))/(4.4791488e7*pow(mkin,23)) + 
   (55*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,24))/(34776.*pow(mkin,23)) - (47848557349149553062713793193.*pow(mcMSmusmc,22))/(1.0790864167779053e31*pow(mkin,21)) + 
   (79809944445306329123.*lnmcmsmkin*pow(mcMSmusmc,22))/(7.242811051244408e21*pow(mkin,21)) + (3078965960711*pow(lnmcmsmkin,2)*pow(mcMSmusmc,22))/(1.24450902576e14*pow(mkin,21)) - 
   (5066*lnmcmusmcms*pow(mcMSmusmc,22))/(4.298427e6*pow(mkin,21)) + (40*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,22))/(3591.*pow(mkin,21)) - 
   (995825*lnmusmkin*pow(mcMSmusmc,22))/(3.120658002e9*pow(mkin,21)) + (250*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,22))/(118503.*pow(mkin,21)) - 
   (308474675741863370556669649.*pow(mcMSmusmc,20))/(5.449931397868209e28*pow(mkin,19)) + (2335763394163431811.*lnmcmsmkin*pow(mcMSmusmc,20))/(1.7558329821198565e20*pow(mkin,19)) + 
   (46137065941*pow(lnmcmsmkin,2)*pow(mcMSmusmc,20))/(1.5084957888e12*pow(mkin,19)) - (5507*lnmcmusmcms*pow(mcMSmusmc,20))/(3.338528e6*pow(mkin,19)) + 
   (9*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,20))/(646.*pow(mkin,19)) - (39163*lnmusmkin*pow(mcMSmusmc,20))/(8.0124672e7*pow(mkin,19)) + (15*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,20))/(5168.*pow(mkin,19)) - 
   (372869427163814153380429.*pow(mcMSmusmc,18))/(5.00576874275692e25*pow(mkin,17)) + (53349597791833939*lnmcmsmkin*pow(mcMSmusmc,18))/(3.2684758005112013e18*pow(mkin,17)) + 
   (1014809351*pow(lnmcmsmkin,2)*pow(mcMSmusmc,18))/(2.615348736e10*pow(mkin,17)) - (7678*lnmcmusmcms*pow(mcMSmusmc,18))/(3.186225e6*pow(mkin,17)) + 
   (32*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,18))/(1785.*pow(mkin,17)) - (16277*lnmusmkin*pow(mcMSmusmc,18))/(2.0646738e7*pow(mkin,17)) + 
   (40*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,18))/(9639.*pow(mkin,17)) - (968011304335065184487.*pow(mcMSmusmc,16))/(9.509568138760686e22*pow(mkin,15)) + 
   (23986410569717*lnmcmsmkin*pow(mcMSmusmc,16))/(1.1780838381312e15*pow(mkin,15)) + (232939907*pow(lnmcmsmkin,2)*pow(mcMSmusmc,16))/(4.576860288e9*pow(mkin,15)) - 
   (283*lnmcmusmcms*pow(mcMSmusmc,16))/(76050.*pow(mkin,15)) + (14*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,16))/(585.*pow(mkin,15)) - (529*lnmusmkin*pow(mcMSmusmc,16))/(389376.*pow(mkin,15)) + 
   (35*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,16))/(5616.*pow(mkin,15)) - (858973162016800669*pow(mcMSmusmc,14))/(5.896309609846656e19*pow(mkin,13)) + 
   (95564969831*lnmcmsmkin*pow(mcMSmusmc,14))/(3.718698984e12*pow(mkin,13)) + (845153*pow(lnmcmsmkin,2)*pow(mcMSmusmc,14))/(1.2108096e7*pow(mkin,13)) - 
   (3166*lnmcmusmcms*pow(mcMSmusmc,14))/(511225.*pow(mkin,13)) + (24*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,14))/(715.*pow(mkin,13)) - (15371*lnmusmkin*pow(mcMSmusmc,14))/(6.012006e6*pow(mkin,13)) + 
   (10*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,14))/(1001.*pow(mkin,13)) - (540989122348457*pow(mcMSmusmc,12))/(2.46471470784e16*pow(mkin,11)) + 
   (199138925387*lnmcmsmkin*pow(mcMSmusmc,12))/(6.22402704e12*pow(mkin,11)) + (15261023*pow(lnmcmsmkin,2)*pow(mcMSmusmc,12))/(1.49688e8*pow(mkin,11)) - 
   (899*lnmcmusmcms*pow(mcMSmusmc,12))/(78408.*pow(mkin,11)) + (5*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,12))/(99.*pow(mkin,11)) - (30725*lnmusmkin*pow(mcMSmusmc,12))/(5.645376e6*pow(mkin,11)) + 
   (125*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,12))/(7128.*pow(mkin,11)) - (174878652127*pow(mcMSmusmc,10))/(5.184974592e12*pow(mkin,9)) + 
   (1395062351*lnmcmsmkin*pow(mcMSmusmc,10))/(4.1150592e10*pow(mkin,9)) + (382589*pow(lnmcmsmkin,2)*pow(mcMSmusmc,10))/(2.3328e6*pow(mkin,9)) - 
   (298*lnmcmusmcms*pow(mcMSmusmc,10))/(11907.*pow(mkin,9)) + (16*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,10))/(189.*pow(mkin,9)) - (997*lnmusmkin*pow(mcMSmusmc,10))/(71442.*pow(mkin,9)) + 
   (20*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,10))/(567.*pow(mkin,9)) - (32093279*pow(mcMSmusmc,8))/(1.185408e9*pow(mkin,7)) - (253711*lnmcmsmkin*pow(mcMSmusmc,8))/(1.27008e7*pow(mkin,7)) + 
   (9451*pow(lnmcmsmkin,2)*pow(mcMSmusmc,8))/(30240.*pow(mkin,7)) - (179*lnmcmusmcms*pow(mcMSmusmc,8))/(2450.*pow(mkin,7)) + (6*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,8))/(35.*pow(mkin,7)) - 
   (463*lnmusmkin*pow(mcMSmusmc,8))/(9408.*pow(mkin,7)) + (5*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,8))/(56.*pow(mkin,7)) + (10893161*pow(mcMSmusmc,6))/(1.7496e7*pow(mkin,5)) - 
   (169129*lnmcmsmkin*pow(mcMSmusmc,6))/(233280.*pow(mkin,5)) + (5329*pow(lnmcmsmkin,2)*pow(mcMSmusmc,6))/(6480.*pow(mkin,5)) - (94*lnmcmusmcms*pow(mcMSmusmc,6))/(225.*pow(mkin,5)) + 
   (8*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,6))/(15.*pow(mkin,5)) - (19*lnmusmkin*pow(mcMSmusmc,6))/(54.*pow(mkin,5)) + (10*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,6))/(27.*pow(mkin,5)) + 
   (10103*pow(mcMSmusmc,4))/(7776.*pow(mkin,3)) - (20*A4*pow(mcMSmusmc,4))/(27.*pow(mkin,3)) - (5*pow(ln2,4)*pow(mcMSmusmc,4))/(162.*pow(mkin,3)) - 
   (4051*lnmcmsmkin*pow(mcMSmusmc,4))/(648.*pow(mkin,3)) + (136*pow(lnmcmsmkin,2)*pow(mcMSmusmc,4))/(27.*pow(mkin,3)) - (7*pow(lnmcmsmkin,3)*pow(mcMSmusmc,4))/(27.*pow(mkin,3)) + 
   (56*lnmcmusmcms*pow(mcMSmusmc,4))/(27.*pow(mkin,3)) - (20*lnmcmsmkin*lnmcmusmcms*pow(mcMSmusmc,4))/(9.*pow(mkin,3)) + (4*pow(lnmcmsmkin,2)*lnmcmusmcms*pow(mcMSmusmc,4))/(3.*pow(mkin,3)) + 
   (3775*lnmusmkin*pow(mcMSmusmc,4))/(1296.*pow(mkin,3)) - (325*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,4))/(108.*pow(mkin,3)) + (25*pow(lnmcmsmkin,2)*lnmusmkin*pow(mcMSmusmc,4))/(18.*pow(mkin,3)) + 
   (86*pow(mcMSmusmc,2))/(9.*mkin) + (4*lnmcmsmkin*pow(mcMSmusmc,2))/mkin + (2*lnmcmusmcms*pow(mcMSmusmc,2))/mkin + (25*lnmusmkin*pow(mcMSmusmc,2))/(6.*mkin) - (8537461*mkin)/93312. + (212*A4*mkin)/27. + 
   (53*pow(ln2,4)*mkin)/162. - (39503*lnmusmkin*mkin)/648. - (12401*pow(lnmusmkin,2)*mkin)/864. - (505*pow(lnmusmkin,3)*mkin)/432. + (7495*mufac)/9. - (3416*lnmufmus*mufac)/9. + 
   (484*pow(lnmufmus,2)*mufac)/9. - (982*lnmusmkin*mufac)/27. + (88*lnmufmus*lnmusmkin*mufac)/9. - (14*pow(lnmusmkin,2)*mufac)/9. + (289*pow(mcMSmusmc,24)*mufac)/(198450.*pow(mkin,24)) - 
   (44*lnmcmsmkin*pow(mcMSmusmc,24)*mufac)/(2835.*pow(mkin,24)) + (62384*pow(mcMSmusmc,22)*mufac)/(3.1843449e7*pow(mkin,22)) - (320*lnmcmsmkin*pow(mcMSmusmc,22)*mufac)/(16929.*pow(mkin,22)) + 
   (1417*pow(mcMSmusmc,20)*mufac)/(520200.*pow(mkin,20)) - (2*lnmcmsmkin*pow(mcMSmusmc,20)*mufac)/(85.*pow(mkin,20)) + (10576*pow(mcMSmusmc,18)*mufac)/(2.679075e6*pow(mkin,18)) - 
   (256*lnmcmsmkin*pow(mcMSmusmc,18)*mufac)/(8505.*pow(mkin,18)) + (661*pow(mcMSmusmc,16)*mufac)/(109512.*pow(mkin,16)) - (14*lnmcmsmkin*pow(mcMSmusmc,16)*mufac)/(351.*pow(mkin,16)) + 
   (13232*pow(mcMSmusmc,14)*mufac)/(1.334025e6*pow(mkin,14)) - (64*lnmcmsmkin*pow(mcMSmusmc,14)*mufac)/(1155.*pow(mkin,14)) + (79*pow(mcMSmusmc,12)*mufac)/(4374.*pow(mkin,12)) - 
   (20*lnmcmsmkin*pow(mcMSmusmc,12)*mufac)/(243.*pow(mkin,12)) + (3824*pow(mcMSmusmc,10)*mufac)/(99225.*pow(mkin,10)) - (128*lnmcmsmkin*pow(mcMSmusmc,10)*mufac)/(945.*pow(mkin,10)) + 
   (49*pow(mcMSmusmc,8)*mufac)/(450.*pow(mkin,8)) - (4*lnmcmsmkin*pow(mcMSmusmc,8)*mufac)/(15.*pow(mkin,8)) + (16*pow(mcMSmusmc,6)*mufac)/(27.*pow(mkin,6)) - 
   (64*lnmcmsmkin*pow(mcMSmusmc,6)*mufac)/(81.*pow(mkin,6)) - (22*pow(mcMSmusmc,4)*mufac)/(9.*pow(mkin,4)) + (8*lnmcmsmkin*pow(mcMSmusmc,4)*mufac)/(3.*pow(mkin,4)) - 
   (16*pow(lnmcmsmkin,2)*pow(mcMSmusmc,4)*mufac)/(9.*pow(mkin,4)) - (16*pow(mcMSmusmc,2)*mufac)/(9.*pow(mkin,2)) + (289*pow(mcMSmusmc,24)*pow(mufac,2))/(529200.*pow(mkin,25)) - 
   (11*lnmcmsmkin*pow(mcMSmusmc,24)*pow(mufac,2))/(1890.*pow(mkin,25)) + (7798*pow(mcMSmusmc,22)*pow(mufac,2))/(1.0614483e7*pow(mkin,23)) - 
   (40*lnmcmsmkin*pow(mcMSmusmc,22)*pow(mufac,2))/(5643.*pow(mkin,23)) + (1417*pow(mcMSmusmc,20)*pow(mufac,2))/(1.3872e6*pow(mkin,21)) - 
   (3*lnmcmsmkin*pow(mcMSmusmc,20)*pow(mufac,2))/(340.*pow(mkin,21)) + (1322*pow(mcMSmusmc,18)*pow(mufac,2))/(893025.*pow(mkin,19)) - 
   (32*lnmcmsmkin*pow(mcMSmusmc,18)*pow(mufac,2))/(2835.*pow(mkin,19)) + (661*pow(mcMSmusmc,16)*pow(mufac,2))/(292032.*pow(mkin,17)) - 
   (7*lnmcmsmkin*pow(mcMSmusmc,16)*pow(mufac,2))/(468.*pow(mkin,17)) + (1654*pow(mcMSmusmc,14)*pow(mufac,2))/(444675.*pow(mkin,15)) - 
   (8*lnmcmsmkin*pow(mcMSmusmc,14)*pow(mufac,2))/(385.*pow(mkin,15)) + (79*pow(mcMSmusmc,12)*pow(mufac,2))/(11664.*pow(mkin,13)) - 
   (5*lnmcmsmkin*pow(mcMSmusmc,12)*pow(mufac,2))/(162.*pow(mkin,13)) + (478*pow(mcMSmusmc,10)*pow(mufac,2))/(33075.*pow(mkin,11)) - 
   (16*lnmcmsmkin*pow(mcMSmusmc,10)*pow(mufac,2))/(315.*pow(mkin,11)) + (49*pow(mcMSmusmc,8)*pow(mufac,2))/(1200.*pow(mkin,9)) - (lnmcmsmkin*pow(mcMSmusmc,8)*pow(mufac,2))/(10.*pow(mkin,9)) + 
   (2*pow(mcMSmusmc,6)*pow(mufac,2))/(9.*pow(mkin,7)) - (8*lnmcmsmkin*pow(mcMSmusmc,6)*pow(mufac,2))/(27.*pow(mkin,7)) - (11*pow(mcMSmusmc,4)*pow(mufac,2))/(12.*pow(mkin,5)) + 
   (lnmcmsmkin*pow(mcMSmusmc,4)*pow(mufac,2))/pow(mkin,5) - (2*pow(lnmcmsmkin,2)*pow(mcMSmusmc,4)*pow(mufac,2))/(3.*pow(mkin,5)) - (2*pow(mcMSmusmc,2)*pow(mufac,2))/(3.*pow(mkin,3)) + 
   (151763*pow(mufac,2))/(648.*mkin) - (733*lnmufmus*pow(mufac,2))/(6.*mkin) + (121*pow(lnmufmus,2)*pow(mufac,2))/(6.*mkin) - (425*lnmusmkin*pow(mufac,2))/(36.*mkin) + 
   (11*lnmufmus*lnmusmkin*pow(mufac,2))/(3.*mkin) - (7*pow(lnmusmkin,2)*pow(mufac,2))/(12.*mkin) - (7801530877413386647*pow(mcMSmusmc,24)*NLMSOS)/(3.7763267488425775e22*pow(mkin,23)) - 
   (73801799*lnmcmsmkin*pow(mcMSmusmc,24)*NLMSOS)/(4.23178780752e12*pow(mkin,23)) + (9727*lnmusmkin*pow(mcMSmusmc,24)*NLMSOS)/(5.598936e8*pow(mkin,23)) - 
   (11*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,24)*NLMSOS)/(86940.*pow(mkin,23)) - (1640519393726677*pow(mcMSmusmc,22)*NLMSOS)/(5.946258455651274e18*pow(mkin,21)) - 
   (14290513*lnmcmsmkin*pow(mcMSmusmc,22)*NLMSOS)/(6.89665418442e11*pow(mkin,21)) + (39833*lnmusmkin*pow(mcMSmusmc,22)*NLMSOS)/(1.560329001e9*pow(mkin,21)) - 
   (20*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,22)*NLMSOS)/(118503.*pow(mkin,21)) - (22189567531163017*pow(mcMSmusmc,20)*NLMSOS)/(5.834712481735737e19*pow(mkin,19)) - 
   (1211963*lnmcmsmkin*pow(mcMSmusmc,20)*NLMSOS)/(5.012799792e10*pow(mkin,19)) + (39163*lnmusmkin*pow(mcMSmusmc,20)*NLMSOS)/(1.0015584e9*pow(mkin,19)) - 
   (3*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,20)*NLMSOS)/(12920.*pow(mkin,19)) - (20555048260909*pow(mcMSmusmc,18)*NLMSOS)/(3.76818092235585e16*pow(mkin,17)) - 
   (197062*lnmcmsmkin*pow(mcMSmusmc,18)*NLMSOS)/(7.381208835e9*pow(mkin,17)) + (16277*lnmusmkin*pow(mcMSmusmc,18)*NLMSOS)/(2.58084225e8*pow(mkin,17)) - 
   (16*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,18)*NLMSOS)/(48195.*pow(mkin,17)) - (214558103603*pow(mcMSmusmc,16)*NLMSOS)/(2.60460712512e14*pow(mkin,15)) - 
   (355*lnmcmsmkin*pow(mcMSmusmc,16)*NLMSOS)/(1.4455584e7*pow(mkin,15)) + (529*lnmusmkin*pow(mcMSmusmc,16)*NLMSOS)/(4.8672e6*pow(mkin,15)) - 
   (7*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,16)*NLMSOS)/(14040.*pow(mkin,15)) - (108352091581*pow(mcMSmusmc,14)*NLMSOS)/(8.1243243081e13*pow(mkin,13)) - 
   (1153*lnmcmsmkin*pow(mcMSmusmc,14)*NLMSOS)/(2.25450225e8*pow(mkin,13)) + (15371*lnmusmkin*pow(mcMSmusmc,14)*NLMSOS)/(7.5150075e7*pow(mkin,13)) - 
   (4*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,14)*NLMSOS)/(5005.*pow(mkin,13)) - (326802499*pow(mcMSmusmc,12)*NLMSOS)/(1.3692859488e11*pow(mkin,11)) + 
   (391*lnmcmsmkin*pow(mcMSmusmc,12)*NLMSOS)/(4.939704e6*pow(mkin,11)) + (1229*lnmusmkin*pow(mcMSmusmc,12)*NLMSOS)/(2.822688e6*pow(mkin,11)) - 
   (5*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,12)*NLMSOS)/(3564.*pow(mkin,11)) - (2816347*pow(mcMSmusmc,10)*NLMSOS)/(5.6260575e8*pow(mkin,9)) + 
   (398*lnmcmsmkin*pow(mcMSmusmc,10)*NLMSOS)/(893025.*pow(mkin,9)) + (997*lnmusmkin*pow(mcMSmusmc,10)*NLMSOS)/(893025.*pow(mkin,9)) - 
   (8*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,10)*NLMSOS)/(2835.*pow(mkin,9)) - (174787*pow(mcMSmusmc,8)*NLMSOS)/(1.2348e7*pow(mkin,7)) + (37*lnmcmsmkin*pow(mcMSmusmc,8)*NLMSOS)/(14700.*pow(mkin,7)) + 
   (463*lnmusmkin*pow(mcMSmusmc,8)*NLMSOS)/(117600.*pow(mkin,7)) - (lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,8)*NLMSOS)/(140.*pow(mkin,7)) - (2729*pow(mcMSmusmc,6)*NLMSOS)/(30375.*pow(mkin,5)) + 
   (2*lnmcmsmkin*pow(mcMSmusmc,6)*NLMSOS)/(75.*pow(mkin,5)) + (19*lnmusmkin*pow(mcMSmusmc,6)*NLMSOS)/(675.*pow(mkin,5)) - (4*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,6)*NLMSOS)/(135.*pow(mkin,5)) + 
   (1423*pow(mcMSmusmc,4)*NLMSOS)/(3888.*pow(mkin,3)) + (lnmcmsmkin*pow(mcMSmusmc,4)*NLMSOS)/(12.*pow(mkin,3)) - (13*pow(lnmcmsmkin,2)*pow(mcMSmusmc,4)*NLMSOS)/(54.*pow(mkin,3)) + 
   (2*pow(lnmcmsmkin,3)*pow(mcMSmusmc,4)*NLMSOS)/(27.*pow(mkin,3)) - (151*lnmusmkin*pow(mcMSmusmc,4)*NLMSOS)/(648.*pow(mkin,3)) + (13*lnmcmsmkin*lnmusmkin*pow(mcMSmusmc,4)*NLMSOS)/(54.*pow(mkin,3)) - 
   (pow(lnmcmsmkin,2)*lnmusmkin*pow(mcMSmusmc,4)*NLMSOS)/(9.*pow(mkin,3)) - (2*pow(mcMSmusmc,2)*NLMSOS)/(9.*mkin) - (lnmusmkin*pow(mcMSmusmc,2)*NLMSOS)/(3.*mkin) + (241937*mkin*NLMSOS)/23328. - 
   (8*A4*mkin*NLMSOS)/27. - (pow(ln2,4)*mkin*NLMSOS)/81. + (2477*lnmusmkin*mkin*NLMSOS)/324. + (911*pow(lnmusmkin,2)*mkin*NLMSOS)/432. + (23*pow(lnmusmkin,3)*mkin*NLMSOS)/108. - (11*mufac*NLMSOS)/27. + 
   (4*lnmusmkin*mufac*NLMSOS)/81. + (4*pow(lnmusmkin,2)*mufac*NLMSOS)/27. - (11*pow(mufac,2)*NLMSOS)/(72.*mkin) + (lnmusmkin*pow(mufac,2)*NLMSOS)/(54.*mkin) + 
   (pow(lnmusmkin,2)*pow(mufac,2)*NLMSOS)/(18.*mkin) - (2353*mkin*pow(NLMSOS,2))/23328. - (89*lnmusmkin*mkin*pow(NLMSOS,2))/648. - (13*pow(lnmusmkin,2)*mkin*pow(NLMSOS,2))/216. - 
   (pow(lnmusmkin,3)*mkin*pow(NLMSOS,2))/108. - (20303*mufac*NLOSKIN)/243. + (3388*lnmufmus*mufac*NLOSKIN)/81. - (176*pow(lnmufmus,2)*mufac*NLOSKIN)/27. + (128*lnmusmkin*mufac*NLOSKIN)/81. - 
   (16*lnmufmus*lnmusmkin*mufac*NLOSKIN)/27. - (14429*pow(mufac,2)*NLOSKIN)/(648.*mkin) + (715*lnmufmus*pow(mufac,2)*NLOSKIN)/(54.*mkin) - (22*pow(lnmufmus,2)*pow(mufac,2)*NLOSKIN)/(9.*mkin) + 
   (13*lnmusmkin*pow(mufac,2)*NLOSKIN)/(27.*mkin) - (2*lnmufmus*lnmusmkin*pow(mufac,2)*NLOSKIN)/(9.*mkin) + (1292*mufac*pow(NLOSKIN,2))/729. - (256*lnmufmus*mufac*pow(NLOSKIN,2))/243. + 
   (16*pow(lnmufmus,2)*mufac*pow(NLOSKIN,2))/81. + (209*pow(mufac,2)*pow(NLOSKIN,2))/(486.*mkin) - (26*lnmufmus*pow(mufac,2)*pow(NLOSKIN,2))/(81.*mkin) + 
   (2*pow(lnmufmus,2)*pow(mufac,2)*pow(NLOSKIN,2))/(27.*mkin) - (30053*mcMSmusmc*pow(Pi,2))/2430. + (1199*ln2*mcMSmusmc*pow(Pi,2))/81. + (31*lnmcmsmkin*mcMSmusmc*pow(Pi,2))/18. - 
   (lnmcmusmcms*mcMSmusmc*pow(Pi,2))/6. - (25*lnmusmkin*mcMSmusmc*pow(Pi,2))/36. - (106559417*pow(mcMSmusmc,25)*pow(Pi,2))/(6.374388636e11*pow(mkin,24)) + 
   (533*lnmcmsmkin*pow(mcMSmusmc,25)*pow(Pi,2))/(277725.*pow(mkin,24)) - (240817793781176357*pow(mcMSmusmc,24)*pow(Pi,2))/(6.28867544642696e20*pow(mkin,23)) + 
   (877591*ln2*pow(mcMSmusmc,24)*pow(Pi,2))/(6.79477248e9*pow(mkin,23)) + (1620816161*lnmcmsmkin*pow(mcMSmusmc,24)*pow(Pi,2))/(8.2216747008e11*pow(mkin,23)) - 
   (171475369*pow(mcMSmusmc,23)*pow(Pi,2))/(7.7816811996e11*pow(mkin,22)) + (445*lnmcmsmkin*pow(mcMSmusmc,23)*pow(Pi,2))/(192717.*pow(mkin,22)) - 
   (174200135864459*pow(mcMSmusmc,22)*pow(Pi,2))/(3.495434721951744e17*pow(mkin,21)) + (38675*ln2*pow(mcMSmusmc,22)*pow(Pi,2))/(2.33570304e8*pow(mkin,21)) + 
   (13926181*lnmcmsmkin*pow(mcMSmusmc,22)*pow(Pi,2))/(5.8392576e9*pow(mkin,21)) - (846435761*pow(mcMSmusmc,21)*pow(Pi,2))/(2.83231951884e12*pow(mkin,20)) + 
   (365*lnmcmsmkin*pow(mcMSmusmc,21)*pow(Pi,2))/(128877.*pow(mkin,20)) - (9007367733163*pow(mcMSmusmc,20)*pow(Pi,2))/(1.34810154565632e16*pow(mkin,19)) + 
   (143*ln2*pow(mcMSmusmc,20)*pow(Pi,2))/(655360.*pow(mkin,19)) + (156353*lnmcmsmkin*pow(mcMSmusmc,20)*pow(Pi,2))/(5.308416e7*pow(mkin,19)) - 
   (86836957*pow(mcMSmusmc,19)*pow(Pi,2))/(2.0687188752e11*pow(mkin,18)) + (293*lnmcmsmkin*pow(mcMSmusmc,19)*pow(Pi,2))/(82365.*pow(mkin,18)) - 
   (26463251891*pow(mcMSmusmc,18)*pow(Pi,2))/(2.845499424768e13*pow(mkin,17)) + (4147*ln2*pow(mcMSmusmc,18)*pow(Pi,2))/(1.3934592e7*pow(mkin,17)) + 
   (52013*lnmcmsmkin*pow(mcMSmusmc,18)*pow(Pi,2))/(1.3934592e7*pow(mkin,17)) - (22757641*pow(mcMSmusmc,17)*pow(Pi,2))/(3.6923796e10*pow(mkin,16)) + 
   (229*lnmcmsmkin*pow(mcMSmusmc,17)*pow(Pi,2))/(49725.*pow(mkin,16)) - (58806560951*pow(mcMSmusmc,16)*pow(Pi,2))/(4.326189170688e13*pow(mkin,15)) + 
   (1001*ln2*pow(mcMSmusmc,16)*pow(Pi,2))/(2.359296e6*pow(mkin,15)) + (565351*lnmcmsmkin*pow(mcMSmusmc,16)*pow(Pi,2))/(1.15605504e8*pow(mkin,15)) - 
   (17346493*pow(mcMSmusmc,15)*pow(Pi,2))/(1.808754948e10*pow(mkin,14)) + (173*lnmcmsmkin*pow(mcMSmusmc,15)*pow(Pi,2))/(27885.*pow(mkin,14)) - 
   (82285201*pow(mcMSmusmc,14)*pow(Pi,2))/(3.87459072e10*pow(mkin,13)) + (23*ln2*pow(mcMSmusmc,14)*pow(Pi,2))/(35840.*pow(mkin,13)) + 
   (2161*lnmcmsmkin*pow(mcMSmusmc,14)*pow(Pi,2))/(322560.*pow(mkin,13)) - (31786481*pow(mcMSmusmc,13)*pow(Pi,2))/(1.967766372e10*pow(mkin,12)) + 
   (125*lnmcmsmkin*pow(mcMSmusmc,13)*pow(Pi,2))/(14157.*pow(mkin,12)) - (84041429*pow(mcMSmusmc,12)*pow(Pi,2))/(2.29920768e10*pow(mkin,11)) + 
   (175*ln2*pow(mcMSmusmc,12)*pow(Pi,2))/(165888.*pow(mkin,11)) + (40553*lnmcmsmkin*pow(mcMSmusmc,12)*pow(Pi,2))/(4.1472e6*pow(mkin,11)) - 
   (1055689*pow(mcMSmusmc,11)*pow(Pi,2))/(3.4577928e8*pow(mkin,10)) + (85*lnmcmsmkin*pow(mcMSmusmc,11)*pow(Pi,2))/(6237.*pow(mkin,10)) - 
   (359801*pow(mcMSmusmc,10)*pow(Pi,2))/(4.89888e7*pow(mkin,9)) + (17*ln2*pow(mcMSmusmc,10)*pow(Pi,2))/(8640.*pow(mkin,9)) + (17*lnmcmsmkin*pow(mcMSmusmc,10)*pow(Pi,2))/(1080.*pow(mkin,9)) - 
   (416029*pow(mcMSmusmc,9)*pow(Pi,2))/(6.001128e7*pow(mkin,8)) + (53*lnmcmsmkin*pow(mcMSmusmc,9)*pow(Pi,2))/(2205.*pow(mkin,8)) - (1789*pow(mcMSmusmc,8)*pow(Pi,2))/(90720.*pow(mkin,7)) + 
   (7*ln2*pow(mcMSmusmc,8)*pow(Pi,2))/(1536.*pow(mkin,7)) + (139*lnmcmsmkin*pow(mcMSmusmc,8)*pow(Pi,2))/(4608.*pow(mkin,7)) - (262769*pow(mcMSmusmc,7)*pow(Pi,2))/(1.1907e7*pow(mkin,6)) + 
   (29*lnmcmsmkin*pow(mcMSmusmc,7)*pow(Pi,2))/(525.*pow(mkin,6)) - (311*pow(mcMSmusmc,6)*pow(Pi,2))/(2592.*pow(mkin,5)) + (11*ln2*pow(mcMSmusmc,6)*pow(Pi,2))/(648.*pow(mkin,5)) + 
   (113*lnmcmsmkin*pow(mcMSmusmc,6)*pow(Pi,2))/(1296.*pow(mkin,5)) - (149*pow(mcMSmusmc,5)*pow(Pi,2))/(1080.*pow(mkin,4)) + (13*lnmcmsmkin*pow(mcMSmusmc,5)*pow(Pi,2))/(45.*pow(mkin,4)) - 
   (pow(mcMSmusmc,4)*pow(Pi,2))/(12.*pow(mkin,3)) - (5*ln2*pow(mcMSmusmc,4)*pow(Pi,2))/(144.*pow(mkin,3)) - (2*pow(ln2,2)*pow(mcMSmusmc,4)*pow(Pi,2))/(81.*pow(mkin,3)) + 
   (67*lnmcmsmkin*pow(mcMSmusmc,4)*pow(Pi,2))/(144.*pow(mkin,3)) + (ln2*lnmcmsmkin*pow(mcMSmusmc,4)*pow(Pi,2))/(9.*pow(mkin,3)) - 
   (pow(lnmcmsmkin,2)*pow(mcMSmusmc,4)*pow(Pi,2))/(4.*pow(mkin,3)) + (2*lnmcmusmcms*pow(mcMSmusmc,4)*pow(Pi,2))/(9.*pow(mkin,3)) + (25*lnmusmkin*pow(mcMSmusmc,4)*pow(Pi,2))/(108.*pow(mkin,3)) - 
   (19507*pow(mcMSmusmc,3)*pow(Pi,2))/(1620.*pow(mkin,2)) + (1199*ln2*pow(mcMSmusmc,3)*pow(Pi,2))/(81.*pow(mkin,2)) + (271*lnmcmsmkin*pow(mcMSmusmc,3)*pow(Pi,2))/(162.*pow(mkin,2)) - 
   (lnmcmusmcms*pow(mcMSmusmc,3)*pow(Pi,2))/(2.*pow(mkin,2)) - (25*lnmusmkin*pow(mcMSmusmc,3)*pow(Pi,2))/(36.*pow(mkin,2)) + (13*pow(mcMSmusmc,2)*pow(Pi,2))/(12.*mkin) + 
   (3*lnmcmsmkin*pow(mcMSmusmc,2)*pow(Pi,2))/(2.*mkin) - (587741*mkin*pow(Pi,2))/38880. + (203*ln2*mkin*pow(Pi,2))/54. + (20*pow(ln2,2)*mkin*pow(Pi,2))/81. - (125*lnmusmkin*mkin*pow(Pi,2))/108. - 
   (25*ln2*lnmusmkin*mkin*pow(Pi,2))/54. - (3154*mufac*pow(Pi,2))/81. - (16*ln2*mufac*pow(Pi,2))/81. + (88*lnmufmus*mufac*pow(Pi,2))/9. + (8*lnmusmkin*mufac*pow(Pi,2))/9. - 
   (8*pow(mcMSmusmc,4)*mufac*pow(Pi,2))/(27.*pow(mkin,4)) + (16*pow(mcMSmusmc,3)*mufac*pow(Pi,2))/(27.*pow(mkin,3)) - (pow(mcMSmusmc,4)*pow(mufac,2)*pow(Pi,2))/(9.*pow(mkin,5)) + 
   (2*pow(mcMSmusmc,3)*pow(mufac,2)*pow(Pi,2))/(9.*pow(mkin,4)) - (1379*pow(mufac,2)*pow(Pi,2))/(108.*mkin) - (2*ln2*pow(mufac,2)*pow(Pi,2))/(27.*mkin) + 
   (11*lnmufmus*pow(mufac,2)*pow(Pi,2))/(3.*mkin) + (lnmusmkin*pow(mufac,2)*pow(Pi,2))/(3.*mkin) + (7*mcMSmusmc*NLMSOS*pow(Pi,2))/27. - (2*ln2*mcMSmusmc*NLMSOS*pow(Pi,2))/9. - 
   (lnmcmsmkin*mcMSmusmc*NLMSOS*pow(Pi,2))/9. + (lnmusmkin*mcMSmusmc*NLMSOS*pow(Pi,2))/18. + (23*pow(mcMSmusmc,25)*NLMSOS*pow(Pi,2))/(415800.*pow(mkin,24)) - 
   (11*pow(mcMSmusmc,24)*NLMSOS*pow(Pi,2))/(260820.*pow(mkin,23)) + (7*pow(mcMSmusmc,23)*NLMSOS*pow(Pi,2))/(96140.*pow(mkin,22)) - (20*pow(mcMSmusmc,22)*NLMSOS*pow(Pi,2))/(355509.*pow(mkin,21)) + 
   (19*pow(mcMSmusmc,21)*NLMSOS*pow(Pi,2))/(192780.*pow(mkin,20)) - (pow(mcMSmusmc,20)*NLMSOS*pow(Pi,2))/(12920.*pow(mkin,19)) + (17*pow(mcMSmusmc,19)*NLMSOS*pow(Pi,2))/(123120.*pow(mkin,18)) - 
   (16*pow(mcMSmusmc,18)*NLMSOS*pow(Pi,2))/(144585.*pow(mkin,17)) + (5*pow(mcMSmusmc,17)*NLMSOS*pow(Pi,2))/(24752.*pow(mkin,16)) - (7*pow(mcMSmusmc,16)*NLMSOS*pow(Pi,2))/(42120.*pow(mkin,15)) + 
   (13*pow(mcMSmusmc,15)*NLMSOS*pow(Pi,2))/(41580.*pow(mkin,14)) - (4*pow(mcMSmusmc,14)*NLMSOS*pow(Pi,2))/(15015.*pow(mkin,13)) + (11*pow(mcMSmusmc,13)*NLMSOS*pow(Pi,2))/(21060.*pow(mkin,12)) - 
   (5*pow(mcMSmusmc,12)*NLMSOS*pow(Pi,2))/(10692.*pow(mkin,11)) + (3*pow(mcMSmusmc,11)*NLMSOS*pow(Pi,2))/(3080.*pow(mkin,10)) - (8*pow(mcMSmusmc,10)*NLMSOS*pow(Pi,2))/(8505.*pow(mkin,9)) + 
   (7*pow(mcMSmusmc,9)*NLMSOS*pow(Pi,2))/(3240.*pow(mkin,8)) - (pow(mcMSmusmc,8)*NLMSOS*pow(Pi,2))/(420.*pow(mkin,7)) + (5*pow(mcMSmusmc,7)*NLMSOS*pow(Pi,2))/(756.*pow(mkin,6)) - 
   (4*pow(mcMSmusmc,6)*NLMSOS*pow(Pi,2))/(405.*pow(mkin,5)) + (pow(mcMSmusmc,5)*NLMSOS*pow(Pi,2))/(20.*pow(mkin,4)) + (13*pow(mcMSmusmc,4)*NLMSOS*pow(Pi,2))/(324.*pow(mkin,3)) - 
   (lnmcmsmkin*pow(mcMSmusmc,4)*NLMSOS*pow(Pi,2))/(27.*pow(mkin,3)) - (lnmusmkin*pow(mcMSmusmc,4)*NLMSOS*pow(Pi,2))/(54.*pow(mkin,3)) + (7*pow(mcMSmusmc,3)*NLMSOS*pow(Pi,2))/(54.*pow(mkin,2)) - 
   (2*ln2*pow(mcMSmusmc,3)*NLMSOS*pow(Pi,2))/(9.*pow(mkin,2)) - (lnmcmsmkin*pow(mcMSmusmc,3)*NLMSOS*pow(Pi,2))/(9.*pow(mkin,2)) + (lnmusmkin*pow(mcMSmusmc,3)*NLMSOS*pow(Pi,2))/(18.*pow(mkin,2)) + 
   (305*mkin*NLMSOS*pow(Pi,2))/216. + (11*ln2*mkin*NLMSOS*pow(Pi,2))/81. - (2*pow(ln2,2)*mkin*NLMSOS*pow(Pi,2))/81. + (35*lnmusmkin*mkin*NLMSOS*pow(Pi,2))/108. + 
   (ln2*lnmusmkin*mkin*NLMSOS*pow(Pi,2))/27. + (8*mufac*NLMSOS*pow(Pi,2))/81. + (pow(mufac,2)*NLMSOS*pow(Pi,2))/(27.*mkin) - (13*mkin*pow(NLMSOS,2)*pow(Pi,2))/324. - 
   (lnmusmkin*mkin*pow(NLMSOS,2)*pow(Pi,2))/54. + (208*mufac*NLOSKIN*pow(Pi,2))/81. - (16*lnmufmus*mufac*NLOSKIN*pow(Pi,2))/27. + (23*pow(mufac,2)*NLOSKIN*pow(Pi,2))/(27.*mkin) - 
   (2*lnmufmus*pow(mufac,2)*NLOSKIN*pow(Pi,2))/(9.*mkin) - (8*mufac*pow(NLOSKIN,2)*pow(Pi,2))/243. - (pow(mufac,2)*pow(NLOSKIN,2)*pow(Pi,2))/(81.*mkin) - (13*mcMSmusmc*pow(Pi,3))/162. + 
   (5698043*pow(mcMSmusmc,25)*pow(Pi,3))/(1.189085184e11*pow(mkin,24)) + (81991*pow(mcMSmusmc,23)*pow(Pi,3))/(1.374683136e9*pow(mkin,22)) + 
   (127699*pow(mcMSmusmc,21)*pow(Pi,3))/(1.684537344e9*pow(mkin,20)) + (99671*pow(mcMSmusmc,19)*pow(Pi,3))/(1.00859904e9*pow(mkin,18)) + 
   (1925*pow(mcMSmusmc,17)*pow(Pi,3))/(1.4483456e7*pow(mkin,16)) + (377*pow(mcMSmusmc,15)*pow(Pi,3))/(2.02752e6*pow(mkin,14)) + (1771*pow(mcMSmusmc,13)*pow(Pi,3))/(6.469632e6*pow(mkin,12)) + 
   (17*pow(mcMSmusmc,11)*pow(Pi,3))/(39424.*pow(mkin,10)) + (77*pow(mcMSmusmc,9)*pow(Pi,3))/(103680.*pow(mkin,8)) + (25*pow(mcMSmusmc,7)*pow(Pi,3))/(18144.*pow(mkin,6)) - 
   (pow(mcMSmusmc,5)*pow(Pi,3))/(240.*pow(mkin,4)) - (7*pow(mcMSmusmc,3)*pow(Pi,3))/(108.*pow(mkin,2)) - (271*pow(mcMSmusmc,4)*pow(Pi,4))/(19440.*pow(mkin,3)) - 
   (lnmcmsmkin*pow(mcMSmusmc,2)*pow(Pi,4))/(12.*mkin) + (451*mkin*pow(Pi,4))/7776. + (2*mufac*pow(Pi,4))/3. + (pow(mufac,2)*pow(Pi,4))/(4.*mkin) - (61*mkin*NLMSOS*pow(Pi,4))/1944. + 
   (2710689767*pow(mcMSmusmc,24)*Zeta3)/(1.64433494016e12*pow(mkin,23)) + (23017987*pow(mcMSmusmc,22)*Zeta3)/(1.16785152e10*pow(mkin,21)) + (254791*pow(mcMSmusmc,20)*Zeta3)/(1.0616832e8*pow(mkin,19)) + 
   (83291*pow(mcMSmusmc,18)*Zeta3)/(2.7869184e7*pow(mkin,17)) + (885457*pow(mcMSmusmc,16)*Zeta3)/(2.31211008e8*pow(mkin,15)) + (3287*pow(mcMSmusmc,14)*Zeta3)/(645120.*pow(mkin,13)) + 
   (59231*pow(mcMSmusmc,12)*Zeta3)/(8.2944e6*pow(mkin,11)) + (187*pow(mcMSmusmc,10)*Zeta3)/(17280.*pow(mkin,9)) + (173*pow(mcMSmusmc,8)*Zeta3)/(9216.*pow(mkin,7)) + 
   (29*pow(mcMSmusmc,6)*Zeta3)/(648.*pow(mkin,5)) + (2309*pow(mcMSmusmc,4)*Zeta3)/(864.*pow(mkin,3)) - (14*lnmcmsmkin*pow(mcMSmusmc,4)*Zeta3)/(9.*pow(mkin,3)) + (11*pow(mcMSmusmc,2)*Zeta3)/(2.*mkin) + 
   (23*mkin*Zeta3)/24. + (85*lnmusmkin*mkin*Zeta3)/36. - (3070*mufac*Zeta3)/27. - (1535*pow(mufac,2)*Zeta3)/(36.*mkin) - (pow(mcMSmusmc,4)*NLMSOS*Zeta3)/(3.*pow(mkin,3)) + (667*mkin*NLMSOS*Zeta3)/216. + 
   (7*lnmusmkin*mkin*NLMSOS*Zeta3)/9. - (7*mkin*pow(NLMSOS,2)*Zeta3)/54. + (140*mufac*NLOSKIN*Zeta3)/27. + (35*pow(mufac,2)*NLOSKIN*Zeta3)/(18.*mkin) - (3*pow(mcMSmusmc,2)*pow(Pi,2)*Zeta3)/(4.*mkin) + 
   (1439*mkin*pow(Pi,2)*Zeta3)/432. - (5*pow(mcMSmusmc,2)*Zeta5)/(2.*mkin) - (1975*mkin*Zeta5)/216.);
    
    
    double mMS = 0.0;
    for(int i = 0; i <= nloops; i++) {
       mMS += ret[i];
    }    
    return mMS;
}

double CRunDec::mkin2mMSD(double mkin, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops) {
    if(nloops<0||nloops>3){
      cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops<<" LOOPS"<<endl;
      RETURN
    }
    double ret[5];
    double lnmusmkin = log(mus*mus/(mkin*mkin));
    double lnmufmus = log(2*mufac/mus);
    double lnmusmcms = 0.0;
    double lnmcmsmkin = 0.0;
    double lnmcmusmcms = 0.0;  
    double mcMSmusmc = 0.0;
    double musmc = 0.0;
    if(mcMSmusmcin > 0.0) {
       lnmusmcms = log(mus*mus/(mcMSmusmcin*mcMSmusmcin));
       lnmcmsmkin = log(mcMSmusmcin/mkin);
       lnmcmusmcms = log(musmcin*musmcin/(mcMSmusmcin*mcMSmusmcin));
       mcMSmusmc = mcMSmusmcin;
       musmc = musmcin;
    }
    
    ret[0] = mkin;
    ret[1] = apinlmus*((-4*mkin)/3. - lnmusmkin*mkin + (16*mufac)/9. + (2*pow(mufac,2))/(3.*mkin));
    ret[2] = apinlmus*apinlmus*((-3019*mkin)/288. - (461*lnmusmkin*mkin)/72. - (23*pow(lnmusmkin,2)*mkin)/24. + (892*mufac)/27. - (88*lnmufmus*mufac)/9. - (16*lnmusmkin*mufac)/9. + (95*pow(mufac,2))/(9.*mkin) - 
   (11*lnmufmus*pow(mufac,2))/(3.*mkin) - (2*lnmusmkin*pow(mufac,2))/(3.*mkin) + (71*mkin*NLMSOS)/144. + (13*lnmusmkin*mkin*NLMSOS)/36. + (pow(lnmusmkin,2)*mkin*NLMSOS)/12. - (128*mufac*NLOSKIN)/81. + 
   (16*lnmufmus*mufac*NLOSKIN)/27. - (13*pow(mufac,2)*NLOSKIN)/(27.*mkin) + (2*lnmufmus*pow(mufac,2)*NLOSKIN)/(9.*mkin) - (mkin*pow(Pi,2))/3. - (ln2*mkin*pow(Pi,2))/9. - (8*mufac*pow(Pi,2))/9. - 
   (pow(mufac,2)*pow(Pi,2))/(3.*mkin) + (mkin*NLMSOS*pow(Pi,2))/18. + (mkin*Zeta3)/6.);
    ret[3] = apinlmus*apinlmus*apinlmus*((-9514621*mkin)/93312. + (220*A4*mkin)/27. + (55*pow(ln2,4)*mkin)/162. - (22273*lnmusmkin*mkin)/324. - (14275*pow(lnmusmkin,2)*mkin)/864. - (601*pow(lnmusmkin,3)*mkin)/432. + (22496*mufac)/27. - 
   (3416*lnmufmus*mufac)/9. + (484*pow(lnmufmus,2)*mufac)/9. - (2950*lnmusmkin*mufac)/81. + (88*lnmufmus*lnmusmkin*mufac)/9. - (46*pow(lnmusmkin,2)*mufac)/27. + (75931*pow(mufac,2))/(324.*mkin) - 
   (733*lnmufmus*pow(mufac,2))/(6.*mkin) + (121*pow(lnmufmus,2)*pow(mufac,2))/(6.*mkin) - (1277*lnmusmkin*pow(mufac,2))/(108.*mkin) + (11*lnmufmus*lnmusmkin*pow(mufac,2))/(3.*mkin) - 
   (23*pow(lnmusmkin,2)*pow(mufac,2))/(36.*mkin) + (246643*mkin*NLMSOS)/23328. - (8*A4*mkin*NLMSOS)/27. - (pow(ln2,4)*mkin*NLMSOS)/81. + (1283*lnmusmkin*mkin*NLMSOS)/162. + 
   (107*pow(lnmusmkin,2)*mkin*NLMSOS)/48. + (25*pow(lnmusmkin,3)*mkin*NLMSOS)/108. - (11*mufac*NLMSOS)/27. + (4*lnmusmkin*mufac*NLMSOS)/81. + (4*pow(lnmusmkin,2)*mufac*NLMSOS)/27. - 
   (11*pow(mufac,2)*NLMSOS)/(72.*mkin) + (lnmusmkin*pow(mufac,2)*NLMSOS)/(54.*mkin) + (pow(lnmusmkin,2)*pow(mufac,2)*NLMSOS)/(18.*mkin) - (2353*mkin*pow(NLMSOS,2))/23328. - 
   (89*lnmusmkin*mkin*pow(NLMSOS,2))/648. - (13*pow(lnmusmkin,2)*mkin*pow(NLMSOS,2))/216. - (pow(lnmusmkin,3)*mkin*pow(NLMSOS,2))/108. - (20303*mufac*NLOSKIN)/243. + (3388*lnmufmus*mufac*NLOSKIN)/81. - 
   (176*pow(lnmufmus,2)*mufac*NLOSKIN)/27. + (128*lnmusmkin*mufac*NLOSKIN)/81. - (16*lnmufmus*lnmusmkin*mufac*NLOSKIN)/27. - (14429*pow(mufac,2)*NLOSKIN)/(648.*mkin) + 
   (715*lnmufmus*pow(mufac,2)*NLOSKIN)/(54.*mkin) - (22*pow(lnmufmus,2)*pow(mufac,2)*NLOSKIN)/(9.*mkin) + (13*lnmusmkin*pow(mufac,2)*NLOSKIN)/(27.*mkin) - 
   (2*lnmufmus*lnmusmkin*pow(mufac,2)*NLOSKIN)/(9.*mkin) + (1292*mufac*pow(NLOSKIN,2))/729. - (256*lnmufmus*mufac*pow(NLOSKIN,2))/243. + (16*pow(lnmufmus,2)*mufac*pow(NLOSKIN,2))/81. + 
   (209*pow(mufac,2)*pow(NLOSKIN,2))/(486.*mkin) - (26*lnmufmus*pow(mufac,2)*pow(NLOSKIN,2))/(81.*mkin) + (2*pow(lnmufmus,2)*pow(mufac,2)*pow(NLOSKIN,2))/(27.*mkin) - 
   (644201*mkin*pow(Pi,2))/38880. + (587*ln2*mkin*pow(Pi,2))/162. + (22*pow(ln2,2)*mkin*pow(Pi,2))/81. - (3*lnmusmkin*mkin*pow(Pi,2))/2. - (ln2*lnmusmkin*mkin*pow(Pi,2))/2. - 
   (1054*mufac*pow(Pi,2))/27. - (16*ln2*mufac*pow(Pi,2))/81. + (88*lnmufmus*mufac*pow(Pi,2))/9. + (8*lnmusmkin*mufac*pow(Pi,2))/9. - (461*pow(mufac,2)*pow(Pi,2))/(36.*mkin) - 
   (2*ln2*pow(mufac,2)*pow(Pi,2))/(27.*mkin) + (11*lnmufmus*pow(mufac,2)*pow(Pi,2))/(3.*mkin) + (lnmusmkin*pow(mufac,2)*pow(Pi,2))/(3.*mkin) + (967*mkin*NLMSOS*pow(Pi,2))/648. + 
   (11*ln2*mkin*NLMSOS*pow(Pi,2))/81. - (2*pow(ln2,2)*mkin*NLMSOS*pow(Pi,2))/81. + (13*lnmusmkin*mkin*NLMSOS*pow(Pi,2))/36. + (ln2*lnmusmkin*mkin*NLMSOS*pow(Pi,2))/27. + (8*mufac*NLMSOS*pow(Pi,2))/81. + 
   (pow(mufac,2)*NLMSOS*pow(Pi,2))/(27.*mkin) - (13*mkin*pow(NLMSOS,2)*pow(Pi,2))/324. - (lnmusmkin*mkin*pow(NLMSOS,2)*pow(Pi,2))/54. + (208*mufac*NLOSKIN*pow(Pi,2))/81. - 
   (16*lnmufmus*mufac*NLOSKIN*pow(Pi,2))/27. + (23*pow(mufac,2)*NLOSKIN*pow(Pi,2))/(27.*mkin) - (2*lnmufmus*pow(mufac,2)*NLOSKIN*pow(Pi,2))/(9.*mkin) - (8*mufac*pow(NLOSKIN,2)*pow(Pi,2))/243. - 
   (pow(mufac,2)*pow(NLOSKIN,2)*pow(Pi,2))/(81.*mkin) + (695*mkin*pow(Pi,4))/7776. + (2*mufac*pow(Pi,4))/3. + (pow(mufac,2)*pow(Pi,4))/(4.*mkin) - (61*mkin*NLMSOS*pow(Pi,4))/1944. - 
   (61*mkin*Zeta3)/27. + (19*lnmusmkin*mkin*Zeta3)/12. - (3070*mufac*Zeta3)/27. - (1535*pow(mufac,2)*Zeta3)/(36.*mkin) + (241*mkin*NLMSOS*Zeta3)/72. + (7*lnmusmkin*mkin*NLMSOS*Zeta3)/9. - 
   (7*mkin*pow(NLMSOS,2)*Zeta3)/54. + (140*mufac*NLOSKIN*Zeta3)/27. + (35*pow(mufac,2)*NLOSKIN*Zeta3)/(18.*mkin) + (1439*mkin*pow(Pi,2)*Zeta3)/432. - (1975*mkin*Zeta5)/216.);
    
    
    double mMS = 0.0;
    for(int i = 0; i <= nloops; i++) {
       mMS += ret[i];
    }    
    return mMS;
}

// Function: double CRunDec::mMS2mKIN(double mMS, std::pair<double,double>* mq,
//                                    double asmus, double mus, double muf, int nl, int nloops, std::string deccase)
double CRunDec::mMS2mKIN(double mMS, std::pair<double,double>* mq,
                         double asmus, double mus, double muf, int nlmsos, int nloskin, int nloops, std::string deccase) {
     if(nloops<0||nloops>3){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops <<" LOOPS"<<endl;
       RETURN  
     }
     double mkin;
     int decc = 0;
     if(deccase == "A") decc = 0;
     else if(deccase == "B") decc = 1;
     else if(deccase == "C") decc = 2;
     else if(deccase == "D") decc = 3;
     else if(deccase == "") decc = 4;
     else {
       cout<<"DECCASE "<< deccase <<" NOT IMPLEMENTED"<<endl;
       RETURN       
     }
     
     switch(decc) {
       case 0:
           mkin = mMS2mkinA(mMS, asmus/Pi, mus, muf, nlmsos, nloskin, mq[0].first, mq[0].second, nloops);
           break;
       case 1:
           mkin = mMS2mkinB(mMS, asmus/Pi, mus, muf, nlmsos, nloskin, mq[0].first, mq[0].second, nloops);
           break;
       case 2:
           mkin = mMS2mkinC(mMS, asmus/Pi, mus, muf, nlmsos, nloskin, mq[0].first, mq[0].second, nloops);
           break;
       case 3:
           mkin = mMS2mkinD(mMS, asmus/Pi, mus, muf, nlmsos, nloskin, mq[0].first, mq[0].second, nloops);
           break;
       case 4:
           mkin = mMS2mkinD(mMS, asmus/Pi, mus, muf, nlmsos, nloskin, 0.0, 0.0, nloops);
           break;
     }
     return mkin;
}

double CRunDec::mMS2mKIN(double mMS, std::pair<double,double>* mq, double asmus, double mus, double muf, int nloops, std::string deccase) {
     if(deccase == "A") return mMS2mKIN(mMS, mq, asmus, mus, muf, 3, 3, nloops, deccase);
     else if(deccase == "B") return mMS2mKIN(mMS, mq, asmus, mus, muf, 3, 3, nloops, deccase);
     else if(deccase == "C") return mMS2mKIN(mMS, mq, asmus, mus, muf, 3, 4, nloops, deccase);
     else if(deccase == "D") return mMS2mKIN(mMS, mq, asmus, mus, muf, 3, 3, nloops, deccase);
     else {
       cout<<"DECCASE "<< deccase <<" NOT IMPLEMENTED"<<endl;
       RETURN       
     }
}

// Function: double CRunDec::mKIN2mMS(double mKIN, std::pair<double,double>* mq,
//                                    double asmus, double mus, double muf, int nlmsos, int nloskin, int nloops, std::string deccase)
double CRunDec::mKIN2mMS(double mKIN, std::pair<double,double>* mq,
                         double asmus, double mus, double muf, int nlmsos, int nloskin, int nloops, std::string deccase) {
     if(nloops<0||nloops>3){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nloops <<" LOOPS"<<endl;
       RETURN  
     }
     double mMS;
     int decc = 0;
     if(deccase == "A") decc = 0;
     else if(deccase == "B") decc = 1;
     else if(deccase == "C") decc = 2;
     else if(deccase == "D") decc = 3;
     else if(deccase == "") decc = 4;
     else {
       cout<<"DECCASE "<< deccase <<" NOT IMPLEMENTED"<<endl;
       RETURN       
     }
     
     switch(decc) {
       case 0:
           mMS = mkin2mMSA(mKIN, asmus/Pi, mus, muf, nlmsos, nloskin, mq[0].first, mq[0].second, nloops);
           break;
       case 1:
           mMS = mkin2mMSB(mKIN, asmus/Pi, mus, muf, nlmsos, nloskin, mq[0].first, mq[0].second, nloops);
           break;
       case 2:
           mMS = mkin2mMSC(mKIN, asmus/Pi, mus, muf, nlmsos, nloskin, mq[0].first, mq[0].second, nloops);
           break;
       case 3:
           mMS = mkin2mMSD(mKIN, asmus/Pi, mus, muf, nlmsos, nloskin, mq[0].first, mq[0].second, nloops);
           break;
       case 4:
           mMS = mkin2mMSD(mKIN, asmus/Pi, mus, muf, nlmsos, nloskin, 0.0, 0.0, nloops);
           break;
     }
     return mMS;
}

double CRunDec::mKIN2mMS(double mKIN, std::pair<double,double>* mq, double asmus, double mus, double muf, int nloops, std::string deccase) {
     if(deccase == "A") return mKIN2mMS(mKIN, mq, asmus, mus, muf, 3, 3, nloops, deccase);
     else if(deccase == "B") return mKIN2mMS(mKIN, mq, asmus, mus, muf, 3, 3, nloops, deccase);
     else if(deccase == "C") return mKIN2mMS(mKIN, mq, asmus, mus, muf, 3, 4, nloops, deccase);
     else if(deccase == "D") return mKIN2mMS(mKIN, mq, asmus, mus, muf, 3, 3, nloops, deccase);
     else {
       cout<<"DECCASE "<< deccase <<" NOT IMPLEMENTED"<<endl;
       RETURN       
     }
}


// Coefficients of eq.(22) of [RunDec]
double CRunDec::fas5to6os(double A, double mass, double mu, double nlq, double nl){
     double lmM  = log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=A*(-lmM)/6.;
     sum[2]=A*A*(-7./24. - (19.*lmM)/24. + (lmM*lmM)/36.);
     sum[3]=A*A*A*(-58933./124416. - (8521.*lmM)/1728. - (131.*lmM*lmM)/576. -
                     (lmM*lmM*lmM)/216. +
       nlq*(2479./31104. + (409.*lmM)/1728. + Zeta2/9.) - (2.*Zeta2)/3. - 
       (2.*Zeta2*log(2))/9. -
       (80507.*Zeta3)/27648.);
     sum[4]=
       A*A*A*A*(-141841753./24494400. - 
       (7693.*lmM*lmM)/1152. - (8371.*lmM*lmM*lmM)/10368. + lmM*lmM*lmM*lmM/1296. - 
       (644201.*Pi*Pi)/116640. - (71102219.*Pi*Pi*Pi*Pi)/195955200. - (49.*Zeta2)/18. - 
       (49*log(2)*Zeta2)/54 + (49*Zeta3)/216 + (587*Pi*Pi*log(2))/486 + 
       (9318467*Pi*Pi*Pi*Pi*log(2))/32659200 - (2913037*Pi*Pi*log(2)*log(2))/1306368 + 
       (340853*Pi*Pi*log(2)*log(2)*log(2))/816480 + (3179149*log(2)*log(2)*log(2)*log(2))/1306368 - 
       (340853*log(2)*log(2)*log(2)*log(2)*log(2))/1360800 + (3179149*A4)/54432 + 
       (340853*A5)/11340 + lmM*(-19696909./746496. - (29*Zeta2)/9 - 
         (29*log(2)*Zeta2)/27 + (59*Zeta3)/108 - (2529743*Zeta3)/165888) + 
       nlq*nlq*(-140825./1492992. - (493*lmM*lmM)/20736 - (13*Pi*Pi)/972 + 
         lmM*(-1679./186624. - Zeta2/27) - (19*Zeta3)/1728) - 
       (2428169183*Zeta3)/87091200 + (1439*Pi*Pi*Zeta3)/1296 + 
       nlq*(1773073./746496. + (6661*lmM*lmM)/10368 + (107*lmM*lmM*lmM)/1728 + 
         (967*Pi*Pi)/1944 - (697709*Pi*Pi*Pi*Pi)/14929920 + (49*Zeta2)/108 + 
         (11*Pi*Pi*log(2))/243 - (1709*Pi*Pi*log(2)*log(2))/124416 + 
         (173*log(2)*log(2)*log(2)*log(2))/124416 + (173*A4)/5184 + 
         (4756441*Zeta3)/995328 + lmM*(1110443./373248. + (41*Zeta2)/54 + 
           (2*log(2)*Zeta2)/27 + (7*Zeta3)/27 + (110779*Zeta3)/82944) + 
         (115*Zeta5)/576) - (40596749*Zeta5)/1451520);

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsDownOS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsDownOS(double als, double massth, double muth, int nl){
       if(nl<1||nl>5){
         cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
         RETURN 
       }
       double erg=(this->fas5to6os(als/Pi, massth, muth, Nf,nl));
       return als*erg;
}

// Coefficients of eq.(25) of [RunDec]
double CRunDec::fas6to5os(double A, double mass, double mu, double nlq, double nl){
     double lmM=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=A*(lmM)/6.;
     sum[2]=A*A*(7./24. + (19.*lmM)/24. + (lmM*lmM)/36.);
     sum[3]=A*A*A*(58933./124416. + (2.*Zeta2)/3. + (2.*Zeta2*log(2))/9. + 
                    (80507.*Zeta3)/27648. + (8941.*lmM)/1728. +
                    (511.*lmM*lmM)/576. + (lmM*lmM*lmM)/216. +
                    nlq*(-2479./31104. - Zeta2/9. - (409.*lmM)/1728));
     sum[4]=A*A*A*A*
       (2313952./382725. + (47039.*lmM*lmM)/3456. + (14149.*lmM*lmM*lmM)/10368. + 
       lmM*lmM*lmM*lmM/1296. + (644201.*Pi*Pi)/116640. + (71102219.*Pi*Pi*Pi*Pi)/195955200. + 
       (49.*Zeta2)/18. + (49.*log(2)*Zeta2)/54. - (49.*Zeta3)/216. - (587.*Pi*Pi*log(2))/486. - 
       (9318467.*Pi*Pi*Pi*Pi*log(2))/32659200. + (2913037.*Pi*Pi*log(2)*log(2))/1306368. - 
       (340853.*Pi*Pi*log(2)*log(2)*log(2))/816480. - (3179149.*log(2)*log(2)*log(2)*log(2))/1306368. + 
       (340853.*log(2)*log(2)*log(2)*log(2)*log(2))/1360800. - (3179149.*A4)/54432. - 
       (340853.*A5)/11340. + nlq*nlq*(140825./1492992. + 
         (493.*lmM*lmM)/20736. + (13.*Pi*Pi)/972. + lmM*(1679./186624. + Zeta2/27.) + 
         (19.*Zeta3)/1728.) + (2428169183.*Zeta3)/87091200. - 
       (1439.*Pi*Pi*Zeta3)/1296. + lmM*(21084715./746496. + (35.*Zeta2)/9. + 
         (35.*log(2)*Zeta2)/27. - (65.*Zeta3)/108. + (3022001.*Zeta3)/165888.) + 
       nlq*(-1773073./746496. - (9115.*lmM*lmM)/10368. - (107.*lmM*lmM*lmM)/1728. - 
         (967.*Pi*Pi)/1944. + (697709.*Pi*Pi*Pi*Pi)/14929920. - (49.*Zeta2)/108. - 
         (11.*Pi*Pi*log(2))/243. + (1709.*Pi*Pi*log(2)*log(2))/124416. - 
         (173.*log(2)*log(2)*log(2)*log(2))/124416. - (173.*A4)/5184. + 
         lmM*(-1140191./373248. - (47.*Zeta2)/54. - (2.*log(2)*Zeta2)/27. - (7.*Zeta3)/27. - 
           (110779.*Zeta3)/82944.) - (4756441.*Zeta3)/995328. - 
         (115.*Zeta5)/576.) + (40596749.*Zeta5)/1451520.);

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsUpOS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsUpOS(double als, double massth, double muth, int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fas6to5os(als/Pi, massth, muth, Nf,nl));
     return als*erg;
}


double CRunDec::fas6to5ms(double A, double mass, double mu, double nlq, double nl){
     double lmm=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=  1.;
     sum[1]=A*(lmm)/6.;
     sum[2]=A*A*(-11./72. + (11.*lmm)/24. + lmm*lmm/36.);
     sum[3]=A*A*A*(-564731./124416. + (2645.*lmm)/1728. + (167.*lmm*lmm)/576. + lmm*lmm*lmm/216. + 
   					(2633./31104. - (67.*lmm)/576. + lmm*lmm/36.)*nlq + (82043.*Zeta3)/27648.);
     sum[4]=A*A*A*A*(-1165152397./24494400. - (3031309*A4)/54432 - (340853*A5)/11340
                     +(1837*lmm*lmm)/1152 + (2909*lmm*lmm*lmm)/10368 + lmm*lmm*lmm*lmm/1296
                     -(3031309*log(2)*log(2)*log(2)*log(2))/1306368 + (340853*log(2)*log(2)*log(2)*log(2)*log(2))/1360800
                     +(3031309*log(2)*log(2)*Pi*Pi)/1306368 - (340853*log(2)*log(2)*log(2)*Pi*Pi)/816480
                     +(76940219*Pi*Pi*Pi*Pi)/195955200 - (9318467*log(2)*Pi*Pi*Pi*Pi)/32659200 + nlq*nlq*(271883./4478976.
                     -(6865*lmm)/186624 + (77*lmm*lmm)/20736 - lmm*lmm*lmm/324 - (167*Zeta3)/5184)
                     +(2362581983*Zeta3)/87091200 + lmm*(-11093717./746496. + (3022001*Zeta3)/165888) + nlq*(4770941./2239488.
                     -(685*A4)/5184 + (277*lmm*lmm)/10368 + (271*lmm*lmm*lmm)/5184 - (685*log(2)*log(2)*log(2)*log(2))/124416
                     +(685*log(2)*log(2)*Pi*Pi)/124416 + (541549*Pi*Pi*Pi*Pi)/14929920 + lmm*(141937./373248.
                     -(110779*Zeta3)/82944) - (3645913*Zeta3)/995328 - (115*Zeta5)/576) + (12057583*Zeta5)/483840);
     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsUpMS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsUpMS(double als, double massth, double muth, int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fas6to5ms(als/Pi, massth, muth, Nf,nl));
     return als*erg;
}

double CRunDec::fas5to6ms(double A, double mass, double mu, double nlq, double nl){
     double lmm=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=  1.;
     sum[1]=A*(-lmm)/6.;
     sum[2]=A*A*(11./72. - (11*lmm)/24 + lmm*lmm/36.);
     sum[3]=A*A*A*(564731./124416. - (955.*lmm)/576. + (53.*lmm*lmm)/576. - lmm*lmm*lmm/216.
                   +(-2633./31104. + (67.*lmm)/576. - lmm*lmm/36.)*nlq - (82043.*Zeta3)/27648.);
     sum[4]=A*A*A*A*(291716893./6123600. + (3031309*A4)/54432 + (340853*A5)/11340
                     +(2177*lmm*lmm)/3456 - (1883*lmm*lmm*lmm)/10368 + lmm*lmm*lmm*lmm/1296
                     +(3031309*log(2)*log(2)*log(2)*log(2))/1306368
                     -(340853*log(2)*log(2)*log(2)*log(2)*log(2))/1360800
                     -(3031309*log(2)*log(2)*Pi*Pi)/1306368 + (340853*log(2)*log(2)*log(2)*Pi*Pi)/816480
                     -(76940219*Pi*Pi*Pi*Pi)/195955200 + (9318467*log(2)*Pi*Pi*Pi*Pi)/32659200
                     +lmm*(7391699./746496. - (2529743*Zeta3)/165888) + nlq*nlq*(-271883./4478976.
                     +(6865*lmm)/186624 - (77*lmm*lmm)/20736 + lmm*lmm*lmm/324 + (167*Zeta3)/5184)
                     -(2362581983*Zeta3)/87091200 + nlq*(-4770941./2239488. + (685*A4)/5184
                     -(1483*lmm*lmm)/10368 - (127*lmm*lmm*lmm)/5184 + (685*log(2)*log(2)*log(2)*log(2))/124416
                     -(685*log(2)*log(2)*Pi*Pi)/124416 - (541549*Pi*Pi*Pi*Pi)/14929920 + (3645913*Zeta3)/995328
                     +lmm*(-110341./373248. + (110779*Zeta3)/82944) + (115*Zeta5)/576) - (12057583*Zeta5)/483840);
     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsDownMS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsDownMS(double als, double massth, double muth, int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fas5to6ms(als/Pi, massth, muth, Nf,nl));
     return als*erg;
}


double CRunDec::fas6to5si(double A, double mass, double mu, double nlq, double nl){
     double lmmu=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0] = 1.0;
     sum[1] = (A*lmmu)/6.;
     sum[2] = A*A*(-11./72. + (19.*lmmu)/24. + lmmu*lmmu/36.);
     sum[3] = A*A*A*(-564731./124416. + (2191.*lmmu)/576. + (511.*lmmu*lmmu)/576. + lmmu*lmmu*lmmu/216. + 
              (2633./31104. - (281.*lmmu)/1728.)*nlq + (82043.*Zeta3)/27648.);
     sum[4] = A*A*A*A*(-1165152397./24494400. - (3031309*A4)/54432 - (340853*A5)/11340 + (33887*lmmu*lmmu)/3456
                       +(14149*lmmu*lmmu*lmmu)/10368 + lmmu*lmmu*lmmu*lmmu/1296
                       -(3031309*log(2)*log(2)*log(2)*log(2))/1306368 + (340853*log(2)*log(2)*log(2)*log(2)*log(2))/1360800
                       +(3031309*log(2)*log(2)*Pi*Pi)/1306368 - (340853*log(2)*log(2)*log(2)*Pi*Pi)/816480
                       +(76940219*Pi*Pi*Pi*Pi)/195955200 - (9318467*log(2)*Pi*Pi*Pi*Pi)/32659200 + nlq*nlq*(271883./4478976.
                       -(8545*lmmu)/186624 + (79*lmmu*lmmu)/6912 - (167*Zeta3)/5184) + (2362581983*Zeta3)/87091200
                       +lmmu*(-1531493./746496. + (2975921*Zeta3)/165888) + nlq*(4770941./2239488. - (685*A4)/5184
                       -(515*lmmu*lmmu)/1152 - (107*lmmu*lmmu*lmmu)/1728 - (685*log(2)*log(2)*log(2)*log(2))/124416
                       +(685*log(2)*log(2)*Pi*Pi)/124416 + (541549*Pi*Pi*Pi*Pi)/14929920 + lmmu*(-158687./373248.
                       -(133819*Zeta3)/82944) - (3645913*Zeta3)/995328 - (115*Zeta5)/576) + (12057583*Zeta5)/483840);
     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsUpSI(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsUpSI(double als, double massth, double muth, int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fas6to5si(als/Pi, massth, muth, Nf,nl));
     return als*erg;
}

double CRunDec::fas5to6si(double A, double mass, double mu, double nlq, double nl){
     double lmmu=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0] = 1;
     sum[1] = -(A*lmmu)/6.;
     sum[2] = A*A*(11./72. - (19.*lmmu)/24. + lmmu*lmmu/36.); 
     sum[3] = A*A*A*(564731./124416. - (6793.*lmmu)/1728. - (131.*lmmu*lmmu)/576. - lmmu*lmmu*lmmu/216. + 
              (-2633./31104. + (281.*lmmu)/1728.)*nlq - (82043.*Zeta3)/27648.);
     sum[4] = A*A*A*A*(291716893./6123600. + (3031309*A4)/54432 + (340853*A5)/11340 - (14023*lmmu*lmmu)/3456
                       -(8371*lmmu*lmmu*lmmu)/10368 + lmmu*lmmu*lmmu*lmmu/1296
                       +(3031309*log(2)*log(2)*log(2)*log(2))/1306368 - (340853*log(2)*log(2)*log(2)*log(2)*log(2))/1360800
                       -(3031309*log(2)*log(2)*Pi*Pi)/1306368 + (340853*log(2)*log(2)*log(2)*Pi*Pi)/816480
                       -(76940219*Pi*Pi*Pi*Pi)/195955200 + (9318467*log(2)*Pi*Pi*Pi*Pi)/32659200 + lmmu*(-2398621./746496.
                       -(2483663*Zeta3)/165888) + nlq*nlq*(-271883./4478976. + (8545*lmmu)/186624 - (79*lmmu*lmmu)/6912
                       +(167*Zeta3)/5184) - (2362581983*Zeta3)/87091200 + nlq*(-4770941./2239488. + (685*A4)/5184
                       +(983*lmmu*lmmu)/3456 + (107*lmmu*lmmu*lmmu)/1728 + (685*log(2)*log(2)*log(2)*log(2))/124416
                       -(685*log(2)*log(2)*Pi*Pi)/124416 - (541549*Pi*Pi*Pi*Pi)/14929920 + (3645913*Zeta3)/995328
                       +lmmu*(190283./373248. + (133819*Zeta3)/82944) + (115*Zeta5)/576) - (12057583*Zeta5)/483840);
     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsDownSI(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsDownSI(double als, double massth, double muth, int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fas5to6si(als/Pi, massth, muth, Nf,nl));
     return als*erg;
}


// Coefficients of eq.(33) of [RunDec]
double CRunDec::fmq6to5os(double A, double mass, double mu, double nlq, double nl){
     double lmM=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=0.;
     sum[2]=-A*A*(89./432. - (5.*lmM)/36. + lmM*lmM/12.);
     sum[3]=-A*A*A*(1871./2916. - (4.*A4)/9. - log(2)*log(2)*log(2)*log(2)/54. + (299.*lmM)/2592.
            + (299.*lmM*lmM)/432. + (35.*lmM*lmM*lmM)/216. + (1327.*nlq)/11664.
            - (53.*lmM*nlq)/432. - (lmM*lmM*lmM*nlq)/108. + (log(2)*log(2)*Zeta2)/9. - (407.*Zeta3)/864.
            - (5.*lmM*Zeta3)/6. - (2.*nlq*Zeta3)/27. + (103.*Zeta4)/72.);
     sum[4]=A*A*A*A*(-122722873709./3292047360. - (6941*lmM*lmM*lmM)/2592 - 
       (1147*lmM*lmM*lmM*lmM)/3456 - (16187201*Pi*Pi*Pi*Pi)/52254720 + (787*Pi*Pi*Pi*Pi*Pi*Pi)/81648 + 
       (5*Zeta2)/9 + (5*log(2)*Zeta2)/27 - (5*Zeta3)/108 + (145*Pi*Pi*Pi*Pi*log(2))/1944 - 
       (1924649*Pi*Pi*log(2)*log(2))/4354560 + (59*Pi*Pi*log(2)*log(2)*log(2))/972 + 
       (1924649*log(2)*log(2)*log(2)*log(2))/4354560 - (59*log(2)*log(2)*log(2)*log(2)*log(2))/1620 - (4*log(2)*log(2)*log(2)*log(2)*log(2)*log(2))/81 + 
       (log(2)*log(2)*log(2)*log(2)*log(2)*log(16))/81 + (1924649*A4)/181440 + 
       (118*A5)/27 + nlq*nlq*(-17671./124416. - (31*lmM*lmM)/1296 - 
         lmM*lmM*lmM*lmM/864 + (7*Pi*Pi*Pi*Pi)/8640 + lmM*(3401./46656. - (7*Zeta3)/108) + 
         (5*Zeta3)/864) + (443509931*Zeta3)/40642560 - 
       (1061*Zeta3*Zeta3)/576 + lmM*lmM*(-34297./3456. + (175*Zeta3)/48) + 
       lmM*(99319./41472. - (481*Pi*Pi*Pi*Pi)/2880 - (2*Zeta2)/3 - (2*log(2)*Zeta2)/9 + 
         Zeta3/18 - (11*Pi*Pi*log(2)*log(2))/72 + (11*log(2)*log(2)*log(2)*log(2))/72 + 
         (11*A4)/3 + (47317*Zeta3)/3072 - (575*Zeta5)/72) + 
       nlq*(2403419./746496. + (10237*lmM*lmM)/10368 + (47*lmM*lmM*lmM)/288 + 
         (17*lmM*lmM*lmM*lmM)/432 + (245*Pi*Pi*Pi*Pi)/62208 - (5*Zeta2)/54 - 
         (49*Pi*Pi*Pi*Pi*log(2))/6480 + (49*Pi*Pi*log(2)*log(2))/2592 - 
         (Pi*Pi*log(2)*log(2)*log(2))/162 - (49*log(2)*log(2)*log(2)*log(2))/2592 + log(2)*log(2)*log(2)*log(2)*log(2)/270 - 
         (49*A4)/108 - (4*A5)/9 + 
         lmM*(-26443./93312. + (163*Pi*Pi*Pi*Pi)/12960 + Zeta2/9 + (Pi*Pi*log(2)*log(2))/108 - 
           log(2)*log(2)*log(2)*log(2)/108 - (2*A4)/9 - (599*Zeta3)/1728) + 
         (1075*Zeta3)/1728 - (497*Zeta5)/288) + (59015*Zeta5)/1728);

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

double CRunDec::fmq6to5ms(double A, double mass, double mu, double nlq, double nl){
     double lmm=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=0.;
     sum[2]=A*A*(-89./432. + (5*lmm)/36 - lmm*lmm/12);
     sum[3]=A*A*A*(-2951./2916. + (16*A4+2./3.*log(2)*log(2)*log(2)*log(2)-2./3.*Pi*Pi*log(2)*log(2)-13./180.*Pi*Pi*Pi*Pi)/36
                   -(155*lmm*lmm)/432 - (35*lmm*lmm*lmm)/216 - Pi*Pi*Pi*Pi/72
                   +nl*(-1327./11664. + (53*lmm)/432 + lmm*lmm*lmm/108 + (2*Zeta3)/27)
                   +lmm*(133./2592. + (5*Zeta3)/6) + (407*Zeta3)/864);
     sum[4]=A*A*A*A*(-131621265869./3292047360. + (1924649*A4)/181440 + (118*A5)/27
                     -(3161*lmm*lmm*lmm)/2592 - (1147*lmm*lmm*lmm*lmm)/3456
                     +(1924649*log(2)*log(2)*log(2)*log(2))/4354560 - (59*log(2)*log(2)*log(2)*log(2)*log(2))/1620
                     +(log(16)*log(2)*log(2)*log(2)*log(2)*log(2))/81 - (4*log(2)*log(2)*log(2)*log(2)*log(2)*log(2))/81
                     -(1924649*log(2)*log(2)*Pi*Pi)/4354560 + (59*log(2)*log(2)*log(2)*Pi*Pi)/972
                     -(16187201*Pi*Pi*Pi*Pi)/52254720 + (145*log(2)*Pi*Pi*Pi*Pi)/1944 + (787*Pi*Pi*Pi*Pi*Pi*Pi)/81648
                     +nlq*nlq*(-17671./124416. - (31*lmm*lmm)/1296 - lmm*lmm*lmm*lmm/864 + (7*Pi*Pi*Pi*Pi)/8640
                     +lmm*(3401./46656. - (7*Zeta3)/108) + (5*Zeta3)/864) + (353193131*Zeta3)/40642560
                     -(1061*Zeta3*Zeta3)/576 + lmm*lmm*(-16193./3456. + (175*Zeta3)/48) + lmm*(279367./41472.
                     +(11*A4)/3 + (11*log(2)*log(2)*log(2)*log(2))/72 - (11*log(2)*log(2)*Pi*Pi)/72
                     -(481*Pi*Pi*Pi*Pi)/2880 + (42197*Zeta3)/3072 - (575*Zeta5)/72) + nlq*(2261435./746496.
                     -(49*A4)/108 - (4*A5)/9 + (8461*lmm*lmm)/10368 + (23*lmm*lmm*lmm)/288 + (17*lmm*lmm*lmm*lmm)/432
                     -(49*log(2)*log(2)*log(2)*log(2))/2592 + log(2)*log(2)*log(2)*log(2)*log(2)/270
                     +(49*log(2)*log(2)*Pi*Pi)/2592 - (log(2)*log(2)*log(2)*Pi*Pi)/162 + (245*Pi*Pi*Pi*Pi)/62208
                     -(49*log(2)*Pi*Pi*Pi*Pi)/6480 + lmm*(-55315./93312. - (2*A4)/9 - log(2)*log(2)*log(2)*log(2)/108
                     +(log(2)*log(2)*Pi*Pi)/108 + (163*Pi*Pi*Pi*Pi)/12960 - (599*Zeta3)/1728)
                     +(1075*Zeta3)/1728 - (497*Zeta5)/288) + (59015*Zeta5)/1728);
     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

double CRunDec::fmq6to5si(double A, double mass, double mu, double nlq, double nl){
     double lmmu=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=0.;
     sum[2]=A*A*(-89./432. + (5*lmmu)/36 - lmmu*lmmu/12);
     sum[3]=A*A*A*(-2951./2916. - (299*lmmu*lmmu)/432 - (35*lmmu*lmmu*lmmu)/216 - Pi*Pi*Pi*Pi/72 + (16*A4
                   +(2*log(2)*log(2)*log(2)*log(2))/3 - (2*log(2)*log(2)*Pi*Pi)/3 - (13*Pi*Pi*Pi*Pi)/180)/36
                   +nlq*(-1327./11664. + (53*lmmu)/432 + lmmu*lmmu*lmmu/108 + (2*Zeta3)/27) + lmmu*(853./2592.
                   +(5*Zeta3)/6) + (407*Zeta3)/864);
     sum[4]=A*A*A*A*(-131621265869./3292047360. + (1924649*A4)/181440 + (118*A5)/27 - (6941*lmmu*lmmu*lmmu)/2592
                     -(1147*lmmu*lmmu*lmmu*lmmu)/3456 + (1924649*log(2)*log(2)*log(2)*log(2))/4354560
                     -(59*log(2)*log(2)*log(2)*log(2)*log(2))/1620 + (log(16)*log(2)*log(2)*log(2)*log(2)*log(2))/81
                     -(4*log(2)*log(2)*log(2)*log(2)*log(2)*log(2))/81 - (1924649*log(2)*log(2)*Pi*Pi)/4354560
                     +(59*log(2)*log(2)*log(2)*Pi*Pi)/972 - (16187201*Pi*Pi*Pi*Pi)/52254720
                     +(145*log(2)*Pi*Pi*Pi*Pi)/1944 + (787*Pi*Pi*Pi*Pi*Pi*Pi)/81648 + nlq*nlq*(-17671./124416.
                     -(31*lmmu*lmmu)/1296 - lmmu*lmmu*lmmu*lmmu/864 + (7*Pi*Pi*Pi*Pi)/8640 + lmmu*(3401./46656.
                     -(7*Zeta3)/108) + (5*Zeta3)/864) + (353193131*Zeta3)/40642560 - (1061*Zeta3*Zeta3)/576
                     +lmmu*lmmu*(-8531./1152. + (175*Zeta3)/48) + lmmu*(330503./41472. + (11*A4)/3
                     +(11*log(2)*log(2)*log(2)*log(2))/72 - (11*log(2)*log(2)*Pi*Pi)/72 - (481*Pi*Pi*Pi*Pi)/2880
                     +(47317*Zeta3)/3072 - (575*Zeta5)/72) + nlq*(2261435./746496. - (49*A4)/108 - (4*A5)/9
                     +(8701*lmmu*lmmu)/10368 + (47*lmmu*lmmu*lmmu)/288 + (17*lmmu*lmmu*lmmu*lmmu)/432
                     -(49*log(2)*log(2)*log(2)*log(2))/2592 + log(2)*log(2)*log(2)*log(2)*log(2)/270
                     +(49*log(2)*log(2)*Pi*Pi)/2592 - (log(2)*log(2)*log(2)*Pi*Pi)/162 + (245*Pi*Pi*Pi*Pi)/62208
                     -(49*log(2)*Pi*Pi*Pi*Pi)/6480 + lmmu*(-36019./93312. - (2*A4)/9 - log(2)*log(2)*log(2)*log(2)/108
                     +(log(2)*log(2)*Pi*Pi)/108 + (163*Pi*Pi*Pi*Pi)/12960 - (599*Zeta3)/1728) + (1075*Zeta3)/1728
                     -(497*Zeta5)/288) + (59015*Zeta5)/1728);
     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecMqUpOS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecMqUpOS(double mq, double als, double massth, double muth, 
                           int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN 
     }
     double erg=(this->fmq6to5os(als/Pi, massth, muth, Nf,nl));
     return mq*erg;
}

// Function double CRunDec::DecMqUpMS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecMqUpMS(double mq, double als, double massth, double muth, 
                           int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN 
     }
     double erg=(this->fmq6to5ms(als/Pi, massth, muth, Nf,nl));
     return mq*erg;
}

// Function double CRunDec::DecMqUpSI(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecMqUpSI(double mq, double als, double massth, double muth, 
                           int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN 
     }
     double erg=(this->fmq6to5si(als/Pi, massth, muth, Nf,nl));
     return mq*erg;
}

// Coefficients of eq.(30) of [RunDec]
double CRunDec::fmq5to6os(double A, double mass, double mu, double nlq, double nl){
     double lmM=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=0.;
     sum[2]=A*A*(89./432. - (5.*lmM)/36. + lmM*lmM/12.);
     sum[3]=A*A*A*(1871./2916. - (4.*A4)/9. - ln2*ln2*ln2*ln2/54. + (121.*lmM)/2592.
                   +(319.*lmM*lmM)/432. + (29.*lmM*lmM*lmM)/216. + (1327.*nlq)/11664.
                   -(53.*lmM*nlq)/432. - (lmM*lmM*lmM*nlq)/108. + (ln2*ln2*Zeta2)/9.
                   -(407.*Zeta3)/864. - (5.*lmM*Zeta3)/6. - (2.*nlq*Zeta3)/27. + (103.*Zeta4)/72.);
     sum[4]=A*A*A*A*(122466970229./3292047360. - (1924649.*A4)/181440. - (118.*A5)/27.
                     -(1924649.*ln2*ln2*ln2*ln2)/4354560. + (59.*ln2*ln2*ln2*ln2*ln2)/1620.
                     -(1126487.*lmM)/373248. - (31.*A4*lmM)/9.-(31.*ln2*ln2*ln2*ln2*lmM)/216.
                     +(104803.*lmM*lmM)/10368. + (1403.*lmM*lmM*lmM)/648. + (305.*lmM*lmM*lmM*lmM)/1152.
                     -(2403419.*nlq)/746496. + (49.*A4*nlq)/108. + (4.*A5*nlq)/9. + (49.*ln2*ln2*ln2*ln2*nlq)/2592.
                     -(ln2*ln2*ln2*ln2*ln2*nlq)/270. + (7045.*lmM*nlq)/31104. + (2.*A4*lmM*nlq)/9.
                     +1./108.*ln2*ln2*ln2*ln2*lmM*nlq - (9601.*lmM*lmM*nlq)/10368. - (47.*lmM*lmM*lmM*nlq)/288.
                     -(5.*lmM*lmM*lmM*lmM*nlq)/144. + (17671.*nlq*nlq)/124416. - (3401.*lmM*nlq*nlq)/46656.
                     +(31.*lmM*lmM*nlq*nlq)/1296. + (lmM*lmM*lmM*lmM*nlq*nlq)/864. - (5.*Zeta2)/9.
                     -(5.*ln2*Zeta2)/27. + (1924649.*ln2*ln2*Zeta2)/725760. - (59.*ln2*ln2*ln2*Zeta2)/162.
                     +(2.*lmM*Zeta2)/3. + (2.*ln2*lmM*Zeta2)/9. + 31./36.*ln2*ln2*lmM*Zeta2 + (5.*nlq*Zeta2)/54.
                     -49./432.*ln2*ln2*nlq*Zeta2 + 1./27.*ln2*ln2*ln2*nlq*Zeta2 - (lmM*nlq*Zeta2)/9.
                     - 1./18.*ln2*ln2*lmM*nlq*Zeta2 - (441628331.*Zeta3)/40642560. - (420877.*lmM*Zeta3)/27648.
                     -(155.*lmM*lmM*Zeta3)/48. - (1075.*nlq*Zeta3)/1728. + (221.*lmM*nlq*Zeta3)/576.
                     -(5.*nlq*nlq*Zeta3)/864. + 7./108.*lmM*nlq*nlq*Zeta3 + (1061.*Zeta3*Zeta3)/576.
                     +(16187201.*Zeta4)/580608. - (725.*ln2*Zeta4)/108. + (4123.*lmM*Zeta4)/288.
                     -(1225.*nlq*Zeta4)/3456. + (49.*ln2*nlq*Zeta4)/72. - (163.*lmM*nlq*Zeta4)/144.
                     -(7.*nlq*nlq*Zeta4)/96. - (59015.*Zeta5)/1728. + (575.*lmM*Zeta5)/72.
                     +(497.*nlq*Zeta5)/288. - (3935.*Zeta6)/432.);

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

double CRunDec::fmq5to6ms(double A, double mass, double mu, double nlq, double nl){
     double lmm=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=0.;
     sum[2]=A*A*(89./432. - (5.*lmm)/36. + lmm*lmm/12.);
     sum[3]=A*A*A*(2951./2916. - (16*A4+2./3.*log(2)*log(2)*log(2)*log(2)-2./3.*Pi*Pi*log(2)*log(2)-13./180.*Pi*Pi*Pi*Pi)/36
                   +(175*lmm*lmm)/432 + (29*lmm*lmm*lmm)/216 + Pi*Pi*Pi*Pi/72
                   +lmm*(-311./2592. - (5*Zeta3)/6)
                   +nlq*(1327./11664. - (53*lmm)/432 - lmm*lmm*lmm/108 - (2*Zeta3)/27) - (407*Zeta3)/864);
     sum[4]=A*A*A*A*(131968227029./3292047360. - (1924649*A4)/181440 - (118*A5)/27 
                     +(301*lmm*lmm*lmm)/324 + (305*lmm*lmm*lmm*lmm)/1152 - (1924649*log(2)*log(2)*log(2)*log(2))/4354560 
                     +(59*log(2)*log(2)*log(2)*log(2)*log(2))/1620 - (log(16)*log(2)*log(2)*log(2)*log(2)*log(2))/81
                     +(4*log(2)*log(2)*log(2)*log(2)*log(2)*log(2))/81 + (1924649*log(2)*log(2)*Pi*Pi)/4354560
                     -(59*log(2)*log(2)*log(2)*Pi*Pi)/972 + (16187201*Pi*Pi*Pi*Pi)/52254720
                     -(145*log(2)*Pi*Pi*Pi*Pi)/1944 - (787*Pi*Pi*Pi*Pi*Pi*Pi)/81648 + lmm*lmm*(51163./10368.
                     -(155*Zeta3)/48) + nlq*nlq*(17671./124416. + (31*lmm*lmm)/1296 + lmm*lmm*lmm*lmm/864
                     -(7*Pi*Pi*Pi*Pi)/8640 + lmm*(-3401./46656. + (7*Zeta3)/108) - (5*Zeta3)/864)
                     -(353193131*Zeta3)/40642560 + (1061*Zeta3*Zeta3)/576 - (59015*Zeta5)/1728
                     +nlq*(-2261435./746496. + (49*A4)/108 + (4*A5)/9 - (7825*lmm*lmm)/10368 - (23*lmm*lmm*lmm)/288
                     -(5*lmm*lmm*lmm*lmm)/144 + (49*log(2)*log(2)*log(2)*log(2))/2592
                     -log(2)*log(2)*log(2)*log(2)*log(2)/270 - (49*log(2)*log(2)*Pi*Pi)/2592 
                     +(log(2)*log(2)*log(2)*Pi*Pi)/162 - (245*Pi*Pi*Pi*Pi)/62208 + (49*log(2)*Pi*Pi*Pi*Pi)/6480
                     +lmm*(16669./31104. + (2*A4)/9 + log(2)*log(2)*log(2)*log(2)/108 - (log(2)*log(2)*Pi*Pi)/108
                     -(163*Pi*Pi*Pi*Pi)/12960 + (221*Zeta3)/576) - (1075*Zeta3)/1728 + (497*Zeta5)/288)
                     +lmm*(-2810855./373248. - (31*A4)/9 - (31*log(2)*log(2)*log(2)*log(2))/216
                     +(31*log(2)*log(2)*Pi*Pi)/216 + (4123*Pi*Pi*Pi*Pi)/25920 - (373261*Zeta3)/27648 + (575*Zeta5)/72));

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

double CRunDec::fmq5to6si(double A, double mass, double mu, double nlq, double nl){
     double lmmu=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=0.;
     sum[2]=A*A*(89./432. - (5*lmmu)/36 + lmmu*lmmu/12);
     sum[3]=A*A*A*(2951./2916. + (319*lmmu*lmmu)/432 + (29*lmmu*lmmu*lmmu)/216 + Pi*Pi*Pi*Pi/72 + (-16*A4
                   -(2*log(2)*log(2)*log(2)*log(2))/3 + (2*log(2)*log(2)*Pi*Pi)/3 + (13*Pi*Pi*Pi*Pi)/180)/36
                   +lmmu*(-1031./2592. - (5*Zeta3)/6) + nlq*(1327./11664. - (53*lmmu)/432 - lmmu*lmmu*lmmu/108
                   -(2*Zeta3)/27) - (407*Zeta3)/864);
     sum[4]=A*A*A*A*(131968227029./3292047360. - (1924649*A4)/181440 - (118*A5)/27 + (1403*lmmu*lmmu*lmmu)/648
                     +(305*lmmu*lmmu*lmmu*lmmu)/1152 - (1924649*log(2)*log(2)*log(2)*log(2))/4354560
                     +(59*log(2)*log(2)*log(2)*log(2)*log(2))/1620 - (log(16)*log(2)*log(2)*log(2)*log(2)*log(2))/81
                     +(4*log(2)*log(2)*log(2)*log(2)*log(2)*log(2))/81 + (1924649*log(2)*log(2)*Pi*Pi)/4354560
                     -(59*log(2)*log(2)*log(2)*Pi*Pi)/972 + (16187201*Pi*Pi*Pi*Pi)/52254720 - (145*log(2)*Pi*Pi*Pi*Pi)/1944
                     -(787*Pi*Pi*Pi*Pi*Pi*Pi)/81648 + lmmu*lmmu*(81763./10368. - (155*Zeta3)/48) + nlq*nlq*(17671./124416.
                     +(31*lmmu*lmmu)/1296 + lmmu*lmmu*lmmu*lmmu/864 - (7*Pi*Pi*Pi*Pi)/8640 + lmmu*(-3401./46656.
                     +(7*Zeta3)/108) - (5*Zeta3)/864) - (353193131*Zeta3)/40642560 + (1061*Zeta3*Zeta3)/576
                     -(59015*Zeta5)/1728 + nlq*(-2261435./746496. + (49*A4)/108 + (4*A5)/9 - (8065*lmmu*lmmu)/10368
                     -(47*lmmu*lmmu*lmmu)/288 - (5*lmmu*lmmu*lmmu*lmmu)/144 + (49*log(2)*log(2)*log(2)*log(2))/2592
                     -log(2)*log(2)*log(2)*log(2)*log(2)/270 - (49*log(2)*log(2)*Pi*Pi)/2592 + (log(2)*log(2)*log(2)*Pi*Pi)/162
                     -(245*Pi*Pi*Pi*Pi)/62208 + (49*log(2)*Pi*Pi*Pi*Pi)/6480 + lmmu*(10237./31104. + (2*A4)/9
                     +log(2)*log(2)*log(2)*log(2)/108 - (log(2)*log(2)*Pi*Pi)/108 - (163*Pi*Pi*Pi*Pi)/12960 + (221*Zeta3)/576)
                     -(1075*Zeta3)/1728 + (497*Zeta5)/288) + lmmu*(-3322343./373248. - (31*A4)/9
                     -(31*log(2)*log(2)*log(2)*log(2))/216 + (31*log(2)*log(2)*Pi*Pi)/216 + (4123*Pi*Pi*Pi*Pi)/25920
                     -(419341*Zeta3)/27648 + (575*Zeta5)/72));

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecMqDownOS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecMqDownOS(double mq, double als, double massth, double muth, 
                           int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fmq5to6os(als/Pi, massth, muth, Nf,nl));
     return mq*erg;
}

// Function double CRunDec::DecMqDownMS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecMqDownMS(double mq, double als, double massth, double muth, 
                           int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fmq5to6ms(als/Pi, massth, muth, Nf,nl));
     return mq*erg;
}

// Function double CRunDec::DecMqDownSI(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecMqDownSI(double mq, double als, double massth, double muth, 
                           int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fmq5to6si(als/Pi, massth, muth, Nf,nl));
     return mq*erg;
}

// Function double CRunDec::DecLambdaDown(double lam, double massth, int nl, 
//                          int nloops)
double CRunDec::DecLambdaDown(double lam, double massth, int nl, 
                           int nloops){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     SetConstants(nl+1);
     double l=log((massth*massth)/(lam*lam));
     double sum[5];
     double b0p2 = Betap[0]*2.;
     double c2 = 11./72.;
     double c3 = (564731./124416. - 2633./31104.*(nl) - (82043.*Zeta3)/27648.);
     double c4 = (291716893./6123600. + (3031309.*A4)/54432. + (340853.*A5)/11340.
                       +(3031309.*log(2)*log(2)*log(2)*log(2))/1306368. - (340853.*log(2)*log(2)*log(2)*log(2)*log(2))/1360800.
                       -(3031309.*log(2)*log(2)*Pi*Pi)/1306368. + (340853.*log(2)*log(2)*log(2)*Pi*Pi)/816480.
                       -(76940219.*Pi*Pi*Pi*Pi)/195955200. + (9318467.*log(2)*Pi*Pi*Pi*Pi)/32659200.
                       +(nl)*(nl)*(-271883./4478976.
                       +(167.*Zeta3)/5184.) - (2362581983.*Zeta3)/87091200. + (nl)*(-4770941./2239488. + (685.*A4)/5184.
                       +(685.*log(2)*log(2)*log(2)*log(2))/124416.
                       -(685.*log(2)*log(2)*Pi*Pi)/124416. - (541549.*Pi*Pi*Pi*Pi)/14929920. + (3645913.*Zeta3)/995328.
                       +(115.*Zeta5)/576.) - (12057583.*Zeta5)/483840.);
     sum[0]= l*(Betap[0]-Beta[0])/b0p2;
     sum[1]= ((Bp[1]-B[1])*log(l) - Bp[1]*log(Betap[0]/Beta[0]))/b0p2;
     sum[2]= (c2 + B[2] - Bp[2] - B[1]*B[1] + Bp[1]*Bp[1] + log(l)*B[1]*(Bp[1]-B[1]))/(b0p2*Beta[0]*l);
     sum[3]= (c3 - (Bp[3]-B[3])/2. + Bp[1]*(Bp[2]-B[2]-c2) - Bp[1]*Bp[1]*Bp[1]/2. - B[1]*(Bp[1]*Bp[1] - Bp[2] + B[2] + c2)*log(l)
              -Bp[1]*B[1]*B[1]*(log(l)*log(l)/2. - log(l) - 1.) + B[1]*B[1]*B[1]*(log(l)*log(l)-1.)/2.)/(b0p2*Beta[0]*Beta[0]*l*l);
     sum[4]= (c4 - (Bp[4]-B[4])/3. - B[1]*B[3]/6. + Bp[1]*(2*Bp[3]/3. - B[3]/2. - c3) - c2*c2 + (Bp[2]-B[2])*(Bp[2]/3. - 2.*B[2]/3. - c2)
              -Bp[1]*Bp[1]*(Bp[2]-B[2]-c2) + Bp[1]*Bp[1]*Bp[1]*Bp[1]/3. + B[1]*log(l)*(Bp[1]*Bp[1]*Bp[1] - 2.*Bp[1]*(Bp[2] - B[2] - c2)
              + Bp[3] - B[3] - 2.*c3) + B[1]*B[1]*(Bp[1]*Bp[1] - Bp[2] + B[2] + c2)*(log(l)*log(l) - log(l) - 1.)
              + Bp[1]*B[1]*B[1]*B[1]*(log(l)*log(l)*log(l)/3.-3.*log(l)*log(l)/2. - log(l) + 1./2.) - B[1]*B[1]*B[1]*B[1]*(
              log(l)*log(l)*log(l)/3. - log(l)*log(l)/2. - log(l) - 1./6.))/(b0p2*Beta[0]*Beta[0]*Beta[0]*l*l*l);
     double erg=0.0;
     for(int i=0;i<nloops;i++){
       erg+=sum[i];
     }

     return (lam*exp(erg));
}

// Function double CRunDec::DecLambdaUp(double lam, double massth, int nl, 
//                          int nloops)
double CRunDec::DecLambdaUp(double lam, double massth, int nl, 
                           int nloops){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     SetConstants(nl+1);
     double l=log((massth*massth)/(lam*lam));
     double sum[5];
     double b02 = Beta[0]*2.;
     double c2 = 11./72.;
     double c3 = (564731./124416. - 2633./31104.*(nl) - (82043.*Zeta3)/27648.);
     double c4 = (291716893./6123600. + (3031309.*A4)/54432. + (340853.*A5)/11340.
                       +(3031309.*log(2)*log(2)*log(2)*log(2))/1306368. - (340853.*log(2)*log(2)*log(2)*log(2)*log(2))/1360800.
                       -(3031309.*log(2)*log(2)*Pi*Pi)/1306368. + (340853.*log(2)*log(2)*log(2)*Pi*Pi)/816480.
                       -(76940219.*Pi*Pi*Pi*Pi)/195955200. + (9318467.*log(2)*Pi*Pi*Pi*Pi)/32659200.
                       +(nl)*(nl)*(-271883./4478976.
                       +(167.*Zeta3)/5184.) - (2362581983.*Zeta3)/87091200. + (nl)*(-4770941./2239488. + (685.*A4)/5184.
                       +(685.*log(2)*log(2)*log(2)*log(2))/124416.
                       -(685.*log(2)*log(2)*Pi*Pi)/124416. - (541549.*Pi*Pi*Pi*Pi)/14929920. + (3645913.*Zeta3)/995328.
                       +(115.*Zeta5)/576.) - (12057583.*Zeta5)/483840.);
     sum[0]= l*(Beta[0]-Betap[0])/b02;
     sum[1]= ((B[1]-Bp[1])*log(l) - B[1]*log(Beta[0]/Betap[0]))/b02;
     sum[2]= (-c2 + Bp[2] - B[2] - Bp[1]*Bp[1] + B[1]*B[1] + log(l)*Bp[1]*(B[1]-Bp[1]))/(b02*Betap[0]*l);
     sum[3]= (-c3 + (Bp[3]-B[3])/2. + B[1]*(B[2]-Bp[2]+c2) - B[1]*B[1]*B[1]/2. - Bp[1]*(B[1]*B[1] - B[2] + Bp[2] - c2)*log(l)
              -B[1]*Bp[1]*Bp[1]*(log(l)*log(l)/2. - log(l) - 1.) + Bp[1]*Bp[1]*Bp[1]*(log(l)*log(l)-1.)/2.)/(b02*Betap[0]*Betap[0]*l*l);
     sum[4]= (-c4 + (Bp[4]-B[4])/3. - Bp[1]*Bp[3]/6. + B[1]*(2*B[3]/3. - Bp[3]/2. + c3) - c2*c2 + (B[2]-Bp[2])*(B[2]/3. - 2.*Bp[2]/3. + c2)
              -B[1]*B[1]*(B[2]-Bp[2]+c2) + B[1]*B[1]*B[1]*B[1]/3. + Bp[1]*log(l)*(B[1]*B[1]*B[1] - 2.*B[1]*(B[2] - Bp[2] + c2)
              + B[3] - Bp[3] + 2.*c3) + Bp[1]*Bp[1]*(B[1]*B[1] - B[2] + Bp[2] - c2)*(log(l)*log(l) - log(l) - 1.)
              + Bp[1]*Bp[1]*Bp[1]*B[1]*(log(l)*log(l)*log(l)/3.-3.*log(l)*log(l)/2. - log(l) + 1./2.) - Bp[1]*Bp[1]*Bp[1]*Bp[1]*(
              log(l)*log(l)*log(l)/3. - log(l)*log(l)/2. - log(l) - 1./6.))/(b02*Betap[0]*Betap[0]*Betap[0]*l*l*l);
     double erg=0.0;
     for(int i=0;i<nloops;i++){
       erg+=sum[i];
     }

     return (lam*exp(erg));
}

// Function double CRunDec::AlL2AlH(double asl, double mu1, TriplenfMmu decpar[],
//                          double mu2, int nl)
double CRunDec::AlL2AlH(double asl, double mu1, TriplenfMmu decpar[], double mu2, 
                         int nl){
     int n=0;
     int help;
     double help2;
     double asini=asl;
     double muini=mu1;
     for(int i=0; i<4; i++){
       if(decpar[i].nf!=0){
         n+=1;
       }
     }
     int zaehler=4;
     while (zaehler!=0){
       for(int i=0; i<zaehler-1; i++){
         if(decpar[i].nf>decpar[i+1].nf){
           help=decpar[i].nf;
           decpar[i].nf=decpar[i+1].nf;
           decpar[i+1].nf=help;
           help2=decpar[i].Mth;
           decpar[i].Mth=decpar[i+1].Mth;
           decpar[i+1].Mth=help2;
           help2=decpar[i].muth;
           decpar[i].muth=decpar[i+1].muth;
           decpar[i+1].muth=help2;
         }
       }
       zaehler-=1;
     }
     for(int i=3-n+2;i<=3;i++){
       if(decpar[i].nf-decpar[i-1].nf!=1){
         cout<<"WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."<<endl;
         RETURN  
       }
     }
     double erg1,erg2;
     int i;
     for(i=3-n+1;i<=3;i++){
       erg1= AlphasExact(asini, muini, decpar[i].muth,decpar[i].nf-1,nl);
       erg2= DecAsUpOS(erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       asini=erg2;
       muini=decpar[i].muth; 
     }
     double alpha= AlphasExact(asini,muini,mu2,decpar[i-1].nf,nl);
     for(int j=0;j<=3;j++){
       decpar[j].nf=0;
       decpar[j].Mth=0.;
       decpar[j].muth=0.;
     }
     return alpha;
}

// Function double CRunDec::AlH2AlL(double ash, double mu1, TriplenfMmu decpar[],
//                          double mu2, int nl)
double CRunDec::AlH2AlL(double ash, double mu1, TriplenfMmu decpar[], double mu2, 
                         int nl){
     int n=0;
     int help;
     double help2;
     double asini=ash;
     double muini=mu1;
     for(int i=0; i<4; i++){
       if(decpar[i].nf!=0){
         n+=1;
       }
     }
     int zaehler=4;
     while (zaehler!=0){
       for(int i=0; i<zaehler-1; i++){
         if(decpar[i].nf<decpar[i+1].nf){
           help=decpar[i].nf;
           decpar[i].nf=decpar[i+1].nf;
           decpar[i+1].nf=help;
           help2=decpar[i].Mth;
           decpar[i].Mth=decpar[i+1].Mth;
           decpar[i+1].Mth=help2;
           help2=decpar[i].muth;
           decpar[i].muth=decpar[i+1].muth;
           decpar[i+1].muth=help2;
         }
       }
       zaehler-=1;
     }
     for(int i=1;i<=n-1;i++){
       if(decpar[i].nf-decpar[i-1].nf!=-1){
         cout<<"WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT.";
         RETURN  
       }
     }
     double erg1,erg2;
     int i;
     for(i=0;i<=n-1;i++){
       erg1= AlphasExact(asini, muini, decpar[i].muth,decpar[i].nf,nl);
       erg2= DecAsDownOS(erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       asini=erg2;
       muini=decpar[i].muth; 
     }
     double alpha= AlphasExact(asini,muini,mu2,decpar[i-1].nf-1,nl);
     for(int j=0;j<=3;j++){
       decpar[j].nf=0;
       decpar[j].Mth=0.;
       decpar[j].muth=0.;
     }
     return alpha;
}

// Function double CRunDec::mL2mH(double mql, double asl, double mu1,
//                         TriplenfMmu decpar[], double mu2, int nl)
double CRunDec::mL2mH(double mql, double asl, double mu1, TriplenfMmu decpar[],
                       double mu2, int nl){
     int n=0;
     int help;
     double help2;
     double asini=asl;
     double muini=mu1;
     double mqini=mql;
     for(int i=0; i<4; i++){
       if(decpar[i].nf!=0){
         n+=1;
       }
     }
     int zaehler=4;
     while (zaehler!=0){
       for(int i=0; i<zaehler-1; i++){
         if(decpar[i].nf>decpar[i+1].nf){
           help=decpar[i].nf;
           decpar[i].nf=decpar[i+1].nf;
           decpar[i+1].nf=help;
           help2=decpar[i].Mth;
           decpar[i].Mth=decpar[i+1].Mth;
           decpar[i+1].Mth=help2;
           help2=decpar[i].muth;
           decpar[i].muth=decpar[i+1].muth;
           decpar[i+1].muth=help2;
         }
       }
       zaehler-=1;
     }
     for(int i=3-n+2;i<=3;i++){
       if(decpar[i].nf-decpar[i-1].nf!=1){
         cout<<"WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT.";
         RETURN  
       }
     }
     
     double erg1,erg2,erg3, erg4, erg5, erg6;
     int i;
     for(i=3-n+1;i<=3;i++){
       erg1= AlphasExact(asini, muini, decpar[i].muth,decpar[i].nf-1,nl);
       erg3= mMS2mMS(mqini,asini,erg1,decpar[i].nf-1,nl);
       erg2= DecAsUpOS(erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       erg4=DecMqUpOS(erg3,erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       asini=erg2;
       mqini=erg4;
       muini=decpar[i].muth; 
     }
     erg5= AlphasExact(asini,muini,mu2,decpar[i-1].nf,nl);
     erg6= mMS2mMS(mqini,asini,erg5,decpar[i-1].nf,nl);
     for(int j=0;j<=3;j++){
       decpar[j].nf=0;
       decpar[j].Mth=0.;
       decpar[j].muth=0.;
     }
     return erg6;
}

// Function double CRunDec::mH2mL(double mqh, double ash, double mu1,
//                         TriplenfMmu decpar[], double mu2, int nl)
double CRunDec::mH2mL(double mqh, double ash, double mu1, TriplenfMmu decpar[],
                      double mu2, int nl){
     int n=0;
     int help;
     double help2;
     double asini=ash;
     double muini=mu1;
     double mqini=mqh;
     for(int i=0; i<4; i++){
       if(decpar[i].nf!=0){
         n+=1;
       }
     }
     int zaehler=4;
     while (zaehler!=0){
       for(int i=0; i<zaehler-1; i++){
         if(decpar[i].nf<decpar[i+1].nf){
           help=decpar[i].nf;
           decpar[i].nf=decpar[i+1].nf;
           decpar[i+1].nf=help;
           help2=decpar[i].Mth;
           decpar[i].Mth=decpar[i+1].Mth;
           decpar[i+1].Mth=help2;
           help2=decpar[i].muth;
           decpar[i].muth=decpar[i+1].muth;
           decpar[i+1].muth=help2;
         }
       }
       zaehler-=1;
     }
     for(int i=1;i<=n-1;i++){
       if(decpar[i].nf-decpar[i-1].nf!=-1){
         cout<<"WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT.";
         RETURN  
       }
     }
     double erg1,erg2,erg3,erg4,erg5,erg6;
     int i;
     for( i=0;i<=n-1;i++){
       erg1= AlphasExact(asini, muini, decpar[i].muth,decpar[i].nf,nl);
       erg3= mMS2mMS(mqini,asini,erg1,decpar[i].nf,nl);
       erg2= DecAsDownOS(erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       erg4=DecMqDownOS(erg3,erg1,decpar[i].Mth,decpar[i].muth,
            decpar[i].nf-1,nl);
       asini=erg2;
       mqini=erg4;
       muini=decpar[i].muth; 
     }
     erg5= AlphasExact(asini,muini,mu2,decpar[i-1].nf-1,nl);
     erg6 =mMS2mMS(mqini,asini,erg5,decpar[i-1].nf-1,nl);
     for(int j=0;j<=3;j++){
       decpar[j].nf=0;
       decpar[j].Mth=0.;
       decpar[j].muth=0.;
     }
     return erg6;
}


// In the following the functions above are overloaded w.r.t. to an
// additional argument, n_f (number of active flavours).
// Use SetConstants(nf)

double CRunDec::LamExpl(double AlphaS, double Mu, int nf, int nl){
     SetConstants(nf);
     return (this->LamExpl(AlphaS, Mu, nl));      
}
  
double CRunDec::LamImpl(double AlphaS, double Mu,int nf,int nl){
     SetConstants(nf);
     return (this->LamImpl(AlphaS, Mu, nl));
}
  
double CRunDec::AlphasLam(double Lambda, double Mu,int nf, int nl){
     SetConstants(nf);
     return (this->AlphasLam(Lambda, Mu, nl));
}
  
double CRunDec::AlphasExact(double AlphaS0, double Mu0, double MuEnd, 
                                int nf,int nl){
     SetConstants(nf);
     return (this->AlphasExact(AlphaS0,Mu0,MuEnd,nl));
}
  
double CRunDec::mMS2mMS(double mmu0, double AlphaSEnd, double AlphaS0,int nf, 
                           int nl){
     SetConstants(nf);
     return (this->mMS2mMS(mmu0, AlphaSEnd, AlphaS0, nl));
}
  
AsmMS CRunDec::AsmMSrunexact(double mMu, double AlphaS0, double Mu0, 
                              double MuEnd,int nf, int nl){
     SetConstants(nf);
     return (this->AsmMSrunexact(mMu, AlphaS0, Mu0, MuEnd, nl));
}

double CRunDec::mMS2mOS(double mMS, std::pair<double,double>* mq, double asmu, double mu,int nf,
                           int nl, double fdelm){
     SetConstants(nf);
     return (this->mMS2mOS(mMS, mq, asmu, mu, nl, fdelm));
} 
  
double CRunDec::mOS2mMS(double mOS, std::pair<double,double>* mq, double asmu, double Mu,int nf,
                           int nl, double fdelm){
     SetConstants(nf);
     return (this->mOS2mMS(mOS, mq, asmu, Mu, nl, fdelm));
}
  
double CRunDec::mMS2mSI(double mMS, double asmu, double mu,int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mSI(mMS, asmu, mu, nl));
}
  
double CRunDec::mRI2mMS(double mRI, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mRI2mMS(mRI, asmu, nl));
}
  
double CRunDec::mMS2mRGI(double mMS, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mRGI(mMS, asmu, nl));
}
  
double CRunDec::mRGI2mMS(double mRGI, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mRGI2mMS(mRGI, asmu, nl));
}
  
double CRunDec::mOS2mSI(double mOS, std::pair<double,double>* mq, double asM,int nf, int nl, double fdelm){
     SetConstants(nf);
     return (this->mOS2mSI(mOS, mq, asM, nl, fdelm));
}

 
double CRunDec::mOS2mMSrun(double mOS, std::pair<double,double>* mq, double asmu, double mu,
                              int nf, int nl){
     SetConstants(nf);
     return (this->mOS2mMSrun(mOS, mq, asmu, mu, nl));
}
   
double CRunDec::mMS2mOSrun(double mMS, std::pair<double,double>* mq, double asmu, double mu,
                              int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mOSrun(mMS, mq, asmu, mu, nl));
}
   
double CRunDec::mMS2mRI(double mMS, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mRI(mMS, asmu, nl));
}
  
double CRunDec::mOS2mMSit(double mOS, std::pair<double,double>* mq, double asmu, double mu,
                             int nf,int nl){
     SetConstants(nf);
     return (this->mOS2mMSit(mOS, mq, asmu, mu, nl));
}
  
double CRunDec::mMS2mRGImod(double mMS, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mRGImod(mMS, asmu, nl));
}       


double CRunDec::mOS2mRS(double mOS, std::pair<double,double>* mq, double asmu, double mu, double nuf, int nl, int nloops) {
     return (this->mOS2mRS(mOS, mq, asmu, mu, nuf, nl, nloops, false));
}

double CRunDec::mOS2mRSp(double mOS, std::pair<double,double>* mq, double asmu, double mu, double nuf, int nl, int nloops) {
     return (this->mOS2mRS(mOS, mq, asmu, mu, nuf, nl, nloops, true));
}

double CRunDec::mMS2mRS(double mMS, std::pair<double,double>* mq, double asmu, double mu, double nuf, int nl, int nloops, double fdelm) {
     return (this->mMS2mRS(mMS, mq, asmu, mu, nuf, nl, nloops, fdelm, false));
}

double CRunDec::mMS2mRSp(double mMS, std::pair<double,double>* mq, double asmu, double mu, double nuf, int nl, int nloops, double fdelm) {
     return (this->mMS2mRS(mMS, mq, asmu, mu, nuf, nl, nloops, fdelm, true));
}

double CRunDec::mRS2mMS(double mRS, std::pair<double,double>* mq, double asmu, double mu, double muf, int nl, int nloops, double fdelm) {
     return (this->mRS2mMS(mRS, mq, asmu, mu, muf, nl, nloops, fdelm, false));
}

double CRunDec::mRSp2mMS(double mRS, std::pair<double,double>* mq, double asmu, double mu, double muf, int nl, int nloops, double fdelm) {
     return (this->mRS2mMS(mRS, mq, asmu, mu, muf, nl, nloops, fdelm, true));
}

double CRunDec::mRS2mSI(double mRS, std::pair<double,double>* mq, double (*as)(double), double muf, int nl, int nloops, double fdelm) {
     return (this->mRS2mSI(mRS, mq, as, muf, nl, nloops, fdelm, false));
}

double CRunDec::mRSp2mSI(double mRS, std::pair<double,double>* mq, double (*as)(double), double muf, int nl, int nloops, double fdelm) {
     return (this->mRS2mSI(mRS, mq, as, muf, nl, nloops, fdelm, true));
}

double CRunDec::DecAsDownOS(double als, double massth, double muth,int nf, 
                               int nl){
     SetConstants(nf);
     return (this->DecAsDownOS(als, massth, muth, nl));
}
  
double CRunDec::DecAsUpOS(double als, double massth, double muth, int nf,
                             int nl){
     SetConstants(nf);
     return (this->DecAsUpOS(als, massth, muth, nl)); 
}

double CRunDec::DecAsUpMS(double als, double massth, double muth, int nf,
                             int nl){
     SetConstants(nf);
     return (this->DecAsUpMS(als, massth, muth, nl)); 
}

double CRunDec::DecAsDownMS(double als, double massth, double muth,int nf, 
                               int nl){
     SetConstants(nf);
     return (this->DecAsDownMS(als, massth, muth, nl));
}

double CRunDec::DecAsUpSI(double als, double massth, double muth, int nf,
                             int nl){
     SetConstants(nf);
     return (this->DecAsUpSI(als, massth, muth, nl)); 
}

double CRunDec::DecAsDownSI(double als, double massth, double muth, int nf,
                             int nl){
     SetConstants(nf);
     return (this->DecAsDownSI(als, massth, muth, nl)); 
}
  
double CRunDec::DecMqUpOS(double mq, double als, double massth, double muth,
                             int nf, int nl){
     SetConstants(nf);
     return (this->DecMqUpOS(mq, als, massth, muth, nl));
}
  
double CRunDec::DecMqDownOS(double mq, double als, double massth, 
                               double muth,int nf, int nl){
     SetConstants(nf);
     return (this->DecMqDownOS(mq, als, massth, muth, nl));
}

double CRunDec::DecMqUpMS(double mq, double als, double massth, double muth,
                             int nf, int nl){
     SetConstants(nf);
     return (this->DecMqUpMS(mq, als, massth, muth, nl));
}
  
double CRunDec::DecMqDownMS(double mq, double als, double massth, 
                               double muth,int nf, int nl){
     SetConstants(nf);
     return (this->DecMqDownMS(mq, als, massth, muth, nl));
}

double CRunDec::DecMqUpSI(double mq, double als, double massth, double muth,
                             int nf, int nl){
     SetConstants(nf);
     return (this->DecMqUpSI(mq, als, massth, muth, nl));
}
  
double CRunDec::DecMqDownSI(double mq, double als, double massth, 
                               double muth,int nf, int nl){
     SetConstants(nf);
     return (this->DecMqDownSI(mq, als, massth, muth, nl));
}

