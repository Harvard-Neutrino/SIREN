/*
  CRunDec.h

  Header file for CRunDec.cpp

  Author: Barbara Schmidt (Jan 2012)
          Florian Herren and Matthias Steinhauser (Jan 2016)
*/

/*
  Minor update: (Sep 2016)
  Remove "static" from definition and initialization of constants for Runge-Kutta procedure.
  The C++11 standard supports the new version, so a C++11 compliant compiler and maybe
  the switch "-std=c++11" should be used.
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


#ifndef CRUNDEC_H_
#define CRUNDEC_H_

#include <utility>

// Default return statement:
#define RETURN return(0);
// The following might be useful for Windows:
//#define RETURN   system("PAUSE"); exit(1);


// Numerical values for input parameters:
static struct RunDec_values {
    double asMz    = 0.1181;
    double asMtau  = 0.332;
    double Mz      = 91.1876;
    double Mh      = 125.09;
    double muc     = 1.279;
    double mc3     = 0.986;
    double mub     = 4.163;
    double Mtau    = 1.77686;
    double Mc      = 1.5;
    double Mb      = 4.8;
    double Mt      = 173.21;
} NumDef;




// Struct for triple {nf, Mth, Muth}:

//! Structure containing: the number of light flavours (nf), the mass of the heavy quark (Mth) and the decoupling scale (muth)
/*!
  This structure is used to pass information on decoupling thresholds to AlH2AlL, AlL2AlH, mH2mL and mL2mH.
*/
struct TriplenfMmu{
    int nf;
    double Mth;
    double muth;
};


//! Structure containing: \f$\alpha_s^{(n_f)}\f$ (Asexact) and \f$m_{MS}^{(n_f)}\f$ (mMSexact)
/*!
  This structure is used to return the results of AsmMsexact.
*/  
struct AsmMS{
    double Asexact;
    double mMSexact;
};

// class declaration of CRunDec:

//! Main class, contains all functions
/*!
  All RunDec functions are accessed via this class.
  For many of the functions versions with and without the parameters nl/nf exist. 
  The ones with the parameters nl/nf call the ones without after calling SetNf(nl)/SetNf(nf). It is
  recommended not to use the ones without nl/nf since some of them change the values of nl/nf during their execution.

  Detailed descriptions are only provided for the the functions
  with parameters nl/nf.    
*/
class CRunDec
{
private: 
    // Aux. constants for implicit Runge-Kutta-Procedure:
    const double a2=0.2, a3=0.3, a4=0.6, a5=1., a6=0.875;
       
    const double b21=0.2, b31=3./40., b32=9./40., b41=0.3, b42=-0.9, 
                 b43=6./5.;
    const double b51=-11./54., b52=2.5, b53=-70./27., b54=35./27.;
    const double b61=1631./55296., b62=175./512., b63=575./13824.;
    const double b64=44275./110592., b65=253./4096.;
       
    const double c1=37./378., c2=0., c3=250./621., c4=125./594., c5=0.;
    const double c6= 512./1771.;
       
    const double dc1=37./378.-2825./27648., dc2=0.-0., 
                 dc3=250./621.-18575./48384.;
    const double dc4=125./594.-13525./55296., dc5=0.-277./14336., 
                 dc6=512./1771.-0.25;
  
    // Coefficients for diff. equations:
    double Beta[5], B[5], Betap[5], Bp[5], Gamma[5], C[5], Nf;

    // Define constants (if not already done with constructor):
    void SetConstants(int n);
  
    // R.h.s. of diff. equations:
    friend double fSetdydx(CRunDec S, double A, int nl);
  
    friend double fSetdydxa1(CRunDec S, double x, double A);
    friend double fSetdydxM1(CRunDec S, double A, double M);

    friend double fSetdydxa2(CRunDec S, double x, double A);
    friend double fSetdydxM2(CRunDec S, double A, double M);

    friend double fSetdydxa3(CRunDec S, double x, double A);
    friend double fSetdydxM3(CRunDec S, double A, double M);
  
    friend double fSetdydxa4(CRunDec S, double x, double A);
    friend double fSetdydxM4(CRunDec S, double A, double M);

    friend double fSetdydxa5(CRunDec S, double x, double A);
    friend double fSetdydxM5(CRunDec S, double A, double M);
  
    // Additional aux. functions:
    int Abbruch(void);
    double fSetAsL(double Lambda, double Mu, int nl, double AlphaS);
    double fSetcx(double x, int nl);
    double fOsFromMs1(double mu, double M);
    double fOsFromMs2(double mu, double M, double nl);
    double fOsFromMs3(double mu, double M, double nl);
    double fOsFromMs4(double mu, double M, double nl, double err);
    double fMsFromOs1(double mu, double M);
    double fMsFromOs2(double mu, double M, double nl);
    double fMsFromOs3(double mu, double M, double nl);
    double fMsFromOs4(double mu, double M, double nl, double err);
    double fZmM(double n);
    double fZmInvM(double n);
    double deltamOS2mMS(double mOS, std::pair<double,double>* mq,
                        double asmu, double mu, int nlq, int nloops); 
    double deltamMS2mOS(double mMS, std::pair<double,double>* mq,
                        double asmu, double mu, int nlq, int nloops); 
    double fMsFromRi1(void);
    double fMsFromRi2(void);
    double fMsFromRi3(void);
    double fMumFromOs1(void);
    double fMumFromOs2(void);
    double fMumFromOs3(void);
    double fMumFromOs4(double err);
    double fRiFromMs(double alpha, double nl);
    double fMsFromRi(double alpha, double nl);
    double fHelpmOS2mMSit(double mMS,double mOS, std::pair<double,double>* mq,
                          double asmu, double mu, int nl);
    double fas5to6os(double alpha, double mass, double mu, double nlq, double nl);
    double fas6to5os(double alpha, double mass, double mu, double nlq, double nl);
    double fas5to6ms(double alpha, double mass, double mu, double nlq, double nl);
    double fas6to5ms(double alpha, double mass, double mu, double nlq, double nl);
    double fas6to5si(double alpha, double mass, double mu, double nlq, double nl);
    double fas5to6si(double alpha, double mass, double mu, double nlq, double nl);
    double fmq5to6os(double A, double mass, double mu, double nlq, double nl);
    double fmq6to5os(double A, double mass, double mu, double nlq, double nl);
    double fmq5to6ms(double A, double mass, double mu, double nlq, double nl);
    double fmq6to5ms(double A, double mass, double mu, double nlq, double nl);
    double fmq5to6si(double A, double mass, double mu, double nlq, double nl);
    double fmq6to5si(double A, double mass, double mu, double nlq, double nl);

    double PSdelta(double asmu, double muf, double mu, int nl, int nloops);
    double E1p(double mOS, double asmu, double mu, int nl, int nloops);
    double exOS2RS(double api, double mmu, double nnuf, int nnl, int nloops);
    double exOS2RSp(double api, double mmu, double nnuf, int nnl, int nloops);
    double mMS2mOSmod(double mMS, std::pair<double,double>* mq,
                      double asmu, double mu, int nf, int nloops, double err);
    double mkin2mMSA(double mkin, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops);
    double mkin2mMSB(double mkin, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops);
    double mkin2mMSC(double mkin, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops);
    double mkin2mMSD(double mkin, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops);
    double mMS2mkinA(double mbMSmus, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops);
    double mMS2mkinB(double mbMSmus, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops);
    double mMS2mkinC(double mbMSmus, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops);
    double mMS2mkinD(double mbMSmus, double apinlmus, double mus, double mufac, int NLMSOS, int NLOSKIN, double mcMSmusmcin, double musmcin, int nloops);
  
    double fRungeKuttaImpl(double &x, double y, double &htry, int nl, 
                           double (*f)(CRunDec, double, int));
    double fRKSchritt(double x, double y, double h, double &yerr,
                      double (*f)(CRunDec, double, double));
    double PolyLog(unsigned int n, double x);

     
public:
    // constructor:
    CRunDec();
    CRunDec(int);
  
    // Arrays and structs to store data:
    std::pair<double,double> mq[4];
    TriplenfMmu nfMmu[4];
    AsmMS AM;
  

    //! GetNf returns the number of light flavours currently in use
    /*!
        \return \f$n_f\f$
    */
    int GetNf();

    //! SetNf sets the number of light flavours
    /*!
        \param nf \f$n_f\f$
    */
    void SetNf(int nf);

    // Functions for the running of alpha_s and m_q:

    //! LamExpl calculates \f$\Lambda^{(n_f)}\f$ from \f$\alpha_s^{(n_f)}\f$ explicitly solving for \f$\Lambda^{(n_f)}\f$
    /*!
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$\Lambda^{(n_f)}\f$
    */
    double LamExpl(double asmu, double mu, int nf, int nloops);

    //! LamImpl calculates \f$\Lambda^{(n_f)}\f$ from \f$\alpha_s^{(n_f)}\f$ implicitly solving for \f$\Lambda^{(n_f)}\f$
    /*!
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$\Lambda^{(n_f)}\f$
    */
    double LamImpl(double asmu, double mu,int nf,int nloops);

    //! AlphasLam calculates \f$\alpha_s^{(n_f)}\f$ from \f$\Lambda^{(n_f)}\f$
    /*!
        \param Lambda \f$\Lambda^{(n_f)}\f$
        \param mu \f$\mu\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$\alpha_s^{(n_f)}(\mu)\f$
    */
    double AlphasLam(double Lambda, double mu,int nf, int nloops);

    //! AlphasExact calculates \f$\alpha_s^{(n_f)}(\mu_1)\f$ from \f$\alpha_s^{(n_f)}(\mu_0)\f$
    /*!
        \param asmu0 \f$\alpha_s^{(n_f)}(\mu_0)\f$
        \param mu0 \f$\mu_0\f$
        \param mu1 \f$\mu_1\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$\alpha_s^{(n_f)}(\mu_1)\f$
    */
    double AlphasExact(double asmu0, double mu0, double mu1, int nf,int nloops);

    //! mMS2mMS calculates \f$m_{MS}^{(n_f}(\mu_1)\f$ from \f$m_{MS}^{(n_f)}(\mu_0)\f$
    /*!
        \param mu0 \f$m_{MS}^{(n_f)}(\mu_0)\f$
        \param asmu0 \f$\alpha_s^{(n_f)}(\mu_0)\f$
        \param asmu1 \f$\alpha_s^{(n_f)}(\mu_1)\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$m_{MS}^{(n_f)}(\mu_1)\f$
    */
    double mMS2mMS(double mu0, double asmu1, double asmu0,int nf, int nloops);

    //! AsmMSrunexact solves simultaneously the differential equations for \f$\alpha_s^{(n_f)}\f$ and \f$m_{MS}^{(n_f)}\f$

    //! Returns a structure of type AsmMS.
    /*!
        \param mmu \f$m_{MS}^{(n_f)}(\mu_0)\f$
        \param asmu0 \f$\alpha_s^{(n_f)}(\mu_0)\f$
        \param mu0 \f$\mu_0\f$
        \param mu1 \f$\mu_1\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$(\alpha_s^{(n_f)}(\mu_1),m_{MS}^{(n_f)}(\mu_1))\f$
    */
    AsmMS AsmMSrunexact(double mmu, double asmu0, double mu0, double mu1,
                        int nf, int nloops);

    // Decoupling relations:
    //! DecLambdaUp calculates \f$\Lambda^{(n_l+1)}\f$ from \f$\Lambda^{(n_l)}\f$
    /*!
        \param lam \f$\Lambda^{(n_l)}\f$
        \param massth quark mass at which the matching is performed (\f$m_{MS}^{(n_l+1)}(m_{MS})\f$)
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$\Lambda^{(n_l+1)}\f$
    */
    double DecLambdaUp(double lam, double massth, int nl, int nloops);

    //! DecLambdaDown calculates \f$\Lambda^{(n_l)}\f$ from \f$\Lambda^{(n_l+1)}\f$
    /*!
        \param lam \f$\Lambda^{(n_l+1)}\f$
        \param massth quark mass at which the matching is performed (\f$m_{MS}^{(n_l+1)}(m_{MS})\f$)
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$\Lambda^{(n_l)}\f$
    */
    double DecLambdaDown(double lam, double massth, int nl, int nloops);


    //! DecAsDownOS calculates \f$\alpha_s^{(n_l)}\f$ from \f$\alpha_s^{(n_l+1)}\f$
    /*!
        \param asmu \f$\alpha_s^{(n_l+1)}(\mu)\f$
        \param massth mass of the heavy quark (\f$M_{OS}\f$)
        \param muth scale \f$\mu\f$ at which the matching is performed
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$\alpha_s^{(n_l)}(\mu)\f$
    */
    double DecAsDownOS(double asmu, double massth, double muth, int nl, int nloops);

    //! DecAsUpOS calculates \f$\alpha_s^{(n_l+1)}\f$ from \f$\alpha_s^{(n_l)}\f$
    /*!
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param massth mass of the heavy quark (\f$M_{OS}\f$)
        \param muth scale \f$\mu\f$ at which the matching is performed
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$\alpha_s^{(n_l+1)}(\mu)\f$
    */
    double DecAsUpOS(double asmu, double massth, double muth, int nl, int nloops);

    //! DecAsDownMS calculates \f$\alpha_s^{(n_l)}\f$ from \f$\alpha_s^{(n_l+1)}\f$
    /*!
        \param asmu \f$\alpha_s^{(n_l+1)}(\mu)\f$
        \param massth mass of the heavy quark (\f$m_{MS}^{(n_l+1)}(\mu)\f$)
        \param muth scale \f$\mu\f$ at which the matching is performed
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$\alpha_s^{(n_l)}(\mu)\f$
    */
    double DecAsDownMS(double asmu, double massth, double muth, int nl, int nloops);

    //! DecAsUpMS calculates \f$\alpha_s^{(n_l+1)}\f$ from \f$\alpha_s^{(n_l)}\f$
    /*!
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param massth mass of the heavy quark (\f$m_{MS}^{(n_l+1)}(\mu)\f$)
        \param muth scale \f$\mu\f$ at which the matching is performed
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$\alpha_s^{(n_l+1)}(\mu)\f$
    */
    double DecAsUpMS(double asmu, double massth, double muth, int nl, int nloops);

    //! DecAsDownSI calculates \f$\alpha_s^{(n_l)}\f$ from \f$\alpha_s^{(n_l+1)}\f$
    /*!
        \param asmu \f$\alpha_s^{(n_l+1)}(\mu)\f$
        \param massth mass of the heavy quark (\f$m_{MS}^{(n_l+1)}(m_{MS})\f$)
        \param muth scale \f$\mu\f$ at which the matching is performed
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$\alpha_s^{(n_l)}(\mu)\f$
    */
    double DecAsDownSI(double asmu, double massth, double muth, int nl, int nloops);

    //! DecAsUpSI calculates \f$\alpha_s^{(n_l+1)}\f$ from \f$\alpha_s^{(n_l)}\f$
    /*!
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param massth mass of the heavy quark (\f$m_{MS}^{(n_l+1)}(m_{MS})\f$)
        \param muth scale \f$\mu\f$ at which the matching is performed
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$\alpha_s^{(n_l+1)}(\mu)\f$
    */
    double DecAsUpSI(double asmu, double massth, double muth, int nl, int nloops);

    //! DecMqUpOS calculates \f$m_{MS}^{(n_l+1)}(\mu)\f$ from \f$m_{MS}^{(n_l)}(\mu)\f$
    /*!
        \param mq \f$m_{MS}^{(n_l)}(\mu)\f$
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param massth mass of the heavy quark (\f$M_{OS}\f$)
        \param muth scale \f$\mu\f$ at which the matching is performed
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$m_{MS}^{(n_l+1)}(\mu)\f$
    */
    double DecMqUpOS(double mq, double asmu, double massth, double muth, int nl,
                     int nloops);

    //! DecMqDownOS calculates \f$m_{MS}^{(n_l)}(\mu)\f$ from \f$m_{MS}^{(n_l+1)}(\mu)\f$
    /*!
        \param mq \f$m_{MS}^{(n_l+1)}(\mu)\f$
        \param asmu \f$\alpha_s^{(n_l+1)}(\mu)\f$
        \param massth mass of the heavy quark (\f$M_{OS}\f$)
        \param muth scale \f$\mu\f$ at which the matching is performed
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$m_{MS}^{(n_l)}(\mu)\f$
    */
    double DecMqDownOS(double mq, double asmu, double massth, double muth, int nl,
                       int nloops);

    //! DecMqUpMS calculates \f$m_{MS}^{(n_l+1)}(\mu)\f$ from \f$m_{MS}^{(n_l)}(\mu)\f$
    /*!
        \param mq \f$m_{MS}^{(n_l)}(\mu)\f$
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param massth mass of the heavy quark (\f$m_{MS}^{(n_l+1)}(\mu)\f$)
        \param muth scale \f$\mu\f$ at which the matching is performed
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$m_{MS}^{(n_l+1)}(\mu)\f$
    */
    double DecMqUpMS(double mq, double asmu, double massth, double muth, int nl,
                     int nloops);

    //! DecMqDownMS calculates \f$m_{MS}^{(n_l)}(\mu)\f$ from \f$m_{MS}^{(n_l+1)}(\mu)\f$
    /*!
        \param mq \f$m_{MS}^{(n_l+1)}(\mu)\f$
        \param asmu \f$\alpha_s^{(n_l+1)}(\mu)\f$
        \param massth mass of the heavy quark (\f$m_{MS}^{(n_l+1)}(\mu)\f$)
        \param muth scale \f$\mu\f$ at which the matching is performed
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$m_{MS}^{(n_l)}(\mu)\f$
    */
    double DecMqDownMS(double mq, double asmu, double massth, double muth, int nl,
                       int nloops);

    //! DecMqUpSI calculates \f$m_{MS}^{(n_l+1)}(\mu)\f$ from \f$m_{MS}^{(n_l)}(\mu)\f$
    /*!
        \param mq \f$m_{MS}^{(n_l)}(\mu)\f$
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param massth mass of the heavy quark (\f$m_{MS}^{(n_l+1)}(m_{MS})\f$)
        \param muth scale \f$\mu\f$ at which the matching is performed
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$m_{MS}^{(n_l+1)}(\mu)\f$
    */
    double DecMqUpSI(double mq, double asmu, double massth, double muth, int nl,
                     int nloops);

    //! DecMqDownSI calculates \f$m_{MS}^{(n_l)}(\mu)\f$ from \f$m_{MS}^{(n_l+1)}(\mu)\f$
    /*!
        \param mq \f$m_{MS}^{(n_l+1)}(\mu)\f$
        \param asmu \f$\alpha_s^{(n_l+1)}(\mu)\f$
        \param massth mass of the heavy quark (\f$m_{MS}^{(n_l+1)}(m_{MS})\f$)
        \param muth scale \f$\mu\f$ at which the matching is performed
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$m_{MS}^{(n_l)}(\mu)\f$
    */
    double DecMqDownSI(double mq, double asmu, double massth, double muth, int nl,
                       int nloops);

    // Running and decoupling:

    //! AlL2AlH calculates \f$\alpha_s(\mu_2)\f$ from \f$\alpha_s(\mu_1)\f$ decoupling at intermediate scales, running from low to high
    /*!
        \param asl \f$\alpha_s(\mu_1)\f$
        \param mu1 \f$\mu_1\f$
        \param decpar arrays of triples indicating the number of flavours, the OS-mass and the scale at which decoupling is performed
        \param mu2 \f$\mu_2\f$
        \param nloops number of loops
        \return \f$\alpha_s(\mu_2)\f$
    */
    double AlL2AlH(double asl, double mu1, TriplenfMmu decpar[], double mu2, int nloops);

    //! AlH2AlL calculates \f$\alpha_s(\mu_2)\f$ from \f$\alpha_s(\mu_1)\f$ decoupling at intermediate scales, running from high to low
    /*!
        \param ash \f$\alpha_s(\mu_1)\f$
        \param mu1 \f$\mu_1\f$
        \param decpar arrays of triples indicating the number of flavours, the OS-mass and the scale at which decoupling is performed
        \param mu2 \f$\mu_2\f$
        \param nloops number of loops
        \return \f$\alpha_s(\mu_2)\f$
    */
    double AlH2AlL(double ash, double mu1, TriplenfMmu decpar[], double mu2, int nloops);

    //! mL2mH calculates \f$m_{MS}(\mu_2)\f$ from \f$m_{MS}(\mu_1)\f$ decoupling at intermediate scales, running from low to high
    /*!
        \param mql \f$m_{MS}(\mu_1)\f$
        \param asl \f$\alpha_s(\mu_1)\f$
        \param mu1 \f$\mu_1\f$
        \param decpar arrays of triples indicating the number of flavours, the OS-mass and the scale at which decoupling is performed
        \param mu2 \f$\mu_2\f$
        \param nloops number of loops
        \return \f$m_{MS}(\mu_2)\f$
    */
    double mL2mH(double mql, double asl, double mu1, TriplenfMmu decpar[], double mu2,
                 int nloops);

    //! mH2mL calculates \f$m_{MS}(\mu_2)\f$ from \f$m_{MS}(\mu_1)\f$ decoupling at intermediate scales, running from high to low
    /*!
        \param mqh \f$m_{MS}(\mu_1)\f$
        \param ash \f$\alpha_s(\mu_1)\f$
        \param mu1 \f$\mu_1\f$
        \param decpar arrays of triples indicating the number of flavours, the OS-mass and the scale at which decoupling is performed
        \param mu2 \f$\mu_2\f$
        \param nloops number of loops
        \return \f$m_{MS}(\mu_2)\f$
    */
    double mH2mL(double mqh, double ash, double mu1, TriplenfMmu decpar[], double mu2,
                 int nloops);
  
    // Mass relations
    
    //! mMS2mOS calculates \f$M_{OS}\f$ from \f$m_{MS}^{(n_f)}(\mu)\f$
    /*!
        \param mMS \f$m_{MS}^{(n_f)}(\mu)\f$
        \param mq pointer to pairs of light quark masses
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$M_{OS}\f$
    */

    double mMS2mOS(double mMS, std::pair<double,double>* mq,
                   double asmu, double mu,int nf, int nloops, double fdelm=1.0);

    //! mOS2mMS calculates \f$m_{MS}^{(n_f)}(\mu)\f$ from \f$M_{OS}\f$
    /*!
        \param mOS \f$M_{OS}\f$
        \param mq pointer to pairs of light quark masses
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m_{MS}^{(n_f)}(\mu)\f$
    */
    double mOS2mMS(double mOS, std::pair<double,double>* mq,
                   double asmu, double mu,int nf, int nloops, double fdelm=1.0);

    //! mMS2mSI calculates \f$m_{MS}^{(n_f)}(m_{MS})\f$ from \f$m_{MS}^{(n_f)}(\mu)\f$
    /*!
        \param mMS \f$m_{MS}^{(n_f)}(\mu)\f$
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$m_{MS}^{(n_f)}(m_{MS})\f$
    */
    double mMS2mSI(double mMS, double asmu, double mu, int nf, int nloops);

    //! mRI2mMS calculates \f$m_{MS}^{(n_f)}(\mu)\f$ from \f$m^{RI}(\mu)\f$
    /*!
        \param mRI \f$m^{RI}(\mu)\f$
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$m_{MS}^{(n_f)}(\mu)\f$
    */
    double mRI2mMS(double mRI, double asmu, int nf, int nloops);

    //! mMS2mRGI calculates \f$m^{RGI}\f$ from \f$m_{MS}^{(n_f)}(\mu)\f$

    //! for the difference between mMS2mRGImod and mMS2mRGI see arXiv:1201:6149
    /*!
        \param mMS \f$m_{MS}^{(n_f)}(\mu)\f$
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$m^{RGI}\f$
    */
    double mMS2mRGI(double mMS, double asmu, int nf, int nloops);

    //! mRGI2mMS calculates \f$m_{MS}^{(n_f)}(\mu)\f$ from \f$m^{RGI}\f$
    /*!
        \param mRGI \f$m^{RGI}\f$
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$m_{MS}^{(n_f)}(\mu)\f$
    */
    double mRGI2mMS(double mRGI, double asmu, int nf, int nloops);

    //! mOS2mSI calculates \f$m_{MS}^{(n_f)}(m_{MS})\f$ from \f$M_{OS}\f$
    /*!
        \param mOS \f$M_{OS}\f$
        \param mq pointer to pairs of light quark masses
        \param asM \f$\alpha_s^{(n_f)}(M_{OS})\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m_{MS}^{(n_f)}(m_{MS})\f$
    */
    double mOS2mSI(double mOS, std::pair<double,double>* mq,
                   double asM, int nf, int nloops, double fdelm=1.0);

    //! mOS2mMSrun calculates \f$m_{MS}^{(n_f)}(\mu)\f$ from \f$M_{OS}\f$

    //! Calculates \f$m_{MS}^{(n_f)}(m_{MS})\f$ in an intermediate step
    /*!
        \param mOS \f$M_{OS}\f$
        \param mq pointer to pairs of light quark masses
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$m_{MS}^{(n_f)}(\mu)\f$
    */ 
    double mOS2mMSrun(double mOS, std::pair<double,double>* mq,
                      double asmu, double mu,int nf, int nloops);

    //! mMS2mOSrun calculates \f$M_{OS}\f$ from \f$m_{MS}^{(n_f)}(\mu)\f$

    //! Calculates \f$m_{MS}^{(n_f)}(m_{MS})\f$ in an intermediate step
    /*!
        \param mMS \f$m_{MS}^{(n_f)}(\mu)\f$
        \param mq pointer to pairs of light quark masses
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$M_{OS}\f$
    */  
    double mMS2mOSrun(double mMS, std::pair<double,double>* mq,
                      double asmu, double mu,int nf, int nloops);

    //! mMS2mRI calculates \f$m^{RI}(\mu)\f$ from \f$m_{MS}^{(n_f)}(\mu)\f$
    /*!
        \param mMS \f$m_{MS}^{(n_f)}(\mu)\f$
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$m^{RI}(\mu)\f$
    */  
    double mMS2mRI(double mMS, double asmu,int nf, int nloops);

    //! mOS2mMSit calculates \f$m_{MS}^{(n_f)}(\mu)\f$ from \f$M_{OS}\f$

    //! DEPRECATED, will be removed in future versions
    /*!
        \param mOS \f$M_{OS}\f$
        \param mq pointer to pairs of light quark masses
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$m_{MS}^{(n_f)}(\mu)\f$
    */  
    double mOS2mMSit(double mOS, std::pair<double,double>* mq,
                     double asmu, double mu, int nf, int nloops);

    //! mMS2mRGImod calculates \f$m^{RGI}\f$ from \f$m_{MS}^{(n_f)}(\mu)\f$

    //! for the difference between mMS2mRGImod and mMS2mRGI see arXiv:1201:6149
    /*!
        \param mMS \f$m_{MS}^{(n_f)}(\mu)\f$
        \param asmu \f$\alpha_s^{(n_f)}(\mu)\f$
        \param nf \f$n_f\f$
        \param nloops number of loops
        \return \f$m^{RGI}\f$
    */
    double mMS2mRGImod(double mMS, double asmu, int nf, int nloops);

    //! mOS2mPS calculates \f$m^{PS}(\mu_f)\f$ from \f$M_{OS}\f$
    /*!
        \param mOS \f$M_{OS}\f$

        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param muf \f$\mu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$m^{PS}(\mu_f)\f$
    */ 
    double mOS2mPS(double mOS, std::pair<double,double>* mq,
                   double asmu, double mu, double muf, int nl, int nloops);

    //! mMS2mPS calculates \f$m^{PS}(\mu_f)\f$ from \f$m_{MS}^{(n_l+1)}(\mu)\f$
    /*!
        \param mMS \f$m_{MS}^{(n_l+1)}(\mu)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param muf \f$\mu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m^{PS}(\mu_f)\f$
    */
    double mMS2mPS(double mMS, std::pair<double,double>* mq,
                   double asmu, double mu, double muf, int nl, int nloops, double fdelm=1.0);

    //! mPS2mMS calculates \f$m_{MS}^{(n_l+1)}(\mu)\f$ from \f$m^{PS}(\mu_f)\f$
    /*!
        \param mPS \f$m^{PS}(\mu_f)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param muf \f$\mu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m_{MS}^{(n_l+1)}(\mu)\f$
    */
    double mPS2mMS(double mPS, std::pair<double,double>* mq,
                   double asmu, double mu, double muf, int nl, int nloops, double fdelm=1.0);

    //! mPS2mSI calculates \f$m_{MS}^{(n_l+1)}(m_{MS})\f$ from \f$m^{PS}(\mu_f)\f$
    /*!
        \param mPS \f$m^{PS}(\mu_f)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param as pointer to a function which calculates \f$\alpha_s^{(n_l)}(\mu)\f$ and takes \f$\mu\f$ as argument
        \param muf \f$\mu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m_{MS}^{(n_l+1)}(m_{MS})\f$
    */
    double mPS2mSI(double mPS, std::pair<double,double>* mq,
                   double (*as)(double), double muf, int nl, int nloops, double fdelm=1.0);

    //! mOS2m1S calculates \f$m^{1S}\f$ from \f$M_{OS}\f$
    /*!
        \param mOS \f$M_{OS}\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$m^{1S}\f$
    */ 
    double mOS2m1S(double mOS, std::pair<double,double>* mq,
                   double asmu, double mu, int nl, int nloops);

    //! mMS2m1S calculates \f$m^{1S}\f$ from \f$m_{MS}^{(n_l+1)}(\mu)\f$
    /*!
        \param mMS \f$m_{MS}^{(n_l+1)}(\mu)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m^{1S}\f$
    */
    double mMS2m1S(double mMS, std::pair<double,double>* mq,
                   double asmu, double mu, int nl, int nloops, double fdelm=1.0);

    //! m1S2mMS calculates \f$m_{MS}^{(n_l+1)}(\mu)\f$ from \f$m^{1S}\f$
    /*!
        \param m1S \f$m^{PS}\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m_{MS}^{(n_l+1)}(\mu)\f$
    */
    double m1S2mMS(double m1S, std::pair<double,double>* mq,
                   double asmu, double mu, int nl, int nloops, double fdelm=1.0);

    //! m1S2mSI calculates \f$m_{MS}^{(n_l+1)}(m_{MS})\f$ from \f$m^{1S}\f$
    /*!
        \param m1S \f$m^{1S}\f$
        \param mq pointer to pairs of light quark masses, not active
        \param as pointer to a function which calculates \f$\alpha_s^{(n_l)}(\mu)\f$ and takes \f$\mu\f$ as argument
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m_{MS}^{(n_l+1)}(m_{MS})\f$
    */
    double m1S2mSI(double m1S, std::pair<double,double>* mq,
                   double (*as)(double), int nl, int nloops, double fdelm=1.0);

    //! mOS2mRS calculates \f$m^{RS}(\nu_f)\f$ or \f$m^{RS'}(\nu_f)\f$ from \f$M_{OS}\f$

    //! This is a helper function, use the version without the parameter prime or mOS2mRSp
    /*!
        \param mOS \f$M_{OS}\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nuf \f$\nu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param prime selects if \f$m^{RS}\f$ or \f$m^{RS'}\f$ should be calculated
        \return \f$m^{RS}(\nu_f)\f$ or \f$m^{RS'}(\nu_f)\f$
    */ 
    double mOS2mRS(double mOS, std::pair<double,double>* mq, double asmu,
                   double mu, double nuf, int nl, int nloops, bool prime);

    //! mMS2mRS calculates \f$m^{RS}(\nu_f)\f$ or \f$m^{RS'}(\nu_f)\f$ from \f$m_{MS}^{(n_l+1)}(\mu)\f$

    //! This is a helper function, use the version without the parameter prime or mMS2mRSp
    /*!
        \param mMS \f$m_{MS}^{(n_l+1)}(\mu)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nuf \f$\nu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \param prime selects if \f$m^{RS}\f$ or \f$m^{RS'}\f$ should be calculated
        \return \f$m^{RS}(\nu_f)\f$ or \f$m^{RS'}(\nu_f)\f$
    */ 
    double mMS2mRS(double mMS, std::pair<double,double>* mq, double asmu,
                   double mu, double nuf, int nl, int nloops, double fdelm, bool prime);

    //! mRS2mMS calculates \f$m_{MS}^{(n_l+1)}(\mu)\f$ from \f$m^{RS}(\nu_f)\f$ or \f$m^{RS'}(\nu_f)\f$

    //! This is a helper function, use the version without the parameter prime or mRSp2mMS
    /*!
        \param mRS \f$m^{RS}(\nu_f)\f$ or \f$m^{RS'}(\nu_f)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nuf \f$\nu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \param prime selects if \f$m^{RS}\f$ or \f$m^{RS'}\f$ is given
        \return \f$m_{MS}^{(n_l+1)}(\mu)\f$
    */ 
    double mRS2mMS(double mRS, std::pair<double,double>* mq, double asmu,
                   double mu, double nuf, int nl, int nloops, double fdelm, bool prime);

    //! mRS2mSI calculates \f$m_{MS}^{(n_l+1)}(m_{MS})\f$ from \f$m^{RS}(\nu_f)\f$ or \f$m^{RS'}(\nu_f)\f$

    //! This is a helper function, use the version without the parameter prime or mRSp2mSI
    /*!
        \param mRS \f$m^{RS}(\nu_f)\f$ or \f$m^{RS'}(\nu_f)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param as pointer to a function which calculates \f$\alpha_s^{(n_l)}(\mu)\f$ and takes \f$\mu\f$ as argument
        \param nuf \f$\nu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \param prime selects if \f$m^{RS}\f$ or \f$m^{RS'}\f$ is given
        \return \f$m_{MS}^{(n_l+1)}(m_{MS})\f$
    */
    double mRS2mSI(double mRS, std::pair<double,double>* mq, double (*as)(double),
                   double nuf, int nl, int nloops, double fdelm, bool prime);


    //! mOS2mRS calculates \f$m^{RS}(\nu_f)\f$ from \f$M_{OS}\f$

    /*!
        \param mOS \f$M_{OS}\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nuf \f$\nu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$m^{RS}(\nu_f)\f$
    */ 
    double mOS2mRS(double mOS, std::pair<double,double>* mq,
                   double asmu, double mu, double nuf, int nl, int nloops);

    //! mMS2mRS calculates \f$m^{RS}(\nu_f)\f$ from \f$m_{MS}^{(n_l+1)}(\mu)\f$

    /*!
        \param mMS \f$m_{MS}^{(n_l+1)}(\mu)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nuf \f$\nu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m^{RS}(\nu_f)\f$
    */ 
    double mMS2mRS(double mMS, std::pair<double,double>* mq,
                   double asmu, double mu, double nuf, int nl, int nloops, double fdelm=1.0);

    //! mRS2mMS calculates \f$m_{MS}^{(n_l+1)}(\mu)\f$ from \f$m^{RS}(\nu_f)\f$

    /*!
        \param mRS \f$m^{RS}(\nu_f)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nuf \f$\nu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m_{MS}^{(n_l+1)}(\mu)\f$
    */ 
    double mRS2mMS(double mRS, std::pair<double,double>* mq,
                   double asmu, double mu, double nuf, int nl, int nloops, double fdelm=1.0);

    //! mRS2mSI calculates \f$m_{MS}^{(n_l+1)}(m_{MS})\f$ from \f$m^{RS}(\nu_f)\f$

    /*!
        \param mRS \f$m^{RS}(\nu_f)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param as pointer to a function which calculates \f$\alpha_s^{(n_l)}(\mu)\f$ and takes \f$\mu\f$ as argument
        \param nuf \f$\nu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m_{MS}^{(n_l+1)}(m_{MS})\f$
    */
    double mRS2mSI(double mRS, std::pair<double,double>* mq,
                   double (*as)(double), double nuf, int nl, int nloops, double fdelm=1.0);

    //! mOS2mRSp calculates \f$m^{RS'}(\nu_f)\f$ from \f$M_{OS}\f$

    /*!
        \param mOS \f$M_{OS}\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nuf \f$\nu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \return \f$m^{RS'}(\nu_f)\f$
    */ 
    double mOS2mRSp(double mOS, std::pair<double,double>* mq,
                    double asmu, double mu, double nuf, int nl, int nloops);

    //! mMS2mRSp calculates \f$m^{RS'}(\nu_f)\f$ from \f$m_{MS}^{(n_l+1)}(\mu)\f$

    /*!
        \param mMS \f$m_{MS}^{(n_l+1)}(\mu)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nuf \f$\nu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m^{RS'}(\nu_f)\f$
    */ 
    double mMS2mRSp(double mMS, std::pair<double,double>* mq,
                    double asmu, double mu, double nuf, int nl, int nloops, double fdelm=1.0);

    //! mRSp2mMS calculates \f$m_{MS}^{(n_l+1)}(\mu)\f$ from \f$m^{RS'}(\nu_f)\f$

    /*!
        \param mRS \f$m^{RS'}(\nu_f)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param asmu \f$\alpha_s^{(n_l)}(\mu)\f$
        \param mu \f$\mu\f$
        \param nuf \f$\nu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m_{MS}^{(n_l+1)}(\mu)\f$
    */ 
    double mRSp2mMS(double mRS, std::pair<double,double>* mq,
                    double asmu, double mu, double nuf, int nl, int nloops, double fdelm=1.0);

    //! mRSp2mSI calculates \f$m_{MS}^{(n_l+1)}(m_{MS})\f$ from \f$m^{RS'}(\nu_f)\f$

    /*!
        \param mRS \f$m^{RS'}(\nu_f)\f$
        \param mq pointer to pairs of light quark masses, not active
        \param as pointer to a function which calculates \f$\alpha_s^{(n_l)}(\mu)\f$ and takes \f$\mu\f$ as argument
        \param nuf \f$\nu_f\f$
        \param nl \f$n_l\f$
        \param nloops number of loops
        \param fdelm factor multiplying the non-logarithmic part of the 4-loop term
        \return \f$m_{MS}^{(n_l+1)}(m_{MS})\f$
    */
    double mRSp2mSI(double mRS, std::pair<double,double>* mq,
                    double (*as)(double), double nuf, int nl, int nloops, double fdelm=1.0);



    //! mMS2mKIN calculates \f$m^{\mathrm{kin}}\f$ from \f$m_{MS}^{(n_f)}(\mu_s)\f$

    /*!
        \param mMS \f$m_{MS}^{(n_f)}(\mu_s)\f$
        \param mq pointer to pairs of light quark masses
        \param asmus \f$\alpha_s^{(n_l)}(\mus)\f$
        \param mus \f$\mu_s\f$
        \param Wilsonian cut off muf \f$\mu\f$
        \param nlmsos \f$n_l\f$ in the MS-OS relation
        \param nloskin \f$n_l\f$ in the OS-KIN relation
        \param nloops number of loops
        \param deccase Specify if there is a light but massive quark
        \return \f$m^{\mathrm{kin}}\f$
    */ 
    double mMS2mKIN(double mMS, std::pair<double,double>* mq,
                    double asmus, double mus, double muf, int nlmsos, int nloskin, int nloops, std::string deccase);

    //! mKIN2mMS calculates \f$m_{MS}^{(n_f)}(\mu_s)\f$ from \f$m^{\mathrm{kin}}\f$

    /*!
        \param mKIN \f$m^{\mathrm{kin}}\f$
        \param mq pointer to pairs of light quark masses
        \param asmus \f$\alpha_s^{(n_l)}(\mus)\f$
        \param mus \f$\mu_s\f$
        \param Wilsonian cut off muf \f$\mu\f$
        \param nlmsos \f$n_l\f$ in the MS-OS relation
        \param nloskin \f$n_l\f$ in the OS-KIN relation
        \param nloops number of loops
        \param deccase Specify if there is a light but massive quark
        \return \f$m_{MS}^{(n_f)}(\mu_s)\f$
    */ 
    double mKIN2mMS(double mKIN, std::pair<double,double>* mq,
                    double asmus, double mus, double muf, int nlmsos, int nloskin, int nloops, std::string deccase);
  
    // Overload functions:
    double LamExpl(double asmu, double mu, int nloops);
    double LamImpl(double asmu, double mu,int nloops);
    double AlphasLam(double Lambda, double mu, int nloops);
    double AlphasExact(double asmu0, double mu0, double mu1, int nloops);
    double mMS2mMS(double mmu0, double asmu0, double asmu1, int nloops);
    AsmMS AsmMSrunexact(double mmu, double asmu0, double mu0, double mu1, 
                      int nloops);

    double DecAsDownOS(double asmu, double massth, double muth, int nloops);
    double DecAsUpOS(double asmu, double massth, double muth, int nloops);
    double DecAsDownMS(double asmu, double massth, double muth, int nloops);
    double DecAsUpMS(double asmu, double massth, double muth, int nloops);
    double DecAsDownSI(double asmu, double massth, double muth, int nloops);
    double DecAsUpSI(double asmu, double massth, double muth, int nloops);

    double DecMqUpOS(double mq, double asmu, double massth, double muth, int nloops);
    double DecMqDownOS(double mq, double asmu, double massth, double muth, int nloops);
    double DecMqUpMS(double mq, double asmu, double massth, double muth, int nloops);
    double DecMqDownMS(double mq, double asmu, double massth, double muth, int nloops);
    double DecMqUpSI(double mq, double asmu, double massth, double muth, int nloops);
    double DecMqDownSI(double mq, double asmu, double massth, double muth, int nloops);    


    double mMS2mOS(double mMS, std::pair<double,double>* mq,
                   double asmu, double mu, int nloops, double fdelm=1.0);
    double mOS2mMS(double mOS, std::pair<double,double>* mq,
                   double asmu, double mu, int nloops, double fdelm=1.0);

    double mMS2mSI(double mMS, double asmu, double mu, int nloops);
    double mRI2mMS(double mRI, double asmu, int nloops);
    double mMS2mRGI(double mMS, double asmu, int nloops);
    double mRGI2mMS(double mRGI, double asmu, int nloops);
    double mOS2mSI(double mOS, std::pair<double,double>* mq,
                   double asM, int nloops, double fdelm=1.0);
    double mOS2mMSrun(double mOS, std::pair<double,double>* mq,
                      double asmu, double mu, int nloops);
    double mMS2mOSrun(double mMS, std::pair<double,double>* mq,
                      double asmu, double mu, int nloops);
    double mMS2mRI(double mMS, double asmu, int nloops);
    double mOS2mMSit(double mOS, std::pair<double,double>* mq,
                     double asmu, double mu, int nloops); 
    double mMS2mRGImod(double mMS, double asmu, int nloops);

    double mMS2mKIN(double mMS, std::pair<double,double>* mq, double asmus, double mus, double muf, int nloops, std::string deccase);
    double mKIN2mMS(double mKIN, std::pair<double,double>* mq, double asmus, double mus, double muf, int nloops, std::string deccase);

};

#endif
