(* ***
 RunDec: a Mathematica package for running and decoupling of the
 strong coupling and quark masses

 First version: K.G. Chetyrkin, J.H. Kuehn and M. Steinhauser (Jan. 2000)

 v2.0: F. Herren and M. Steinhauser (Jan. 2016)
 v2.1: F. Herren and M. Steinhauser (Jun. 2016)
 v2.2: F. Herren and M. Steinhauser (Jan. 2017)
 v3.0: F. Herren and M. Steinhauser (Feb. 2017)
 v3.1: F. Herren and M. Steinhauser (Jul. 2020)
         
 Improvements and extensions for v3.1: 
 - Coefficients ctil[3,{3,4,5}] updated.
   Thanks to Christopher Lepenik for providing the new values.
 - Relation between kinetic and MSbar mass implemented.

 Improvements and extensions for v3.0:
 - five-loop beta function implemented (arXiv:1606.08659)
   Updated: AlpahsExact, mMS2mMS,     
            Dec{As,Mq}Up{OS,MS,SI}, Dec{As,Mq}Down{OS,MS,SI}, 
            AlphasLam, LamImpl, LamExpl, DecLambdaUp, DecLambdaDown
 - remaining four-loop decoupling terms of \zeta_\alpha_s and \zeta_m implemented
 - several minor improvements and bug fixes
 - improved version of mMS2mSI[]: mMS2mSInew[] -> mMS2mSI[]
   and old version renamed to mMS2mSIold[]
 - AsmMSrunexact[] implemented
 - five-loop result for \gamma_m implemented (arXiv:1402.6611)
   preparations for implementation of five-loop beta function 
 - four-loop decoupling terms of \zeta_\alpha_s and \zeta_m implemented
   (only for on-shell heavy quark masses), 
   (hep-ph/0512058, hep-ph/0512060,  arXiv:1502.04719)
 - mMS2OS, mOS2MS, mOS2SI extedned to four loops (arXiv:1502.01030)
 - relations between PS, 1S, RS, RS' and OS, MS, SI masses implemented to N3LO
   (arXiv:1502.01030)

References:

  K.G. Chetyrkin, J.H. Kuhn and M. Steinhauser,
  "RunDec: A Mathematica package for running and decoupling of the strong
  coupling and quark masses"
  Comput. Phys. Commun.  133 (2000) 43,
  arXiv:hep-ph/0004189

  "CRunDec: a C++ package for running and decoupling of the
  strong coupling and quark masses"
  B. Schmidt, M. Steinhauser
  Comput.Phys.Commun. 183 (2012) 1845-1848,
  arXiv:1201.6149, SFB/CPP-12-03, TTP12-02

  "Version 3 of RunDec and CRunDec"
  F. Herren, M. Steinhauser
  arXiv:1703.03751, TTP17-011

*** *)

(* ************************************************************ *)
(* ************************************************************ *)

(* v1.0: The examples shown in the paper can be found in the the file
   RunDec_ex.m *)

(* ************************************************************ *)

Print["RunDec: a Mathematica package for running and decoupling of the"];
Print["        strong coupling and quark masses"];
Print["by K.G. Chetyrkin, J.H. K\\\"uhn and M. Steinhauser (January 2000)"];
Print["by F. Herren and M. Steinhauser (July 2020, v3.1)"];

(* ************************************************************ *)
(* ************************************************************ *)

(*
 default values of the numerical constants
*)

NumDef = {
    asMz   -> 0.1181,
    asMtau -> 0.332,
    Mz     -> 91.1876,
    Mh     -> 125.09,
    muc    -> 1.279,
    mc3    -> 0.986,
    mub    -> 4.163,
    Mtau   -> 1.77686,
    Mc     -> 1.5,
    Mb     -> 4.8,
    Mt     -> 173.21
    };

(* ************************************************************ *)
(* ************************************************************ *)

BeginPackage["RunDec`"];

 fdelm::usage = 
    "Factor in front of 4-loop non-logarithmic coefficient; introduced
in mMS2mOS[], mOS2mMS[] and mOS2mSI[]. This option is also available
in the relations between the MS or SI and thresold masses."

 LamExpl::usage =
    "LamExpl[als, mu, nf, l] computes \\Lambda^(nf) with l-loop accuracy
using the explicite formulae \\ln(\\mu/\\Lambda) = f(als,mu,nf).";

 LamImpl::usage =
    "LamImpl[als, mu, nf, l] computes \\Lambda^(nf) with l-loop accuracy
using the implicite formulae \\alpha_s = f(mu,\\Lambda,nf).";

 AlphasLam::usage = 
    "AlphasLam[lam ,mu ,nf ,l] computes \\alpha_s^(nf)(mu) to l-loop accuracy
using the formulae \\alpha_s = f(mu,\\Lambda,nf).";

 AlphasExact::usage = 
    "AlphasExact[als0 ,mu0 ,mu ,nf ,l] computes \\alpha_s^(nf)(mu) integrating
the \\beta function numerically. The initial condition \\alpha_s(mu0)=als0
is used.";

 mOS2mMS::usage =
    "mOS2mMS[mOS, {mq}, asmu, mu, nf, l] computes the MS-bar mass at the
scale mu. {mq} represents a set of light quark masses in the on-shell 
scheme leading to corrections of O(as^2 mq/mOS).
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.
mOS2mMS[mOS, nf, l] computes the MS-bar mass at the the scale mOS. 
\\alpha_s^(nf)(M) is computed from \\alpha_s^(5)(Mz) as defined in NumDef.";

 mMS2mOS::usage =
    "mMS2mOS[mMS, {mq}, asmu, mu, nf, l] computes the on-shell mass.
{mq} represents a set of light quark masses in the on-shell scheme
leading to corrections of O(as^2 mq/mMS).
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.
mMS2mOS[mum, nf, l] computes the on-shell mass where mum is 
defined through mum=m^(nf)(mum).
\\alpha_s^(nf)(mum) is computed from \\alpha_s^(5)(Mz) as defined in NumDef.";

 mOS2mMSrun::usage =
    "mOS2mMSrun[mOS, {mq}, asmu, mu, nf, l] computes the MS-bar mass at the
scale mu.
In contrast to mOS2mMS[..] in a first step mMS(mMS) is computed
and only then the conversion to mOS is performed.
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mOSrun::usage =
    "mMS2mOSrun[mMS, {mq}, asmu, mu, nf, l] computes the on-shell mass.
In contrast to mMS2mOS[..] in a first step mMS(mMS) is computed
and afterwards mMS(mu) is evaluated.
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mOS2mMSit::usage =
    "mOS2mMSit[mOS, {mq}, asmu, mu, nf, l] computes the MS-bar mass at the
scale mu. However, in contrast to mOS2mMS[..], the relation
M_OS = m_MS * (1 + ...) is solved iteratively. This has the advantage
that the light quark masses can be evaluated in the MS-bar scheme.";

 mOS2mSI::usage =
    "mOS2mSI[mOS, {mq}, asM, nf, l] computes the mass
\\mu_m = m_MS(\\mu_m). {mq} represents a set of light quark masses in 
the on-shell scheme leading to corrections of O(as^2 mq/mOS).
asM = \\alpha_s^(nf)(M_OS), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mMS::usage =
    "mMS2mMS[mmu0, asmu0, asmu, nf, l] computes m(mu) from the knowledge
of m(mu0). asmu = \\alpha_s(\\mu), asmu0 = \\alpha_s(\\mu0), nf is the number 
of active flavours and l represents the number of loops.";

 mMS2mSI::usage =
    "mMS2mSI[mMS, asmu, mu, nf, l] computes the scale
invarant mass mMS(mMS). mMS is the MS-bar mass at the scale mu,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mSIold::usage = 
    "mMS2mSIold[mMS, asmu, mu, nf, l] old version of mMS2mSI[];
uses AlphasLam[]";

 mMS2mRI::usage =
    "mMS2mRI[mMS, asmu, nf, l] computes the regularization invariant
mass. mMS is the MS-bar mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mRI2mMS::usage =
    "mRI2mMS[mRI, asmu, nf, l] computes the MS-bar quark mass.
mRI is the regularization invariant mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mRGI::usage =
    "mMS2mRGI[mMS, asmu, nf, l] computes the renormalization group
invariant mass: mRGI = mMS/c(alphas/Pi). mMS is the MS-bar mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mRGImod::usage =
    "mMS2mRGImod[mMS, asmu, nf, l] computes the renormalization group
invariant mass defined through: mRGI = mMS/c(2*beta0*alphas/Pi). mMS is the MS-bar mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.
(slightly modified version of mMS2mRGI[], see factor '2 beta_0')";

 mRGI2mMS::usage =
    "mRGI2mMS[mRGI, asmu, nf, l] computes the MS-bar quark mass.
mRGI is the renormalization group invariant mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mOS2mPS::usage =
    "mOS2mPS[mOS, mq, asmu, mu, muf, nl, loops] computes the PS mass mPS(muf) at 
the scale mu. {mq} represents a set of light quark masses in the on-shell 
scheme leading to corrections of O(as^2 mq/mOS) (not yet implemented).
asmu = \\alpha_s^(nl)(\\mu) and loops represents the number of loops.";

 mMS2mPS::usage =
    "mMS2mPS[mMS, mq, asmu, mu, muf, nl, loops, fdelm] 
computes the PS mass mPS(muf) at 
the scale mu. mMS is the MSbar mass at the scale mu. 
asmu = \\alpha_s^(nl)(\\mu), nf=nl+1 is the number
of active flavours and loops represents the number of loops.
fdelm (=1) can be omitted.";

 mPS2mMS::usage =
    "mPS2mMS[mPS, mq, asnlmu, mu, muf, nl, loops, fdelm] computes
mMS(mu) from mPS(muf). asnmu = \\alpha_s^(nl)(\\mu).
fdelm (=1) can be omitted.";

 mPS2mSI::usage =
    "mPS2mSI[mPS_, mq_, asfct_, muf_, nl_, loops_, fdelm_] computes
mMS(mMS) from mPS(muf). asfct is the head of the function which computes \\alpha_s^{(nl)}.
fdelm (=1) can be omitted.";

 mOS2mRS::usage =
    "mOS2mRS[mOS, mq, asmu, mu, nuf, nl, loops]  computes the RS
mass mRS(nuf) at the scale mu. {mq} represents a set of light quark masses in 
the on-shell scheme leading to corrections of O(as^2 mq/mOS) (not yet implemented).
asmu = \\alpha_s^(nl)(\\mu) and loops represents the number of loops.";

 mMS2mRS::usage =
    "mMS2mRS[mMS, mq, asmu, mu, nuf, nl, loops, fdelm]
computes mRS(nuf) from mMS(mu). asmu = \\alpha_s^(nl)(mu).
fdelm (=1) can be omitted.";

 mRS2mMS::usage =
    "mRS2mMS[mRS, mq, asnlmu, mu, nuf, nl, loops, fdelm]
computes mMS(mu) from mRS(nuf). asnlmu = \\alpha_s^(nl)(mu).
fdelm (=1) can be omitted.";

 mRS2mSI::usage =
    "mRS2mSI[mRS, mq, asfct, nuf, nl, loops, fdelm]
computes mMS(mMS) from mRS(nuf). asfct is the head of the function which computes \\alpha_s^{(nl)}.
fdelm (=1) can be omitted.";

 mOS2mRSp::usage =
    "As mOS2mRS[] but for RS' mass.";

 mMS2mRSp::usage =
    "As mMS2mRS[] but for RS' mass.";

 mRSp2mMS::usage =
    "As mRS2mMS[] but for RS' mass."; 

 mRSp2mSI::usage =
    "As mRS2mSI[] but for RS' mass.";

 mOS2m1S::usage =
    "mOS2m1S[mOS, mq, asmu, mu, nl, loops] computes the 1S 
mass at the scale mu. {mq} represents a set of light quark masses in 
the on-shell scheme leading to corrections of O(as^2 mq/mOS) (not yet implemented).
asmu = \\alpha_s^(nl)(\\mu) and loops represents the number of loops.";

 mMS2m1S::usage =
    "mMS2m1S[mMS, mq, asmu, mu, nl, loops, fdelm]
computes m1S from mMS(mu). asmu = \\alpha_s^(nl)(mu).
fdelm (=1) can be omitted.";

 m1S2mMS::usage =
    "m1S2mMS[m1S, mq, asnlmu, mu, nl, loops, fdelm]
compues mMS(mu) from m1S. asnlmu = \\alpha_s^(nl)(mu).
fdelm (=1) can be omitted.";

 m1S2mSI::usage =
    "m1S2mSI[m1S, mq, asfct, nl, loops, fdelm]
computes mMS(mMS) from m1S. asfct is the head of the function which computes \\alpha_s^{(nl)}.
fdelm (=1) can be omitted.";

 mKIN2mMS::usage =
    "mKIN2mMS[mbKIN, {mcMS,musmc}, asmus, mus, muf, nlMSOS, nlOSKIN, nloops, CASE, KINverbose]
computes mMS(mus) from mKIN(muf). The routines are developed for bottom and charm quarks.
However, they can also be applied to other qurk flavours.
A theory with nf active flavours is assumed.
nlMSOS and nlOSKIN are the number of massless quarks in the MS-OS and OS-kinetic relation,
respectively. nloops=1,2 or 3 is the number of loops.
CASE defines the scheme which has to be adapted for the lighter quark (charm) mass.
{mcMS,musmc} specifies the lighter quark mass in the MSbar scheme: mcMS(musmc).
{zmcMS,zmusmc} = {0,0}: mass relations fro the light quark are computed; CASE is ignored.
Four pre-defined cases are implemetned, see ??mKIN2mMS.";

 mMS2mKIN::usage =
    "mMS2mKIN[mbMS, {mcMS,musmc}, asmus, mus, muf, nlMSOS, nlOSKIN, nloops, CASE, KINverbose]
computes mKIN(muf) from mMS(mus). See description of mKIN2mMS[] for more details.
Four pre-defined cases are implemetned, see ??mKIN2mMS.";


 DecAsUpOS::usage = 
    "DecAsUpOS[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl+1)(muth)
from the knowledge of als=\\alpha_s^(nl)(muth). massth is the on-shell
mass value of the heavy quark. l specifies the number of loops.
It coincides with the parameter used for the running, i.e. for l=1
one has \\alpha_s^(nl+1)(muth)=\\alpha_s^(nl)(muth), for l=2 the
one-loop decoupling formula is used, ...";

 DecAsDownOS::usage = 
    "DecAsDownOS[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl)(muth)
from the knowledge of als=\\alpha_s^(nl+1)(muth). For the other parameters see
DecAsUpOS[als ,massth ,muth ,nl ,l ]";

 DecAsUpMS::usage = 
    "DecAsUpMS[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl+1)(muth)
from the knowledge of als=\\alpha_s^(nl)(muth). massth is the MS-bar
mass value of the heavy quark. l specifies the number of loops.
It coincides with the parameter used for the running, i.e. for l=1
one has \\alpha_s^(nl+1)(muth)=\\alpha_s^(nl)(muth), for l=2 the
one-loop decoupling formula is used, ...";

 DecAsDownMS::usage = 
    "DecAsDownMS[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl)(muth)
from the knowledge of als=\\alpha_s^(nl+1)(muth). For the other parameters see
DecAsUpMS[als ,massth ,muth ,nl ,l ]";

 DecAsUpSI::usage = 
    "DecAsUpSI[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl+1)(muth)
from the knowledge of als=\\alpha_s^(nl)(muth). massth is the
scale invariant
mass value of the heavy quark. l specifies the number of loops.
It coincides with the parameter used for the running, i.e. for l=1
one has \\alpha_s^(nl+1)(muth)=\\alpha_s^(nl)(muth), for l=2 the
one-loop decoupling formula is used, ...";

 DecAsDownSI::usage = 
    "DecAsDownSI[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl)(muth)
from the knowledge of als=\\alpha_s^(nl+1)(muth). For the other parameters see
DecAsUpSI[als ,massth ,muth ,nl ,l ]";

 DecMqUpOS::usage = 
    "DecMqUpOS[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl+1)(muth)
from the knowledge of mq=m_q^(nl)(muth). als is \\alpha_s^(nl)(muth).
massth is the on-shell mass value
of the heavy quark. l specifies the number of loops. It coincides with the
parameter used for the running, i.e. for l=1 one has 
m_q^(nl+1)(muth)=m_q^(nl)(muth), for l=2 the one-loop 
decoupling formula is used, ...";

 DecMqDownOS::usage = 
    "DecMqDownOS[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl)(muth)
from the knowledge of mq=m_q^(nl+1)(muth). als is \\alpha_s^(nl+1)(muth).
For the other parameters see DecMqUpOS[mq, als ,massth ,muth ,nl ,l ]";

 DecMqUpMS::usage = 
    "DecMqUpMS[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl+1)(muth)
from the knowledge of mq=m_q^(nl)(muth). als is \\alpha_s^(nl)(muth).
massth is the MS-bar mass value
of the heavy quark. l specifies the number of loops. It coincides with the
parameter used for the running, i.e. for l=1 one has 
m_q^(nl+1)(muth)=m_q^(nl)(muth), for l=2 the one-loop 
decoupling formula is used, ...";

 DecMqDownMS::usage = 
    "DecMqDownMS[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl)(muth)
from the knowledge of mq=m_q^(nl+1)(muth). als is \\alpha_s^(nl+1)(muth).
For the other parameters see DecMqUpMS[mq, als ,massth ,muth ,nl ,l ]";

 DecMqUpSI::usage = 
    "DecMqUpSI[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl+1)(muth)
from the knowledge of mq=m_q^(nl)(muth). als is \\alpha_s^(nl)(muth).
massth is the scale invariant mass value
of the heavy quark. l specifies the number of loops. It coincides with the
parameter used for the running, i.e. for l=1 one has 
m_q^(nl+1)(muth)=m_q^(nl)(muth), for l=2 the one-loop 
decoupling formula is used, ...";

 DecMqDownSI::usage = 
    "DecMqDownSI[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl)(muth)
from the knowledge of mq=m_q^(nl+1)(muth). als is \\alpha_s^(nl+1)(muth).
For the other parameters see DecMqUpSI[mq, als ,massth ,muth ,nl ,l ]";

 DecLambdaUp::usage = 
    "DecLambdaUp[lam, massth, nl ,l ] computes \\Lambda^(nl+1) 
from the knowledge of lam=\\Lambda^(nl).
massth is the scale invariant mass value
of the heavy quark.
l specifies the number of loops. It coincides with the
parameter used for the running.";

 DecLambdaDown::usage = 
    "DecLambdaDown[lam, massth, nl ,l ] computes \\Lambda^(nl) 
from the knowledge of lam=\\Lambda^(nl+1).
massth is the scale invariant mass value
of the heavy quark.
l specifies the number of loops. It coincides with the
parameter used for the running.";

 AlL2AlH::usage = 
    "AlL2AlH[als, mu1, decpar, mu2, nloops] computes alphas(mu2) with h active
flavours from the knowledge of als=alphas(mu1) with l active flavours.
h > l must be fulfilled. decpar is a set which 
may contain several triples  indicating the
number of flavours, the heavy (on-shell)
quark mass and the scale at which the
decoupling is performed. The running is performed with
AlphasExact[..] and the decoupling with DecAsUpOS[..].
    ";

 AlH2AlL::usage = 
    "AlH2AlL[als, mu1, decpar, mu2, nloops] computes alphas(mu2) with l active
flavours from the knowledge of als=alphas(mu1) with h active flavours.
h > l must be fulfilled. decpar is a set which 
may contain several triples  indicating the
number of flavours, the heavy (on-shell)
quark mass and the scale at which the
decoupling is performed. The running is performed with
AlphasExact[..] and the decoupling with DecAsDownOS[..].
    ";

 mL2mH::usage = 
    "mL2mH[mql, asl, mu1, decpar, mu2, l] computes mq(mu2) with h active
flavours from the knowledge of mql=mq(mu1) with l active flavours.
asl=\\alpha_s(mu1) with l active flavours.
h > l must be fulfilled. decpar is a set which 
may contain several triples  indicating the
number of flavours, the heavy (on-shell)
quark mass and the scale at which the
decoupling is performed. The running of the coupling constant and 
the quark mass is performed with
AlphasExact[..] and mMS2mMS[..], respectively, and 
the decoupling with DecMqUpOS[..] and DecAsUpOS[..].
    ";

 mH2mL::usage = 
    "mH2L[mqh ,ash ,mu1 , decpar, mu2, l] computes mq(mu2) with l active
flavours from the knowledge of mqh=mq(mu1) with h active flavours.
ash=\\alpha_s(mu1) with h active flavours.
h > l must be fulfilled. decpar is a set which 
may contain several triples  indicating the
number of flavours, the heavy (on-shell)
quark mass and the scale at which the
decoupling is performed. The running of the coupling constant and 
the quark mass is performed with
AlphasExact[..] and mMS2mMS[..], respectively, and 
the decoupling with DecMqDownOS[..] and DecAsDownOS[..].
    ";

 AsRunDec::usage = 
    "AsRunDec[als ,mu0 ,mu , l] computes \\alpha_s(mu) from the knowledge of
als=\\alpha_s(mu0) The number of active flavours is determined automatically
from the choices of mu and mu0, respectively. The decoupling is performed at 
the pole mass of the respective heavy quark.";

 AsmMSrunexact::usage = 
    "AsmMSrunexact[mmu0, asmu0, mu0, mu, nf, l] solves simultanesously
the renormlaization group functions for \\alpha_s and mq.
Computes \\alpha_s(mu) and mq(mu) using asmu0=\\alpha_s(mu0) and 
mmu0 = mq(mu0) as starting values.";

 Mc5Mzfrommuc4::usage = 
    "Mc5Mzfrommuc4[asMz_,muc4_,Mb_,mub_,Mz_,loops_] computes mc^(5)(Mz)
using \\alpha_s^(5), mc(mc), mb^OS, mu^b_dec, M_Z and the loop order as input.";

 AsMbfromAsMz::usage = 
    "AsMbfromAsMz[asMz_,Mb_,loops_] or AsMbfromAsMz[asMz_,loops_] compare
running of \\alpha_s from mu=M_Z to mu=mb^OS for 5 active flavours using
AlphasExact[] and AlphasLam[].";

 setbeta::usage = 
    "Coefficients of the QCD beta function."

 setgamma::usage = 
    "Coefficients of the quark mass anomalous dimension in QCD."

 setdec::usage = 
    "Decoupling constants."

(* ************************************************************ *)

Begin["`Modules`"];

cut[x_, n_Integer] := (x^p_Integer /; p > n -> 0);

(*
 argument used in N[] 
*)

$NumPrec=20;
numprec=$NumPrec;

(*
 minimal allowed ratio of \mu/\Lambda for which no warning is printed
*)

rmulam = 1.5;

(*
 numerical expression of some symbols ... 
*)

num1 = {
    z2 -> Zeta[2],
    z3 -> Zeta[3],
    z4 -> Zeta[4],
    z5 -> Zeta[5],
    log2 -> Log[2],
    B4->16*PolyLog[4,1/2]+2/3*Log[2]^4-2/3*Pi^2*Log[2]^2-13/180*Pi^4
};

(*
 numerical expression for colour factors (QCD)
*)

cf2num = {
    cf -> 4/3,
    ca -> 3,
    tr -> 1/2
};

(*
 coefficients of \beta function 
*)

setbeta = {
    b0 -> 11/4 - nf/6, 
    b1 -> 51/8 - (19*nf)/24,
    b2 -> 2857/128 - (5033*nf)/1152 + (325*nf^2)/3456,
    b3 -> 149753/1536 - (1078361*nf)/41472 + (50065*nf^2)/41472 +
	(1093*nf^3)/186624 + (891*Zeta[3])/64 - (1627*nf*Zeta[3])/1728 +
	    (809*nf^2*Zeta[3])/2592,
    b4 -> 8157455/16384 - (336460813*nf)/1990656 + (25960913*nf^2)/1990656 - 
 (630559*nf^3)/5971968 + (1205*nf^4)/2985984 + (621885*z3)/2048 - 
 (1202791*nf*z3)/20736 + (698531*nf^2*z3)/82944 - (24361*nf^3*z3)/124416 - 
 (19*nf^4*z3)/10368 - (88209*z4)/2048 + (33935*nf*z4)/6144 - 
 (5263*nf^2*z4)/4608 + (809*nf^3*z4)/13824 - (144045*z5)/512 + 
 (1358995*nf*z5)/27648 - (5965*nf^2*z5)/1296 + (115*nf^3*z5)/2304
	    } /. num1;

(*
 coefficients of \gamma_m function 
*)

setgamma = {
    g0 -> 1,
    g1 -> (101/2 - (5*nf)/3)/12,
    g2 -> (3747/4 - (nf*(1108/9 + (70*nf)/27 + 80*Zeta[3]))/2)/48,
    g3 -> 4603055/41472 - (91723*nf)/6912 + (2621*nf^2)/31104 - 
	(83*nf^3)/15552 +
     (11*nf*Pi^4)/288 - (nf^2*Pi^4)/432 + (530*Zeta[3])/27 -
	 (2137*nf*Zeta[3])/144 + (25*nf^2*Zeta[3])/72 + (nf^3*Zeta[3])/108 -
	     (275*Zeta[5])/8 + (575*nf*Zeta[5])/72,
    g4 ->
          (1/4)^5*(99512327/162 + 46402466/243*Zeta[3] + 96800*Zeta[3]^2-698126/9*Zeta[4]
            -231757160/243*Zeta[5] + 242000*Zeta[6] + 412720*Zeta[7]
          +nf*(-150736283/1458 - 12538016/81*Zeta[3] - 75680/9*Zeta[3]^2 + 2038742/27*Zeta[4] 
               +49876180/243*Zeta[5] - 638000/9*Zeta[6] - 1820000/27*Zeta[7])
          +nf^2*(1320742/729 + 2010824/243*Zeta[3] + 46400/27*Zeta[3]^2 - 166300/27*Zeta[4] -264040/81*Zeta[5] + 92000/27*Zeta[6])
          +nf^3*(91865/1458 + 12848/81*Zeta[3] +448/9*Zeta[4] - 5120/27*Zeta[5]) 
          +nf^4*(-260/243 - 320/243*Zeta[3] + 64/27*Zeta[4]))    
    };

(*
 \alpha_s as function of \mu and \Lambda expanded in 1/\Lambda\beta_0
*)

setasL = {
    as0 -> 1/L/b0,
    as1 -> 1/L/b0 - 1/(b0*L)^2*b1/b0*Log[L],
    as2 -> ( 1/L/b0 - 1/(b0*L)^2*b1/b0*Log[L] +
	    1/(b0*L)^3*(b1^2/b0^2*(Log[L]^2-Log[L]-1)+b2/b0)
	    ),
    as3 -> (
	    1/(b0*L) - (b1*Log[L])/(b0^3*L^2) + 
	    (b2/b0^4 + (b1^2*(-1 - Log[L] + Log[L]^2))/b0^5)/L^3 + 
	    (b3/(2*b0^5) - (3*b1*b2*Log[L])/b0^6 + 
	     (b1^3*(-1/2 + 2*Log[L] + (5*Log[L]^2)/2 - Log[L]^3))/b0^7)/L^4
	    ),
    as4 -> (
	    1/(b0*L) - (b1*Log[L])/(b0^3*L^2) + 
	    (b2/b0^4 + (b1^2*(-1 - Log[L] + Log[L]^2))/b0^5)/L^3 + 
	    (b3/(2*b0^5) - (3*b1*b2*Log[L])/b0^6 + 
	     (b1^3*(-1/2 + 2*Log[L] + (5*Log[L]^2)/2 - Log[L]^3))/b0^7)/L^4
	    +(b4/3/b0^6 + b2^2*5/3/b0^7 - b1*b3*(2*Log[L]+1/6)/b0^7
	      + 3*b1^2*b2*(2*Log[L]^2-Log[L]-1)/b0^8
	      + b1^4*(Log[L]^4-13/3*Log[L]^3-3/2*Log[L]^2+4*Log[L]+7/6)/b0^9 
	      )/L^5
	    )
    };


(*
 \zeta_g, 1/\zeta_g, ...
 as6to5:  api^(nf) = api^(nf-1)*(1 + O(api^(nf-1)) + ...)
*)

sa4 = PolyLog[4,1/2];
sa5 = PolyLog[5,1/2];

setdec = {

    as6to5os ->
  1 + (7*api^2)/24 + (58933*api^3)/124416 + (api*lmm)/6 + (19*api^2*lmm)/24 +
   (8941*api^3*lmm)/1728 + (api^2*lmm^2)/36 + (511*api^3*lmm^2)/576 +
   (api^3*lmm^3)/216 - (2479*api^3*nl)/31104 - (409*api^3*lmm*nl)/1728 +
   (2*api^3*z2)/3 - (api^3*nl*z2)/9 + (2*api^3*z2*Log[2])/9 +
   (80507*api^3*Zeta[3])/27648
   +      api^4*(2313952/382725 + (47039*lmm^2)/3456 + (14149*lmm^3)/10368 + 
       lmm^4/1296 + (644201*Pi^2)/116640 + (71102219*Pi^4)/195955200 + 
       (49*z2)/18 + (49*log2*z2)/54 - (49*z3)/216 - (587*Pi^2*Log[2])/486 - 
       (9318467*Pi^4*Log[2])/32659200 + (2913037*Pi^2*Log[2]^2)/1306368 - 
       (340853*Pi^2*Log[2]^3)/816480 - (3179149*Log[2]^4)/1306368 + 
       (340853*Log[2]^5)/1360800 - (3179149*PolyLog[4, 1/2])/54432 - 
       (340853*PolyLog[5, 1/2])/11340 + nl^2*(140825/1492992 + 
         (493*lmm^2)/20736 + (13*Pi^2)/972 + lmm*(1679/186624 + z2/27) + 
         (19*Zeta[3])/1728) + (2428169183*Zeta[3])/87091200 - 
       (1439*Pi^2*Zeta[3])/1296 + lmm*(21084715/746496 + (35*z2)/9 + 
         (35*log2*z2)/27 - (65*z3)/108 + (3022001*Zeta[3])/165888) + 
       nl*(-1773073/746496 - (9115*lmm^2)/10368 - (107*lmm^3)/1728 - 
         (967*Pi^2)/1944 + (697709*Pi^4)/14929920 - (49*z2)/108 - 
         (11*Pi^2*Log[2])/243 + (1709*Pi^2*Log[2]^2)/124416 - 
         (173*Log[2]^4)/124416 - (173*PolyLog[4, 1/2])/5184 + 
         lmm*(-1140191/373248 - (47*z2)/54 - (2*log2*z2)/27 - (7*z3)/27 - 
           (110779*Zeta[3])/82944) - (4756441*Zeta[3])/995328 - 
         (115*Zeta[5])/576) + (40596749*Zeta[5])/1451520)
,

    as6to5ms ->
  1 + (api*lmm)/6 + api^2*(-11/72 + (11*lmm)/24 + lmm^2/36) + 
   api^3*(-564731/124416 + (2645*lmm)/1728 + (167*lmm^2)/576 + lmm^3/216 + 
   (2633/31104 - (67*lmm)/576 + lmm^2/36)*nl + (82043*Zeta[3])/27648)
   +  api^4*(-1165152397/24494400 + (1837*lmm^2)/1152 + (2909*lmm^3)/10368 + 
       lmm^4/1296 + (76940219*Pi^4)/195955200 - (9318467*Pi^4*Log[2])/
        32659200 + (3031309*Pi^2*Log[2]^2)/1306368 - (340853*Pi^2*Log[2]^3)/
        816480 - (3031309*Log[2]^4)/1306368 + (340853*Log[2]^5)/1360800 - 
       (3031309*PolyLog[4, 1/2])/54432 - (340853*PolyLog[5, 1/2])/11340 + 
       nl^2*(271883/4478976 - (6865*lmm)/186624 + (77*lmm^2)/20736 - 
         lmm^3/324 - (167*Zeta[3])/5184) + (2362581983*Zeta[3])/87091200 + 
       lmm*(-11093717/746496 + (3022001*Zeta[3])/165888) + 
       nl*(4770941/2239488 + (277*lmm^2)/10368 + (271*lmm^3)/5184 + 
         (541549*Pi^4)/14929920 + (685*Pi^2*Log[2]^2)/124416 - 
         (685*Log[2]^4)/124416 - (685*PolyLog[4, 1/2])/5184 + 
         lmm*(141937/373248 - (110779*Zeta[3])/82944) - 
         (3645913*Zeta[3])/995328 - (115*Zeta[5])/576) + 
       (12057583*Zeta[5])/483840)
,

    as6to5si ->
  1 + (api*lmmu)/6 + api^2*(-11/72 + (19*lmmu)/24 + lmmu^2/36) + 
   api^3*(-564731/124416 + (2191*lmmu)/576 + (511*lmmu^2)/576 + lmmu^3/216 + 
   (2633/31104 - (281*lmmu)/1728)*nl + (82043*Zeta[3])/27648)
   + api^4*(-1165152397/24494400 - 
       (1531493*lmmu)/746496 + (33887*lmmu^2)/3456 + (14149*lmmu^3)/10368 + 
       lmmu^4/1296 + (2*Pi^2)/9 + (76940219*Pi^4)/195955200 - (4*z2)/3 - 
       (4*log2*z2)/9 + z3/9 - (5*lmmu*z3)/18 + (2*Pi^2*Log[2])/27 - 
       (9318467*Pi^4*Log[2])/32659200 + (3031309*Pi^2*Log[2]^2)/1306368 - 
       (340853*Pi^2*Log[2]^3)/816480 - (3031309*Log[2]^4)/1306368 + 
       (340853*Log[2]^5)/1360800 - (3031309*PolyLog[4, 1/2])/54432 - 
       (340853*PolyLog[5, 1/2])/11340 + nl^2*(271883/4478976 - 
         (8545*lmmu)/186624 + (79*lmmu^2)/6912 - (167*Zeta[3])/5184) + 
       (2352905183*Zeta[3])/87091200 + (3022001*lmmu*Zeta[3])/165888 + 
       nl*(4770941/2239488 - (158687*lmmu)/373248 - (515*lmmu^2)/1152 - 
         (107*lmmu^3)/1728 - Pi^2/27 + (541549*Pi^4)/14929920 + (2*z2)/9 - 
         (5*lmmu*z3)/18 + (685*Pi^2*Log[2]^2)/124416 - 
         (685*Log[2]^4)/124416 - (685*PolyLog[4, 1/2])/5184 - 
         (3645913*Zeta[3])/995328 - (110779*lmmu*Zeta[3])/82944 - 
         (115*Zeta[5])/576) + (12057583*Zeta[5])/483840)
 ,

    as5to6os ->
  1 - (api*lmm)/6 + api^2*(-7/24 - (19*lmm)/24 + lmm^2/36) +
   api^3*(-58933/124416 - (8521*lmm)/1728 - (131*lmm^2)/576 - lmm^3/216 +
      nl*(2479/31104 + (409*lmm)/1728 + z2/9) - (2*z2)/3 - (2*z2*Log[2])/9 -
      (80507*Zeta[3])/27648)
+ api^4*(-141841753/24494400 - 
       (7693*lmm^2)/1152 - (8371*lmm^3)/10368 + lmm^4/1296 - 
       (644201*Pi^2)/116640 - (71102219*Pi^4)/195955200 - (49*z2)/18 - 
       (49*log2*z2)/54 + (49*z3)/216 + (587*Pi^2*Log[2])/486 + 
       (9318467*Pi^4*Log[2])/32659200 - (2913037*Pi^2*Log[2]^2)/1306368 + 
       (340853*Pi^2*Log[2]^3)/816480 + (3179149*Log[2]^4)/1306368 - 
       (340853*Log[2]^5)/1360800 + (3179149*PolyLog[4, 1/2])/54432 + 
       (340853*PolyLog[5, 1/2])/11340 + lmm*(-19696909/746496 - (29*z2)/9 - 
         (29*log2*z2)/27 + (59*z3)/108 - (2529743*Zeta[3])/165888) + 
       nl^2*(-140825/1492992 - (493*lmm^2)/20736 - (13*Pi^2)/972 + 
         lmm*(-1679/186624 - z2/27) - (19*Zeta[3])/1728) - 
       (2428169183*Zeta[3])/87091200 + (1439*Pi^2*Zeta[3])/1296 + 
       nl*(1773073/746496 + (6661*lmm^2)/10368 + (107*lmm^3)/1728 + 
         (967*Pi^2)/1944 - (697709*Pi^4)/14929920 + (49*z2)/108 + 
         (11*Pi^2*Log[2])/243 - (1709*Pi^2*Log[2]^2)/124416 + 
         (173*Log[2]^4)/124416 + (173*PolyLog[4, 1/2])/5184 + 
         (4756441*Zeta[3])/995328 + lmm*(1110443/373248 + (41*z2)/54 + 
           (2*log2*z2)/27 + (7*z3)/27 + (110779*Zeta[3])/82944) + 
         (115*Zeta[5])/576) - (40596749*Zeta[5])/1451520)
,

    as5to6ms ->
  1 - (api*lmm)/6 + api^2*(11/72 - (11*lmm)/24 + lmm^2/36) +
   api^3*(564731/124416 - (955*lmm)/576 + (53*lmm^2)/576 - lmm^3/216 +
	(-2633/31104 + (67*lmm)/576 - lmm^2/36)*nl - (82043*Zeta[3])/27648)
    + api^4*(291716893/6123600 + (2177*lmm^2)/3456 - (1883*lmm^3)/10368 + 
       lmm^4/1296 - (76940219*Pi^4)/195955200 + (9318467*Pi^4*Log[2])/
        32659200 - (3031309*Pi^2*Log[2]^2)/1306368 + (340853*Pi^2*Log[2]^3)/
        816480 + (3031309*Log[2]^4)/1306368 - (340853*Log[2]^5)/1360800 + 
       (3031309*PolyLog[4, 1/2])/54432 + (340853*PolyLog[5, 1/2])/11340 + 
       lmm*(7391699/746496 - (2529743*Zeta[3])/165888) + 
       nl^2*(-271883/4478976 + (6865*lmm)/186624 - (77*lmm^2)/20736 + 
         lmm^3/324 + (167*Zeta[3])/5184) - (2362581983*Zeta[3])/87091200 + 
       nl*(-4770941/2239488 - (1483*lmm^2)/10368 - (127*lmm^3)/5184 - 
         (541549*Pi^4)/14929920 - (685*Pi^2*Log[2]^2)/124416 + 
         (685*Log[2]^4)/124416 + (685*PolyLog[4, 1/2])/5184 + 
         (3645913*Zeta[3])/995328 + lmm*(-110341/373248 + (110779*Zeta[3])/
            82944) + (115*Zeta[5])/576) - (12057583*Zeta[5])/483840)
,

    as5to6si ->
  1 - (api*lmmu)/6 + api^2*(11/72 - (19*lmmu)/24 + lmmu^2/36) + 
   api^3*(564731/124416 - (6793*lmmu)/1728 - (131*lmmu^2)/576 - lmmu^3/216 + 
   (-2633/31104 + (281*lmmu)/1728)*nl - (82043*Zeta[3])/27648)
  + api^4*(291716893/6123600 - 
       (2398621*lmmu)/746496 - (14023*lmmu^2)/3456 - (8371*lmmu^3)/10368 + 
       lmmu^4/1296 - (2*Pi^2)/9 - (76940219*Pi^4)/195955200 + (4*z2)/3 + 
       (4*log2*z2)/9 - z3/9 + (5*lmmu*z3)/18 - (2*Pi^2*Log[2])/27 + 
       (9318467*Pi^4*Log[2])/32659200 - (3031309*Pi^2*Log[2]^2)/1306368 + 
       (340853*Pi^2*Log[2]^3)/816480 + (3031309*Log[2]^4)/1306368 - 
       (340853*Log[2]^5)/1360800 + (3031309*PolyLog[4, 1/2])/54432 + 
       (340853*PolyLog[5, 1/2])/11340 + nl^2*(-271883/4478976 + 
         (8545*lmmu)/186624 - (79*lmmu^2)/6912 + (167*Zeta[3])/5184) - 
       (2352905183*Zeta[3])/87091200 - (2529743*lmmu*Zeta[3])/165888 + 
       nl*(-4770941/2239488 + (190283*lmmu)/373248 + (983*lmmu^2)/3456 + 
         (107*lmmu^3)/1728 + Pi^2/27 - (541549*Pi^4)/14929920 - (2*z2)/9 + 
         (5*lmmu*z3)/18 - (685*Pi^2*Log[2]^2)/124416 + 
         (685*Log[2]^4)/124416 + (685*PolyLog[4, 1/2])/5184 + 
         (3645913*Zeta[3])/995328 + (110779*lmmu*Zeta[3])/82944 + 
         (115*Zeta[5])/576) - (12057583*Zeta[5])/483840)
,

    mq5to6os ->
  1 + api^2*(89/432 - (5*lmm)/36 + lmm^2/12) +
   api^3*(1871/2916 - B4/36 + (121*lmm)/2592 + (319*lmm^2)/432 +
      (29*lmm^3)/216 + nl*(1327/11664 - (53*lmm)/432 - lmm^3/108 -
         (2*z3)/27) - (407*z3)/864 - (5*lmm*z3)/6 + (5*z4)/4)

     +api^4*(122466970229/3292047360 + (1403*lmm^3)/648 + (305*lmm^4)/1152 + 
       (16187201*Pi^4)/52254720 - (787*Pi^6)/81648 - (5*z2)/9 - 
       (5*log2*z2)/27 + (5*z3)/108 - (145*Pi^4*Log[2])/1944 + 
       (1924649*Pi^2*Log[2]^2)/4354560 - (59*Pi^2*Log[2]^3)/972 - 
       (1924649*Log[2]^4)/4354560 + (59*Log[2]^5)/1620 + (4*Log[2]^6)/81 - 
       (Log[2]^5*Log[16])/81 - (1924649*PolyLog[4, 1/2])/181440 - 
       (118*PolyLog[5, 1/2])/27 + lmm^2*(104803/10368 - (155*Zeta[3])/48) + 
       nl^2*(17671/124416 + (31*lmm^2)/1296 + lmm^4/864 - (7*Pi^4)/8640 + 
         lmm*(-3401/46656 + (7*Zeta[3])/108) - (5*Zeta[3])/864) - 
       (443509931*Zeta[3])/40642560 + (1061*Zeta[3]^2)/576 - 
       (59015*Zeta[5])/1728 + nl*(-2403419/746496 - (9601*lmm^2)/10368 - 
         (47*lmm^3)/288 - (5*lmm^4)/144 - (245*Pi^4)/62208 + (5*z2)/54 + 
         (49*Pi^4*Log[2])/6480 - (49*Pi^2*Log[2]^2)/2592 + 
         (Pi^2*Log[2]^3)/162 + (49*Log[2]^4)/2592 - Log[2]^5/270 + 
         (49*PolyLog[4, 1/2])/108 + (4*PolyLog[5, 1/2])/9 + 
         lmm*(7045/31104 - (163*Pi^4)/12960 - z2/9 - (Pi^2*Log[2]^2)/108 + 
           Log[2]^4/108 + (2*PolyLog[4, 1/2])/9 + (221*Zeta[3])/576) - 
         (1075*Zeta[3])/1728 + (497*Zeta[5])/288) + 
       lmm*(-1126487/373248 + (4123*Pi^4)/25920 + (2*z2)/3 + (2*log2*z2)/9 - 
         z3/18 + (31*Pi^2*Log[2]^2)/216 - (31*Log[2]^4)/216 - 
         (31*PolyLog[4, 1/2])/9 - (419341*Zeta[3])/27648 + (575*Zeta[5])/72))
,
    mq5to6ms ->
  1 + api^2*(89/432 - (5*lmm)/36 + lmm^2/12) + 
   api^3*(2951/2916 - B4/36 + (175*lmm^2)/432 + (29*lmm^3)/216 + 
   lmm*(-311/2592 - (5*z3)/6) + nl*(1327/11664 - (53*lmm)/432 - lmm^3/108 - 
     (2*z3)/27) - (407*z3)/864 + (5*z4)/4)
   +  api^4*(131968227029/3292047360 + (301*lmm^3)/324 + (305*lmm^4)/1152 + 
       (16187201*Pi^4)/52254720 - (787*Pi^6)/81648 - (145*Pi^4*Log[2])/1944 + 
       (1924649*Pi^2*Log[2]^2)/4354560 - (59*Pi^2*Log[2]^3)/972 - 
       (1924649*Log[2]^4)/4354560 + (59*Log[2]^5)/1620 + (4*Log[2]^6)/81 - 
       (Log[2]^5*Log[16])/81 - (1924649*PolyLog[4, 1/2])/181440 - 
       (118*PolyLog[5, 1/2])/27 + lmm^2*(51163/10368 - (155*Zeta[3])/48) + 
       nl^2*(17671/124416 + (31*lmm^2)/1296 + lmm^4/864 - (7*Pi^4)/8640 + 
         lmm*(-3401/46656 + (7*Zeta[3])/108) - (5*Zeta[3])/864) - 
       (353193131*Zeta[3])/40642560 + (1061*Zeta[3]^2)/576 - 
       (59015*Zeta[5])/1728 + nl*(-2261435/746496 - (7825*lmm^2)/10368 - 
         (23*lmm^3)/288 - (5*lmm^4)/144 - (245*Pi^4)/62208 + 
         (49*Pi^4*Log[2])/6480 - (49*Pi^2*Log[2]^2)/2592 + 
         (Pi^2*Log[2]^3)/162 + (49*Log[2]^4)/2592 - Log[2]^5/270 + 
         (49*PolyLog[4, 1/2])/108 + (4*PolyLog[5, 1/2])/9 + 
         lmm*(16669/31104 - (163*Pi^4)/12960 - (Pi^2*Log[2]^2)/108 + 
           Log[2]^4/108 + (2*PolyLog[4, 1/2])/9 + (221*Zeta[3])/576) - 
         (1075*Zeta[3])/1728 + (497*Zeta[5])/288) + 
       lmm*(-2810855/373248 + (4123*Pi^4)/25920 + (31*Pi^2*Log[2]^2)/216 - 
         (31*Log[2]^4)/216 - (31*PolyLog[4, 1/2])/9 - (373261*Zeta[3])/
          27648 + (575*Zeta[5])/72))

  ,

    mq5to6si ->
  1 + api^2*(89/432 - (5*lmmu)/36 + lmmu^2/12) + 
   api^3*(2951/2916 - B4/36 - (1031*lmmu)/2592 + (319*lmmu^2)/432 + 
   (29*lmmu^3)/216 + nl*(1327/11664 - (53*lmmu)/432 - lmmu^3/108 - 
     (2*z3)/27) - (407*z3)/864 - (5*lmmu*z3)/6 + (5*z4)/4)
   +      api^4*(131968227029/3292047360 - (3322343*lmmu)/373248 + 
       (81763*lmmu^2)/10368 + (1403*lmmu^3)/648 + (305*lmmu^4)/1152 + 
       (16187201*Pi^4)/52254720 + (4123*lmmu*Pi^4)/25920 - (787*Pi^6)/81648 - 
       (145*Pi^4*Log[2])/1944 + (1924649*Pi^2*Log[2]^2)/4354560 + 
       (31*lmmu*Pi^2*Log[2]^2)/216 - (59*Pi^2*Log[2]^3)/972 - 
       (1924649*Log[2]^4)/4354560 - (31*lmmu*Log[2]^4)/216 + 
       (59*Log[2]^5)/1620 + (4*Log[2]^6)/81 - (Log[2]^5*Log[16])/81 - 
       (1924649*PolyLog[4, 1/2])/181440 - (31*lmmu*PolyLog[4, 1/2])/9 - 
       (118*PolyLog[5, 1/2])/27 - (353193131*Zeta[3])/40642560 - 
       (419341*lmmu*Zeta[3])/27648 - (155*lmmu^2*Zeta[3])/48 + 
       (1061*Zeta[3]^2)/576 + nl^2*(17671/124416 - (3401*lmmu)/46656 + 
         (31*lmmu^2)/1296 + lmmu^4/864 - (7*Pi^4)/8640 - (5*Zeta[3])/864 + 
         (7*lmmu*Zeta[3])/108) - (59015*Zeta[5])/1728 + 
       (575*lmmu*Zeta[5])/72 + nl*(-2261435/746496 + (10237*lmmu)/31104 - 
         (8065*lmmu^2)/10368 - (47*lmmu^3)/288 - (5*lmmu^4)/144 - 
         (245*Pi^4)/62208 - (163*lmmu*Pi^4)/12960 + (49*Pi^4*Log[2])/6480 - 
         (49*Pi^2*Log[2]^2)/2592 - (lmmu*Pi^2*Log[2]^2)/108 + 
         (Pi^2*Log[2]^3)/162 + (49*Log[2]^4)/2592 + (lmmu*Log[2]^4)/108 - 
         Log[2]^5/270 + (49*PolyLog[4, 1/2])/108 + (2*lmmu*PolyLog[4, 1/2])/
          9 + (4*PolyLog[5, 1/2])/9 - (1075*Zeta[3])/1728 + 
         (221*lmmu*Zeta[3])/576 + (497*Zeta[5])/288))

,

    mq6to5os ->
  1 + api^2*(-89/432 + (5*lmm)/36 - lmm^2/12) +
   api^3*(-1871/2916 + B4/36 - (299*lmm)/2592 - (299*lmm^2)/432 -
      (35*lmm^3)/216 - (1327*nl)/11664 + (53*lmm*nl)/432 + (lmm^3*nl)/108 +
      (407*z3)/864 + (5*lmm*z3)/6 + (2*nl*z3)/27 - (5*z4)/4)
     +      api^4*(-122722873709/3292047360 - (6941*lmm^3)/2592 - 
       (1147*lmm^4)/3456 - (16187201*Pi^4)/52254720 + (787*Pi^6)/81648 + 
       (5*z2)/9 + (5*log2*z2)/27 - (5*z3)/108 + (145*Pi^4*Log[2])/1944 - 
       (1924649*Pi^2*Log[2]^2)/4354560 + (59*Pi^2*Log[2]^3)/972 + 
       (1924649*Log[2]^4)/4354560 - (59*Log[2]^5)/1620 - (4*Log[2]^6)/81 + 
       (Log[2]^5*Log[16])/81 + (1924649*PolyLog[4, 1/2])/181440 + 
       (118*PolyLog[5, 1/2])/27 + nl^2*(-17671/124416 - (31*lmm^2)/1296 - 
         lmm^4/864 + (7*Pi^4)/8640 + lmm*(3401/46656 - (7*Zeta[3])/108) + 
         (5*Zeta[3])/864) + (443509931*Zeta[3])/40642560 - 
       (1061*Zeta[3]^2)/576 + lmm^2*(-34297/3456 + (175*Zeta[3])/48) + 
       lmm*(99319/41472 - (481*Pi^4)/2880 - (2*z2)/3 - (2*log2*z2)/9 + 
         z3/18 - (11*Pi^2*Log[2]^2)/72 + (11*Log[2]^4)/72 + 
         (11*PolyLog[4, 1/2])/3 + (47317*Zeta[3])/3072 - (575*Zeta[5])/72) + 
       nl*(2403419/746496 + (10237*lmm^2)/10368 + (47*lmm^3)/288 + 
         (17*lmm^4)/432 + (245*Pi^4)/62208 - (5*z2)/54 - 
         (49*Pi^4*Log[2])/6480 + (49*Pi^2*Log[2]^2)/2592 - 
         (Pi^2*Log[2]^3)/162 - (49*Log[2]^4)/2592 + Log[2]^5/270 - 
         (49*PolyLog[4, 1/2])/108 - (4*PolyLog[5, 1/2])/9 + 
         lmm*(-26443/93312 + (163*Pi^4)/12960 + z2/9 + (Pi^2*Log[2]^2)/108 - 
           Log[2]^4/108 - (2*PolyLog[4, 1/2])/9 - (599*Zeta[3])/1728) + 
         (1075*Zeta[3])/1728 - (497*Zeta[5])/288) + (59015*Zeta[5])/1728)

,

    mq6to5ms ->
  1 + api^2*(-89/432 + (5*lmm)/36 - lmm^2/12) + 
   api^3*(-2951/2916 + B4/36 - (155*lmm^2)/432 - (35*lmm^3)/216 + 
   nl*(-1327/11664 + (53*lmm)/432 + lmm^3/108 + (2*z3)/27) + 
   lmm*(133/2592 + (5*z3)/6) + (407*z3)/864 - (5*z4)/4)
  +     api^4*(-131621265869/3292047360 - (3161*lmm^3)/2592 - 
       (1147*lmm^4)/3456 - (16187201*Pi^4)/52254720 + (787*Pi^6)/81648 + 
       (145*Pi^4*Log[2])/1944 - (1924649*Pi^2*Log[2]^2)/4354560 + 
       (59*Pi^2*Log[2]^3)/972 + (1924649*Log[2]^4)/4354560 - 
       (59*Log[2]^5)/1620 - (4*Log[2]^6)/81 + (Log[2]^5*Log[16])/81 + 
       (1924649*PolyLog[4, 1/2])/181440 + (118*PolyLog[5, 1/2])/27 + 
       nl^2*(-17671/124416 - (31*lmm^2)/1296 - lmm^4/864 + (7*Pi^4)/8640 + 
         lmm*(3401/46656 - (7*Zeta[3])/108) + (5*Zeta[3])/864) + 
       (353193131*Zeta[3])/40642560 - (1061*Zeta[3]^2)/576 + 
       lmm^2*(-16193/3456 + (175*Zeta[3])/48) + 
       lmm*(279367/41472 - (481*Pi^4)/2880 - (11*Pi^2*Log[2]^2)/72 + 
         (11*Log[2]^4)/72 + (11*PolyLog[4, 1/2])/3 + (42197*Zeta[3])/3072 - 
         (575*Zeta[5])/72) + nl*(2261435/746496 + (8461*lmm^2)/10368 + 
         (23*lmm^3)/288 + (17*lmm^4)/432 + (245*Pi^4)/62208 - 
         (49*Pi^4*Log[2])/6480 + (49*Pi^2*Log[2]^2)/2592 - 
         (Pi^2*Log[2]^3)/162 - (49*Log[2]^4)/2592 + Log[2]^5/270 - 
         (49*PolyLog[4, 1/2])/108 - (4*PolyLog[5, 1/2])/9 + 
         lmm*(-55315/93312 + (163*Pi^4)/12960 + (Pi^2*Log[2]^2)/108 - 
           Log[2]^4/108 - (2*PolyLog[4, 1/2])/9 - (599*Zeta[3])/1728) + 
         (1075*Zeta[3])/1728 - (497*Zeta[5])/288) + (59015*Zeta[5])/1728)
,

    mq6to5si ->
  1 + api^2*(-89/432 + (5*lmmu)/36 - lmmu^2/12) + 
   api^3*(-2951/2916 + B4/36 + (853*lmmu)/2592 - (299*lmmu^2)/432 - 
   (35*lmmu^3)/216 + nl*(-1327/11664 + (53*lmmu)/432 + lmmu^3/108 + 
     (2*z3)/27) + (407*z3)/864 + (5*lmmu*z3)/6 - (5*z4)/4)
   +     api^4*(-131621265869/3292047360 + (330503*lmmu)/41472 - 
       (8531*lmmu^2)/1152 - (6941*lmmu^3)/2592 - (1147*lmmu^4)/3456 - 
       (16187201*Pi^4)/52254720 - (481*lmmu*Pi^4)/2880 + (787*Pi^6)/81648 + 
       (145*Pi^4*Log[2])/1944 - (1924649*Pi^2*Log[2]^2)/4354560 - 
       (11*lmmu*Pi^2*Log[2]^2)/72 + (59*Pi^2*Log[2]^3)/972 + 
       (1924649*Log[2]^4)/4354560 + (11*lmmu*Log[2]^4)/72 - 
       (59*Log[2]^5)/1620 - (4*Log[2]^6)/81 + (Log[2]^5*Log[16])/81 + 
       (1924649*PolyLog[4, 1/2])/181440 + (11*lmmu*PolyLog[4, 1/2])/3 + 
       (118*PolyLog[5, 1/2])/27 + (353193131*Zeta[3])/40642560 + 
       (47317*lmmu*Zeta[3])/3072 + (175*lmmu^2*Zeta[3])/48 - 
       (1061*Zeta[3]^2)/576 + nl^2*(-17671/124416 + (3401*lmmu)/46656 - 
         (31*lmmu^2)/1296 - lmmu^4/864 + (7*Pi^4)/8640 + (5*Zeta[3])/864 - 
         (7*lmmu*Zeta[3])/108) + nl*(2261435/746496 - (36019*lmmu)/93312 + 
         (8701*lmmu^2)/10368 + (47*lmmu^3)/288 + (17*lmmu^4)/432 + 
         (245*Pi^4)/62208 + (163*lmmu*Pi^4)/12960 - (49*Pi^4*Log[2])/6480 + 
         (49*Pi^2*Log[2]^2)/2592 + (lmmu*Pi^2*Log[2]^2)/108 - 
         (Pi^2*Log[2]^3)/162 - (49*Log[2]^4)/2592 - (lmmu*Log[2]^4)/108 + 
         Log[2]^5/270 - (49*PolyLog[4, 1/2])/108 - (2*lmmu*PolyLog[4, 1/2])/
          9 - (4*PolyLog[5, 1/2])/9 + (1075*Zeta[3])/1728 - 
         (599*lmmu*Zeta[3])/1728 - (497*Zeta[5])/288) + 
       (59015*Zeta[5])/1728 - (575*lmmu*Zeta[5])/72)

   };

(ddd4 = b1b^4*g0b/4 - 3*b1b^2*b2b*g0b/4 + b2b^2*g0b/4 
 + b1b*b3b*g0b/2 - b4b*g0b/4 - b1b^3*g1b/4 + b1b*b2b*g1b/2
 - b3b*g1b/4 + b1b^2*g2b/4 - b2b*g2b/4 - b1b*g3b/4 + g4b/4
 /. {g0b->g0/b0, g1b->g1/b0, g2b->g2/b0, g3b->g3/b0, g4b->g4/b0,
     b1b->b1/b0, b2b->b2/b0, b3b->b3/b0, b4b->b4/b0 (*Global`b4num/b0*)});

ddd1 = g1/b0 - b1*g0/b0^2;
ddd2 = 1/2*(+ g2/b0 + b1^2*g0/b0^3 - b1*g1/b0^2 - b2*g0/b0^2);
ddd3 = 1/3*( g3/b0 - b1^3*g0/b0^4 + 2*b1*b2*g0/b0^3
	    - b3*g0/b0^2 + b1^2*g1/b0^3 - b2*g1/b0^2 - b1*g2/b0^2);

setMS2OS = {

  ms2msc ->  1 
      + ddd1 * x
      + ( 1/2 * ddd1^2 + ddd2 ) * x^2
      + ( 1/6*ddd1^3 + ddd1*ddd2 + ddd3) * x^3
      + ( ddd1^4/24 + ddd1^2*ddd2/2 + ddd2^2/2 + ddd1*ddd3 + ddd4) * x^4 ,

    msfromos -> 1 + 
	api*(-cf - (3*cf*lmM)/4) + api^2*((-1111*ca*cf)/384 + (7*cf^2)/128 -
   (185*ca*cf*lmM)/96 + (21*cf^2*lmM)/32 - (11*ca*cf*lmM^2)/32 +
   (9*cf^2*lmM^2)/32 + (143*cf*tr)/96 + (13*cf*lmM*tr)/24 + (cf*lmM^2*tr)/8 +
   (71*cf*nl*tr)/96 + (13*cf*lmM*nl*tr)/24 + (cf*lmM^2*nl*tr)/8 +
   (ca*cf*z2)/2 - (15*cf^2*z2)/8 - (3*ca*cf*log2*z2)/2 + 3*cf^2*log2*z2 -
   cf*tr*z2 + (cf*nl*tr*z2)/2 + (3*ca*cf*z3)/8 - (3*cf^2*z3)/4) +
 api^3*((lmM^2*(-2341*ca^2*cf + 1962*ca*cf^2 - 243*cf^3 + 1492*ca*cf*tr -
      468*cf^2*tr + 1492*ca*cf*nl*tr - 468*cf^2*nl*tr - 208*cf*tr^2 -
      416*cf*nl*tr^2 - 208*cf*nl^2*tr^2))/1152 +
   (lmM^3*(-242*ca^2*cf + 297*ca*cf^2 - 81*cf^3 + 176*ca*cf*tr -
      108*cf^2*tr + 176*ca*cf*nl*tr - 108*cf^2*nl*tr - 32*cf*tr^2 -
      64*cf*nl*tr^2 - 32*cf*nl^2*tr^2))/1152 +
   (lmM*(-105944*ca^2*cf + 52317*ca*cf^2 - 13203*cf^3 + 74624*ca*cf*tr -
      5436*cf^2*tr + 55616*ca*cf*nl*tr + 2340*cf^2*nl*tr - 12608*cf*tr^2 -
      18304*cf*nl*tr^2 - 5696*cf*nl^2*tr^2 + 12672*ca^2*cf*z2 -
      52704*ca*cf^2*z2 + 19440*cf^3*z2 - 38016*ca^2*cf*log2*z2 +
      91584*ca*cf^2*log2*z2 - 31104*cf^3*log2*z2 - 29952*ca*cf*tr*z2 +
      27648*cf^2*tr*z2 + 13824*ca*cf*log2*tr*z2 - 27648*cf^2*log2*tr*z2 +
      8064*ca*cf*nl*tr*z2 + 12096*cf^2*nl*tr*z2 + 13824*ca*cf*log2*nl*tr*z2 -
      27648*cf^2*log2*nl*tr*z2 + 9216*cf*tr^2*z2 + 4608*cf*nl*tr^2*z2 -
      4608*cf*nl^2*tr^2*z2 + 9504*ca^2*cf*z3 - 22896*ca*cf^2*z3 +
      7776*cf^3*z3 + 6912*ca*cf*tr*z3 - 3456*cf^2*tr*z3 +
      6912*ca*cf*nl*tr*z3 - 3456*cf^2*nl*tr*z3))/13824
	),

    msfromos3set[0] -> -202,
    msfromos3set[1] -> -176,
    msfromos3set[2] -> -150,
    msfromos3set[3] -> -126,
    msfromos3set[4] -> -103,
    msfromos3set[5] -> -82,

    msfromos3 -> -9478333/93312 - (644201*Pi^2)/38880 + (695*Pi^4)/7776 + 
 (587*Pi^2*Log[2])/162 + (22*Pi^2*Log[2]^2)/81 + (55*Log[2]^4)/162 + 
 (220*PolyLog[4, 1/2])/27 + nl^2*(-2353/23328 - (13*Pi^2)/324 - 
   (7*Zeta[3])/54) - (61*Zeta[3])/27 + (1439*Pi^2*Zeta[3])/432 + 
 nl*(246643/23328 + (967*Pi^2)/648 - (61*Pi^4)/1944 + (11*Pi^2*Log[2])/81 - 
   (2*Pi^2*Log[2]^2)/81 - Log[2]^4/81 - (8*PolyLog[4, 1/2])/27 + 
   (241*Zeta[3])/72) - (1975*Zeta[5])/216,


    osfromms -> 1 + 
	api*(cf + (3*cf*lmm)/4) + api^2*((1111*ca*cf)/384 - (71*cf^2)/128 -
   (143*cf*tr)/96 - (71*cf*nl*tr)/96 +
   lmm*((185*ca*cf)/96 - (9*cf^2)/32 - (13*cf*tr)/24 - (13*cf*nl*tr)/24) +
   lmm^2*((11*ca*cf)/32 + (9*cf^2)/32 - (cf*tr)/8 - (cf*nl*tr)/8) -
   (ca*cf*z2)/2 + (15*cf^2*z2)/8 + (3*ca*cf*log2*z2)/2 - 3*cf^2*log2*z2 +
   cf*tr*z2 - (cf*nl*tr*z2)/2 - (3*ca*cf*z3)/8 + (3*cf^2*z3)/4) +
   api^3*(lmm^3*((121*ca^2*cf)/576 + (33*ca*cf^2)/128 +
     (9*cf^3)/128 - (11*ca*cf*tr)/72 - (3*cf^2*tr)/32 - (11*ca*cf*nl*tr)/72 -
     (3*cf^2*nl*tr)/32 + (cf*tr^2)/36 + (cf*nl*tr^2)/18 +
     (cf*nl^2*tr^2)/36) + lmm^2*((2341*ca^2*cf)/1152 + (21*ca*cf^2)/64 -
     (63*cf^3)/128 - (373*ca*cf*tr)/288 - (3*cf^2*tr)/32 -
     (373*ca*cf*nl*tr)/288 - (3*cf^2*nl*tr)/32 + (13*cf*tr^2)/72 +
     (13*cf*nl*tr^2)/36 + (13*cf*nl^2*tr^2)/72) +
   lmm*((13243*ca^2*cf)/1728 - (4219*ca*cf^2)/1536 + (495*cf^3)/512 -
     (583*ca*cf*tr)/108 - (307*cf^2*tr)/384 - (869*ca*cf*nl*tr)/216 -
     (91*cf^2*nl*tr)/384 + (197*cf*tr^2)/216 + (143*cf*nl*tr^2)/108 +
     (89*cf*nl^2*tr^2)/216 - (11*ca^2*cf*z2)/12 + (49*ca*cf^2*z2)/16 +
     (45*cf^3*z2)/32 + (11*ca^2*cf*log2*z2)/4 - (35*ca*cf^2*log2*z2)/8 -
     (9*cf^3*log2*z2)/4 + (13*ca*cf*tr*z2)/6 - (cf^2*tr*z2)/2 -
     ca*cf*log2*tr*z2 + 2*cf^2*log2*tr*z2 - (7*ca*cf*nl*tr*z2)/12 -
     (13*cf^2*nl*tr*z2)/8 - ca*cf*log2*nl*tr*z2 + 2*cf^2*log2*nl*tr*z2 -
     (2*cf*tr^2*z2)/3 - (cf*nl*tr^2*z2)/3 + (cf*nl^2*tr^2*z2)/3 -
     (11*ca^2*cf*z3)/16 + (35*ca*cf^2*z3)/32 + (9*cf^3*z3)/16 -
     (ca*cf*tr*z3)/2 + (cf^2*tr*z3)/4 - (ca*cf*nl*tr*z3)/2 +
     (cf^2*nl*tr*z3)/4)),

    osfrommsset -> 1 + 
	api*(cf + (3*cf*lmm)/4) + api^2*((1111*ca*cf)/384 - (71*cf^2)/128 -
   (143*cf*tr)/96 - (71*cf*nl*tr)/96 +
   lmm*((185*ca*cf)/96 - (9*cf^2)/32 - (13*cf*tr)/24 - (13*cf*nl*tr)/24) +
   lmm^2*((11*ca*cf)/32 + (9*cf^2)/32 - (cf*tr)/8 - (cf*nl*tr)/8) -
   (ca*cf*z2)/2 + (15*cf^2*z2)/8 + (3*ca*cf*log2*z2)/2 - 3*cf^2*log2*z2 +
   cf*tr*z2 - (cf*nl*tr*z2)/2 - (3*ca*cf*z3)/8 + (3*cf^2*z3)/4) +
   api^3*(lmm^3*((121*ca^2*cf)/576 + (33*ca*cf^2)/128 +
     (9*cf^3)/128 - (11*ca*cf*tr)/72 - (3*cf^2*tr)/32 - (11*ca*cf*nl*tr)/72 -
     (3*cf^2*nl*tr)/32 + (cf*tr^2)/36 + (cf*nl*tr^2)/18 +
     (cf*nl^2*tr^2)/36) + lmm^2*((2341*ca^2*cf)/1152 + (21*ca*cf^2)/64 -
     (63*cf^3)/128 - (373*ca*cf*tr)/288 - (3*cf^2*tr)/32 -
     (373*ca*cf*nl*tr)/288 - (3*cf^2*nl*tr)/32 + (13*cf*tr^2)/72 +
     (13*cf*nl*tr^2)/36 + (13*cf*nl^2*tr^2)/72) - (ca*cf^2*z2)/4 +
   (15*cf^3*z2)/16 + (3*ca*cf^2*log2*z2)/4 - (3*cf^3*log2*z2)/2 +
   (cf^2*tr*z2)/2 - (cf^2*nl*tr*z2)/4 - (3*ca*cf^2*z3)/16 + (3*cf^3*z3)/8 +
   lmm*((13243*ca^2*cf)/1728 - (4219*ca*cf^2)/1536 + (495*cf^3)/512 -
     (583*ca*cf*tr)/108 - (307*cf^2*tr)/384 - (869*ca*cf*nl*tr)/216 -
     (91*cf^2*nl*tr)/384 + (197*cf*tr^2)/216 + (143*cf*nl*tr^2)/108 +
     (89*cf*nl^2*tr^2)/216 - (11*ca^2*cf*z2)/12 + (49*ca*cf^2*z2)/16 +
     (45*cf^3*z2)/32 + (11*ca^2*cf*log2*z2)/4 - (35*ca*cf^2*log2*z2)/8 -
     (9*cf^3*log2*z2)/4 + (13*ca*cf*tr*z2)/6 - (cf^2*tr*z2)/2 -
     ca*cf*log2*tr*z2 + 2*cf^2*log2*tr*z2 - (7*ca*cf*nl*tr*z2)/12 -
     (13*cf^2*nl*tr*z2)/8 - ca*cf*log2*nl*tr*z2 + 2*cf^2*log2*nl*tr*z2 -
     (2*cf*tr^2*z2)/3 - (cf*nl*tr^2*z2)/3 + (cf*nl^2*tr^2*z2)/3 -
     (11*ca^2*cf*z3)/16 + (35*ca*cf^2*z3)/32 + (9*cf^3*z3)/16 -
     (ca*cf*tr*z3)/2 + (cf^2*tr*z3)/4 - (ca*cf*nl*tr*z3)/2 +
     (cf^2*nl*tr*z3)/4)),

    osfromms3set[0] -> 194,
    osfromms3set[1] -> 168,
    osfromms3set[2] -> 143,
    osfromms3set[3] -> 119,
    osfromms3set[4] -> 96,
    osfromms3set[5] -> 75,

    osfromms3 -> 8481925/93312 + 
 (137*nl)/216 + (652841*Pi^2)/38880 - (nl*Pi^2)/27 - 
 (695*Pi^4)/7776 - (575*Pi^2*Log[2])/162 - (22*Pi^2*Log[2]^2)/81 - 
 (55*Log[2]^4)/162 - (220*PolyLog[4, 1/2])/27 - 
 nl^2*(-2353/23328 - (13*Pi^2)/324 - (7*Zeta[3])/54) + (58*Zeta[3])/27 - 
 (1439*Pi^2*Zeta[3])/432 - nl*(246643/23328 + (967*Pi^2)/648 - 
   (61*Pi^4)/1944 + (11*Pi^2*Log[2])/81 - (2*Pi^2*Log[2]^2)/81 - 
   Log[2]^4/81 - (8*PolyLog[4, 1/2])/27 + (241*Zeta[3])/72) + 
 (1975*Zeta[5])/216,

    mumfromos -> 1 
   - apiM*cf + apiM^2*((-1111*ca*cf)/384 + (199*cf^2)/128 + (143*cf*tr)/96 +
   (71*cf*nl*tr)/96 + (ca*cf*z2)/2 - (15*cf^2*z2)/8 - (3*ca*cf*log2*z2)/2 +
   3*cf^2*log2*z2 - cf*tr*z2 + (cf*nl*tr*z2)/2 + (3*ca*cf*z3)/8 -
   (3*cf^2*z3)/4),

    mumfromos3set[0] -> -170,
    mumfromos3set[1] -> -146,
    mumfromos3set[2] -> -123,
    mumfromos3set[3] -> -101,
    mumfromos3set[4] -> -81,
    mumfromos3set[5] -> -62,

    mumfromos3 -> -7172965/93312 - 
	(293*nl)/216 - (618281*Pi^2)/38880 - (nl*Pi^2)/9 + 
 (695*Pi^4)/7776 + (623*Pi^2*Log[2])/162 + (22*Pi^2*Log[2]^2)/81 + 
 (55*Log[2]^4)/162 + (220*PolyLog[4, 1/2])/27 + 
 nl^2*(-2353/23328 - (13*Pi^2)/324 - (7*Zeta[3])/54) - (70*Zeta[3])/27 + 
 (1439*Pi^2*Zeta[3])/432 + nl*(246643/23328 + (967*Pi^2)/648 - 
   (61*Pi^4)/1944 + (11*Pi^2*Log[2])/81 - (2*Pi^2*Log[2]^2)/81 - 
   Log[2]^4/81 - (8*PolyLog[4, 1/2])/27 + (241*Zeta[3])/72) - 
 (1975*Zeta[5])/216,

    msfromri -> 1 - (4*api)/3 + api^2*(-995/72 + (89*nf)/144 + (19*z3)/6) + 
 api^3*(-6663911/41472 + (118325*nf)/7776 - (4459*nf^2)/23328 + 
   (5*nf*z4)/12 - (185*z5)/36 + (408007*z3)/6912 - 
   (617*nf*z3)/216 - (nf^2*z3)/54),

    rifromms -> 1 + (4*api)/3 + api^2*(1123/72 - (89*nf)/144 - (19*z3)/6) + 
 api^3*(6663911/41472 - (118325*nf)/7776 + (4459*nf^2)/23328 + 
   (4*(1123/72 - (89*nf)/144 - (19*z3)/6))/3 - (408007*z3)/6912 + 
   (617*nf*z3)/216 + (nf^2*z3)/54 - (4*(-995/72 + (89*nf)/144 + (19*z3)/6))/
    3 - (5*nf*z4)/12 + (185*z5)/36)

    };

(* ************************************************************ *)

(* Expressions for mKIN <-> mMS relation: *)


RunDecmMS2mkinnl = mbMSmus + apinlmus*((-16*mufac)/9 - 
      (2*mufac^2)/(3*mbMSmus) + mbMSmus*(4/3 + Log[mus^2/mbMSmus^2])) + 
    apinlmus^2*((mcMSmusmc*Pi^2*TAGnc)/6 + (mcMSmusmc^3*Pi^2*TAGnc)/
       (6*mbMSmus^2) + ((19*mcMSmusmc^6*TAGnc)/225 - 
        (4*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mbMSmus])/45)/mbMSmus^5 + 
      ((463*mcMSmusmc^8*TAGnc)/39200 - (3*mcMSmusmc^8*TAGnc*
          Log[mcMSmusmc/mbMSmus])/140)/mbMSmus^7 + 
      ((997*mcMSmusmc^10*TAGnc)/297675 - (8*mcMSmusmc^10*TAGnc*
          Log[mcMSmusmc/mbMSmus])/945)/mbMSmus^9 + 
      ((1229*mcMSmusmc^12*TAGnc)/940896 - (5*mcMSmusmc^12*TAGnc*
          Log[mcMSmusmc/mbMSmus])/1188)/mbMSmus^11 + 
      ((15371*mcMSmusmc^14*TAGnc)/25050025 - (12*mcMSmusmc^14*TAGnc*
          Log[mcMSmusmc/mbMSmus])/5005)/mbMSmus^13 + 
      ((529*mcMSmusmc^16*TAGnc)/1622400 - (7*mcMSmusmc^16*TAGnc*
          Log[mcMSmusmc/mbMSmus])/4680)/mbMSmus^15 + 
      ((16277*mcMSmusmc^18*TAGnc)/86028075 - (16*mcMSmusmc^18*TAGnc*
          Log[mcMSmusmc/mbMSmus])/16065)/mbMSmus^17 + 
      ((39163*mcMSmusmc^20*TAGnc)/333852800 - 
        (9*mcMSmusmc^20*TAGnc*Log[mcMSmusmc/mbMSmus])/12920)/mbMSmus^19 + 
      ((39833*mcMSmusmc^22*TAGnc)/520109667 - 
        (20*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mbMSmus])/39501)/mbMSmus^21 + 
      ((9727*mcMSmusmc^24*TAGnc)/186631200 - (11*mcMSmusmc^24*TAGnc*
          Log[mcMSmusmc/mbMSmus])/28980)/mbMSmus^23 + 
      ((-151*mcMSmusmc^4*TAGnc)/216 - (mcMSmusmc^4*Pi^2*TAGnc)/18 + 
        (13*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mbMSmus])/18 - 
        (mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mbMSmus]^2)/3)/mbMSmus^3 + 
      mufac*(-860/27 + (128*NLxOSKIN)/81 + (8*Pi^2)/9 + 
        (88/9 - (16*NLxOSKIN)/27)*Log[(2*mufac)/mus]) + 
      (-(mcMSmusmc^2*TAGnc) + mufac^2*(-83/9 + (13*NLxOSKIN)/27 + Pi^2/3 + 
          (11/3 - (2*NLxOSKIN)/9)*Log[(2*mufac)/mus] + 
          (2*Log[mus^2/mbMSmus^2])/3))/mbMSmus + 
      mbMSmus*(307/32 + Pi^2/3 - (71*TAGnc)/144 - (Pi^2*TAGnc)/18 + 
        (Pi^2*Log[2])/9 + (47/24 - TAGnc/12)*Log[mus^2/mbMSmus^2]^2 + 
        NLxMSOS*(-71/144 - Pi^2/18 - (13*Log[mus^2/mbMSmus^2])/36 - 
          Log[mus^2/mbMSmus^2]^2/12) + (2*TAGdecmcMSOS*
          Log[mus^2/mcMSmusmc^2])/9 + Log[mus^2/mbMSmus^2]*
         (509/72 - (13*TAGnc)/36 + (TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/6) - 
        Zeta[3]/6)) + apinlmus^3*((5971*mcMSmusmc*Pi^2*TAGnc)/486 + 
      (13*mcMSmusmc*Pi^3*TAGnc)/162 - (4*mcMSmusmc*Pi^2*TAGnc^2)/135 - 
      (1199*mcMSmusmc*Pi^2*TAGnc*Log[2])/81 + 
      ((-11*mcMSmusmc*Pi^2*TAGnc)/6 + (mcMSmusmc*Pi^2*TAGnc^2)/9)*
       Log[mcMSmusmc/mbMSmus] + 
      (mufac^2*((9727*mcMSmusmc^24*TAGnc)/279946800 - 
         (11*mcMSmusmc^24*TAGnc*Log[mcMSmusmc/mbMSmus])/43470))/mbMSmus^25 + 
      ((149*mcMSmusmc^5*Pi^2*TAGnc)/1080 + (mcMSmusmc^3*mufac^2*Pi^2*TAGnc)/
         9 - (mcMSmusmc^5*NLxMSOS*Pi^2*TAGnc)/20 + (mcMSmusmc^5*Pi^3*TAGnc)/
         240 - (13*mcMSmusmc^5*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/45)/
       mbMSmus^4 + ((262769*mcMSmusmc^7*Pi^2*TAGnc)/11907000 - 
        (5*mcMSmusmc^7*NLxMSOS*Pi^2*TAGnc)/756 - (25*mcMSmusmc^7*Pi^3*TAGnc)/
         18144 - (29*mcMSmusmc^7*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/525)/
       mbMSmus^6 + ((416029*mcMSmusmc^9*Pi^2*TAGnc)/60011280 - 
        (7*mcMSmusmc^9*NLxMSOS*Pi^2*TAGnc)/3240 - (77*mcMSmusmc^9*Pi^3*TAGnc)/
         103680 - (53*mcMSmusmc^9*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/2205)/
       mbMSmus^8 + ((1055689*mcMSmusmc^11*Pi^2*TAGnc)/345779280 - 
        (3*mcMSmusmc^11*NLxMSOS*Pi^2*TAGnc)/3080 - 
        (17*mcMSmusmc^11*Pi^3*TAGnc)/39424 - (85*mcMSmusmc^11*Pi^2*TAGnc*
          Log[mcMSmusmc/mbMSmus])/6237)/mbMSmus^10 + 
      ((31786481*mcMSmusmc^13*Pi^2*TAGnc)/19677663720 - 
        (11*mcMSmusmc^13*NLxMSOS*Pi^2*TAGnc)/21060 - 
        (1771*mcMSmusmc^13*Pi^3*TAGnc)/6469632 - 
        (125*mcMSmusmc^13*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/14157)/
       mbMSmus^12 + ((17346493*mcMSmusmc^15*Pi^2*TAGnc)/18087549480 - 
        (13*mcMSmusmc^15*NLxMSOS*Pi^2*TAGnc)/41580 - 
        (377*mcMSmusmc^15*Pi^3*TAGnc)/2027520 - 
        (173*mcMSmusmc^15*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/27885)/
       mbMSmus^14 + ((22757641*mcMSmusmc^17*Pi^2*TAGnc)/36923796000 - 
        (5*mcMSmusmc^17*NLxMSOS*Pi^2*TAGnc)/24752 - 
        (1925*mcMSmusmc^17*Pi^3*TAGnc)/14483456 - 
        (229*mcMSmusmc^17*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/49725)/
       mbMSmus^16 + ((86836957*mcMSmusmc^19*Pi^2*TAGnc)/206871887520 - 
        (17*mcMSmusmc^19*NLxMSOS*Pi^2*TAGnc)/123120 - 
        (99671*mcMSmusmc^19*Pi^3*TAGnc)/1008599040 - 
        (293*mcMSmusmc^19*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/82365)/
       mbMSmus^18 + ((846435761*mcMSmusmc^21*Pi^2*TAGnc)/2832319518840 - 
        (19*mcMSmusmc^21*NLxMSOS*Pi^2*TAGnc)/192780 - 
        (127699*mcMSmusmc^21*Pi^3*TAGnc)/1684537344 - 
        (365*mcMSmusmc^21*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/128877)/
       mbMSmus^20 + ((171475369*mcMSmusmc^23*Pi^2*TAGnc)/778168119960 - 
        (7*mcMSmusmc^23*NLxMSOS*Pi^2*TAGnc)/96140 - 
        (81991*mcMSmusmc^23*Pi^3*TAGnc)/1374683136 - 
        (445*mcMSmusmc^23*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/192717)/
       mbMSmus^22 + ((106559417*mcMSmusmc^25*Pi^2*TAGnc)/637438863600 - 
        (23*mcMSmusmc^25*NLxMSOS*Pi^2*TAGnc)/415800 - 
        (5698043*mcMSmusmc^25*Pi^3*TAGnc)/118908518400 - 
        (533*mcMSmusmc^25*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/277725)/
       mbMSmus^24 + ((11*mcMSmusmc*Pi^2*TAGnc)/12 - (mcMSmusmc*Pi^2*TAGnc^2)/
         18)*Log[mus^2/mbMSmus^2] + NLxMSOS*((-7*mcMSmusmc*Pi^2*TAGnc)/27 + 
        (2*mcMSmusmc*Pi^2*TAGnc*Log[2])/9 + (mcMSmusmc*Pi^2*TAGnc*
          Log[mcMSmusmc/mbMSmus])/9 - (mcMSmusmc*Pi^2*TAGnc*
          Log[mus^2/mbMSmus^2])/18) + (mcMSmusmc*Pi^2*TAGdecmcMSOS*TAGnc*
        Log[mus^2/mcMSmusmc^2])/18 + (mcMSmusmc*Pi^2*TAGnc*
        Log[musmc^2/mcMSmusmc^2])/6 + ((18607*mcMSmusmc^3*Pi^2*TAGnc)/1620 + 
        (mcMSmusmc*mufac^2*Pi^2*TAGnc)/9 + (7*mcMSmusmc^3*Pi^3*TAGnc)/108 - 
        (1199*mcMSmusmc^3*Pi^2*TAGnc*Log[2])/81 + 
        ((-289*mcMSmusmc^3*Pi^2*TAGnc)/162 + (mcMSmusmc^3*Pi^2*TAGnc^2)/9)*
         Log[mcMSmusmc/mbMSmus] + ((7*mcMSmusmc^3*Pi^2*TAGnc)/12 - 
          (mcMSmusmc^3*Pi^2*TAGnc^2)/18)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((-7*mcMSmusmc^3*Pi^2*TAGnc)/54 + 
          (2*mcMSmusmc^3*Pi^2*TAGnc*Log[2])/9 + (mcMSmusmc^3*Pi^2*TAGnc*
            Log[mcMSmusmc/mbMSmus])/9 - (mcMSmusmc^3*Pi^2*TAGnc*
            Log[mus^2/mbMSmus^2])/18) + (mcMSmusmc^3*Pi^2*TAGdecmcMSOS*TAGnc*
          Log[mus^2/mcMSmusmc^2])/18 + (mcMSmusmc^3*Pi^2*TAGnc*
          Log[musmc^2/mcMSmusmc^2])/2)/mbMSmus^2 + 
      mufac*(-130867/162 + (20047*NLxOSKIN)/243 - (1292*NLxOSKIN^2)/729 + 
        (1022*Pi^2)/27 - (208*NLxOSKIN*Pi^2)/81 + (8*NLxOSKIN^2*Pi^2)/243 - 
        (2*Pi^4)/3 + (10072/27 - (3356*NLxOSKIN)/81 + (256*NLxOSKIN^2)/243 - 
          (88*Pi^2)/9 + (16*NLxOSKIN*Pi^2)/27)*Log[(2*mufac)/mus] + 
        (-484/9 + (176*NLxOSKIN)/27 - (16*NLxOSKIN^2)/81)*
         Log[(2*mufac)/mus]^2 + 114*Zeta[3] - (140*NLxOSKIN*Zeta[3])/27) + 
      ((-3529813*mcMSmusmc^6*TAGnc)/3499200 + (1427*mcMSmusmc^6*Pi^2*TAGnc)/
         12960 - (3481*mcMSmusmc^6*TAGnc^2)/30375 + 
        (4*mcMSmusmc^6*Pi^2*TAGnc^2)/405 - (11*mcMSmusmc^6*Pi^2*TAGnc*Log[2])/
         648 - (5329*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mbMSmus]^2)/6480 + 
        mufac^2*((-151*mcMSmusmc^4*TAGnc)/324 - (mcMSmusmc^4*Pi^2*TAGnc)/27 + 
          (13*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mbMSmus])/27 - 
          (2*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mbMSmus]^2)/9) + 
        ((59*mcMSmusmc^6*TAGnc)/450 - (19*mcMSmusmc^6*TAGnc^2)/675)*
         Log[mus^2/mbMSmus^2] + NLxMSOS*((2729*mcMSmusmc^6*TAGnc)/30375 + 
          (4*mcMSmusmc^6*Pi^2*TAGnc)/405 - (19*mcMSmusmc^6*TAGnc*
            Log[mus^2/mbMSmus^2])/675 + Log[mcMSmusmc/mbMSmus]*
           ((-2*mcMSmusmc^6*TAGnc)/75 + (4*mcMSmusmc^6*TAGnc*
              Log[mus^2/mbMSmus^2])/135)) + (19*mcMSmusmc^6*TAGdecmcMSOS*
          TAGnc*Log[mus^2/mcMSmusmc^2])/675 + 
        (94*mcMSmusmc^6*TAGnc*Log[musmc^2/mcMSmusmc^2])/225 + 
        Log[mcMSmusmc/mbMSmus]*((1481549*mcMSmusmc^6*TAGnc)/1166400 - 
          (113*mcMSmusmc^6*Pi^2*TAGnc)/1296 + (8*mcMSmusmc^6*TAGnc^2)/75 + 
          ((-2*mcMSmusmc^6*TAGnc)/45 + (4*mcMSmusmc^6*TAGnc^2)/135)*
           Log[mus^2/mbMSmus^2] - (4*mcMSmusmc^6*TAGdecmcMSOS*TAGnc*
            Log[mus^2/mcMSmusmc^2])/135 - (8*mcMSmusmc^6*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/15) - (29*mcMSmusmc^6*TAGnc*Zeta[3])/
         648)/mbMSmus^5 + ((-33561121*mcMSmusmc^8*TAGnc)/1185408000 + 
        (1789*mcMSmusmc^8*Pi^2*TAGnc)/90720 - (241*mcMSmusmc^8*TAGnc^2)/
         7056 - (7*mcMSmusmc^8*Pi^2*TAGnc*Log[2])/1536 + 
        ((-9019*mcMSmusmc^8*TAGnc)/30240 - (mcMSmusmc^8*TAGnc^2)/70)*
         Log[mcMSmusmc/mbMSmus]^2 + mufac^2*((38*mcMSmusmc^6*TAGnc)/675 - 
          (8*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mbMSmus])/135) + 
        ((291*mcMSmusmc^8*TAGnc)/78400 - (463*mcMSmusmc^8*TAGnc^2)/117600)*
         Log[mus^2/mbMSmus^2] + NLxMSOS*((174787*mcMSmusmc^8*TAGnc)/
           12348000 + (mcMSmusmc^8*Pi^2*TAGnc)/420 - 
          (463*mcMSmusmc^8*TAGnc*Log[mus^2/mbMSmus^2])/117600 + 
          Log[mcMSmusmc/mbMSmus]*((-37*mcMSmusmc^8*TAGnc)/14700 + 
            (mcMSmusmc^8*TAGnc*Log[mus^2/mbMSmus^2])/140)) + 
        (463*mcMSmusmc^8*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/117600 + 
        (179*mcMSmusmc^8*TAGnc*Log[musmc^2/mcMSmusmc^2])/2450 + 
        Log[mcMSmusmc/mbMSmus]*((2308303*mcMSmusmc^8*TAGnc)/12700800 - 
          (139*mcMSmusmc^8*Pi^2*TAGnc)/4608 + (193*mcMSmusmc^8*TAGnc^2)/
           3675 + ((9*mcMSmusmc^8*TAGnc)/280 + (mcMSmusmc^8*TAGnc^2)/140)*
           Log[mus^2/mbMSmus^2] - (mcMSmusmc^8*TAGdecmcMSOS*TAGnc*
            Log[mus^2/mcMSmusmc^2])/140 - (6*mcMSmusmc^8*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/35) - (173*mcMSmusmc^8*TAGnc*Zeta[3])/
         9216)/mbMSmus^7 + ((47157278431*mcMSmusmc^10*TAGnc)/5184974592000 + 
        (359801*mcMSmusmc^10*Pi^2*TAGnc)/48988800 - 
        (1778*mcMSmusmc^10*TAGnc^2)/273375 - (17*mcMSmusmc^10*Pi^2*TAGnc*
          Log[2])/8640 + ((-2585963*mcMSmusmc^10*TAGnc)/16329600 - 
          (16*mcMSmusmc^10*TAGnc^2)/2835)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((463*mcMSmusmc^8*TAGnc)/58800 - 
          (mcMSmusmc^8*TAGnc*Log[mcMSmusmc/mbMSmus])/70) + 
        ((-277*mcMSmusmc^10*TAGnc)/85050 - (997*mcMSmusmc^10*TAGnc^2)/893025)*
         Log[mus^2/mbMSmus^2] + NLxMSOS*((2816347*mcMSmusmc^10*TAGnc)/
           562605750 + (8*mcMSmusmc^10*Pi^2*TAGnc)/8505 - 
          (997*mcMSmusmc^10*TAGnc*Log[mus^2/mbMSmus^2])/893025 + 
          Log[mcMSmusmc/mbMSmus]*((-398*mcMSmusmc^10*TAGnc)/893025 + 
            (8*mcMSmusmc^10*TAGnc*Log[mus^2/mbMSmus^2])/2835)) + 
        (997*mcMSmusmc^10*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/893025 + 
        (298*mcMSmusmc^10*TAGnc*Log[musmc^2/mcMSmusmc^2])/11907 + 
        Log[mcMSmusmc/mbMSmus]*((2370226609*mcMSmusmc^10*TAGnc)/41150592000 - 
          (17*mcMSmusmc^10*Pi^2*TAGnc)/1080 + (14048*mcMSmusmc^10*TAGnc^2)/
           893025 + ((4*mcMSmusmc^10*TAGnc)/135 + (8*mcMSmusmc^10*TAGnc^2)/
             2835)*Log[mus^2/mbMSmus^2] - (8*mcMSmusmc^10*TAGdecmcMSOS*TAGnc*
            Log[mus^2/mcMSmusmc^2])/2835 - (16*mcMSmusmc^10*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/189) - (187*mcMSmusmc^10*TAGnc*Zeta[3])/
         17280)/mbMSmus^9 + ((237903978448457*mcMSmusmc^12*TAGnc)/
         24647147078400000 + (84041429*mcMSmusmc^12*Pi^2*TAGnc)/22992076800 - 
        (3989*mcMSmusmc^12*TAGnc^2)/1881792 - (175*mcMSmusmc^12*Pi^2*TAGnc*
          Log[2])/165888 + ((-14841023*mcMSmusmc^12*TAGnc)/149688000 - 
          (5*mcMSmusmc^12*TAGnc^2)/1782)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((1994*mcMSmusmc^10*TAGnc)/893025 - 
          (16*mcMSmusmc^10*TAGnc*Log[mcMSmusmc/mbMSmus])/2835) + 
        ((-509*mcMSmusmc^12*TAGnc)/171072 - (1229*mcMSmusmc^12*TAGnc^2)/
           2822688)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((326802499*mcMSmusmc^12*TAGnc)/136928594880 + 
          (5*mcMSmusmc^12*Pi^2*TAGnc)/10692 - (1229*mcMSmusmc^12*TAGnc*
            Log[mus^2/mbMSmus^2])/2822688 + Log[mcMSmusmc/mbMSmus]*
           ((-391*mcMSmusmc^12*TAGnc)/4939704 + (5*mcMSmusmc^12*TAGnc*
              Log[mus^2/mbMSmus^2])/3564)) + (1229*mcMSmusmc^12*TAGdecmcMSOS*
          TAGnc*Log[mus^2/mcMSmusmc^2])/2822688 + 
        (899*mcMSmusmc^12*TAGnc*Log[musmc^2/mcMSmusmc^2])/78408 + 
        Log[mcMSmusmc/mbMSmus]*((160430424613*mcMSmusmc^12*TAGnc)/
           6224027040000 - (40553*mcMSmusmc^12*Pi^2*TAGnc)/4147200 + 
          (9545*mcMSmusmc^12*TAGnc^2)/1411344 + ((5*mcMSmusmc^12*TAGnc)/216 + 
            (5*mcMSmusmc^12*TAGnc^2)/3564)*Log[mus^2/mbMSmus^2] - 
          (5*mcMSmusmc^12*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/3564 - 
          (5*mcMSmusmc^12*TAGnc*Log[musmc^2/mcMSmusmc^2])/99) - 
        (59231*mcMSmusmc^12*TAGnc*Zeta[3])/8294400)/mbMSmus^11 + 
      ((448931017732868509*mcMSmusmc^14*TAGnc)/58963096098466560000 + 
        (82285201*mcMSmusmc^14*Pi^2*TAGnc)/38745907200 - 
        (233201*mcMSmusmc^14*TAGnc^2)/260851500 - 
        (23*mcMSmusmc^14*Pi^2*TAGnc*Log[2])/35840 + 
        ((-4128997*mcMSmusmc^14*TAGnc)/60540480 - (8*mcMSmusmc^14*TAGnc^2)/
           5005)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((1229*mcMSmusmc^12*TAGnc)/1411344 - 
          (5*mcMSmusmc^12*TAGnc*Log[mcMSmusmc/mbMSmus])/1782) + 
        ((-22089*mcMSmusmc^14*TAGnc)/10020010 - (15371*mcMSmusmc^14*TAGnc^2)/
           75150075)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((108352091581*mcMSmusmc^14*TAGnc)/81243243081000 + 
          (4*mcMSmusmc^14*Pi^2*TAGnc)/15015 - (15371*mcMSmusmc^14*TAGnc*
            Log[mus^2/mbMSmus^2])/75150075 + Log[mcMSmusmc/mbMSmus]*
           ((1153*mcMSmusmc^14*TAGnc)/225450225 + (4*mcMSmusmc^14*TAGnc*
              Log[mus^2/mbMSmus^2])/5005)) + (15371*mcMSmusmc^14*TAGdecmcMSOS*
          TAGnc*Log[mus^2/mcMSmusmc^2])/75150075 + 
        (3166*mcMSmusmc^14*TAGnc*Log[musmc^2/mcMSmusmc^2])/511225 + 
        Log[mcMSmusmc/mbMSmus]*((570692369939*mcMSmusmc^14*TAGnc)/
           40905688824000 - (2161*mcMSmusmc^14*Pi^2*TAGnc)/322560 + 
          (263546*mcMSmusmc^14*TAGnc^2)/75150075 + 
          ((18*mcMSmusmc^14*TAGnc)/1001 + (4*mcMSmusmc^14*TAGnc^2)/5005)*
           Log[mus^2/mbMSmus^2] - (4*mcMSmusmc^14*TAGdecmcMSOS*TAGnc*
            Log[mus^2/mcMSmusmc^2])/5005 - (24*mcMSmusmc^14*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/715) - (3287*mcMSmusmc^14*TAGnc*
          Zeta[3])/645120)/mbMSmus^13 + 
      ((558828635339341201127*mcMSmusmc^16*TAGnc)/95095681387606867968000 + 
        (58806560951*mcMSmusmc^16*Pi^2*TAGnc)/43261891706880 - 
        (1373*mcMSmusmc^16*TAGnc^2)/3110400 - (1001*mcMSmusmc^16*Pi^2*TAGnc*
          Log[2])/2359296 + ((-1141880431*mcMSmusmc^16*TAGnc)/22884301440 - 
          (7*mcMSmusmc^16*TAGnc^2)/7020)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((30742*mcMSmusmc^14*TAGnc)/75150075 - 
          (8*mcMSmusmc^14*TAGnc*Log[mcMSmusmc/mbMSmus])/5005) + 
        ((-15593*mcMSmusmc^16*TAGnc)/9734400 - (529*mcMSmusmc^16*TAGnc^2)/
           4867200)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((214558103603*mcMSmusmc^16*TAGnc)/260460712512000 + 
          (7*mcMSmusmc^16*Pi^2*TAGnc)/42120 - (529*mcMSmusmc^16*TAGnc*
            Log[mus^2/mbMSmus^2])/4867200 + Log[mcMSmusmc/mbMSmus]*
           ((355*mcMSmusmc^16*TAGnc)/14455584 + (7*mcMSmusmc^16*TAGnc*
              Log[mus^2/mbMSmus^2])/14040)) + (529*mcMSmusmc^16*TAGdecmcMSOS*
          TAGnc*Log[mus^2/mcMSmusmc^2])/4867200 + 
        (283*mcMSmusmc^16*TAGnc*Log[musmc^2/mcMSmusmc^2])/76050 + 
        Log[mcMSmusmc/mbMSmus]*((10019690708971*mcMSmusmc^16*TAGnc)/
           1178083838131200 - (565351*mcMSmusmc^16*Pi^2*TAGnc)/115605504 + 
          (7469*mcMSmusmc^16*TAGnc^2)/3650400 + 
          ((133*mcMSmusmc^16*TAGnc)/9360 + (7*mcMSmusmc^16*TAGnc^2)/14040)*
           Log[mus^2/mbMSmus^2] - (7*mcMSmusmc^16*TAGdecmcMSOS*TAGnc*
            Log[mus^2/mcMSmusmc^2])/14040 - (14*mcMSmusmc^16*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/585) - (885457*mcMSmusmc^16*TAGnc*
          Zeta[3])/231211008)/mbMSmus^15 + 
      ((230500151130233896919629*mcMSmusmc^18*TAGnc)/
         50057687427569200963584000 + (26463251891*mcMSmusmc^18*Pi^2*TAGnc)/
         28454994247680 - (157882*mcMSmusmc^18*TAGnc^2)/650372247 - 
        (4147*mcMSmusmc^18*Pi^2*TAGnc*Log[2])/13934592 + 
        ((-3391310411*mcMSmusmc^18*TAGnc)/88921857024 - 
          (32*mcMSmusmc^18*TAGnc^2)/48195)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((529*mcMSmusmc^16*TAGnc)/2433600 - 
          (7*mcMSmusmc^16*TAGnc*Log[mcMSmusmc/mbMSmus])/7020) + 
        ((-203011*mcMSmusmc^18*TAGnc)/172056150 - 
          (16277*mcMSmusmc^18*TAGnc^2)/258084225)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((20555048260909*mcMSmusmc^18*TAGnc)/37681809223558500 + 
          (16*mcMSmusmc^18*Pi^2*TAGnc)/144585 - (16277*mcMSmusmc^18*TAGnc*
            Log[mus^2/mbMSmus^2])/258084225 + Log[mcMSmusmc/mbMSmus]*
           ((197062*mcMSmusmc^18*TAGnc)/7381208835 + 
            (16*mcMSmusmc^18*TAGnc*Log[mus^2/mbMSmus^2])/48195)) + 
        (16277*mcMSmusmc^18*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/
         258084225 + (7678*mcMSmusmc^18*TAGnc*Log[musmc^2/mcMSmusmc^2])/
         3186225 + Log[mcMSmusmc/mbMSmus]*((91862575617799009*mcMSmusmc^18*
            TAGnc)/16342379002556006400 - (52013*mcMSmusmc^18*Pi^2*TAGnc)/
           13934592 + (334304*mcMSmusmc^18*TAGnc^2)/258084225 + 
          ((184*mcMSmusmc^18*TAGnc)/16065 + (16*mcMSmusmc^18*TAGnc^2)/48195)*
           Log[mus^2/mbMSmus^2] - (16*mcMSmusmc^18*TAGdecmcMSOS*TAGnc*
            Log[mus^2/mcMSmusmc^2])/48195 - (32*mcMSmusmc^18*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/1785) - (83291*mcMSmusmc^18*TAGnc*
          Zeta[3])/27869184)/mbMSmus^17 + 
      ((200744102824691620675013329*mcMSmusmc^20*TAGnc)/
         54499313978682087207813120000 + (9007367733163*mcMSmusmc^20*Pi^2*
          TAGnc)/13481015456563200 - (367909*mcMSmusmc^20*TAGnc^2)/
         2547216000 - (143*mcMSmusmc^20*Pi^2*TAGnc*Log[2])/655360 + 
        ((-45436526101*mcMSmusmc^20*TAGnc)/1508495788800 - 
          (3*mcMSmusmc^20*TAGnc^2)/6460)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((32554*mcMSmusmc^18*TAGnc)/258084225 - 
          (32*mcMSmusmc^18*TAGnc*Log[mcMSmusmc/mbMSmus])/48195) + 
        ((-592281*mcMSmusmc^20*TAGnc)/667705600 - 
          (39163*mcMSmusmc^20*TAGnc^2)/1001558400)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((22189567531163017*mcMSmusmc^20*TAGnc)/
           58347124817357376000 + (mcMSmusmc^20*Pi^2*TAGnc)/12920 - 
          (39163*mcMSmusmc^20*TAGnc*Log[mus^2/mbMSmus^2])/1001558400 + 
          Log[mcMSmusmc/mbMSmus]*((1211963*mcMSmusmc^20*TAGnc)/50127997920 + 
            (3*mcMSmusmc^20*TAGnc*Log[mus^2/mbMSmus^2])/12920)) + 
        (39163*mcMSmusmc^20*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/
         1001558400 + (5507*mcMSmusmc^20*TAGnc*Log[musmc^2/mcMSmusmc^2])/
         3338528 + Log[mcMSmusmc/mbMSmus]*((691405824368106269*mcMSmusmc^20*
            TAGnc)/175583298211985664000 - (156353*mcMSmusmc^20*Pi^2*TAGnc)/
           53084160 + (1017529*mcMSmusmc^20*TAGnc^2)/1168484800 + 
          ((243*mcMSmusmc^20*TAGnc)/25840 + (3*mcMSmusmc^20*TAGnc^2)/12920)*
           Log[mus^2/mbMSmus^2] - (3*mcMSmusmc^20*TAGdecmcMSOS*TAGnc*
            Log[mus^2/mcMSmusmc^2])/12920 - (9*mcMSmusmc^20*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/646) - (254791*mcMSmusmc^20*TAGnc*
          Zeta[3])/106168320)/mbMSmus^19 + 
      ((32427296704361728583901633193*mcMSmusmc^22*TAGnc)/
         10790864167779053267146997760000 + (174200135864459*mcMSmusmc^22*
          Pi^2*TAGnc)/349543472195174400 - (624853*mcMSmusmc^22*TAGnc^2)/
         6846429744 - (38675*mcMSmusmc^22*Pi^2*TAGnc*Log[2])/233570304 + 
        ((-3036958280711*mcMSmusmc^22*TAGnc)/124450902576000 - 
          (40*mcMSmusmc^22*TAGnc^2)/118503)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((39163*mcMSmusmc^20*TAGnc)/500779200 - 
          (3*mcMSmusmc^20*TAGnc*Log[mcMSmusmc/mbMSmus])/6460) + 
        ((-708143*mcMSmusmc^22*TAGnc)/1040219334 - 
          (39833*mcMSmusmc^22*TAGnc^2)/1560329001)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((1640519393726677*mcMSmusmc^22*TAGnc)/5946258455651273760 + 
          (20*mcMSmusmc^22*Pi^2*TAGnc)/355509 - (39833*mcMSmusmc^22*TAGnc*
            Log[mus^2/mbMSmus^2])/1560329001 + Log[mcMSmusmc/mbMSmus]*
           ((14290513*mcMSmusmc^22*TAGnc)/689665418442 + 
            (20*mcMSmusmc^22*TAGnc*Log[mus^2/mbMSmus^2])/118503)) + 
        (39833*mcMSmusmc^22*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/
         1560329001 + (5066*mcMSmusmc^22*TAGnc*Log[musmc^2/mcMSmusmc^2])/
         4298427 + Log[mcMSmusmc/mbMSmus]*((20875025059584870877*mcMSmusmc^22*
            TAGnc)/7242811051244408640000 - (13926181*mcMSmusmc^22*Pi^2*
            TAGnc)/5839257600 + (318845*mcMSmusmc^22*TAGnc^2)/520109667 + 
          ((310*mcMSmusmc^22*TAGnc)/39501 + (20*mcMSmusmc^22*TAGnc^2)/118503)*
           Log[mus^2/mbMSmus^2] - (20*mcMSmusmc^22*TAGdecmcMSOS*TAGnc*
            Log[mus^2/mcMSmusmc^2])/118503 - (40*mcMSmusmc^22*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/3591) - (23017987*mcMSmusmc^22*TAGnc*
          Zeta[3])/11678515200)/mbMSmus^21 + 
      ((5912372795358292426094681002539449*mcMSmusmc^24*TAGnc)/
         2369996943791663839368455777157120000 + 
        (240817793781176357*mcMSmusmc^24*Pi^2*TAGnc)/628867544642695987200 - 
        (279881*mcMSmusmc^24*TAGnc^2)/4627692000 - 
        (877591*mcMSmusmc^24*Pi^2*TAGnc*Log[2])/6794772480 + 
        ((-941278536953*mcMSmusmc^24*TAGnc)/46646042002560 - 
          (11*mcMSmusmc^24*TAGnc^2)/43470)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((79666*mcMSmusmc^22*TAGnc)/1560329001 - 
          (40*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mbMSmus])/118503) + 
        ((-631*mcMSmusmc^24*TAGnc)/1184960 - (9727*mcMSmusmc^24*TAGnc^2)/
           559893600)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((7801530877413386647*mcMSmusmc^24*TAGnc)/
           37763267488425775296000 + (11*mcMSmusmc^24*Pi^2*TAGnc)/260820 - 
          (9727*mcMSmusmc^24*TAGnc*Log[mus^2/mbMSmus^2])/559893600 + 
          Log[mcMSmusmc/mbMSmus]*((73801799*mcMSmusmc^24*TAGnc)/
             4231787807520 + (11*mcMSmusmc^24*TAGnc*Log[mus^2/mbMSmus^2])/
             86940)) + (9727*mcMSmusmc^24*TAGdecmcMSOS*TAGnc*
          Log[mus^2/mcMSmusmc^2])/559893600 + (10163*mcMSmusmc^24*TAGnc*
          Log[musmc^2/mcMSmusmc^2])/11664450 + Log[mcMSmusmc/mbMSmus]*
         ((43809980868157642153*mcMSmusmc^24*TAGnc)/20069484527233911369600 - 
          (1620816161*mcMSmusmc^24*Pi^2*TAGnc)/822167470080 + 
          (281941*mcMSmusmc^24*TAGnc^2)/629880300 + 
          ((11*mcMSmusmc^24*TAGnc)/1656 + (11*mcMSmusmc^24*TAGnc^2)/86940)*
           Log[mus^2/mbMSmus^2] - (11*mcMSmusmc^24*TAGdecmcMSOS*TAGnc*
            Log[mus^2/mcMSmusmc^2])/86940 - (22*mcMSmusmc^24*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/2415) - (2710689767*mcMSmusmc^24*TAGnc*
          Zeta[3])/1644334940160)/mbMSmus^23 + 
      ((2749*mcMSmusmc^4*TAGnc)/2592 - (2*mcMSmusmc^2*mufac^2*TAGnc)/3 + 
        (31*mcMSmusmc^4*Pi^2*TAGnc)/81 + (271*mcMSmusmc^4*Pi^4*TAGnc)/19440 - 
        (235*mcMSmusmc^4*TAGnc^2)/3888 - (13*mcMSmusmc^4*Pi^2*TAGnc^2)/324 + 
        (2*mcMSmusmc^4*Pi^2*TAGnc*Log[2]^2)/81 + 
        (5*mcMSmusmc^4*TAGnc*Log[2]^4)/162 + ((mcMSmusmc^4*TAGnc)/3 - 
          (2*mcMSmusmc^4*TAGnc^2)/27)*Log[mcMSmusmc/mbMSmus]^3 + 
        Log[2]*((5*mcMSmusmc^4*Pi^2*TAGnc)/144 - (mcMSmusmc^4*Pi^2*TAGnc*
            Log[mcMSmusmc/mbMSmus])/9) + ((-1067*mcMSmusmc^4*TAGnc)/432 - 
          (5*mcMSmusmc^4*Pi^2*TAGnc)/36 + (151*mcMSmusmc^4*TAGnc^2)/648 + 
          (mcMSmusmc^4*Pi^2*TAGnc^2)/54)*Log[mus^2/mbMSmus^2] + 
        ((-151*mcMSmusmc^4*TAGdecmcMSOS*TAGnc)/648 - 
          (mcMSmusmc^4*Pi^2*TAGdecmcMSOS*TAGnc)/54)*Log[mus^2/mcMSmusmc^2] + 
        ((-56*mcMSmusmc^4*TAGnc)/27 - (2*mcMSmusmc^4*Pi^2*TAGnc)/9)*
         Log[musmc^2/mcMSmusmc^2] + Log[mcMSmusmc/mbMSmus]^2*
         ((-67*mcMSmusmc^4*TAGnc)/18 + (mcMSmusmc^4*Pi^2*TAGnc)/4 + 
          (13*mcMSmusmc^4*TAGnc^2)/54 + ((-5*mcMSmusmc^4*TAGnc)/6 + 
            (mcMSmusmc^4*TAGnc^2)/9)*Log[mus^2/mbMSmus^2] - 
          (mcMSmusmc^4*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/9 - 
          (4*mcMSmusmc^4*TAGnc*Log[musmc^2/mcMSmusmc^2])/3) + 
        (20*mcMSmusmc^4*TAGnc*PolyLog[4, 1/2])/27 - 
        (2309*mcMSmusmc^4*TAGnc*Zeta[3])/864 + 
        NLxMSOS*((-1423*mcMSmusmc^4*TAGnc)/3888 - (13*mcMSmusmc^4*Pi^2*TAGnc)/
           324 - (2*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mbMSmus]^3)/27 + 
          ((151*mcMSmusmc^4*TAGnc)/648 + (mcMSmusmc^4*Pi^2*TAGnc)/54)*
           Log[mus^2/mbMSmus^2] + Log[mcMSmusmc/mbMSmus]*
           (-(mcMSmusmc^4*TAGnc)/12 + (mcMSmusmc^4*Pi^2*TAGnc)/27 - 
            (13*mcMSmusmc^4*TAGnc*Log[mus^2/mbMSmus^2])/54) + 
          Log[mcMSmusmc/mbMSmus]^2*((13*mcMSmusmc^4*TAGnc)/54 + 
            (mcMSmusmc^4*TAGnc*Log[mus^2/mbMSmus^2])/9) + 
          (mcMSmusmc^4*TAGnc*Zeta[3])/3) + Log[mcMSmusmc/mbMSmus]*
         ((2713*mcMSmusmc^4*TAGnc)/648 - (217*mcMSmusmc^4*Pi^2*TAGnc)/432 - 
          (5*mcMSmusmc^4*TAGnc^2)/12 + (mcMSmusmc^4*Pi^2*TAGnc^2)/27 + 
          ((89*mcMSmusmc^4*TAGnc)/36 - (13*mcMSmusmc^4*TAGnc^2)/54)*
           Log[mus^2/mbMSmus^2] + (13*mcMSmusmc^4*TAGdecmcMSOS*TAGnc*
            Log[mus^2/mcMSmusmc^2])/54 + (20*mcMSmusmc^4*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/9 + (14*mcMSmusmc^4*TAGnc*Zeta[3])/9))/
       mbMSmus^3 + mbMSmus*(8462917/93312 + (652841*Pi^2)/38880 - 
        (695*Pi^4)/7776 - (11*TAGdecmcMSOS)/54 - (231847*TAGnc)/23328 - 
        (991*Pi^2*TAGnc)/648 + (61*Pi^4*TAGnc)/1944 + (2353*TAGnc^2)/23328 + 
        (13*Pi^2*TAGnc^2)/324 + ((-22*Pi^2)/81 + (2*Pi^2*TAGnc)/81)*
         Log[2]^2 + (-55/162 + TAGnc/81)*Log[2]^4 + 
        (1861/432 - (43*TAGnc)/108 + TAGnc^2/108)*Log[mus^2/mbMSmus^2]^3 + 
        (TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]^2)/27 + 
        Log[2]*((-575*Pi^2)/162 - (11*Pi^2*TAGnc)/81 + 
          ((13*Pi^2)/18 - (Pi^2*TAGnc)/27)*Log[mus^2/mbMSmus^2] + 
          (Pi^2*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/27) + 
        Log[mus^2/mbMSmus^2]^2*(21715/864 - (385*TAGnc)/144 + 
          (13*TAGnc^2)/216 + ((47*TAGdecmcMSOS)/72 - (TAGdecmcMSOS*TAGnc)/36)*
           Log[mus^2/mcMSmusmc^2]) - (4*TAGdecmcMSOS*
          Log[musmc^2/mcMSmusmc^2])/9 - (220*PolyLog[4, 1/2])/27 + 
        (8*TAGnc*PolyLog[4, 1/2])/27 + NLxMSOS^2*(2353/23328 + 
          (13*Pi^2)/324 + (89/648 + Pi^2/54)*Log[mus^2/mbMSmus^2] + 
          (13*Log[mus^2/mbMSmus^2]^2)/216 + Log[mus^2/mbMSmus^2]^3/108 + 
          (7*Zeta[3])/54) + (58*Zeta[3])/27 - (1439*Pi^2*Zeta[3])/432 - 
        (241*TAGnc*Zeta[3])/72 + (7*TAGnc^2*Zeta[3])/54 + 
        Log[mus^2/mcMSmusmc^2]*((1225*TAGdecmcMSOS)/288 + 
          (Pi^2*TAGdecmcMSOS)/9 - (71*TAGdecmcMSOS*TAGnc)/432 - 
          (Pi^2*TAGdecmcMSOS*TAGnc)/54 - (TAGdecmcMSOS*Zeta[3])/18) + 
        Log[mus^2/mbMSmus^2]*(93391/1296 + (13*Pi^2)/6 - (11*TAGdecmcMSOS)/
           72 - (5171*TAGnc)/648 - (17*Pi^2*TAGnc)/36 + (89*TAGnc^2)/648 + 
          (Pi^2*TAGnc^2)/54 + ((85*TAGdecmcMSOS)/27 - (13*TAGdecmcMSOS*TAGnc)/
             108)*Log[mus^2/mcMSmusmc^2] + 
          (TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]^2)/36 - 
          (TAGdecmcMSOS*Log[musmc^2/mcMSmusmc^2])/3 - (23*Zeta[3])/12 - 
          (7*TAGnc*Zeta[3])/9) + NLxMSOS*(-231847/23328 - (991*Pi^2)/648 + 
          (61*Pi^4)/1944 + (2353*TAGnc)/11664 + (13*Pi^2*TAGnc)/162 + 
          (2*Pi^2*Log[2]^2)/81 + Log[2]^4/81 + (-43/108 + TAGnc/54)*
           Log[mus^2/mbMSmus^2]^3 + Log[2]*((-11*Pi^2)/81 - 
            (Pi^2*Log[mus^2/mbMSmus^2])/27) + ((-71*TAGdecmcMSOS)/432 - 
            (Pi^2*TAGdecmcMSOS)/54)*Log[mus^2/mcMSmusmc^2] + 
          Log[mus^2/mbMSmus^2]^2*(-385/144 + (13*TAGnc)/108 - 
            (TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/36) + (8*PolyLog[4, 1/2])/
           27 + Log[mus^2/mbMSmus^2]*(-5171/648 - (17*Pi^2)/36 + 
            (89*TAGnc)/324 + (Pi^2*TAGnc)/27 - (13*TAGdecmcMSOS*
              Log[mus^2/mcMSmusmc^2])/108 - (7*Zeta[3])/9) - 
          (241*Zeta[3])/72 + (7*TAGnc*Zeta[3])/27) + (1975*Zeta[5])/216) + 
      ((-70*mcMSmusmc^2*TAGnc)/9 - (13*mcMSmusmc^2*Pi^2*TAGnc)/12 + 
        (2*mcMSmusmc^2*TAGnc^2)/9 + (-4*mcMSmusmc^2*TAGnc - 
          (3*mcMSmusmc^2*Pi^2*TAGnc)/2 + (mcMSmusmc^2*Pi^4*TAGnc)/12)*
         Log[mcMSmusmc/mbMSmus] + ((-9*mcMSmusmc^2*TAGnc)/2 + 
          (mcMSmusmc^2*TAGnc^2)/3)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((2*mcMSmusmc^2*TAGnc)/9 + (mcMSmusmc^2*TAGnc*
            Log[mus^2/mbMSmus^2])/3) - (mcMSmusmc^2*TAGdecmcMSOS*TAGnc*
          Log[mus^2/mcMSmusmc^2])/3 - 2*mcMSmusmc^2*TAGnc*
         Log[musmc^2/mcMSmusmc^2] - (11*mcMSmusmc^2*TAGnc*Zeta[3])/2 + 
        (3*mcMSmusmc^2*Pi^2*TAGnc*Zeta[3])/4 + 
        mufac^2*(-22055/108 + (13805*NLxOSKIN)/648 - (209*NLxOSKIN^2)/486 + 
          (437*Pi^2)/36 - (23*NLxOSKIN*Pi^2)/27 + (NLxOSKIN^2*Pi^2)/81 - 
          Pi^4/4 - (71*TAGnc)/216 - (Pi^2*TAGnc)/27 + (2*Pi^2*Log[2])/27 + 
          (-121/6 + (22*NLxOSKIN)/9 - (2*NLxOSKIN^2)/27)*Log[(2*mufac)/mus]^
            2 + (23/36 - TAGnc/18)*Log[mus^2/mbMSmus^2]^2 + 
          Log[(2*mufac)/mus]*(689/6 - (691*NLxOSKIN)/54 + (26*NLxOSKIN^2)/
             81 - (11*Pi^2)/3 + (2*NLxOSKIN*Pi^2)/9 + 
            (-11/3 + (2*NLxOSKIN)/9)*Log[mus^2/mbMSmus^2]) + 
          NLxMSOS*(-71/216 - Pi^2/27 - (13*Log[mus^2/mbMSmus^2])/54 - 
            Log[mus^2/mbMSmus^2]^2/18) + (4*TAGdecmcMSOS*
            Log[mus^2/mcMSmusmc^2])/27 + Log[mus^2/mbMSmus^2]*
           (1409/108 - (13*NLxOSKIN)/27 - Pi^2/3 - (13*TAGnc)/54 + 
            (TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/9) + (1535*Zeta[3])/36 - 
          (35*NLxOSKIN*Zeta[3])/18) + (5*mcMSmusmc^2*TAGnc*Zeta[5])/2)/
       mbMSmus)


RunDecmMS2mkinmc = mbMSmus + apinl1mus*((-16*mufac)/9 - 
      (2*mufac^2)/(3*mbMSmus) + mbMSmus*(4/3 + Log[mus^2/mbMSmus^2])) + 
    apinl1mus^2*((mcMSmusmc*Pi^2*TAGnc)/6 + (mcMSmusmc^3*Pi^2*TAGnc)/
       (6*mbMSmus^2) + ((19*mcMSmusmc^6*TAGnc)/225 - 
        (4*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mbMSmus])/45)/mbMSmus^5 + 
      ((463*mcMSmusmc^8*TAGnc)/39200 - (3*mcMSmusmc^8*TAGnc*
          Log[mcMSmusmc/mbMSmus])/140)/mbMSmus^7 + 
      ((997*mcMSmusmc^10*TAGnc)/297675 - (8*mcMSmusmc^10*TAGnc*
          Log[mcMSmusmc/mbMSmus])/945)/mbMSmus^9 + 
      ((1229*mcMSmusmc^12*TAGnc)/940896 - (5*mcMSmusmc^12*TAGnc*
          Log[mcMSmusmc/mbMSmus])/1188)/mbMSmus^11 + 
      ((15371*mcMSmusmc^14*TAGnc)/25050025 - (12*mcMSmusmc^14*TAGnc*
          Log[mcMSmusmc/mbMSmus])/5005)/mbMSmus^13 + 
      ((529*mcMSmusmc^16*TAGnc)/1622400 - (7*mcMSmusmc^16*TAGnc*
          Log[mcMSmusmc/mbMSmus])/4680)/mbMSmus^15 + 
      ((16277*mcMSmusmc^18*TAGnc)/86028075 - (16*mcMSmusmc^18*TAGnc*
          Log[mcMSmusmc/mbMSmus])/16065)/mbMSmus^17 + 
      ((39163*mcMSmusmc^20*TAGnc)/333852800 - 
        (9*mcMSmusmc^20*TAGnc*Log[mcMSmusmc/mbMSmus])/12920)/mbMSmus^19 + 
      ((39833*mcMSmusmc^22*TAGnc)/520109667 - 
        (20*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mbMSmus])/39501)/mbMSmus^21 + 
      ((9727*mcMSmusmc^24*TAGnc)/186631200 - (11*mcMSmusmc^24*TAGnc*
          Log[mcMSmusmc/mbMSmus])/28980)/mbMSmus^23 + 
      ((-151*mcMSmusmc^4*TAGnc)/216 - (mcMSmusmc^4*Pi^2*TAGnc)/18 + 
        (13*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mbMSmus])/18 - 
        (mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mbMSmus]^2)/3)/mbMSmus^3 + 
      mufac*(-860/27 + (128*NLxOSKIN)/81 + (8*Pi^2)/9 + 
        (88/9 - (16*NLxOSKIN)/27)*Log[(2*mufac)/mus] + 
        (8*TAGkinmc*Log[mus^2/mcMSmusmc^2])/27) + 
      (-(mcMSmusmc^2*TAGnc) + mufac^2*(-83/9 + (13*NLxOSKIN)/27 + Pi^2/3 + 
          (11/3 - (2*NLxOSKIN)/9)*Log[(2*mufac)/mus] + 
          (2*Log[mus^2/mbMSmus^2])/3 + (TAGkinmc*Log[mus^2/mcMSmusmc^2])/9))/
       mbMSmus + mbMSmus*(307/32 + Pi^2/3 - (71*TAGnc)/144 - 
        (Pi^2*TAGnc)/18 + (Pi^2*Log[2])/9 + (509/72 - (13*TAGnc)/36)*
         Log[mus^2/mbMSmus^2] + (47/24 - TAGnc/12)*Log[mus^2/mbMSmus^2]^2 + 
        NLxMSOS*(-71/144 - Pi^2/18 - (13*Log[mus^2/mbMSmus^2])/36 - 
          Log[mus^2/mbMSmus^2]^2/12) - Zeta[3]/6)) + 
    apinl1mus^3*((5971*mcMSmusmc*Pi^2*TAGnc)/486 + (13*mcMSmusmc*Pi^3*TAGnc)/
       162 - (4*mcMSmusmc*Pi^2*TAGnc^2)/135 - 
      (1199*mcMSmusmc*Pi^2*TAGnc*Log[2])/81 + 
      ((-11*mcMSmusmc*Pi^2*TAGnc)/6 + (mcMSmusmc*Pi^2*TAGnc^2)/9)*
       Log[mcMSmusmc/mbMSmus] + 
      (mufac^2*((9727*mcMSmusmc^24*TAGnc)/279946800 - 
         (11*mcMSmusmc^24*TAGnc*Log[mcMSmusmc/mbMSmus])/43470))/mbMSmus^25 + 
      ((149*mcMSmusmc^5*Pi^2*TAGnc)/1080 + (mcMSmusmc^3*mufac^2*Pi^2*TAGnc)/
         9 - (mcMSmusmc^5*NLxMSOS*Pi^2*TAGnc)/20 + (mcMSmusmc^5*Pi^3*TAGnc)/
         240 - (13*mcMSmusmc^5*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/45)/
       mbMSmus^4 + ((262769*mcMSmusmc^7*Pi^2*TAGnc)/11907000 - 
        (5*mcMSmusmc^7*NLxMSOS*Pi^2*TAGnc)/756 - (25*mcMSmusmc^7*Pi^3*TAGnc)/
         18144 - (29*mcMSmusmc^7*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/525)/
       mbMSmus^6 + ((416029*mcMSmusmc^9*Pi^2*TAGnc)/60011280 - 
        (7*mcMSmusmc^9*NLxMSOS*Pi^2*TAGnc)/3240 - (77*mcMSmusmc^9*Pi^3*TAGnc)/
         103680 - (53*mcMSmusmc^9*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/2205)/
       mbMSmus^8 + ((1055689*mcMSmusmc^11*Pi^2*TAGnc)/345779280 - 
        (3*mcMSmusmc^11*NLxMSOS*Pi^2*TAGnc)/3080 - 
        (17*mcMSmusmc^11*Pi^3*TAGnc)/39424 - (85*mcMSmusmc^11*Pi^2*TAGnc*
          Log[mcMSmusmc/mbMSmus])/6237)/mbMSmus^10 + 
      ((31786481*mcMSmusmc^13*Pi^2*TAGnc)/19677663720 - 
        (11*mcMSmusmc^13*NLxMSOS*Pi^2*TAGnc)/21060 - 
        (1771*mcMSmusmc^13*Pi^3*TAGnc)/6469632 - 
        (125*mcMSmusmc^13*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/14157)/
       mbMSmus^12 + ((17346493*mcMSmusmc^15*Pi^2*TAGnc)/18087549480 - 
        (13*mcMSmusmc^15*NLxMSOS*Pi^2*TAGnc)/41580 - 
        (377*mcMSmusmc^15*Pi^3*TAGnc)/2027520 - 
        (173*mcMSmusmc^15*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/27885)/
       mbMSmus^14 + ((22757641*mcMSmusmc^17*Pi^2*TAGnc)/36923796000 - 
        (5*mcMSmusmc^17*NLxMSOS*Pi^2*TAGnc)/24752 - 
        (1925*mcMSmusmc^17*Pi^3*TAGnc)/14483456 - 
        (229*mcMSmusmc^17*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/49725)/
       mbMSmus^16 + ((86836957*mcMSmusmc^19*Pi^2*TAGnc)/206871887520 - 
        (17*mcMSmusmc^19*NLxMSOS*Pi^2*TAGnc)/123120 - 
        (99671*mcMSmusmc^19*Pi^3*TAGnc)/1008599040 - 
        (293*mcMSmusmc^19*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/82365)/
       mbMSmus^18 + ((846435761*mcMSmusmc^21*Pi^2*TAGnc)/2832319518840 - 
        (19*mcMSmusmc^21*NLxMSOS*Pi^2*TAGnc)/192780 - 
        (127699*mcMSmusmc^21*Pi^3*TAGnc)/1684537344 - 
        (365*mcMSmusmc^21*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/128877)/
       mbMSmus^20 + ((171475369*mcMSmusmc^23*Pi^2*TAGnc)/778168119960 - 
        (7*mcMSmusmc^23*NLxMSOS*Pi^2*TAGnc)/96140 - 
        (81991*mcMSmusmc^23*Pi^3*TAGnc)/1374683136 - 
        (445*mcMSmusmc^23*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/192717)/
       mbMSmus^22 + ((106559417*mcMSmusmc^25*Pi^2*TAGnc)/637438863600 - 
        (23*mcMSmusmc^25*NLxMSOS*Pi^2*TAGnc)/415800 - 
        (5698043*mcMSmusmc^25*Pi^3*TAGnc)/118908518400 - 
        (533*mcMSmusmc^25*Pi^2*TAGnc*Log[mcMSmusmc/mbMSmus])/277725)/
       mbMSmus^24 + ((11*mcMSmusmc*Pi^2*TAGnc)/12 - (mcMSmusmc*Pi^2*TAGnc^2)/
         18)*Log[mus^2/mbMSmus^2] + NLxMSOS*((-7*mcMSmusmc*Pi^2*TAGnc)/27 + 
        (2*mcMSmusmc*Pi^2*TAGnc*Log[2])/9 + (mcMSmusmc*Pi^2*TAGnc*
          Log[mcMSmusmc/mbMSmus])/9 - (mcMSmusmc*Pi^2*TAGnc*
          Log[mus^2/mbMSmus^2])/18) + (mcMSmusmc*Pi^2*TAGnc*
        Log[musmc^2/mcMSmusmc^2])/6 + ((18607*mcMSmusmc^3*Pi^2*TAGnc)/1620 + 
        (mcMSmusmc*mufac^2*Pi^2*TAGnc)/9 + (7*mcMSmusmc^3*Pi^3*TAGnc)/108 - 
        (1199*mcMSmusmc^3*Pi^2*TAGnc*Log[2])/81 + 
        ((-289*mcMSmusmc^3*Pi^2*TAGnc)/162 + (mcMSmusmc^3*Pi^2*TAGnc^2)/9)*
         Log[mcMSmusmc/mbMSmus] + ((7*mcMSmusmc^3*Pi^2*TAGnc)/12 - 
          (mcMSmusmc^3*Pi^2*TAGnc^2)/18)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((-7*mcMSmusmc^3*Pi^2*TAGnc)/54 + 
          (2*mcMSmusmc^3*Pi^2*TAGnc*Log[2])/9 + (mcMSmusmc^3*Pi^2*TAGnc*
            Log[mcMSmusmc/mbMSmus])/9 - (mcMSmusmc^3*Pi^2*TAGnc*
            Log[mus^2/mbMSmus^2])/18) + (mcMSmusmc^3*Pi^2*TAGnc*
          Log[musmc^2/mcMSmusmc^2])/2)/mbMSmus^2 + 
      mufac*(-130867/162 + (20047*NLxOSKIN)/243 - (1292*NLxOSKIN^2)/729 + 
        (1022*Pi^2)/27 - (208*NLxOSKIN*Pi^2)/81 + (8*NLxOSKIN^2*Pi^2)/243 - 
        (2*Pi^4)/3 - (22*TAGkinmc)/81 + (-484/9 + (176*NLxOSKIN)/27 - 
          (16*NLxOSKIN^2)/81)*Log[(2*mufac)/mus]^2 + 
        ((926*TAGkinmc)/81 - (128*NLxOSKIN*TAGkinmc)/243 - 
          (8*Pi^2*TAGkinmc)/27)*Log[mus^2/mcMSmusmc^2] - 
        (4*TAGkinmc*Log[mus^2/mcMSmusmc^2]^2)/81 + Log[(2*mufac)/mus]*
         (10072/27 - (3356*NLxOSKIN)/81 + (256*NLxOSKIN^2)/243 - 
          (88*Pi^2)/9 + (16*NLxOSKIN*Pi^2)/27 + ((-88*TAGkinmc)/27 + 
            (16*NLxOSKIN*TAGkinmc)/81)*Log[mus^2/mcMSmusmc^2]) + 
        (16*TAGkinmc*Log[mus^2/musmc^2])/27 + 114*Zeta[3] - 
        (140*NLxOSKIN*Zeta[3])/27) + ((-3529813*mcMSmusmc^6*TAGnc)/3499200 + 
        (1427*mcMSmusmc^6*Pi^2*TAGnc)/12960 - (3481*mcMSmusmc^6*TAGnc^2)/
         30375 + (4*mcMSmusmc^6*Pi^2*TAGnc^2)/405 - 
        (11*mcMSmusmc^6*Pi^2*TAGnc*Log[2])/648 - 
        (5329*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mbMSmus]^2)/6480 + 
        mufac^2*((-151*mcMSmusmc^4*TAGnc)/324 - (mcMSmusmc^4*Pi^2*TAGnc)/27 + 
          (13*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mbMSmus])/27 - 
          (2*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mbMSmus]^2)/9) + 
        ((59*mcMSmusmc^6*TAGnc)/450 - (19*mcMSmusmc^6*TAGnc^2)/675)*
         Log[mus^2/mbMSmus^2] + NLxMSOS*((2729*mcMSmusmc^6*TAGnc)/30375 + 
          (4*mcMSmusmc^6*Pi^2*TAGnc)/405 - (19*mcMSmusmc^6*TAGnc*
            Log[mus^2/mbMSmus^2])/675 + Log[mcMSmusmc/mbMSmus]*
           ((-2*mcMSmusmc^6*TAGnc)/75 + (4*mcMSmusmc^6*TAGnc*
              Log[mus^2/mbMSmus^2])/135)) + (94*mcMSmusmc^6*TAGnc*
          Log[musmc^2/mcMSmusmc^2])/225 + Log[mcMSmusmc/mbMSmus]*
         ((1481549*mcMSmusmc^6*TAGnc)/1166400 - (113*mcMSmusmc^6*Pi^2*TAGnc)/
           1296 + (8*mcMSmusmc^6*TAGnc^2)/75 + ((-2*mcMSmusmc^6*TAGnc)/45 + 
            (4*mcMSmusmc^6*TAGnc^2)/135)*Log[mus^2/mbMSmus^2] - 
          (8*mcMSmusmc^6*TAGnc*Log[musmc^2/mcMSmusmc^2])/15) - 
        (29*mcMSmusmc^6*TAGnc*Zeta[3])/648)/mbMSmus^5 + 
      ((-33561121*mcMSmusmc^8*TAGnc)/1185408000 + 
        (1789*mcMSmusmc^8*Pi^2*TAGnc)/90720 - (241*mcMSmusmc^8*TAGnc^2)/
         7056 - (7*mcMSmusmc^8*Pi^2*TAGnc*Log[2])/1536 + 
        ((-9019*mcMSmusmc^8*TAGnc)/30240 - (mcMSmusmc^8*TAGnc^2)/70)*
         Log[mcMSmusmc/mbMSmus]^2 + mufac^2*((38*mcMSmusmc^6*TAGnc)/675 - 
          (8*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mbMSmus])/135) + 
        ((291*mcMSmusmc^8*TAGnc)/78400 - (463*mcMSmusmc^8*TAGnc^2)/117600)*
         Log[mus^2/mbMSmus^2] + NLxMSOS*((174787*mcMSmusmc^8*TAGnc)/
           12348000 + (mcMSmusmc^8*Pi^2*TAGnc)/420 - 
          (463*mcMSmusmc^8*TAGnc*Log[mus^2/mbMSmus^2])/117600 + 
          Log[mcMSmusmc/mbMSmus]*((-37*mcMSmusmc^8*TAGnc)/14700 + 
            (mcMSmusmc^8*TAGnc*Log[mus^2/mbMSmus^2])/140)) + 
        (179*mcMSmusmc^8*TAGnc*Log[musmc^2/mcMSmusmc^2])/2450 + 
        Log[mcMSmusmc/mbMSmus]*((2308303*mcMSmusmc^8*TAGnc)/12700800 - 
          (139*mcMSmusmc^8*Pi^2*TAGnc)/4608 + (193*mcMSmusmc^8*TAGnc^2)/
           3675 + ((9*mcMSmusmc^8*TAGnc)/280 + (mcMSmusmc^8*TAGnc^2)/140)*
           Log[mus^2/mbMSmus^2] - (6*mcMSmusmc^8*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/35) - (173*mcMSmusmc^8*TAGnc*Zeta[3])/
         9216)/mbMSmus^7 + ((47157278431*mcMSmusmc^10*TAGnc)/5184974592000 + 
        (359801*mcMSmusmc^10*Pi^2*TAGnc)/48988800 - 
        (1778*mcMSmusmc^10*TAGnc^2)/273375 - (17*mcMSmusmc^10*Pi^2*TAGnc*
          Log[2])/8640 + ((-2585963*mcMSmusmc^10*TAGnc)/16329600 - 
          (16*mcMSmusmc^10*TAGnc^2)/2835)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((463*mcMSmusmc^8*TAGnc)/58800 - 
          (mcMSmusmc^8*TAGnc*Log[mcMSmusmc/mbMSmus])/70) + 
        ((-277*mcMSmusmc^10*TAGnc)/85050 - (997*mcMSmusmc^10*TAGnc^2)/893025)*
         Log[mus^2/mbMSmus^2] + NLxMSOS*((2816347*mcMSmusmc^10*TAGnc)/
           562605750 + (8*mcMSmusmc^10*Pi^2*TAGnc)/8505 - 
          (997*mcMSmusmc^10*TAGnc*Log[mus^2/mbMSmus^2])/893025 + 
          Log[mcMSmusmc/mbMSmus]*((-398*mcMSmusmc^10*TAGnc)/893025 + 
            (8*mcMSmusmc^10*TAGnc*Log[mus^2/mbMSmus^2])/2835)) + 
        (298*mcMSmusmc^10*TAGnc*Log[musmc^2/mcMSmusmc^2])/11907 + 
        Log[mcMSmusmc/mbMSmus]*((2370226609*mcMSmusmc^10*TAGnc)/41150592000 - 
          (17*mcMSmusmc^10*Pi^2*TAGnc)/1080 + (14048*mcMSmusmc^10*TAGnc^2)/
           893025 + ((4*mcMSmusmc^10*TAGnc)/135 + (8*mcMSmusmc^10*TAGnc^2)/
             2835)*Log[mus^2/mbMSmus^2] - (16*mcMSmusmc^10*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/189) - (187*mcMSmusmc^10*TAGnc*Zeta[3])/
         17280)/mbMSmus^9 + ((237903978448457*mcMSmusmc^12*TAGnc)/
         24647147078400000 + (84041429*mcMSmusmc^12*Pi^2*TAGnc)/22992076800 - 
        (3989*mcMSmusmc^12*TAGnc^2)/1881792 - (175*mcMSmusmc^12*Pi^2*TAGnc*
          Log[2])/165888 + ((-14841023*mcMSmusmc^12*TAGnc)/149688000 - 
          (5*mcMSmusmc^12*TAGnc^2)/1782)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((1994*mcMSmusmc^10*TAGnc)/893025 - 
          (16*mcMSmusmc^10*TAGnc*Log[mcMSmusmc/mbMSmus])/2835) + 
        ((-509*mcMSmusmc^12*TAGnc)/171072 - (1229*mcMSmusmc^12*TAGnc^2)/
           2822688)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((326802499*mcMSmusmc^12*TAGnc)/136928594880 + 
          (5*mcMSmusmc^12*Pi^2*TAGnc)/10692 - (1229*mcMSmusmc^12*TAGnc*
            Log[mus^2/mbMSmus^2])/2822688 + Log[mcMSmusmc/mbMSmus]*
           ((-391*mcMSmusmc^12*TAGnc)/4939704 + (5*mcMSmusmc^12*TAGnc*
              Log[mus^2/mbMSmus^2])/3564)) + (899*mcMSmusmc^12*TAGnc*
          Log[musmc^2/mcMSmusmc^2])/78408 + Log[mcMSmusmc/mbMSmus]*
         ((160430424613*mcMSmusmc^12*TAGnc)/6224027040000 - 
          (40553*mcMSmusmc^12*Pi^2*TAGnc)/4147200 + 
          (9545*mcMSmusmc^12*TAGnc^2)/1411344 + ((5*mcMSmusmc^12*TAGnc)/216 + 
            (5*mcMSmusmc^12*TAGnc^2)/3564)*Log[mus^2/mbMSmus^2] - 
          (5*mcMSmusmc^12*TAGnc*Log[musmc^2/mcMSmusmc^2])/99) - 
        (59231*mcMSmusmc^12*TAGnc*Zeta[3])/8294400)/mbMSmus^11 + 
      ((448931017732868509*mcMSmusmc^14*TAGnc)/58963096098466560000 + 
        (82285201*mcMSmusmc^14*Pi^2*TAGnc)/38745907200 - 
        (233201*mcMSmusmc^14*TAGnc^2)/260851500 - 
        (23*mcMSmusmc^14*Pi^2*TAGnc*Log[2])/35840 + 
        ((-4128997*mcMSmusmc^14*TAGnc)/60540480 - (8*mcMSmusmc^14*TAGnc^2)/
           5005)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((1229*mcMSmusmc^12*TAGnc)/1411344 - 
          (5*mcMSmusmc^12*TAGnc*Log[mcMSmusmc/mbMSmus])/1782) + 
        ((-22089*mcMSmusmc^14*TAGnc)/10020010 - (15371*mcMSmusmc^14*TAGnc^2)/
           75150075)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((108352091581*mcMSmusmc^14*TAGnc)/81243243081000 + 
          (4*mcMSmusmc^14*Pi^2*TAGnc)/15015 - (15371*mcMSmusmc^14*TAGnc*
            Log[mus^2/mbMSmus^2])/75150075 + Log[mcMSmusmc/mbMSmus]*
           ((1153*mcMSmusmc^14*TAGnc)/225450225 + (4*mcMSmusmc^14*TAGnc*
              Log[mus^2/mbMSmus^2])/5005)) + (3166*mcMSmusmc^14*TAGnc*
          Log[musmc^2/mcMSmusmc^2])/511225 + Log[mcMSmusmc/mbMSmus]*
         ((570692369939*mcMSmusmc^14*TAGnc)/40905688824000 - 
          (2161*mcMSmusmc^14*Pi^2*TAGnc)/322560 + 
          (263546*mcMSmusmc^14*TAGnc^2)/75150075 + 
          ((18*mcMSmusmc^14*TAGnc)/1001 + (4*mcMSmusmc^14*TAGnc^2)/5005)*
           Log[mus^2/mbMSmus^2] - (24*mcMSmusmc^14*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/715) - (3287*mcMSmusmc^14*TAGnc*
          Zeta[3])/645120)/mbMSmus^13 + 
      ((558828635339341201127*mcMSmusmc^16*TAGnc)/95095681387606867968000 + 
        (58806560951*mcMSmusmc^16*Pi^2*TAGnc)/43261891706880 - 
        (1373*mcMSmusmc^16*TAGnc^2)/3110400 - (1001*mcMSmusmc^16*Pi^2*TAGnc*
          Log[2])/2359296 + ((-1141880431*mcMSmusmc^16*TAGnc)/22884301440 - 
          (7*mcMSmusmc^16*TAGnc^2)/7020)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((30742*mcMSmusmc^14*TAGnc)/75150075 - 
          (8*mcMSmusmc^14*TAGnc*Log[mcMSmusmc/mbMSmus])/5005) + 
        ((-15593*mcMSmusmc^16*TAGnc)/9734400 - (529*mcMSmusmc^16*TAGnc^2)/
           4867200)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((214558103603*mcMSmusmc^16*TAGnc)/260460712512000 + 
          (7*mcMSmusmc^16*Pi^2*TAGnc)/42120 - (529*mcMSmusmc^16*TAGnc*
            Log[mus^2/mbMSmus^2])/4867200 + Log[mcMSmusmc/mbMSmus]*
           ((355*mcMSmusmc^16*TAGnc)/14455584 + (7*mcMSmusmc^16*TAGnc*
              Log[mus^2/mbMSmus^2])/14040)) + 
        (283*mcMSmusmc^16*TAGnc*Log[musmc^2/mcMSmusmc^2])/76050 + 
        Log[mcMSmusmc/mbMSmus]*((10019690708971*mcMSmusmc^16*TAGnc)/
           1178083838131200 - (565351*mcMSmusmc^16*Pi^2*TAGnc)/115605504 + 
          (7469*mcMSmusmc^16*TAGnc^2)/3650400 + 
          ((133*mcMSmusmc^16*TAGnc)/9360 + (7*mcMSmusmc^16*TAGnc^2)/14040)*
           Log[mus^2/mbMSmus^2] - (14*mcMSmusmc^16*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/585) - (885457*mcMSmusmc^16*TAGnc*
          Zeta[3])/231211008)/mbMSmus^15 + 
      ((230500151130233896919629*mcMSmusmc^18*TAGnc)/
         50057687427569200963584000 + (26463251891*mcMSmusmc^18*Pi^2*TAGnc)/
         28454994247680 - (157882*mcMSmusmc^18*TAGnc^2)/650372247 - 
        (4147*mcMSmusmc^18*Pi^2*TAGnc*Log[2])/13934592 + 
        ((-3391310411*mcMSmusmc^18*TAGnc)/88921857024 - 
          (32*mcMSmusmc^18*TAGnc^2)/48195)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((529*mcMSmusmc^16*TAGnc)/2433600 - 
          (7*mcMSmusmc^16*TAGnc*Log[mcMSmusmc/mbMSmus])/7020) + 
        ((-203011*mcMSmusmc^18*TAGnc)/172056150 - 
          (16277*mcMSmusmc^18*TAGnc^2)/258084225)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((20555048260909*mcMSmusmc^18*TAGnc)/37681809223558500 + 
          (16*mcMSmusmc^18*Pi^2*TAGnc)/144585 - (16277*mcMSmusmc^18*TAGnc*
            Log[mus^2/mbMSmus^2])/258084225 + Log[mcMSmusmc/mbMSmus]*
           ((197062*mcMSmusmc^18*TAGnc)/7381208835 + 
            (16*mcMSmusmc^18*TAGnc*Log[mus^2/mbMSmus^2])/48195)) + 
        (7678*mcMSmusmc^18*TAGnc*Log[musmc^2/mcMSmusmc^2])/3186225 + 
        Log[mcMSmusmc/mbMSmus]*((91862575617799009*mcMSmusmc^18*TAGnc)/
           16342379002556006400 - (52013*mcMSmusmc^18*Pi^2*TAGnc)/13934592 + 
          (334304*mcMSmusmc^18*TAGnc^2)/258084225 + 
          ((184*mcMSmusmc^18*TAGnc)/16065 + (16*mcMSmusmc^18*TAGnc^2)/48195)*
           Log[mus^2/mbMSmus^2] - (32*mcMSmusmc^18*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/1785) - (83291*mcMSmusmc^18*TAGnc*
          Zeta[3])/27869184)/mbMSmus^17 + 
      ((200744102824691620675013329*mcMSmusmc^20*TAGnc)/
         54499313978682087207813120000 + (9007367733163*mcMSmusmc^20*Pi^2*
          TAGnc)/13481015456563200 - (367909*mcMSmusmc^20*TAGnc^2)/
         2547216000 - (143*mcMSmusmc^20*Pi^2*TAGnc*Log[2])/655360 + 
        ((-45436526101*mcMSmusmc^20*TAGnc)/1508495788800 - 
          (3*mcMSmusmc^20*TAGnc^2)/6460)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((32554*mcMSmusmc^18*TAGnc)/258084225 - 
          (32*mcMSmusmc^18*TAGnc*Log[mcMSmusmc/mbMSmus])/48195) + 
        ((-592281*mcMSmusmc^20*TAGnc)/667705600 - 
          (39163*mcMSmusmc^20*TAGnc^2)/1001558400)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((22189567531163017*mcMSmusmc^20*TAGnc)/
           58347124817357376000 + (mcMSmusmc^20*Pi^2*TAGnc)/12920 - 
          (39163*mcMSmusmc^20*TAGnc*Log[mus^2/mbMSmus^2])/1001558400 + 
          Log[mcMSmusmc/mbMSmus]*((1211963*mcMSmusmc^20*TAGnc)/50127997920 + 
            (3*mcMSmusmc^20*TAGnc*Log[mus^2/mbMSmus^2])/12920)) + 
        (5507*mcMSmusmc^20*TAGnc*Log[musmc^2/mcMSmusmc^2])/3338528 + 
        Log[mcMSmusmc/mbMSmus]*((691405824368106269*mcMSmusmc^20*TAGnc)/
           175583298211985664000 - (156353*mcMSmusmc^20*Pi^2*TAGnc)/
           53084160 + (1017529*mcMSmusmc^20*TAGnc^2)/1168484800 + 
          ((243*mcMSmusmc^20*TAGnc)/25840 + (3*mcMSmusmc^20*TAGnc^2)/12920)*
           Log[mus^2/mbMSmus^2] - (9*mcMSmusmc^20*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/646) - (254791*mcMSmusmc^20*TAGnc*
          Zeta[3])/106168320)/mbMSmus^19 + 
      ((32427296704361728583901633193*mcMSmusmc^22*TAGnc)/
         10790864167779053267146997760000 + (174200135864459*mcMSmusmc^22*
          Pi^2*TAGnc)/349543472195174400 - (624853*mcMSmusmc^22*TAGnc^2)/
         6846429744 - (38675*mcMSmusmc^22*Pi^2*TAGnc*Log[2])/233570304 + 
        ((-3036958280711*mcMSmusmc^22*TAGnc)/124450902576000 - 
          (40*mcMSmusmc^22*TAGnc^2)/118503)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((39163*mcMSmusmc^20*TAGnc)/500779200 - 
          (3*mcMSmusmc^20*TAGnc*Log[mcMSmusmc/mbMSmus])/6460) + 
        ((-708143*mcMSmusmc^22*TAGnc)/1040219334 - 
          (39833*mcMSmusmc^22*TAGnc^2)/1560329001)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((1640519393726677*mcMSmusmc^22*TAGnc)/5946258455651273760 + 
          (20*mcMSmusmc^22*Pi^2*TAGnc)/355509 - (39833*mcMSmusmc^22*TAGnc*
            Log[mus^2/mbMSmus^2])/1560329001 + Log[mcMSmusmc/mbMSmus]*
           ((14290513*mcMSmusmc^22*TAGnc)/689665418442 + 
            (20*mcMSmusmc^22*TAGnc*Log[mus^2/mbMSmus^2])/118503)) + 
        (5066*mcMSmusmc^22*TAGnc*Log[musmc^2/mcMSmusmc^2])/4298427 + 
        Log[mcMSmusmc/mbMSmus]*((20875025059584870877*mcMSmusmc^22*TAGnc)/
           7242811051244408640000 - (13926181*mcMSmusmc^22*Pi^2*TAGnc)/
           5839257600 + (318845*mcMSmusmc^22*TAGnc^2)/520109667 + 
          ((310*mcMSmusmc^22*TAGnc)/39501 + (20*mcMSmusmc^22*TAGnc^2)/118503)*
           Log[mus^2/mbMSmus^2] - (40*mcMSmusmc^22*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/3591) - (23017987*mcMSmusmc^22*TAGnc*
          Zeta[3])/11678515200)/mbMSmus^21 + 
      ((5912372795358292426094681002539449*mcMSmusmc^24*TAGnc)/
         2369996943791663839368455777157120000 + 
        (240817793781176357*mcMSmusmc^24*Pi^2*TAGnc)/628867544642695987200 - 
        (279881*mcMSmusmc^24*TAGnc^2)/4627692000 - 
        (877591*mcMSmusmc^24*Pi^2*TAGnc*Log[2])/6794772480 + 
        ((-941278536953*mcMSmusmc^24*TAGnc)/46646042002560 - 
          (11*mcMSmusmc^24*TAGnc^2)/43470)*Log[mcMSmusmc/mbMSmus]^2 + 
        mufac^2*((79666*mcMSmusmc^22*TAGnc)/1560329001 - 
          (40*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mbMSmus])/118503) + 
        ((-631*mcMSmusmc^24*TAGnc)/1184960 - (9727*mcMSmusmc^24*TAGnc^2)/
           559893600)*Log[mus^2/mbMSmus^2] + 
        NLxMSOS*((7801530877413386647*mcMSmusmc^24*TAGnc)/
           37763267488425775296000 + (11*mcMSmusmc^24*Pi^2*TAGnc)/260820 - 
          (9727*mcMSmusmc^24*TAGnc*Log[mus^2/mbMSmus^2])/559893600 + 
          Log[mcMSmusmc/mbMSmus]*((73801799*mcMSmusmc^24*TAGnc)/
             4231787807520 + (11*mcMSmusmc^24*TAGnc*Log[mus^2/mbMSmus^2])/
             86940)) + (10163*mcMSmusmc^24*TAGnc*Log[musmc^2/mcMSmusmc^2])/
         11664450 + Log[mcMSmusmc/mbMSmus]*
         ((43809980868157642153*mcMSmusmc^24*TAGnc)/20069484527233911369600 - 
          (1620816161*mcMSmusmc^24*Pi^2*TAGnc)/822167470080 + 
          (281941*mcMSmusmc^24*TAGnc^2)/629880300 + 
          ((11*mcMSmusmc^24*TAGnc)/1656 + (11*mcMSmusmc^24*TAGnc^2)/86940)*
           Log[mus^2/mbMSmus^2] - (22*mcMSmusmc^24*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/2415) - (2710689767*mcMSmusmc^24*TAGnc*
          Zeta[3])/1644334940160)/mbMSmus^23 + 
      ((2749*mcMSmusmc^4*TAGnc)/2592 - (2*mcMSmusmc^2*mufac^2*TAGnc)/3 + 
        (31*mcMSmusmc^4*Pi^2*TAGnc)/81 + (271*mcMSmusmc^4*Pi^4*TAGnc)/19440 - 
        (235*mcMSmusmc^4*TAGnc^2)/3888 - (13*mcMSmusmc^4*Pi^2*TAGnc^2)/324 + 
        (2*mcMSmusmc^4*Pi^2*TAGnc*Log[2]^2)/81 + 
        (5*mcMSmusmc^4*TAGnc*Log[2]^4)/162 + ((mcMSmusmc^4*TAGnc)/3 - 
          (2*mcMSmusmc^4*TAGnc^2)/27)*Log[mcMSmusmc/mbMSmus]^3 + 
        Log[2]*((5*mcMSmusmc^4*Pi^2*TAGnc)/144 - (mcMSmusmc^4*Pi^2*TAGnc*
            Log[mcMSmusmc/mbMSmus])/9) + ((-1067*mcMSmusmc^4*TAGnc)/432 - 
          (5*mcMSmusmc^4*Pi^2*TAGnc)/36 + (151*mcMSmusmc^4*TAGnc^2)/648 + 
          (mcMSmusmc^4*Pi^2*TAGnc^2)/54)*Log[mus^2/mbMSmus^2] + 
        ((-56*mcMSmusmc^4*TAGnc)/27 - (2*mcMSmusmc^4*Pi^2*TAGnc)/9)*
         Log[musmc^2/mcMSmusmc^2] + Log[mcMSmusmc/mbMSmus]^2*
         ((-67*mcMSmusmc^4*TAGnc)/18 + (mcMSmusmc^4*Pi^2*TAGnc)/4 + 
          (13*mcMSmusmc^4*TAGnc^2)/54 + ((-5*mcMSmusmc^4*TAGnc)/6 + 
            (mcMSmusmc^4*TAGnc^2)/9)*Log[mus^2/mbMSmus^2] - 
          (4*mcMSmusmc^4*TAGnc*Log[musmc^2/mcMSmusmc^2])/3) + 
        (20*mcMSmusmc^4*TAGnc*PolyLog[4, 1/2])/27 - 
        (2309*mcMSmusmc^4*TAGnc*Zeta[3])/864 + 
        NLxMSOS*((-1423*mcMSmusmc^4*TAGnc)/3888 - (13*mcMSmusmc^4*Pi^2*TAGnc)/
           324 - (2*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mbMSmus]^3)/27 + 
          ((151*mcMSmusmc^4*TAGnc)/648 + (mcMSmusmc^4*Pi^2*TAGnc)/54)*
           Log[mus^2/mbMSmus^2] + Log[mcMSmusmc/mbMSmus]*
           (-(mcMSmusmc^4*TAGnc)/12 + (mcMSmusmc^4*Pi^2*TAGnc)/27 - 
            (13*mcMSmusmc^4*TAGnc*Log[mus^2/mbMSmus^2])/54) + 
          Log[mcMSmusmc/mbMSmus]^2*((13*mcMSmusmc^4*TAGnc)/54 + 
            (mcMSmusmc^4*TAGnc*Log[mus^2/mbMSmus^2])/9) + 
          (mcMSmusmc^4*TAGnc*Zeta[3])/3) + Log[mcMSmusmc/mbMSmus]*
         ((2713*mcMSmusmc^4*TAGnc)/648 - (217*mcMSmusmc^4*Pi^2*TAGnc)/432 - 
          (5*mcMSmusmc^4*TAGnc^2)/12 + (mcMSmusmc^4*Pi^2*TAGnc^2)/27 + 
          ((89*mcMSmusmc^4*TAGnc)/36 - (13*mcMSmusmc^4*TAGnc^2)/54)*
           Log[mus^2/mbMSmus^2] + (20*mcMSmusmc^4*TAGnc*
            Log[musmc^2/mcMSmusmc^2])/9 + (14*mcMSmusmc^4*TAGnc*Zeta[3])/9))/
       mbMSmus^3 + mbMSmus*(8462917/93312 + (652841*Pi^2)/38880 - 
        (695*Pi^4)/7776 - (231847*TAGnc)/23328 - (991*Pi^2*TAGnc)/648 + 
        (61*Pi^4*TAGnc)/1944 + (2353*TAGnc^2)/23328 + (13*Pi^2*TAGnc^2)/324 + 
        ((-22*Pi^2)/81 + (2*Pi^2*TAGnc)/81)*Log[2]^2 + (-55/162 + TAGnc/81)*
         Log[2]^4 + (21715/864 - (385*TAGnc)/144 + (13*TAGnc^2)/216)*
         Log[mus^2/mbMSmus^2]^2 + (1861/432 - (43*TAGnc)/108 + TAGnc^2/108)*
         Log[mus^2/mbMSmus^2]^3 + Log[2]*((-575*Pi^2)/162 - 
          (11*Pi^2*TAGnc)/81 + ((13*Pi^2)/18 - (Pi^2*TAGnc)/27)*
           Log[mus^2/mbMSmus^2]) + (4*TAGkinmc*Log[mus^2/mcMSmusmc^2])/9 + 
        ((2*TAGkinmc)/27 - (2*TAGkinmc^2)/27)*Log[mus^2/mcMSmusmc^2]^2 - 
        (4*TAGkinmc*Log[mus^2/musmc^2])/9 - 
        (4*TAGkinmc*Log[musmc^2/mcMSmusmc^2])/9 - (220*PolyLog[4, 1/2])/27 + 
        (8*TAGnc*PolyLog[4, 1/2])/27 + NLxMSOS^2*(2353/23328 + 
          (13*Pi^2)/324 + (89/648 + Pi^2/54)*Log[mus^2/mbMSmus^2] + 
          (13*Log[mus^2/mbMSmus^2]^2)/216 + Log[mus^2/mbMSmus^2]^3/108 + 
          (7*Zeta[3])/54) + (58*Zeta[3])/27 - (1439*Pi^2*Zeta[3])/432 - 
        (241*TAGnc*Zeta[3])/72 + (7*TAGnc^2*Zeta[3])/54 + 
        Log[mus^2/mbMSmus^2]*(93391/1296 + (13*Pi^2)/6 - (5171*TAGnc)/648 - 
          (17*Pi^2*TAGnc)/36 + (89*TAGnc^2)/648 + (Pi^2*TAGnc^2)/54 + 
          (TAGkinmc*Log[mus^2/mcMSmusmc^2])/3 + (TAGkinmc/18 - TAGkinmc^2/18)*
           Log[mus^2/mcMSmusmc^2]^2 - (TAGkinmc*Log[mus^2/musmc^2])/3 - 
          (TAGkinmc*Log[musmc^2/mcMSmusmc^2])/3 - (23*Zeta[3])/12 - 
          (7*TAGnc*Zeta[3])/9) + NLxMSOS*(-231847/23328 - (991*Pi^2)/648 + 
          (61*Pi^4)/1944 + (2353*TAGnc)/11664 + (13*Pi^2*TAGnc)/162 + 
          (2*Pi^2*Log[2]^2)/81 + Log[2]^4/81 + (-385/144 + (13*TAGnc)/108)*
           Log[mus^2/mbMSmus^2]^2 + (-43/108 + TAGnc/54)*Log[mus^2/mbMSmus^2]^
            3 + Log[2]*((-11*Pi^2)/81 - (Pi^2*Log[mus^2/mbMSmus^2])/27) + 
          (8*PolyLog[4, 1/2])/27 + Log[mus^2/mbMSmus^2]*(-5171/648 - 
            (17*Pi^2)/36 + (89*TAGnc)/324 + (Pi^2*TAGnc)/27 - 
            (7*Zeta[3])/9) - (241*Zeta[3])/72 + (7*TAGnc*Zeta[3])/27) + 
        (1975*Zeta[5])/216) + ((-70*mcMSmusmc^2*TAGnc)/9 - 
        (13*mcMSmusmc^2*Pi^2*TAGnc)/12 + (2*mcMSmusmc^2*TAGnc^2)/9 + 
        (-4*mcMSmusmc^2*TAGnc - (3*mcMSmusmc^2*Pi^2*TAGnc)/2 + 
          (mcMSmusmc^2*Pi^4*TAGnc)/12)*Log[mcMSmusmc/mbMSmus] + 
        ((-9*mcMSmusmc^2*TAGnc)/2 + (mcMSmusmc^2*TAGnc^2)/3)*
         Log[mus^2/mbMSmus^2] + NLxMSOS*((2*mcMSmusmc^2*TAGnc)/9 + 
          (mcMSmusmc^2*TAGnc*Log[mus^2/mbMSmus^2])/3) - 
        2*mcMSmusmc^2*TAGnc*Log[musmc^2/mcMSmusmc^2] - 
        (11*mcMSmusmc^2*TAGnc*Zeta[3])/2 + (3*mcMSmusmc^2*Pi^2*TAGnc*Zeta[3])/
         4 + mufac^2*(-22055/108 + (13805*NLxOSKIN)/648 - 
          (209*NLxOSKIN^2)/486 + (437*Pi^2)/36 - (23*NLxOSKIN*Pi^2)/27 + 
          (NLxOSKIN^2*Pi^2)/81 - Pi^4/4 - (11*TAGkinmc)/108 - 
          (71*TAGnc)/216 - (Pi^2*TAGnc)/27 + (2*Pi^2*Log[2])/27 + 
          (-121/6 + (22*NLxOSKIN)/9 - (2*NLxOSKIN^2)/27)*Log[(2*mufac)/mus]^
            2 + (23/36 - TAGnc/18)*Log[mus^2/mbMSmus^2]^2 + 
          NLxMSOS*(-71/216 - Pi^2/27 - (13*Log[mus^2/mbMSmus^2])/54 - 
            Log[mus^2/mbMSmus^2]^2/18) + ((127*TAGkinmc)/36 - 
            (13*NLxOSKIN*TAGkinmc)/81 - (Pi^2*TAGkinmc)/9)*
           Log[mus^2/mcMSmusmc^2] - (TAGkinmc*Log[mus^2/mcMSmusmc^2]^2)/54 + 
          Log[mus^2/mbMSmus^2]*(1409/108 - (13*NLxOSKIN)/27 - Pi^2/3 - 
            (13*TAGnc)/54 - (TAGkinmc*Log[mus^2/mcMSmusmc^2])/9) + 
          Log[(2*mufac)/mus]*(689/6 - (691*NLxOSKIN)/54 + (26*NLxOSKIN^2)/
             81 - (11*Pi^2)/3 + (2*NLxOSKIN*Pi^2)/9 + 
            (-11/3 + (2*NLxOSKIN)/9)*Log[mus^2/mbMSmus^2] + 
            ((-11*TAGkinmc)/9 + (2*NLxOSKIN*TAGkinmc)/27)*
             Log[mus^2/mcMSmusmc^2]) + (2*TAGkinmc*Log[mus^2/musmc^2])/9 + 
          (1535*Zeta[3])/36 - (35*NLxOSKIN*Zeta[3])/18) + 
        (5*mcMSmusmc^2*TAGnc*Zeta[5])/2)/mbMSmus)


RunDecmkin2mMSnl = mkin - (4*apinlmus*mkin)/3 - (3019*apinlmus^2*mkin)/288 - 
    (9514621*apinlmus^3*mkin)/93312 + (16*apinlmus*mufac)/9 + 
    (892*apinlmus^2*mufac)/27 + (22496*apinlmus^3*mufac)/27 + 
    (2*apinlmus*mufac^2)/(3*mkin) + (95*apinlmus^2*mufac^2)/(9*mkin) + 
    (75931*apinlmus^3*mufac^2)/(324*mkin) + (71*apinlmus^2*mkin*NLxMSOS)/
     144 + (246643*apinlmus^3*mkin*NLxMSOS)/23328 - 
    (11*apinlmus^3*mufac*NLxMSOS)/27 - (11*apinlmus^3*mufac^2*NLxMSOS)/
     (72*mkin) - (2353*apinlmus^3*mkin*NLxMSOS^2)/23328 - 
    (128*apinlmus^2*mufac*NLxOSKIN)/81 - (20303*apinlmus^3*mufac*NLxOSKIN)/
     243 - (13*apinlmus^2*mufac^2*NLxOSKIN)/(27*mkin) - 
    (14429*apinlmus^3*mufac^2*NLxOSKIN)/(648*mkin) + 
    (1292*apinlmus^3*mufac*NLxOSKIN^2)/729 + 
    (209*apinlmus^3*mufac^2*NLxOSKIN^2)/(486*mkin) - 
    (apinlmus^2*mkin*Pi^2)/3 - (644201*apinlmus^3*mkin*Pi^2)/38880 - 
    (8*apinlmus^2*mufac*Pi^2)/9 - (1054*apinlmus^3*mufac*Pi^2)/27 - 
    (apinlmus^2*mufac^2*Pi^2)/(3*mkin) - (461*apinlmus^3*mufac^2*Pi^2)/
     (36*mkin) + (apinlmus^2*mkin*NLxMSOS*Pi^2)/18 + 
    (967*apinlmus^3*mkin*NLxMSOS*Pi^2)/648 + 
    (8*apinlmus^3*mufac*NLxMSOS*Pi^2)/81 + (apinlmus^3*mufac^2*NLxMSOS*Pi^2)/
     (27*mkin) - (13*apinlmus^3*mkin*NLxMSOS^2*Pi^2)/324 + 
    (208*apinlmus^3*mufac*NLxOSKIN*Pi^2)/81 + 
    (23*apinlmus^3*mufac^2*NLxOSKIN*Pi^2)/(27*mkin) - 
    (8*apinlmus^3*mufac*NLxOSKIN^2*Pi^2)/243 - 
    (apinlmus^3*mufac^2*NLxOSKIN^2*Pi^2)/(81*mkin) + 
    (695*apinlmus^3*mkin*Pi^4)/7776 + (2*apinlmus^3*mufac*Pi^4)/3 + 
    (apinlmus^3*mufac^2*Pi^4)/(4*mkin) - (61*apinlmus^3*mkin*NLxMSOS*Pi^4)/
     1944 + (11*apinlmus^3*mkin*TAGdecmcMSOS)/54 - 
    (9727*apinlmus^2*mcMSmusmc^24*TAGnc)/(186631200*mkin^23) - 
    (8583265857182497000029852773317049*apinlmus^3*mcMSmusmc^24*TAGnc)/
     (2369996943791663839368455777157120000*mkin^23) - 
    (39833*apinlmus^2*mcMSmusmc^22*TAGnc)/(520109667*mkin^21) - 
    (48833406986106441810878913193*apinlmus^3*mcMSmusmc^22*TAGnc)/
     (10790864167779053267146997760000*mkin^21) - 
    (39163*apinlmus^2*mcMSmusmc^20*TAGnc)/(333852800*mkin^19) - 
    (316346323888931757384696529*apinlmus^3*mcMSmusmc^20*TAGnc)/
     (54499313978682087207813120000*mkin^19) - 
    (16277*apinlmus^2*mcMSmusmc^18*TAGnc)/(86028075*mkin^17) - 
    (385021249205567851684429*apinlmus^3*mcMSmusmc^18*TAGnc)/
     (50057687427569200963584000*mkin^17) - 
    (529*apinlmus^2*mcMSmusmc^16*TAGnc)/(1622400*mkin^15) - 
    (1009988661120425340647*apinlmus^3*mcMSmusmc^16*TAGnc)/
     (95095681387606867968000*mkin^15) - 
    (15371*apinlmus^2*mcMSmusmc^14*TAGnc)/(25050025*mkin^13) - 
    (911686115452983709*apinlmus^3*mcMSmusmc^14*TAGnc)/
     (58963096098466560000*mkin^13) - (1229*apinlmus^2*mcMSmusmc^12*TAGnc)/
     (940896*mkin^11) - (593235847648457*apinlmus^3*mcMSmusmc^12*TAGnc)/
     (24647147078400000*mkin^11) - (997*apinlmus^2*mcMSmusmc^10*TAGnc)/
     (297675*mkin^9) - (208601138911*apinlmus^3*mcMSmusmc^10*TAGnc)/
     (5184974592000*mkin^9) - (463*apinlmus^2*mcMSmusmc^8*TAGnc)/
     (39200*mkin^7) - (72581279*apinlmus^3*mcMSmusmc^8*TAGnc)/
     (1185408000*mkin^7) - (19*apinlmus^2*mcMSmusmc^6*TAGnc)/(225*mkin^5) + 
    (1777621*apinlmus^3*mcMSmusmc^6*TAGnc)/(3499200*mkin^5) + 
    (151*apinlmus^2*mcMSmusmc^4*TAGnc)/(216*mkin^3) + 
    (3211*apinlmus^3*mcMSmusmc^4*TAGnc)/(2592*mkin^3) + 
    (apinlmus^2*mcMSmusmc^2*TAGnc)/mkin + (88*apinlmus^3*mcMSmusmc^2*TAGnc)/
     (9*mkin) + (71*apinlmus^2*mkin*TAGnc)/144 + 
    (246643*apinlmus^3*mkin*TAGnc)/23328 - (11*apinlmus^3*mufac*TAGnc)/27 + 
    (289*apinlmus^3*mcMSmusmc^24*mufac*TAGnc)/(198450*mkin^24) + 
    (62384*apinlmus^3*mcMSmusmc^22*mufac*TAGnc)/(31843449*mkin^22) + 
    (1417*apinlmus^3*mcMSmusmc^20*mufac*TAGnc)/(520200*mkin^20) + 
    (10576*apinlmus^3*mcMSmusmc^18*mufac*TAGnc)/(2679075*mkin^18) + 
    (661*apinlmus^3*mcMSmusmc^16*mufac*TAGnc)/(109512*mkin^16) + 
    (13232*apinlmus^3*mcMSmusmc^14*mufac*TAGnc)/(1334025*mkin^14) + 
    (79*apinlmus^3*mcMSmusmc^12*mufac*TAGnc)/(4374*mkin^12) + 
    (3824*apinlmus^3*mcMSmusmc^10*mufac*TAGnc)/(99225*mkin^10) + 
    (49*apinlmus^3*mcMSmusmc^8*mufac*TAGnc)/(450*mkin^8) + 
    (16*apinlmus^3*mcMSmusmc^6*mufac*TAGnc)/(27*mkin^6) - 
    (22*apinlmus^3*mcMSmusmc^4*mufac*TAGnc)/(9*mkin^4) - 
    (16*apinlmus^3*mcMSmusmc^2*mufac*TAGnc)/(9*mkin^2) + 
    (289*apinlmus^3*mcMSmusmc^24*mufac^2*TAGnc)/(529200*mkin^25) + 
    (7798*apinlmus^3*mcMSmusmc^22*mufac^2*TAGnc)/(10614483*mkin^23) + 
    (1417*apinlmus^3*mcMSmusmc^20*mufac^2*TAGnc)/(1387200*mkin^21) + 
    (1322*apinlmus^3*mcMSmusmc^18*mufac^2*TAGnc)/(893025*mkin^19) + 
    (661*apinlmus^3*mcMSmusmc^16*mufac^2*TAGnc)/(292032*mkin^17) + 
    (1654*apinlmus^3*mcMSmusmc^14*mufac^2*TAGnc)/(444675*mkin^15) + 
    (79*apinlmus^3*mcMSmusmc^12*mufac^2*TAGnc)/(11664*mkin^13) + 
    (478*apinlmus^3*mcMSmusmc^10*mufac^2*TAGnc)/(33075*mkin^11) + 
    (49*apinlmus^3*mcMSmusmc^8*mufac^2*TAGnc)/(1200*mkin^9) + 
    (2*apinlmus^3*mcMSmusmc^6*mufac^2*TAGnc)/(9*mkin^7) - 
    (11*apinlmus^3*mcMSmusmc^4*mufac^2*TAGnc)/(12*mkin^5) - 
    (2*apinlmus^3*mcMSmusmc^2*mufac^2*TAGnc)/(3*mkin^3) - 
    (11*apinlmus^3*mufac^2*TAGnc)/(72*mkin) - 
    (7801530877413386647*apinlmus^3*mcMSmusmc^24*NLxMSOS*TAGnc)/
     (37763267488425775296000*mkin^23) - 
    (1640519393726677*apinlmus^3*mcMSmusmc^22*NLxMSOS*TAGnc)/
     (5946258455651273760*mkin^21) - (22189567531163017*apinlmus^3*
      mcMSmusmc^20*NLxMSOS*TAGnc)/(58347124817357376000*mkin^19) - 
    (20555048260909*apinlmus^3*mcMSmusmc^18*NLxMSOS*TAGnc)/
     (37681809223558500*mkin^17) - (214558103603*apinlmus^3*mcMSmusmc^16*
      NLxMSOS*TAGnc)/(260460712512000*mkin^15) - 
    (108352091581*apinlmus^3*mcMSmusmc^14*NLxMSOS*TAGnc)/
     (81243243081000*mkin^13) - (326802499*apinlmus^3*mcMSmusmc^12*NLxMSOS*
      TAGnc)/(136928594880*mkin^11) - (2816347*apinlmus^3*mcMSmusmc^10*
      NLxMSOS*TAGnc)/(562605750*mkin^9) - 
    (174787*apinlmus^3*mcMSmusmc^8*NLxMSOS*TAGnc)/(12348000*mkin^7) - 
    (2729*apinlmus^3*mcMSmusmc^6*NLxMSOS*TAGnc)/(30375*mkin^5) + 
    (1423*apinlmus^3*mcMSmusmc^4*NLxMSOS*TAGnc)/(3888*mkin^3) - 
    (2*apinlmus^3*mcMSmusmc^2*NLxMSOS*TAGnc)/(9*mkin) - 
    (2353*apinlmus^3*mkin*NLxMSOS*TAGnc)/11664 - 
    (apinlmus^2*mcMSmusmc*Pi^2*TAGnc)/6 - 
    (6025*apinlmus^3*mcMSmusmc*Pi^2*TAGnc)/486 - 
    (106559417*apinlmus^3*mcMSmusmc^25*Pi^2*TAGnc)/(637438863600*mkin^24) - 
    (240817793781176357*apinlmus^3*mcMSmusmc^24*Pi^2*TAGnc)/
     (628867544642695987200*mkin^23) - (171475369*apinlmus^3*mcMSmusmc^23*
      Pi^2*TAGnc)/(778168119960*mkin^22) - 
    (174200135864459*apinlmus^3*mcMSmusmc^22*Pi^2*TAGnc)/
     (349543472195174400*mkin^21) - (846435761*apinlmus^3*mcMSmusmc^21*Pi^2*
      TAGnc)/(2832319518840*mkin^20) - (9007367733163*apinlmus^3*mcMSmusmc^20*
      Pi^2*TAGnc)/(13481015456563200*mkin^19) - 
    (86836957*apinlmus^3*mcMSmusmc^19*Pi^2*TAGnc)/(206871887520*mkin^18) - 
    (26463251891*apinlmus^3*mcMSmusmc^18*Pi^2*TAGnc)/
     (28454994247680*mkin^17) - (22757641*apinlmus^3*mcMSmusmc^17*Pi^2*TAGnc)/
     (36923796000*mkin^16) - (58806560951*apinlmus^3*mcMSmusmc^16*Pi^2*TAGnc)/
     (43261891706880*mkin^15) - (17346493*apinlmus^3*mcMSmusmc^15*Pi^2*TAGnc)/
     (18087549480*mkin^14) - (82285201*apinlmus^3*mcMSmusmc^14*Pi^2*TAGnc)/
     (38745907200*mkin^13) - (31786481*apinlmus^3*mcMSmusmc^13*Pi^2*TAGnc)/
     (19677663720*mkin^12) - (84041429*apinlmus^3*mcMSmusmc^12*Pi^2*TAGnc)/
     (22992076800*mkin^11) - (1055689*apinlmus^3*mcMSmusmc^11*Pi^2*TAGnc)/
     (345779280*mkin^10) - (359801*apinlmus^3*mcMSmusmc^10*Pi^2*TAGnc)/
     (48988800*mkin^9) - (416029*apinlmus^3*mcMSmusmc^9*Pi^2*TAGnc)/
     (60011280*mkin^8) - (1789*apinlmus^3*mcMSmusmc^8*Pi^2*TAGnc)/
     (90720*mkin^7) - (262769*apinlmus^3*mcMSmusmc^7*Pi^2*TAGnc)/
     (11907000*mkin^6) - (1427*apinlmus^3*mcMSmusmc^6*Pi^2*TAGnc)/
     (12960*mkin^5) - (149*apinlmus^3*mcMSmusmc^5*Pi^2*TAGnc)/(1080*mkin^4) + 
    (apinlmus^2*mcMSmusmc^4*Pi^2*TAGnc)/(18*mkin^3) - 
    (10*apinlmus^3*mcMSmusmc^4*Pi^2*TAGnc)/(81*mkin^3) - 
    (apinlmus^2*mcMSmusmc^3*Pi^2*TAGnc)/(6*mkin^2) - 
    (19507*apinlmus^3*mcMSmusmc^3*Pi^2*TAGnc)/(1620*mkin^2) + 
    (13*apinlmus^3*mcMSmusmc^2*Pi^2*TAGnc)/(12*mkin) + 
    (apinlmus^2*mkin*Pi^2*TAGnc)/18 + (967*apinlmus^3*mkin*Pi^2*TAGnc)/648 + 
    (8*apinlmus^3*mufac*Pi^2*TAGnc)/81 - 
    (8*apinlmus^3*mcMSmusmc^4*mufac*Pi^2*TAGnc)/(27*mkin^4) + 
    (16*apinlmus^3*mcMSmusmc^3*mufac*Pi^2*TAGnc)/(27*mkin^3) - 
    (apinlmus^3*mcMSmusmc^4*mufac^2*Pi^2*TAGnc)/(9*mkin^5) + 
    (2*apinlmus^3*mcMSmusmc^3*mufac^2*Pi^2*TAGnc)/(9*mkin^4) + 
    (apinlmus^3*mufac^2*Pi^2*TAGnc)/(27*mkin) + 
    (7*apinlmus^3*mcMSmusmc*NLxMSOS*Pi^2*TAGnc)/27 + 
    (23*apinlmus^3*mcMSmusmc^25*NLxMSOS*Pi^2*TAGnc)/(415800*mkin^24) - 
    (11*apinlmus^3*mcMSmusmc^24*NLxMSOS*Pi^2*TAGnc)/(260820*mkin^23) + 
    (7*apinlmus^3*mcMSmusmc^23*NLxMSOS*Pi^2*TAGnc)/(96140*mkin^22) - 
    (20*apinlmus^3*mcMSmusmc^22*NLxMSOS*Pi^2*TAGnc)/(355509*mkin^21) + 
    (19*apinlmus^3*mcMSmusmc^21*NLxMSOS*Pi^2*TAGnc)/(192780*mkin^20) - 
    (apinlmus^3*mcMSmusmc^20*NLxMSOS*Pi^2*TAGnc)/(12920*mkin^19) + 
    (17*apinlmus^3*mcMSmusmc^19*NLxMSOS*Pi^2*TAGnc)/(123120*mkin^18) - 
    (16*apinlmus^3*mcMSmusmc^18*NLxMSOS*Pi^2*TAGnc)/(144585*mkin^17) + 
    (5*apinlmus^3*mcMSmusmc^17*NLxMSOS*Pi^2*TAGnc)/(24752*mkin^16) - 
    (7*apinlmus^3*mcMSmusmc^16*NLxMSOS*Pi^2*TAGnc)/(42120*mkin^15) + 
    (13*apinlmus^3*mcMSmusmc^15*NLxMSOS*Pi^2*TAGnc)/(41580*mkin^14) - 
    (4*apinlmus^3*mcMSmusmc^14*NLxMSOS*Pi^2*TAGnc)/(15015*mkin^13) + 
    (11*apinlmus^3*mcMSmusmc^13*NLxMSOS*Pi^2*TAGnc)/(21060*mkin^12) - 
    (5*apinlmus^3*mcMSmusmc^12*NLxMSOS*Pi^2*TAGnc)/(10692*mkin^11) + 
    (3*apinlmus^3*mcMSmusmc^11*NLxMSOS*Pi^2*TAGnc)/(3080*mkin^10) - 
    (8*apinlmus^3*mcMSmusmc^10*NLxMSOS*Pi^2*TAGnc)/(8505*mkin^9) + 
    (7*apinlmus^3*mcMSmusmc^9*NLxMSOS*Pi^2*TAGnc)/(3240*mkin^8) - 
    (apinlmus^3*mcMSmusmc^8*NLxMSOS*Pi^2*TAGnc)/(420*mkin^7) + 
    (5*apinlmus^3*mcMSmusmc^7*NLxMSOS*Pi^2*TAGnc)/(756*mkin^6) - 
    (4*apinlmus^3*mcMSmusmc^6*NLxMSOS*Pi^2*TAGnc)/(405*mkin^5) + 
    (apinlmus^3*mcMSmusmc^5*NLxMSOS*Pi^2*TAGnc)/(20*mkin^4) + 
    (13*apinlmus^3*mcMSmusmc^4*NLxMSOS*Pi^2*TAGnc)/(324*mkin^3) + 
    (7*apinlmus^3*mcMSmusmc^3*NLxMSOS*Pi^2*TAGnc)/(54*mkin^2) - 
    (13*apinlmus^3*mkin*NLxMSOS*Pi^2*TAGnc)/162 - 
    (13*apinlmus^3*mcMSmusmc*Pi^3*TAGnc)/162 + 
    (5698043*apinlmus^3*mcMSmusmc^25*Pi^3*TAGnc)/(118908518400*mkin^24) + 
    (81991*apinlmus^3*mcMSmusmc^23*Pi^3*TAGnc)/(1374683136*mkin^22) + 
    (127699*apinlmus^3*mcMSmusmc^21*Pi^3*TAGnc)/(1684537344*mkin^20) + 
    (99671*apinlmus^3*mcMSmusmc^19*Pi^3*TAGnc)/(1008599040*mkin^18) + 
    (1925*apinlmus^3*mcMSmusmc^17*Pi^3*TAGnc)/(14483456*mkin^16) + 
    (377*apinlmus^3*mcMSmusmc^15*Pi^3*TAGnc)/(2027520*mkin^14) + 
    (1771*apinlmus^3*mcMSmusmc^13*Pi^3*TAGnc)/(6469632*mkin^12) + 
    (17*apinlmus^3*mcMSmusmc^11*Pi^3*TAGnc)/(39424*mkin^10) + 
    (77*apinlmus^3*mcMSmusmc^9*Pi^3*TAGnc)/(103680*mkin^8) + 
    (25*apinlmus^3*mcMSmusmc^7*Pi^3*TAGnc)/(18144*mkin^6) - 
    (apinlmus^3*mcMSmusmc^5*Pi^3*TAGnc)/(240*mkin^4) - 
    (7*apinlmus^3*mcMSmusmc^3*Pi^3*TAGnc)/(108*mkin^2) - 
    (271*apinlmus^3*mcMSmusmc^4*Pi^4*TAGnc)/(19440*mkin^3) - 
    (61*apinlmus^3*mkin*Pi^4*TAGnc)/1944 + 
    (279881*apinlmus^3*mcMSmusmc^24*TAGnc^2)/(4627692000*mkin^23) + 
    (624853*apinlmus^3*mcMSmusmc^22*TAGnc^2)/(6846429744*mkin^21) + 
    (367909*apinlmus^3*mcMSmusmc^20*TAGnc^2)/(2547216000*mkin^19) + 
    (157882*apinlmus^3*mcMSmusmc^18*TAGnc^2)/(650372247*mkin^17) + 
    (1373*apinlmus^3*mcMSmusmc^16*TAGnc^2)/(3110400*mkin^15) + 
    (233201*apinlmus^3*mcMSmusmc^14*TAGnc^2)/(260851500*mkin^13) + 
    (3989*apinlmus^3*mcMSmusmc^12*TAGnc^2)/(1881792*mkin^11) + 
    (1778*apinlmus^3*mcMSmusmc^10*TAGnc^2)/(273375*mkin^9) + 
    (241*apinlmus^3*mcMSmusmc^8*TAGnc^2)/(7056*mkin^7) + 
    (3481*apinlmus^3*mcMSmusmc^6*TAGnc^2)/(30375*mkin^5) + 
    (235*apinlmus^3*mcMSmusmc^4*TAGnc^2)/(3888*mkin^3) - 
    (2*apinlmus^3*mcMSmusmc^2*TAGnc^2)/(9*mkin) - 
    (2353*apinlmus^3*mkin*TAGnc^2)/23328 + 
    (4*apinlmus^3*mcMSmusmc*Pi^2*TAGnc^2)/135 - 
    (4*apinlmus^3*mcMSmusmc^6*Pi^2*TAGnc^2)/(405*mkin^5) + 
    (13*apinlmus^3*mcMSmusmc^4*Pi^2*TAGnc^2)/(324*mkin^3) - 
    (13*apinlmus^3*mkin*Pi^2*TAGnc^2)/324 - (apinlmus^2*mkin*Pi^2*Log[2])/9 + 
    (587*apinlmus^3*mkin*Pi^2*Log[2])/162 - (16*apinlmus^3*mufac*Pi^2*Log[2])/
     81 - (2*apinlmus^3*mufac^2*Pi^2*Log[2])/(27*mkin) + 
    (11*apinlmus^3*mkin*NLxMSOS*Pi^2*Log[2])/81 + 
    (1199*apinlmus^3*mcMSmusmc*Pi^2*TAGnc*Log[2])/81 + 
    (877591*apinlmus^3*mcMSmusmc^24*Pi^2*TAGnc*Log[2])/(6794772480*mkin^23) + 
    (38675*apinlmus^3*mcMSmusmc^22*Pi^2*TAGnc*Log[2])/(233570304*mkin^21) + 
    (143*apinlmus^3*mcMSmusmc^20*Pi^2*TAGnc*Log[2])/(655360*mkin^19) + 
    (4147*apinlmus^3*mcMSmusmc^18*Pi^2*TAGnc*Log[2])/(13934592*mkin^17) + 
    (1001*apinlmus^3*mcMSmusmc^16*Pi^2*TAGnc*Log[2])/(2359296*mkin^15) + 
    (23*apinlmus^3*mcMSmusmc^14*Pi^2*TAGnc*Log[2])/(35840*mkin^13) + 
    (175*apinlmus^3*mcMSmusmc^12*Pi^2*TAGnc*Log[2])/(165888*mkin^11) + 
    (17*apinlmus^3*mcMSmusmc^10*Pi^2*TAGnc*Log[2])/(8640*mkin^9) + 
    (7*apinlmus^3*mcMSmusmc^8*Pi^2*TAGnc*Log[2])/(1536*mkin^7) + 
    (11*apinlmus^3*mcMSmusmc^6*Pi^2*TAGnc*Log[2])/(648*mkin^5) - 
    (5*apinlmus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[2])/(144*mkin^3) + 
    (1199*apinlmus^3*mcMSmusmc^3*Pi^2*TAGnc*Log[2])/(81*mkin^2) + 
    (11*apinlmus^3*mkin*Pi^2*TAGnc*Log[2])/81 - 
    (2*apinlmus^3*mcMSmusmc*NLxMSOS*Pi^2*TAGnc*Log[2])/9 - 
    (2*apinlmus^3*mcMSmusmc^3*NLxMSOS*Pi^2*TAGnc*Log[2])/(9*mkin^2) + 
    (22*apinlmus^3*mkin*Pi^2*Log[2]^2)/81 - 
    (2*apinlmus^3*mkin*NLxMSOS*Pi^2*Log[2]^2)/81 - 
    (2*apinlmus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[2]^2)/(81*mkin^3) - 
    (2*apinlmus^3*mkin*Pi^2*TAGnc*Log[2]^2)/81 + 
    (55*apinlmus^3*mkin*Log[2]^4)/162 - (apinlmus^3*mkin*NLxMSOS*Log[2]^4)/
     81 - (5*apinlmus^3*mcMSmusmc^4*TAGnc*Log[2]^4)/(162*mkin^3) - 
    (apinlmus^3*mkin*TAGnc*Log[2]^4)/81 + 
    (11*apinlmus^2*mcMSmusmc^24*TAGnc*Log[mcMSmusmc/mkin])/(28980*mkin^23) + 
    (194881611047644800407*apinlmus^3*mcMSmusmc^24*TAGnc*Log[mcMSmusmc/mkin])/
     (20069484527233911369600*mkin^23) + 
    (20*apinlmus^2*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mkin])/(39501*mkin^21) + 
    (84250035134168729123*apinlmus^3*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mkin])/
     (7242811051244408640000*mkin^21) + 
    (9*apinlmus^2*mcMSmusmc^20*TAGnc*Log[mcMSmusmc/mkin])/(12920*mkin^19) + 
    (2488663198975906531*apinlmus^3*mcMSmusmc^20*TAGnc*Log[mcMSmusmc/mkin])/
     (175583298211985664000*mkin^19) + (16*apinlmus^2*mcMSmusmc^18*TAGnc*
      Log[mcMSmusmc/mkin])/(16065*mkin^17) + 
    (287916747607903391*apinlmus^3*mcMSmusmc^18*TAGnc*Log[mcMSmusmc/mkin])/
     (16342379002556006400*mkin^17) + (7*apinlmus^2*mcMSmusmc^16*TAGnc*
      Log[mcMSmusmc/mkin])/(4680*mkin^15) + 
    (26396860982549*apinlmus^3*mcMSmusmc^16*TAGnc*Log[mcMSmusmc/mkin])/
     (1178083838131200*mkin^15) + (12*apinlmus^2*mcMSmusmc^14*TAGnc*
      Log[mcMSmusmc/mkin])/(5005*mkin^13) + 
    (1194668026861*apinlmus^3*mcMSmusmc^14*TAGnc*Log[mcMSmusmc/mkin])/
     (40905688824000*mkin^13) + (5*apinlmus^2*mcMSmusmc^12*TAGnc*
      Log[mcMSmusmc/mkin])/(1188*mkin^11) + 
    (241232375387*apinlmus^3*mcMSmusmc^12*TAGnc*Log[mcMSmusmc/mkin])/
     (6224027040000*mkin^11) + (8*apinlmus^2*mcMSmusmc^10*TAGnc*
      Log[mcMSmusmc/mkin])/(945*mkin^9) + 
    (2042394191*apinlmus^3*mcMSmusmc^10*TAGnc*Log[mcMSmusmc/mkin])/
     (41150592000*mkin^9) + (3*apinlmus^2*mcMSmusmc^8*TAGnc*
      Log[mcMSmusmc/mkin])/(140*mkin^7) + 
    (413297*apinlmus^3*mcMSmusmc^8*TAGnc*Log[mcMSmusmc/mkin])/
     (12700800*mkin^7) + (4*apinlmus^2*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mkin])/
     (45*mkin^5) - (721229*apinlmus^3*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mkin])/
     (1166400*mkin^5) - (13*apinlmus^2*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin])/
     (18*mkin^3) - (4321*apinlmus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin])/
     (648*mkin^3) + (4*apinlmus^3*mcMSmusmc^2*TAGnc*Log[mcMSmusmc/mkin])/
     mkin - (44*apinlmus^3*mcMSmusmc^24*mufac*TAGnc*Log[mcMSmusmc/mkin])/
     (2835*mkin^24) - (320*apinlmus^3*mcMSmusmc^22*mufac*TAGnc*
      Log[mcMSmusmc/mkin])/(16929*mkin^22) - 
    (2*apinlmus^3*mcMSmusmc^20*mufac*TAGnc*Log[mcMSmusmc/mkin])/
     (85*mkin^20) - (256*apinlmus^3*mcMSmusmc^18*mufac*TAGnc*
      Log[mcMSmusmc/mkin])/(8505*mkin^18) - 
    (14*apinlmus^3*mcMSmusmc^16*mufac*TAGnc*Log[mcMSmusmc/mkin])/
     (351*mkin^16) - (64*apinlmus^3*mcMSmusmc^14*mufac*TAGnc*
      Log[mcMSmusmc/mkin])/(1155*mkin^14) - 
    (20*apinlmus^3*mcMSmusmc^12*mufac*TAGnc*Log[mcMSmusmc/mkin])/
     (243*mkin^12) - (128*apinlmus^3*mcMSmusmc^10*mufac*TAGnc*
      Log[mcMSmusmc/mkin])/(945*mkin^10) - 
    (4*apinlmus^3*mcMSmusmc^8*mufac*TAGnc*Log[mcMSmusmc/mkin])/(15*mkin^8) - 
    (64*apinlmus^3*mcMSmusmc^6*mufac*TAGnc*Log[mcMSmusmc/mkin])/(81*mkin^6) + 
    (8*apinlmus^3*mcMSmusmc^4*mufac*TAGnc*Log[mcMSmusmc/mkin])/(3*mkin^4) - 
    (11*apinlmus^3*mcMSmusmc^24*mufac^2*TAGnc*Log[mcMSmusmc/mkin])/
     (1890*mkin^25) - (40*apinlmus^3*mcMSmusmc^22*mufac^2*TAGnc*
      Log[mcMSmusmc/mkin])/(5643*mkin^23) - 
    (3*apinlmus^3*mcMSmusmc^20*mufac^2*TAGnc*Log[mcMSmusmc/mkin])/
     (340*mkin^21) - (32*apinlmus^3*mcMSmusmc^18*mufac^2*TAGnc*
      Log[mcMSmusmc/mkin])/(2835*mkin^19) - 
    (7*apinlmus^3*mcMSmusmc^16*mufac^2*TAGnc*Log[mcMSmusmc/mkin])/
     (468*mkin^17) - (8*apinlmus^3*mcMSmusmc^14*mufac^2*TAGnc*
      Log[mcMSmusmc/mkin])/(385*mkin^15) - 
    (5*apinlmus^3*mcMSmusmc^12*mufac^2*TAGnc*Log[mcMSmusmc/mkin])/
     (162*mkin^13) - (16*apinlmus^3*mcMSmusmc^10*mufac^2*TAGnc*
      Log[mcMSmusmc/mkin])/(315*mkin^11) - 
    (apinlmus^3*mcMSmusmc^8*mufac^2*TAGnc*Log[mcMSmusmc/mkin])/(10*mkin^9) - 
    (8*apinlmus^3*mcMSmusmc^6*mufac^2*TAGnc*Log[mcMSmusmc/mkin])/
     (27*mkin^7) + (apinlmus^3*mcMSmusmc^4*mufac^2*TAGnc*Log[mcMSmusmc/mkin])/
     mkin^5 - (73801799*apinlmus^3*mcMSmusmc^24*NLxMSOS*TAGnc*
      Log[mcMSmusmc/mkin])/(4231787807520*mkin^23) - 
    (14290513*apinlmus^3*mcMSmusmc^22*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin])/
     (689665418442*mkin^21) - (1211963*apinlmus^3*mcMSmusmc^20*NLxMSOS*TAGnc*
      Log[mcMSmusmc/mkin])/(50127997920*mkin^19) - 
    (197062*apinlmus^3*mcMSmusmc^18*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin])/
     (7381208835*mkin^17) - (355*apinlmus^3*mcMSmusmc^16*NLxMSOS*TAGnc*
      Log[mcMSmusmc/mkin])/(14455584*mkin^15) - 
    (1153*apinlmus^3*mcMSmusmc^14*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin])/
     (225450225*mkin^13) + (391*apinlmus^3*mcMSmusmc^12*NLxMSOS*TAGnc*
      Log[mcMSmusmc/mkin])/(4939704*mkin^11) + 
    (398*apinlmus^3*mcMSmusmc^10*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin])/
     (893025*mkin^9) + (37*apinlmus^3*mcMSmusmc^8*NLxMSOS*TAGnc*
      Log[mcMSmusmc/mkin])/(14700*mkin^7) + 
    (2*apinlmus^3*mcMSmusmc^6*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin])/
     (75*mkin^5) + (apinlmus^3*mcMSmusmc^4*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin])/
     (12*mkin^3) + (11*apinlmus^3*mcMSmusmc*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     6 + (533*apinlmus^3*mcMSmusmc^25*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (277725*mkin^24) + (1620816161*apinlmus^3*mcMSmusmc^24*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(822167470080*mkin^23) + 
    (445*apinlmus^3*mcMSmusmc^23*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (192717*mkin^22) + (13926181*apinlmus^3*mcMSmusmc^22*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(5839257600*mkin^21) + 
    (365*apinlmus^3*mcMSmusmc^21*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (128877*mkin^20) + (156353*apinlmus^3*mcMSmusmc^20*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(53084160*mkin^19) + 
    (293*apinlmus^3*mcMSmusmc^19*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (82365*mkin^18) + (52013*apinlmus^3*mcMSmusmc^18*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(13934592*mkin^17) + 
    (229*apinlmus^3*mcMSmusmc^17*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (49725*mkin^16) + (565351*apinlmus^3*mcMSmusmc^16*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(115605504*mkin^15) + 
    (173*apinlmus^3*mcMSmusmc^15*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (27885*mkin^14) + (2161*apinlmus^3*mcMSmusmc^14*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(322560*mkin^13) + 
    (125*apinlmus^3*mcMSmusmc^13*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (14157*mkin^12) + (40553*apinlmus^3*mcMSmusmc^12*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(4147200*mkin^11) + 
    (85*apinlmus^3*mcMSmusmc^11*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (6237*mkin^10) + (17*apinlmus^3*mcMSmusmc^10*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(1080*mkin^9) + 
    (53*apinlmus^3*mcMSmusmc^9*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (2205*mkin^8) + (139*apinlmus^3*mcMSmusmc^8*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(4608*mkin^7) + 
    (29*apinlmus^3*mcMSmusmc^7*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/(525*mkin^6) + 
    (113*apinlmus^3*mcMSmusmc^6*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (1296*mkin^5) + (13*apinlmus^3*mcMSmusmc^5*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(45*mkin^4) + 
    (217*apinlmus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (432*mkin^3) + (289*apinlmus^3*mcMSmusmc^3*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(162*mkin^2) + 
    (3*apinlmus^3*mcMSmusmc^2*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/(2*mkin) - 
    (apinlmus^3*mcMSmusmc*NLxMSOS*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/9 - 
    (apinlmus^3*mcMSmusmc^4*NLxMSOS*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (27*mkin^3) - (apinlmus^3*mcMSmusmc^3*NLxMSOS*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(9*mkin^2) - 
    (apinlmus^3*mcMSmusmc^2*Pi^4*TAGnc*Log[mcMSmusmc/mkin])/(12*mkin) - 
    (281941*apinlmus^3*mcMSmusmc^24*TAGnc^2*Log[mcMSmusmc/mkin])/
     (629880300*mkin^23) - (318845*apinlmus^3*mcMSmusmc^22*TAGnc^2*
      Log[mcMSmusmc/mkin])/(520109667*mkin^21) - 
    (1017529*apinlmus^3*mcMSmusmc^20*TAGnc^2*Log[mcMSmusmc/mkin])/
     (1168484800*mkin^19) - (334304*apinlmus^3*mcMSmusmc^18*TAGnc^2*
      Log[mcMSmusmc/mkin])/(258084225*mkin^17) - 
    (7469*apinlmus^3*mcMSmusmc^16*TAGnc^2*Log[mcMSmusmc/mkin])/
     (3650400*mkin^15) - (263546*apinlmus^3*mcMSmusmc^14*TAGnc^2*
      Log[mcMSmusmc/mkin])/(75150075*mkin^13) - 
    (9545*apinlmus^3*mcMSmusmc^12*TAGnc^2*Log[mcMSmusmc/mkin])/
     (1411344*mkin^11) - (14048*apinlmus^3*mcMSmusmc^10*TAGnc^2*
      Log[mcMSmusmc/mkin])/(893025*mkin^9) - 
    (193*apinlmus^3*mcMSmusmc^8*TAGnc^2*Log[mcMSmusmc/mkin])/(3675*mkin^7) - 
    (8*apinlmus^3*mcMSmusmc^6*TAGnc^2*Log[mcMSmusmc/mkin])/(75*mkin^5) + 
    (5*apinlmus^3*mcMSmusmc^4*TAGnc^2*Log[mcMSmusmc/mkin])/(12*mkin^3) - 
    (apinlmus^3*mcMSmusmc*Pi^2*TAGnc^2*Log[mcMSmusmc/mkin])/9 - 
    (apinlmus^3*mcMSmusmc^4*Pi^2*TAGnc^2*Log[mcMSmusmc/mkin])/(27*mkin^3) - 
    (apinlmus^3*mcMSmusmc^3*Pi^2*TAGnc^2*Log[mcMSmusmc/mkin])/(9*mkin^2) + 
    (apinlmus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[2]*Log[mcMSmusmc/mkin])/
     (9*mkin^3) + (941278536953*apinlmus^3*mcMSmusmc^24*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(46646042002560*mkin^23) + 
    (3036958280711*apinlmus^3*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (124450902576000*mkin^21) + (45436526101*apinlmus^3*mcMSmusmc^20*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(1508495788800*mkin^19) + 
    (3391310411*apinlmus^3*mcMSmusmc^18*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (88921857024*mkin^17) + (1141880431*apinlmus^3*mcMSmusmc^16*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(22884301440*mkin^15) + 
    (4128997*apinlmus^3*mcMSmusmc^14*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (60540480*mkin^13) + (14841023*apinlmus^3*mcMSmusmc^12*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(149688000*mkin^11) + 
    (2585963*apinlmus^3*mcMSmusmc^10*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (16329600*mkin^9) + (9019*apinlmus^3*mcMSmusmc^8*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(30240*mkin^7) + 
    (5329*apinlmus^3*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mkin]^2)/(6480*mkin^5) + 
    (apinlmus^2*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]^2)/(3*mkin^3) + 
    (95*apinlmus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]^2)/(18*mkin^3) - 
    (16*apinlmus^3*mcMSmusmc^4*mufac*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (9*mkin^4) - (2*apinlmus^3*mcMSmusmc^4*mufac^2*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(3*mkin^5) - 
    (13*apinlmus^3*mcMSmusmc^4*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (54*mkin^3) - (apinlmus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (4*mkin^3) + (11*apinlmus^3*mcMSmusmc^24*TAGnc^2*Log[mcMSmusmc/mkin]^2)/
     (43470*mkin^23) + (40*apinlmus^3*mcMSmusmc^22*TAGnc^2*
      Log[mcMSmusmc/mkin]^2)/(118503*mkin^21) + 
    (3*apinlmus^3*mcMSmusmc^20*TAGnc^2*Log[mcMSmusmc/mkin]^2)/
     (6460*mkin^19) + (32*apinlmus^3*mcMSmusmc^18*TAGnc^2*
      Log[mcMSmusmc/mkin]^2)/(48195*mkin^17) + 
    (7*apinlmus^3*mcMSmusmc^16*TAGnc^2*Log[mcMSmusmc/mkin]^2)/
     (7020*mkin^15) + (8*apinlmus^3*mcMSmusmc^14*TAGnc^2*
      Log[mcMSmusmc/mkin]^2)/(5005*mkin^13) + 
    (5*apinlmus^3*mcMSmusmc^12*TAGnc^2*Log[mcMSmusmc/mkin]^2)/
     (1782*mkin^11) + (16*apinlmus^3*mcMSmusmc^10*TAGnc^2*
      Log[mcMSmusmc/mkin]^2)/(2835*mkin^9) + 
    (apinlmus^3*mcMSmusmc^8*TAGnc^2*Log[mcMSmusmc/mkin]^2)/(70*mkin^7) - 
    (13*apinlmus^3*mcMSmusmc^4*TAGnc^2*Log[mcMSmusmc/mkin]^2)/(54*mkin^3) - 
    (apinlmus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]^3)/(3*mkin^3) + 
    (2*apinlmus^3*mcMSmusmc^4*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]^3)/
     (27*mkin^3) + (2*apinlmus^3*mcMSmusmc^4*TAGnc^2*Log[mcMSmusmc/mkin]^3)/
     (27*mkin^3) - (88*apinlmus^2*mufac*Log[(2*mufac)/mus])/9 - 
    (3416*apinlmus^3*mufac*Log[(2*mufac)/mus])/9 - 
    (11*apinlmus^2*mufac^2*Log[(2*mufac)/mus])/(3*mkin) - 
    (733*apinlmus^3*mufac^2*Log[(2*mufac)/mus])/(6*mkin) + 
    (16*apinlmus^2*mufac*NLxOSKIN*Log[(2*mufac)/mus])/27 + 
    (3388*apinlmus^3*mufac*NLxOSKIN*Log[(2*mufac)/mus])/81 + 
    (2*apinlmus^2*mufac^2*NLxOSKIN*Log[(2*mufac)/mus])/(9*mkin) + 
    (715*apinlmus^3*mufac^2*NLxOSKIN*Log[(2*mufac)/mus])/(54*mkin) - 
    (256*apinlmus^3*mufac*NLxOSKIN^2*Log[(2*mufac)/mus])/243 - 
    (26*apinlmus^3*mufac^2*NLxOSKIN^2*Log[(2*mufac)/mus])/(81*mkin) + 
    (88*apinlmus^3*mufac*Pi^2*Log[(2*mufac)/mus])/9 + 
    (11*apinlmus^3*mufac^2*Pi^2*Log[(2*mufac)/mus])/(3*mkin) - 
    (16*apinlmus^3*mufac*NLxOSKIN*Pi^2*Log[(2*mufac)/mus])/27 - 
    (2*apinlmus^3*mufac^2*NLxOSKIN*Pi^2*Log[(2*mufac)/mus])/(9*mkin) + 
    (484*apinlmus^3*mufac*Log[(2*mufac)/mus]^2)/9 + 
    (121*apinlmus^3*mufac^2*Log[(2*mufac)/mus]^2)/(6*mkin) - 
    (176*apinlmus^3*mufac*NLxOSKIN*Log[(2*mufac)/mus]^2)/27 - 
    (22*apinlmus^3*mufac^2*NLxOSKIN*Log[(2*mufac)/mus]^2)/(9*mkin) + 
    (16*apinlmus^3*mufac*NLxOSKIN^2*Log[(2*mufac)/mus]^2)/81 + 
    (2*apinlmus^3*mufac^2*NLxOSKIN^2*Log[(2*mufac)/mus]^2)/(27*mkin) - 
    (2*apinlmus^2*mkin*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/9 - 
    (3931*apinlmus^3*mkin*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/864 + 
    (16*apinlmus^3*mufac*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/81 + 
    (2*apinlmus^3*mufac^2*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/(27*mkin) + 
    (71*apinlmus^3*mkin*NLxMSOS*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/432 - 
    (apinlmus^3*mkin*Pi^2*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/9 + 
    (apinlmus^3*mkin*NLxMSOS*Pi^2*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/54 - 
    (9727*apinlmus^3*mcMSmusmc^24*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/
     (559893600*mkin^23) - (39833*apinlmus^3*mcMSmusmc^22*TAGdecmcMSOS*TAGnc*
      Log[mus^2/mcMSmusmc^2])/(1560329001*mkin^21) - 
    (39163*apinlmus^3*mcMSmusmc^20*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/
     (1001558400*mkin^19) - (16277*apinlmus^3*mcMSmusmc^18*TAGdecmcMSOS*TAGnc*
      Log[mus^2/mcMSmusmc^2])/(258084225*mkin^17) - 
    (529*apinlmus^3*mcMSmusmc^16*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/
     (4867200*mkin^15) - (15371*apinlmus^3*mcMSmusmc^14*TAGdecmcMSOS*TAGnc*
      Log[mus^2/mcMSmusmc^2])/(75150075*mkin^13) - 
    (1229*apinlmus^3*mcMSmusmc^12*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/
     (2822688*mkin^11) - (997*apinlmus^3*mcMSmusmc^10*TAGdecmcMSOS*TAGnc*
      Log[mus^2/mcMSmusmc^2])/(893025*mkin^9) - 
    (463*apinlmus^3*mcMSmusmc^8*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/
     (117600*mkin^7) - (19*apinlmus^3*mcMSmusmc^6*TAGdecmcMSOS*TAGnc*
      Log[mus^2/mcMSmusmc^2])/(675*mkin^5) + 
    (151*apinlmus^3*mcMSmusmc^4*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/
     (648*mkin^3) + (apinlmus^3*mcMSmusmc^2*TAGdecmcMSOS*TAGnc*
      Log[mus^2/mcMSmusmc^2])/(3*mkin) + 
    (71*apinlmus^3*mkin*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/432 - 
    (apinlmus^3*mcMSmusmc*Pi^2*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/
     18 + (apinlmus^3*mcMSmusmc^4*Pi^2*TAGdecmcMSOS*TAGnc*
      Log[mus^2/mcMSmusmc^2])/(54*mkin^3) - 
    (apinlmus^3*mcMSmusmc^3*Pi^2*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2])/
     (18*mkin^2) + (apinlmus^3*mkin*Pi^2*TAGdecmcMSOS*TAGnc*
      Log[mus^2/mcMSmusmc^2])/54 - (apinlmus^3*mkin*Pi^2*TAGdecmcMSOS*Log[2]*
      Log[mus^2/mcMSmusmc^2])/27 + (11*apinlmus^3*mcMSmusmc^24*TAGdecmcMSOS*
      TAGnc*Log[mcMSmusmc/mkin]*Log[mus^2/mcMSmusmc^2])/(86940*mkin^23) + 
    (20*apinlmus^3*mcMSmusmc^22*TAGdecmcMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mcMSmusmc^2])/(118503*mkin^21) + 
    (3*apinlmus^3*mcMSmusmc^20*TAGdecmcMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mcMSmusmc^2])/(12920*mkin^19) + 
    (16*apinlmus^3*mcMSmusmc^18*TAGdecmcMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mcMSmusmc^2])/(48195*mkin^17) + 
    (7*apinlmus^3*mcMSmusmc^16*TAGdecmcMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mcMSmusmc^2])/(14040*mkin^15) + 
    (4*apinlmus^3*mcMSmusmc^14*TAGdecmcMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mcMSmusmc^2])/(5005*mkin^13) + 
    (5*apinlmus^3*mcMSmusmc^12*TAGdecmcMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mcMSmusmc^2])/(3564*mkin^11) + 
    (8*apinlmus^3*mcMSmusmc^10*TAGdecmcMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mcMSmusmc^2])/(2835*mkin^9) + 
    (apinlmus^3*mcMSmusmc^8*TAGdecmcMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mcMSmusmc^2])/(140*mkin^7) + 
    (4*apinlmus^3*mcMSmusmc^6*TAGdecmcMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mcMSmusmc^2])/(135*mkin^5) - 
    (13*apinlmus^3*mcMSmusmc^4*TAGdecmcMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mcMSmusmc^2])/(54*mkin^3) + 
    (apinlmus^3*mcMSmusmc^4*TAGdecmcMSOS*TAGnc*Log[mcMSmusmc/mkin]^2*
      Log[mus^2/mcMSmusmc^2])/(9*mkin^3) - 
    (apinlmus^3*mkin*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]^2)/27 - 
    apinlmus*mkin*Log[mus^2/mkin^2] - (461*apinlmus^2*mkin*Log[mus^2/mkin^2])/
     72 - (22273*apinlmus^3*mkin*Log[mus^2/mkin^2])/324 - 
    (16*apinlmus^2*mufac*Log[mus^2/mkin^2])/9 - 
    (2950*apinlmus^3*mufac*Log[mus^2/mkin^2])/81 - 
    (2*apinlmus^2*mufac^2*Log[mus^2/mkin^2])/(3*mkin) - 
    (1277*apinlmus^3*mufac^2*Log[mus^2/mkin^2])/(108*mkin) + 
    (13*apinlmus^2*mkin*NLxMSOS*Log[mus^2/mkin^2])/36 + 
    (1283*apinlmus^3*mkin*NLxMSOS*Log[mus^2/mkin^2])/162 + 
    (4*apinlmus^3*mufac*NLxMSOS*Log[mus^2/mkin^2])/81 + 
    (apinlmus^3*mufac^2*NLxMSOS*Log[mus^2/mkin^2])/(54*mkin) - 
    (89*apinlmus^3*mkin*NLxMSOS^2*Log[mus^2/mkin^2])/648 + 
    (128*apinlmus^3*mufac*NLxOSKIN*Log[mus^2/mkin^2])/81 + 
    (13*apinlmus^3*mufac^2*NLxOSKIN*Log[mus^2/mkin^2])/(27*mkin) - 
    (3*apinlmus^3*mkin*Pi^2*Log[mus^2/mkin^2])/2 + 
    (8*apinlmus^3*mufac*Pi^2*Log[mus^2/mkin^2])/9 + 
    (apinlmus^3*mufac^2*Pi^2*Log[mus^2/mkin^2])/(3*mkin) + 
    (13*apinlmus^3*mkin*NLxMSOS*Pi^2*Log[mus^2/mkin^2])/36 - 
    (apinlmus^3*mkin*NLxMSOS^2*Pi^2*Log[mus^2/mkin^2])/54 + 
    (11*apinlmus^3*mkin*TAGdecmcMSOS*Log[mus^2/mkin^2])/72 - 
    (9727*apinlmus^3*mcMSmusmc^24*TAGnc*Log[mus^2/mkin^2])/
     (41473600*mkin^23) - (39833*apinlmus^3*mcMSmusmc^22*TAGnc*
      Log[mus^2/mkin^2])/(115579926*mkin^21) - 
    (352467*apinlmus^3*mcMSmusmc^20*TAGnc*Log[mus^2/mkin^2])/
     (667705600*mkin^19) - (16277*apinlmus^3*mcMSmusmc^18*TAGnc*
      Log[mus^2/mkin^2])/(19117350*mkin^17) - 
    (1587*apinlmus^3*mcMSmusmc^16*TAGnc*Log[mus^2/mkin^2])/
     (1081600*mkin^15) - (138339*apinlmus^3*mcMSmusmc^14*TAGnc*
      Log[mus^2/mkin^2])/(50100050*mkin^13) - 
    (1229*apinlmus^3*mcMSmusmc^12*TAGnc*Log[mus^2/mkin^2])/(209088*mkin^11) - 
    (997*apinlmus^3*mcMSmusmc^10*TAGnc*Log[mus^2/mkin^2])/(66150*mkin^9) - 
    (4167*apinlmus^3*mcMSmusmc^8*TAGnc*Log[mus^2/mkin^2])/(78400*mkin^7) - 
    (19*apinlmus^3*mcMSmusmc^6*TAGnc*Log[mus^2/mkin^2])/(50*mkin^5) + 
    (151*apinlmus^3*mcMSmusmc^4*TAGnc*Log[mus^2/mkin^2])/(48*mkin^3) + 
    (9*apinlmus^3*mcMSmusmc^2*TAGnc*Log[mus^2/mkin^2])/(2*mkin) + 
    (13*apinlmus^2*mkin*TAGnc*Log[mus^2/mkin^2])/36 + 
    (1283*apinlmus^3*mkin*TAGnc*Log[mus^2/mkin^2])/162 + 
    (4*apinlmus^3*mufac*TAGnc*Log[mus^2/mkin^2])/81 + 
    (apinlmus^3*mufac^2*TAGnc*Log[mus^2/mkin^2])/(54*mkin) + 
    (9727*apinlmus^3*mcMSmusmc^24*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/
     (559893600*mkin^23) + (39833*apinlmus^3*mcMSmusmc^22*NLxMSOS*TAGnc*
      Log[mus^2/mkin^2])/(1560329001*mkin^21) + 
    (39163*apinlmus^3*mcMSmusmc^20*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/
     (1001558400*mkin^19) + (16277*apinlmus^3*mcMSmusmc^18*NLxMSOS*TAGnc*
      Log[mus^2/mkin^2])/(258084225*mkin^17) + 
    (529*apinlmus^3*mcMSmusmc^16*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/
     (4867200*mkin^15) + (15371*apinlmus^3*mcMSmusmc^14*NLxMSOS*TAGnc*
      Log[mus^2/mkin^2])/(75150075*mkin^13) + 
    (1229*apinlmus^3*mcMSmusmc^12*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/
     (2822688*mkin^11) + (997*apinlmus^3*mcMSmusmc^10*NLxMSOS*TAGnc*
      Log[mus^2/mkin^2])/(893025*mkin^9) + 
    (463*apinlmus^3*mcMSmusmc^8*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/
     (117600*mkin^7) + (19*apinlmus^3*mcMSmusmc^6*NLxMSOS*TAGnc*
      Log[mus^2/mkin^2])/(675*mkin^5) - 
    (151*apinlmus^3*mcMSmusmc^4*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/
     (648*mkin^3) - (apinlmus^3*mcMSmusmc^2*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/
     (3*mkin) - (89*apinlmus^3*mkin*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/324 - 
    (3*apinlmus^3*mcMSmusmc*Pi^2*TAGnc*Log[mus^2/mkin^2])/4 + 
    (apinlmus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[mus^2/mkin^2])/(4*mkin^3) - 
    (3*apinlmus^3*mcMSmusmc^3*Pi^2*TAGnc*Log[mus^2/mkin^2])/(4*mkin^2) + 
    (13*apinlmus^3*mkin*Pi^2*TAGnc*Log[mus^2/mkin^2])/36 + 
    (apinlmus^3*mcMSmusmc*NLxMSOS*Pi^2*TAGnc*Log[mus^2/mkin^2])/18 - 
    (apinlmus^3*mcMSmusmc^4*NLxMSOS*Pi^2*TAGnc*Log[mus^2/mkin^2])/
     (54*mkin^3) + (apinlmus^3*mcMSmusmc^3*NLxMSOS*Pi^2*TAGnc*
      Log[mus^2/mkin^2])/(18*mkin^2) - (apinlmus^3*mkin*NLxMSOS*Pi^2*TAGnc*
      Log[mus^2/mkin^2])/27 + (9727*apinlmus^3*mcMSmusmc^24*TAGnc^2*
      Log[mus^2/mkin^2])/(559893600*mkin^23) + 
    (39833*apinlmus^3*mcMSmusmc^22*TAGnc^2*Log[mus^2/mkin^2])/
     (1560329001*mkin^21) + (39163*apinlmus^3*mcMSmusmc^20*TAGnc^2*
      Log[mus^2/mkin^2])/(1001558400*mkin^19) + 
    (16277*apinlmus^3*mcMSmusmc^18*TAGnc^2*Log[mus^2/mkin^2])/
     (258084225*mkin^17) + (529*apinlmus^3*mcMSmusmc^16*TAGnc^2*
      Log[mus^2/mkin^2])/(4867200*mkin^15) + 
    (15371*apinlmus^3*mcMSmusmc^14*TAGnc^2*Log[mus^2/mkin^2])/
     (75150075*mkin^13) + (1229*apinlmus^3*mcMSmusmc^12*TAGnc^2*
      Log[mus^2/mkin^2])/(2822688*mkin^11) + 
    (997*apinlmus^3*mcMSmusmc^10*TAGnc^2*Log[mus^2/mkin^2])/(893025*mkin^9) + 
    (463*apinlmus^3*mcMSmusmc^8*TAGnc^2*Log[mus^2/mkin^2])/(117600*mkin^7) + 
    (19*apinlmus^3*mcMSmusmc^6*TAGnc^2*Log[mus^2/mkin^2])/(675*mkin^5) - 
    (151*apinlmus^3*mcMSmusmc^4*TAGnc^2*Log[mus^2/mkin^2])/(648*mkin^3) - 
    (apinlmus^3*mcMSmusmc^2*TAGnc^2*Log[mus^2/mkin^2])/(3*mkin) - 
    (89*apinlmus^3*mkin*TAGnc^2*Log[mus^2/mkin^2])/648 + 
    (apinlmus^3*mcMSmusmc*Pi^2*TAGnc^2*Log[mus^2/mkin^2])/18 - 
    (apinlmus^3*mcMSmusmc^4*Pi^2*TAGnc^2*Log[mus^2/mkin^2])/(54*mkin^3) + 
    (apinlmus^3*mcMSmusmc^3*Pi^2*TAGnc^2*Log[mus^2/mkin^2])/(18*mkin^2) - 
    (apinlmus^3*mkin*Pi^2*TAGnc^2*Log[mus^2/mkin^2])/54 - 
    (apinlmus^3*mkin*Pi^2*Log[2]*Log[mus^2/mkin^2])/2 + 
    (apinlmus^3*mkin*NLxMSOS*Pi^2*Log[2]*Log[mus^2/mkin^2])/27 + 
    (apinlmus^3*mkin*Pi^2*TAGnc*Log[2]*Log[mus^2/mkin^2])/27 + 
    (11*apinlmus^3*mcMSmusmc^24*TAGnc*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (6440*mkin^23) + (10*apinlmus^3*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(4389*mkin^21) + 
    (81*apinlmus^3*mcMSmusmc^20*TAGnc*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (25840*mkin^19) + (8*apinlmus^3*mcMSmusmc^18*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(1785*mkin^17) + 
    (7*apinlmus^3*mcMSmusmc^16*TAGnc*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (1040*mkin^15) + (54*apinlmus^3*mcMSmusmc^14*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(5005*mkin^13) + 
    (5*apinlmus^3*mcMSmusmc^12*TAGnc*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (264*mkin^11) + (4*apinlmus^3*mcMSmusmc^10*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(105*mkin^9) + 
    (27*apinlmus^3*mcMSmusmc^8*TAGnc*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (280*mkin^7) + (2*apinlmus^3*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(5*mkin^5) - (13*apinlmus^3*mcMSmusmc^4*TAGnc*
      Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/(4*mkin^3) - 
    (11*apinlmus^3*mcMSmusmc^24*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(86940*mkin^23) - 
    (20*apinlmus^3*mcMSmusmc^22*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(118503*mkin^21) - 
    (3*apinlmus^3*mcMSmusmc^20*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(12920*mkin^19) - 
    (16*apinlmus^3*mcMSmusmc^18*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(48195*mkin^17) - 
    (7*apinlmus^3*mcMSmusmc^16*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(14040*mkin^15) - 
    (4*apinlmus^3*mcMSmusmc^14*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(5005*mkin^13) - 
    (5*apinlmus^3*mcMSmusmc^12*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(3564*mkin^11) - 
    (8*apinlmus^3*mcMSmusmc^10*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(2835*mkin^9) - 
    (apinlmus^3*mcMSmusmc^8*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(140*mkin^7) - 
    (4*apinlmus^3*mcMSmusmc^6*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(135*mkin^5) + 
    (13*apinlmus^3*mcMSmusmc^4*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(54*mkin^3) - (11*apinlmus^3*mcMSmusmc^24*TAGnc^2*
      Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/(86940*mkin^23) - 
    (20*apinlmus^3*mcMSmusmc^22*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(118503*mkin^21) - 
    (3*apinlmus^3*mcMSmusmc^20*TAGnc^2*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (12920*mkin^19) - (16*apinlmus^3*mcMSmusmc^18*TAGnc^2*
      Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/(48195*mkin^17) - 
    (7*apinlmus^3*mcMSmusmc^16*TAGnc^2*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (14040*mkin^15) - (4*apinlmus^3*mcMSmusmc^14*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(5005*mkin^13) - 
    (5*apinlmus^3*mcMSmusmc^12*TAGnc^2*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (3564*mkin^11) - (8*apinlmus^3*mcMSmusmc^10*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(2835*mkin^9) - 
    (apinlmus^3*mcMSmusmc^8*TAGnc^2*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (140*mkin^7) - (4*apinlmus^3*mcMSmusmc^6*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(135*mkin^5) + 
    (13*apinlmus^3*mcMSmusmc^4*TAGnc^2*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (54*mkin^3) + (3*apinlmus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]^2*
      Log[mus^2/mkin^2])/(2*mkin^3) - (apinlmus^3*mcMSmusmc^4*NLxMSOS*TAGnc*
      Log[mcMSmusmc/mkin]^2*Log[mus^2/mkin^2])/(9*mkin^3) - 
    (apinlmus^3*mcMSmusmc^4*TAGnc^2*Log[mcMSmusmc/mkin]^2*Log[mus^2/mkin^2])/
     (9*mkin^3) + (88*apinlmus^3*mufac*Log[(2*mufac)/mus]*Log[mus^2/mkin^2])/
     9 + (11*apinlmus^3*mufac^2*Log[(2*mufac)/mus]*Log[mus^2/mkin^2])/
     (3*mkin) - (16*apinlmus^3*mufac*NLxOSKIN*Log[(2*mufac)/mus]*
      Log[mus^2/mkin^2])/27 - (2*apinlmus^3*mufac^2*NLxOSKIN*
      Log[(2*mufac)/mus]*Log[mus^2/mkin^2])/(9*mkin) - 
    (apinlmus^2*mkin*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]*Log[mus^2/mkin^2])/
     6 - (79*apinlmus^3*mkin*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]*
      Log[mus^2/mkin^2])/27 - (8*apinlmus^3*mufac*TAGdecmcMSOS*
      Log[mus^2/mcMSmusmc^2]*Log[mus^2/mkin^2])/27 - 
    (apinlmus^3*mufac^2*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]*
      Log[mus^2/mkin^2])/(9*mkin) + (13*apinlmus^3*mkin*NLxMSOS*TAGdecmcMSOS*
      Log[mus^2/mcMSmusmc^2]*Log[mus^2/mkin^2])/108 + 
    (13*apinlmus^3*mkin*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2]*
      Log[mus^2/mkin^2])/108 - (apinlmus^3*mkin*TAGdecmcMSOS*
      Log[mus^2/mcMSmusmc^2]^2*Log[mus^2/mkin^2])/36 - 
    (23*apinlmus^2*mkin*Log[mus^2/mkin^2]^2)/24 - 
    (14275*apinlmus^3*mkin*Log[mus^2/mkin^2]^2)/864 - 
    (46*apinlmus^3*mufac*Log[mus^2/mkin^2]^2)/27 - 
    (23*apinlmus^3*mufac^2*Log[mus^2/mkin^2]^2)/(36*mkin) + 
    (apinlmus^2*mkin*NLxMSOS*Log[mus^2/mkin^2]^2)/12 + 
    (107*apinlmus^3*mkin*NLxMSOS*Log[mus^2/mkin^2]^2)/48 + 
    (4*apinlmus^3*mufac*NLxMSOS*Log[mus^2/mkin^2]^2)/27 + 
    (apinlmus^3*mufac^2*NLxMSOS*Log[mus^2/mkin^2]^2)/(18*mkin) - 
    (13*apinlmus^3*mkin*NLxMSOS^2*Log[mus^2/mkin^2]^2)/216 + 
    (apinlmus^2*mkin*TAGnc*Log[mus^2/mkin^2]^2)/12 + 
    (107*apinlmus^3*mkin*TAGnc*Log[mus^2/mkin^2]^2)/48 + 
    (4*apinlmus^3*mufac*TAGnc*Log[mus^2/mkin^2]^2)/27 + 
    (apinlmus^3*mufac^2*TAGnc*Log[mus^2/mkin^2]^2)/(18*mkin) - 
    (13*apinlmus^3*mkin*NLxMSOS*TAGnc*Log[mus^2/mkin^2]^2)/108 - 
    (13*apinlmus^3*mkin*TAGnc^2*Log[mus^2/mkin^2]^2)/216 - 
    (23*apinlmus^3*mkin*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]*
      Log[mus^2/mkin^2]^2)/72 + (apinlmus^3*mkin*NLxMSOS*TAGdecmcMSOS*
      Log[mus^2/mcMSmusmc^2]*Log[mus^2/mkin^2]^2)/36 + 
    (apinlmus^3*mkin*TAGdecmcMSOS*TAGnc*Log[mus^2/mcMSmusmc^2]*
      Log[mus^2/mkin^2]^2)/36 - (601*apinlmus^3*mkin*Log[mus^2/mkin^2]^3)/
     432 + (25*apinlmus^3*mkin*NLxMSOS*Log[mus^2/mkin^2]^3)/108 - 
    (apinlmus^3*mkin*NLxMSOS^2*Log[mus^2/mkin^2]^3)/108 + 
    (25*apinlmus^3*mkin*TAGnc*Log[mus^2/mkin^2]^3)/108 - 
    (apinlmus^3*mkin*NLxMSOS*TAGnc*Log[mus^2/mkin^2]^3)/54 - 
    (apinlmus^3*mkin*TAGnc^2*Log[mus^2/mkin^2]^3)/108 + 
    (4*apinlmus^3*mkin*TAGdecmcMSOS*Log[musmc^2/mcMSmusmc^2])/9 - 
    (10163*apinlmus^3*mcMSmusmc^24*TAGnc*Log[musmc^2/mcMSmusmc^2])/
     (11664450*mkin^23) - (5066*apinlmus^3*mcMSmusmc^22*TAGnc*
      Log[musmc^2/mcMSmusmc^2])/(4298427*mkin^21) - 
    (5507*apinlmus^3*mcMSmusmc^20*TAGnc*Log[musmc^2/mcMSmusmc^2])/
     (3338528*mkin^19) - (7678*apinlmus^3*mcMSmusmc^18*TAGnc*
      Log[musmc^2/mcMSmusmc^2])/(3186225*mkin^17) - 
    (283*apinlmus^3*mcMSmusmc^16*TAGnc*Log[musmc^2/mcMSmusmc^2])/
     (76050*mkin^15) - (3166*apinlmus^3*mcMSmusmc^14*TAGnc*
      Log[musmc^2/mcMSmusmc^2])/(511225*mkin^13) - 
    (899*apinlmus^3*mcMSmusmc^12*TAGnc*Log[musmc^2/mcMSmusmc^2])/
     (78408*mkin^11) - (298*apinlmus^3*mcMSmusmc^10*TAGnc*
      Log[musmc^2/mcMSmusmc^2])/(11907*mkin^9) - 
    (179*apinlmus^3*mcMSmusmc^8*TAGnc*Log[musmc^2/mcMSmusmc^2])/
     (2450*mkin^7) - (94*apinlmus^3*mcMSmusmc^6*TAGnc*
      Log[musmc^2/mcMSmusmc^2])/(225*mkin^5) + 
    (56*apinlmus^3*mcMSmusmc^4*TAGnc*Log[musmc^2/mcMSmusmc^2])/(27*mkin^3) + 
    (2*apinlmus^3*mcMSmusmc^2*TAGnc*Log[musmc^2/mcMSmusmc^2])/mkin - 
    (apinlmus^3*mcMSmusmc*Pi^2*TAGnc*Log[musmc^2/mcMSmusmc^2])/6 + 
    (2*apinlmus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[musmc^2/mcMSmusmc^2])/
     (9*mkin^3) - (apinlmus^3*mcMSmusmc^3*Pi^2*TAGnc*
      Log[musmc^2/mcMSmusmc^2])/(2*mkin^2) + 
    (22*apinlmus^3*mcMSmusmc^24*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(2415*mkin^23) + 
    (40*apinlmus^3*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(3591*mkin^21) + 
    (9*apinlmus^3*mcMSmusmc^20*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(646*mkin^19) + 
    (32*apinlmus^3*mcMSmusmc^18*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(1785*mkin^17) + 
    (14*apinlmus^3*mcMSmusmc^16*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(585*mkin^15) + 
    (24*apinlmus^3*mcMSmusmc^14*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(715*mkin^13) + 
    (5*apinlmus^3*mcMSmusmc^12*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(99*mkin^11) + 
    (16*apinlmus^3*mcMSmusmc^10*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(189*mkin^9) + 
    (6*apinlmus^3*mcMSmusmc^8*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(35*mkin^7) + 
    (8*apinlmus^3*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(15*mkin^5) - 
    (20*apinlmus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(9*mkin^3) + 
    (4*apinlmus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]^2*
      Log[musmc^2/mcMSmusmc^2])/(3*mkin^3) + 
    (apinlmus^3*mkin*TAGdecmcMSOS*Log[mus^2/mkin^2]*Log[musmc^2/mcMSmusmc^2])/
     3 + (220*apinlmus^3*mkin*PolyLog[4, 1/2])/27 - 
    (8*apinlmus^3*mkin*NLxMSOS*PolyLog[4, 1/2])/27 - 
    (20*apinlmus^3*mcMSmusmc^4*TAGnc*PolyLog[4, 1/2])/(27*mkin^3) - 
    (8*apinlmus^3*mkin*TAGnc*PolyLog[4, 1/2])/27 + 
    (apinlmus^2*mkin*Zeta[3])/6 - (61*apinlmus^3*mkin*Zeta[3])/27 - 
    (3070*apinlmus^3*mufac*Zeta[3])/27 - (1535*apinlmus^3*mufac^2*Zeta[3])/
     (36*mkin) + (241*apinlmus^3*mkin*NLxMSOS*Zeta[3])/72 - 
    (7*apinlmus^3*mkin*NLxMSOS^2*Zeta[3])/54 + 
    (140*apinlmus^3*mufac*NLxOSKIN*Zeta[3])/27 + 
    (35*apinlmus^3*mufac^2*NLxOSKIN*Zeta[3])/(18*mkin) + 
    (1439*apinlmus^3*mkin*Pi^2*Zeta[3])/432 + 
    (2710689767*apinlmus^3*mcMSmusmc^24*TAGnc*Zeta[3])/
     (1644334940160*mkin^23) + (23017987*apinlmus^3*mcMSmusmc^22*TAGnc*
      Zeta[3])/(11678515200*mkin^21) + (254791*apinlmus^3*mcMSmusmc^20*TAGnc*
      Zeta[3])/(106168320*mkin^19) + (83291*apinlmus^3*mcMSmusmc^18*TAGnc*
      Zeta[3])/(27869184*mkin^17) + (885457*apinlmus^3*mcMSmusmc^16*TAGnc*
      Zeta[3])/(231211008*mkin^15) + (3287*apinlmus^3*mcMSmusmc^14*TAGnc*
      Zeta[3])/(645120*mkin^13) + (59231*apinlmus^3*mcMSmusmc^12*TAGnc*
      Zeta[3])/(8294400*mkin^11) + (187*apinlmus^3*mcMSmusmc^10*TAGnc*
      Zeta[3])/(17280*mkin^9) + (173*apinlmus^3*mcMSmusmc^8*TAGnc*Zeta[3])/
     (9216*mkin^7) + (29*apinlmus^3*mcMSmusmc^6*TAGnc*Zeta[3])/(648*mkin^5) + 
    (2309*apinlmus^3*mcMSmusmc^4*TAGnc*Zeta[3])/(864*mkin^3) + 
    (11*apinlmus^3*mcMSmusmc^2*TAGnc*Zeta[3])/(2*mkin) + 
    (241*apinlmus^3*mkin*TAGnc*Zeta[3])/72 - 
    (apinlmus^3*mcMSmusmc^4*NLxMSOS*TAGnc*Zeta[3])/(3*mkin^3) - 
    (7*apinlmus^3*mkin*NLxMSOS*TAGnc*Zeta[3])/27 - 
    (3*apinlmus^3*mcMSmusmc^2*Pi^2*TAGnc*Zeta[3])/(4*mkin) - 
    (7*apinlmus^3*mkin*TAGnc^2*Zeta[3])/54 - 
    (14*apinlmus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]*Zeta[3])/
     (9*mkin^3) + (apinlmus^3*mkin*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]*
      Zeta[3])/18 + (19*apinlmus^3*mkin*Log[mus^2/mkin^2]*Zeta[3])/12 + 
    (7*apinlmus^3*mkin*NLxMSOS*Log[mus^2/mkin^2]*Zeta[3])/9 + 
    (7*apinlmus^3*mkin*TAGnc*Log[mus^2/mkin^2]*Zeta[3])/9 - 
    (1975*apinlmus^3*mkin*Zeta[5])/216 - 
    (5*apinlmus^3*mcMSmusmc^2*TAGnc*Zeta[5])/(2*mkin)

RunDecmkin2mMSmc = mkin - (4*apinl1mus*mkin)/3 - (3019*apinl1mus^2*mkin)/
     288 - (9514621*apinl1mus^3*mkin)/93312 + (16*apinl1mus*mufac)/9 + 
    (892*apinl1mus^2*mufac)/27 + (22496*apinl1mus^3*mufac)/27 + 
    (2*apinl1mus*mufac^2)/(3*mkin) + (95*apinl1mus^2*mufac^2)/(9*mkin) + 
    (75931*apinl1mus^3*mufac^2)/(324*mkin) + (71*apinl1mus^2*mkin*NLxMSOS)/
     144 + (246643*apinl1mus^3*mkin*NLxMSOS)/23328 - 
    (11*apinl1mus^3*mufac*NLxMSOS)/27 - (11*apinl1mus^3*mufac^2*NLxMSOS)/
     (72*mkin) - (2353*apinl1mus^3*mkin*NLxMSOS^2)/23328 - 
    (128*apinl1mus^2*mufac*NLxOSKIN)/81 - (20303*apinl1mus^3*mufac*NLxOSKIN)/
     243 - (13*apinl1mus^2*mufac^2*NLxOSKIN)/(27*mkin) - 
    (14429*apinl1mus^3*mufac^2*NLxOSKIN)/(648*mkin) + 
    (1292*apinl1mus^3*mufac*NLxOSKIN^2)/729 + 
    (209*apinl1mus^3*mufac^2*NLxOSKIN^2)/(486*mkin) - 
    (apinl1mus^2*mkin*Pi^2)/3 - (644201*apinl1mus^3*mkin*Pi^2)/38880 - 
    (8*apinl1mus^2*mufac*Pi^2)/9 - (1054*apinl1mus^3*mufac*Pi^2)/27 - 
    (apinl1mus^2*mufac^2*Pi^2)/(3*mkin) - (461*apinl1mus^3*mufac^2*Pi^2)/
     (36*mkin) + (apinl1mus^2*mkin*NLxMSOS*Pi^2)/18 + 
    (967*apinl1mus^3*mkin*NLxMSOS*Pi^2)/648 + 
    (8*apinl1mus^3*mufac*NLxMSOS*Pi^2)/81 + 
    (apinl1mus^3*mufac^2*NLxMSOS*Pi^2)/(27*mkin) - 
    (13*apinl1mus^3*mkin*NLxMSOS^2*Pi^2)/324 + 
    (208*apinl1mus^3*mufac*NLxOSKIN*Pi^2)/81 + 
    (23*apinl1mus^3*mufac^2*NLxOSKIN*Pi^2)/(27*mkin) - 
    (8*apinl1mus^3*mufac*NLxOSKIN^2*Pi^2)/243 - 
    (apinl1mus^3*mufac^2*NLxOSKIN^2*Pi^2)/(81*mkin) + 
    (695*apinl1mus^3*mkin*Pi^4)/7776 + (2*apinl1mus^3*mufac*Pi^4)/3 + 
    (apinl1mus^3*mufac^2*Pi^4)/(4*mkin) - (61*apinl1mus^3*mkin*NLxMSOS*Pi^4)/
     1944 + (22*apinl1mus^3*mufac*TAGkinmc)/81 + 
    (11*apinl1mus^3*mufac^2*TAGkinmc)/(108*mkin) - 
    (9727*apinl1mus^2*mcMSmusmc^24*TAGnc)/(186631200*mkin^23) - 
    (8583265857182497000029852773317049*apinl1mus^3*mcMSmusmc^24*TAGnc)/
     (2369996943791663839368455777157120000*mkin^23) - 
    (39833*apinl1mus^2*mcMSmusmc^22*TAGnc)/(520109667*mkin^21) - 
    (48833406986106441810878913193*apinl1mus^3*mcMSmusmc^22*TAGnc)/
     (10790864167779053267146997760000*mkin^21) - 
    (39163*apinl1mus^2*mcMSmusmc^20*TAGnc)/(333852800*mkin^19) - 
    (316346323888931757384696529*apinl1mus^3*mcMSmusmc^20*TAGnc)/
     (54499313978682087207813120000*mkin^19) - 
    (16277*apinl1mus^2*mcMSmusmc^18*TAGnc)/(86028075*mkin^17) - 
    (385021249205567851684429*apinl1mus^3*mcMSmusmc^18*TAGnc)/
     (50057687427569200963584000*mkin^17) - 
    (529*apinl1mus^2*mcMSmusmc^16*TAGnc)/(1622400*mkin^15) - 
    (1009988661120425340647*apinl1mus^3*mcMSmusmc^16*TAGnc)/
     (95095681387606867968000*mkin^15) - 
    (15371*apinl1mus^2*mcMSmusmc^14*TAGnc)/(25050025*mkin^13) - 
    (911686115452983709*apinl1mus^3*mcMSmusmc^14*TAGnc)/
     (58963096098466560000*mkin^13) - (1229*apinl1mus^2*mcMSmusmc^12*TAGnc)/
     (940896*mkin^11) - (593235847648457*apinl1mus^3*mcMSmusmc^12*TAGnc)/
     (24647147078400000*mkin^11) - (997*apinl1mus^2*mcMSmusmc^10*TAGnc)/
     (297675*mkin^9) - (208601138911*apinl1mus^3*mcMSmusmc^10*TAGnc)/
     (5184974592000*mkin^9) - (463*apinl1mus^2*mcMSmusmc^8*TAGnc)/
     (39200*mkin^7) - (72581279*apinl1mus^3*mcMSmusmc^8*TAGnc)/
     (1185408000*mkin^7) - (19*apinl1mus^2*mcMSmusmc^6*TAGnc)/(225*mkin^5) + 
    (1777621*apinl1mus^3*mcMSmusmc^6*TAGnc)/(3499200*mkin^5) + 
    (151*apinl1mus^2*mcMSmusmc^4*TAGnc)/(216*mkin^3) + 
    (3211*apinl1mus^3*mcMSmusmc^4*TAGnc)/(2592*mkin^3) + 
    (apinl1mus^2*mcMSmusmc^2*TAGnc)/mkin + (88*apinl1mus^3*mcMSmusmc^2*TAGnc)/
     (9*mkin) + (71*apinl1mus^2*mkin*TAGnc)/144 + 
    (246643*apinl1mus^3*mkin*TAGnc)/23328 - (11*apinl1mus^3*mufac*TAGnc)/27 + 
    (289*apinl1mus^3*mcMSmusmc^24*mufac*TAGnc)/(198450*mkin^24) + 
    (62384*apinl1mus^3*mcMSmusmc^22*mufac*TAGnc)/(31843449*mkin^22) + 
    (1417*apinl1mus^3*mcMSmusmc^20*mufac*TAGnc)/(520200*mkin^20) + 
    (10576*apinl1mus^3*mcMSmusmc^18*mufac*TAGnc)/(2679075*mkin^18) + 
    (661*apinl1mus^3*mcMSmusmc^16*mufac*TAGnc)/(109512*mkin^16) + 
    (13232*apinl1mus^3*mcMSmusmc^14*mufac*TAGnc)/(1334025*mkin^14) + 
    (79*apinl1mus^3*mcMSmusmc^12*mufac*TAGnc)/(4374*mkin^12) + 
    (3824*apinl1mus^3*mcMSmusmc^10*mufac*TAGnc)/(99225*mkin^10) + 
    (49*apinl1mus^3*mcMSmusmc^8*mufac*TAGnc)/(450*mkin^8) + 
    (16*apinl1mus^3*mcMSmusmc^6*mufac*TAGnc)/(27*mkin^6) - 
    (22*apinl1mus^3*mcMSmusmc^4*mufac*TAGnc)/(9*mkin^4) - 
    (16*apinl1mus^3*mcMSmusmc^2*mufac*TAGnc)/(9*mkin^2) + 
    (289*apinl1mus^3*mcMSmusmc^24*mufac^2*TAGnc)/(529200*mkin^25) + 
    (7798*apinl1mus^3*mcMSmusmc^22*mufac^2*TAGnc)/(10614483*mkin^23) + 
    (1417*apinl1mus^3*mcMSmusmc^20*mufac^2*TAGnc)/(1387200*mkin^21) + 
    (1322*apinl1mus^3*mcMSmusmc^18*mufac^2*TAGnc)/(893025*mkin^19) + 
    (661*apinl1mus^3*mcMSmusmc^16*mufac^2*TAGnc)/(292032*mkin^17) + 
    (1654*apinl1mus^3*mcMSmusmc^14*mufac^2*TAGnc)/(444675*mkin^15) + 
    (79*apinl1mus^3*mcMSmusmc^12*mufac^2*TAGnc)/(11664*mkin^13) + 
    (478*apinl1mus^3*mcMSmusmc^10*mufac^2*TAGnc)/(33075*mkin^11) + 
    (49*apinl1mus^3*mcMSmusmc^8*mufac^2*TAGnc)/(1200*mkin^9) + 
    (2*apinl1mus^3*mcMSmusmc^6*mufac^2*TAGnc)/(9*mkin^7) - 
    (11*apinl1mus^3*mcMSmusmc^4*mufac^2*TAGnc)/(12*mkin^5) - 
    (2*apinl1mus^3*mcMSmusmc^2*mufac^2*TAGnc)/(3*mkin^3) - 
    (11*apinl1mus^3*mufac^2*TAGnc)/(72*mkin) - 
    (7801530877413386647*apinl1mus^3*mcMSmusmc^24*NLxMSOS*TAGnc)/
     (37763267488425775296000*mkin^23) - 
    (1640519393726677*apinl1mus^3*mcMSmusmc^22*NLxMSOS*TAGnc)/
     (5946258455651273760*mkin^21) - (22189567531163017*apinl1mus^3*
      mcMSmusmc^20*NLxMSOS*TAGnc)/(58347124817357376000*mkin^19) - 
    (20555048260909*apinl1mus^3*mcMSmusmc^18*NLxMSOS*TAGnc)/
     (37681809223558500*mkin^17) - (214558103603*apinl1mus^3*mcMSmusmc^16*
      NLxMSOS*TAGnc)/(260460712512000*mkin^15) - 
    (108352091581*apinl1mus^3*mcMSmusmc^14*NLxMSOS*TAGnc)/
     (81243243081000*mkin^13) - (326802499*apinl1mus^3*mcMSmusmc^12*NLxMSOS*
      TAGnc)/(136928594880*mkin^11) - (2816347*apinl1mus^3*mcMSmusmc^10*
      NLxMSOS*TAGnc)/(562605750*mkin^9) - 
    (174787*apinl1mus^3*mcMSmusmc^8*NLxMSOS*TAGnc)/(12348000*mkin^7) - 
    (2729*apinl1mus^3*mcMSmusmc^6*NLxMSOS*TAGnc)/(30375*mkin^5) + 
    (1423*apinl1mus^3*mcMSmusmc^4*NLxMSOS*TAGnc)/(3888*mkin^3) - 
    (2*apinl1mus^3*mcMSmusmc^2*NLxMSOS*TAGnc)/(9*mkin) - 
    (2353*apinl1mus^3*mkin*NLxMSOS*TAGnc)/11664 - 
    (apinl1mus^2*mcMSmusmc*Pi^2*TAGnc)/6 - 
    (6025*apinl1mus^3*mcMSmusmc*Pi^2*TAGnc)/486 - 
    (106559417*apinl1mus^3*mcMSmusmc^25*Pi^2*TAGnc)/(637438863600*mkin^24) - 
    (240817793781176357*apinl1mus^3*mcMSmusmc^24*Pi^2*TAGnc)/
     (628867544642695987200*mkin^23) - (171475369*apinl1mus^3*mcMSmusmc^23*
      Pi^2*TAGnc)/(778168119960*mkin^22) - 
    (174200135864459*apinl1mus^3*mcMSmusmc^22*Pi^2*TAGnc)/
     (349543472195174400*mkin^21) - (846435761*apinl1mus^3*mcMSmusmc^21*Pi^2*
      TAGnc)/(2832319518840*mkin^20) - 
    (9007367733163*apinl1mus^3*mcMSmusmc^20*Pi^2*TAGnc)/
     (13481015456563200*mkin^19) - (86836957*apinl1mus^3*mcMSmusmc^19*Pi^2*
      TAGnc)/(206871887520*mkin^18) - (26463251891*apinl1mus^3*mcMSmusmc^18*
      Pi^2*TAGnc)/(28454994247680*mkin^17) - 
    (22757641*apinl1mus^3*mcMSmusmc^17*Pi^2*TAGnc)/(36923796000*mkin^16) - 
    (58806560951*apinl1mus^3*mcMSmusmc^16*Pi^2*TAGnc)/
     (43261891706880*mkin^15) - (17346493*apinl1mus^3*mcMSmusmc^15*Pi^2*
      TAGnc)/(18087549480*mkin^14) - (82285201*apinl1mus^3*mcMSmusmc^14*Pi^2*
      TAGnc)/(38745907200*mkin^13) - (31786481*apinl1mus^3*mcMSmusmc^13*Pi^2*
      TAGnc)/(19677663720*mkin^12) - (84041429*apinl1mus^3*mcMSmusmc^12*Pi^2*
      TAGnc)/(22992076800*mkin^11) - (1055689*apinl1mus^3*mcMSmusmc^11*Pi^2*
      TAGnc)/(345779280*mkin^10) - (359801*apinl1mus^3*mcMSmusmc^10*Pi^2*
      TAGnc)/(48988800*mkin^9) - (416029*apinl1mus^3*mcMSmusmc^9*Pi^2*TAGnc)/
     (60011280*mkin^8) - (1789*apinl1mus^3*mcMSmusmc^8*Pi^2*TAGnc)/
     (90720*mkin^7) - (262769*apinl1mus^3*mcMSmusmc^7*Pi^2*TAGnc)/
     (11907000*mkin^6) - (1427*apinl1mus^3*mcMSmusmc^6*Pi^2*TAGnc)/
     (12960*mkin^5) - (149*apinl1mus^3*mcMSmusmc^5*Pi^2*TAGnc)/
     (1080*mkin^4) + (apinl1mus^2*mcMSmusmc^4*Pi^2*TAGnc)/(18*mkin^3) - 
    (10*apinl1mus^3*mcMSmusmc^4*Pi^2*TAGnc)/(81*mkin^3) - 
    (apinl1mus^2*mcMSmusmc^3*Pi^2*TAGnc)/(6*mkin^2) - 
    (19507*apinl1mus^3*mcMSmusmc^3*Pi^2*TAGnc)/(1620*mkin^2) + 
    (13*apinl1mus^3*mcMSmusmc^2*Pi^2*TAGnc)/(12*mkin) + 
    (apinl1mus^2*mkin*Pi^2*TAGnc)/18 + (967*apinl1mus^3*mkin*Pi^2*TAGnc)/
     648 + (8*apinl1mus^3*mufac*Pi^2*TAGnc)/81 - 
    (8*apinl1mus^3*mcMSmusmc^4*mufac*Pi^2*TAGnc)/(27*mkin^4) + 
    (16*apinl1mus^3*mcMSmusmc^3*mufac*Pi^2*TAGnc)/(27*mkin^3) - 
    (apinl1mus^3*mcMSmusmc^4*mufac^2*Pi^2*TAGnc)/(9*mkin^5) + 
    (2*apinl1mus^3*mcMSmusmc^3*mufac^2*Pi^2*TAGnc)/(9*mkin^4) + 
    (apinl1mus^3*mufac^2*Pi^2*TAGnc)/(27*mkin) + 
    (7*apinl1mus^3*mcMSmusmc*NLxMSOS*Pi^2*TAGnc)/27 + 
    (23*apinl1mus^3*mcMSmusmc^25*NLxMSOS*Pi^2*TAGnc)/(415800*mkin^24) - 
    (11*apinl1mus^3*mcMSmusmc^24*NLxMSOS*Pi^2*TAGnc)/(260820*mkin^23) + 
    (7*apinl1mus^3*mcMSmusmc^23*NLxMSOS*Pi^2*TAGnc)/(96140*mkin^22) - 
    (20*apinl1mus^3*mcMSmusmc^22*NLxMSOS*Pi^2*TAGnc)/(355509*mkin^21) + 
    (19*apinl1mus^3*mcMSmusmc^21*NLxMSOS*Pi^2*TAGnc)/(192780*mkin^20) - 
    (apinl1mus^3*mcMSmusmc^20*NLxMSOS*Pi^2*TAGnc)/(12920*mkin^19) + 
    (17*apinl1mus^3*mcMSmusmc^19*NLxMSOS*Pi^2*TAGnc)/(123120*mkin^18) - 
    (16*apinl1mus^3*mcMSmusmc^18*NLxMSOS*Pi^2*TAGnc)/(144585*mkin^17) + 
    (5*apinl1mus^3*mcMSmusmc^17*NLxMSOS*Pi^2*TAGnc)/(24752*mkin^16) - 
    (7*apinl1mus^3*mcMSmusmc^16*NLxMSOS*Pi^2*TAGnc)/(42120*mkin^15) + 
    (13*apinl1mus^3*mcMSmusmc^15*NLxMSOS*Pi^2*TAGnc)/(41580*mkin^14) - 
    (4*apinl1mus^3*mcMSmusmc^14*NLxMSOS*Pi^2*TAGnc)/(15015*mkin^13) + 
    (11*apinl1mus^3*mcMSmusmc^13*NLxMSOS*Pi^2*TAGnc)/(21060*mkin^12) - 
    (5*apinl1mus^3*mcMSmusmc^12*NLxMSOS*Pi^2*TAGnc)/(10692*mkin^11) + 
    (3*apinl1mus^3*mcMSmusmc^11*NLxMSOS*Pi^2*TAGnc)/(3080*mkin^10) - 
    (8*apinl1mus^3*mcMSmusmc^10*NLxMSOS*Pi^2*TAGnc)/(8505*mkin^9) + 
    (7*apinl1mus^3*mcMSmusmc^9*NLxMSOS*Pi^2*TAGnc)/(3240*mkin^8) - 
    (apinl1mus^3*mcMSmusmc^8*NLxMSOS*Pi^2*TAGnc)/(420*mkin^7) + 
    (5*apinl1mus^3*mcMSmusmc^7*NLxMSOS*Pi^2*TAGnc)/(756*mkin^6) - 
    (4*apinl1mus^3*mcMSmusmc^6*NLxMSOS*Pi^2*TAGnc)/(405*mkin^5) + 
    (apinl1mus^3*mcMSmusmc^5*NLxMSOS*Pi^2*TAGnc)/(20*mkin^4) + 
    (13*apinl1mus^3*mcMSmusmc^4*NLxMSOS*Pi^2*TAGnc)/(324*mkin^3) + 
    (7*apinl1mus^3*mcMSmusmc^3*NLxMSOS*Pi^2*TAGnc)/(54*mkin^2) - 
    (13*apinl1mus^3*mkin*NLxMSOS*Pi^2*TAGnc)/162 - 
    (13*apinl1mus^3*mcMSmusmc*Pi^3*TAGnc)/162 + 
    (5698043*apinl1mus^3*mcMSmusmc^25*Pi^3*TAGnc)/(118908518400*mkin^24) + 
    (81991*apinl1mus^3*mcMSmusmc^23*Pi^3*TAGnc)/(1374683136*mkin^22) + 
    (127699*apinl1mus^3*mcMSmusmc^21*Pi^3*TAGnc)/(1684537344*mkin^20) + 
    (99671*apinl1mus^3*mcMSmusmc^19*Pi^3*TAGnc)/(1008599040*mkin^18) + 
    (1925*apinl1mus^3*mcMSmusmc^17*Pi^3*TAGnc)/(14483456*mkin^16) + 
    (377*apinl1mus^3*mcMSmusmc^15*Pi^3*TAGnc)/(2027520*mkin^14) + 
    (1771*apinl1mus^3*mcMSmusmc^13*Pi^3*TAGnc)/(6469632*mkin^12) + 
    (17*apinl1mus^3*mcMSmusmc^11*Pi^3*TAGnc)/(39424*mkin^10) + 
    (77*apinl1mus^3*mcMSmusmc^9*Pi^3*TAGnc)/(103680*mkin^8) + 
    (25*apinl1mus^3*mcMSmusmc^7*Pi^3*TAGnc)/(18144*mkin^6) - 
    (apinl1mus^3*mcMSmusmc^5*Pi^3*TAGnc)/(240*mkin^4) - 
    (7*apinl1mus^3*mcMSmusmc^3*Pi^3*TAGnc)/(108*mkin^2) - 
    (271*apinl1mus^3*mcMSmusmc^4*Pi^4*TAGnc)/(19440*mkin^3) - 
    (61*apinl1mus^3*mkin*Pi^4*TAGnc)/1944 + 
    (279881*apinl1mus^3*mcMSmusmc^24*TAGnc^2)/(4627692000*mkin^23) + 
    (624853*apinl1mus^3*mcMSmusmc^22*TAGnc^2)/(6846429744*mkin^21) + 
    (367909*apinl1mus^3*mcMSmusmc^20*TAGnc^2)/(2547216000*mkin^19) + 
    (157882*apinl1mus^3*mcMSmusmc^18*TAGnc^2)/(650372247*mkin^17) + 
    (1373*apinl1mus^3*mcMSmusmc^16*TAGnc^2)/(3110400*mkin^15) + 
    (233201*apinl1mus^3*mcMSmusmc^14*TAGnc^2)/(260851500*mkin^13) + 
    (3989*apinl1mus^3*mcMSmusmc^12*TAGnc^2)/(1881792*mkin^11) + 
    (1778*apinl1mus^3*mcMSmusmc^10*TAGnc^2)/(273375*mkin^9) + 
    (241*apinl1mus^3*mcMSmusmc^8*TAGnc^2)/(7056*mkin^7) + 
    (3481*apinl1mus^3*mcMSmusmc^6*TAGnc^2)/(30375*mkin^5) + 
    (235*apinl1mus^3*mcMSmusmc^4*TAGnc^2)/(3888*mkin^3) - 
    (2*apinl1mus^3*mcMSmusmc^2*TAGnc^2)/(9*mkin) - 
    (2353*apinl1mus^3*mkin*TAGnc^2)/23328 + 
    (4*apinl1mus^3*mcMSmusmc*Pi^2*TAGnc^2)/135 - 
    (4*apinl1mus^3*mcMSmusmc^6*Pi^2*TAGnc^2)/(405*mkin^5) + 
    (13*apinl1mus^3*mcMSmusmc^4*Pi^2*TAGnc^2)/(324*mkin^3) - 
    (13*apinl1mus^3*mkin*Pi^2*TAGnc^2)/324 - (apinl1mus^2*mkin*Pi^2*Log[2])/
     9 + (587*apinl1mus^3*mkin*Pi^2*Log[2])/162 - 
    (16*apinl1mus^3*mufac*Pi^2*Log[2])/81 - 
    (2*apinl1mus^3*mufac^2*Pi^2*Log[2])/(27*mkin) + 
    (11*apinl1mus^3*mkin*NLxMSOS*Pi^2*Log[2])/81 + 
    (1199*apinl1mus^3*mcMSmusmc*Pi^2*TAGnc*Log[2])/81 + 
    (877591*apinl1mus^3*mcMSmusmc^24*Pi^2*TAGnc*Log[2])/
     (6794772480*mkin^23) + (38675*apinl1mus^3*mcMSmusmc^22*Pi^2*TAGnc*
      Log[2])/(233570304*mkin^21) + (143*apinl1mus^3*mcMSmusmc^20*Pi^2*TAGnc*
      Log[2])/(655360*mkin^19) + (4147*apinl1mus^3*mcMSmusmc^18*Pi^2*TAGnc*
      Log[2])/(13934592*mkin^17) + (1001*apinl1mus^3*mcMSmusmc^16*Pi^2*TAGnc*
      Log[2])/(2359296*mkin^15) + (23*apinl1mus^3*mcMSmusmc^14*Pi^2*TAGnc*
      Log[2])/(35840*mkin^13) + (175*apinl1mus^3*mcMSmusmc^12*Pi^2*TAGnc*
      Log[2])/(165888*mkin^11) + (17*apinl1mus^3*mcMSmusmc^10*Pi^2*TAGnc*
      Log[2])/(8640*mkin^9) + (7*apinl1mus^3*mcMSmusmc^8*Pi^2*TAGnc*Log[2])/
     (1536*mkin^7) + (11*apinl1mus^3*mcMSmusmc^6*Pi^2*TAGnc*Log[2])/
     (648*mkin^5) - (5*apinl1mus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[2])/
     (144*mkin^3) + (1199*apinl1mus^3*mcMSmusmc^3*Pi^2*TAGnc*Log[2])/
     (81*mkin^2) + (11*apinl1mus^3*mkin*Pi^2*TAGnc*Log[2])/81 - 
    (2*apinl1mus^3*mcMSmusmc*NLxMSOS*Pi^2*TAGnc*Log[2])/9 - 
    (2*apinl1mus^3*mcMSmusmc^3*NLxMSOS*Pi^2*TAGnc*Log[2])/(9*mkin^2) + 
    (22*apinl1mus^3*mkin*Pi^2*Log[2]^2)/81 - 
    (2*apinl1mus^3*mkin*NLxMSOS*Pi^2*Log[2]^2)/81 - 
    (2*apinl1mus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[2]^2)/(81*mkin^3) - 
    (2*apinl1mus^3*mkin*Pi^2*TAGnc*Log[2]^2)/81 + 
    (55*apinl1mus^3*mkin*Log[2]^4)/162 - (apinl1mus^3*mkin*NLxMSOS*Log[2]^4)/
     81 - (5*apinl1mus^3*mcMSmusmc^4*TAGnc*Log[2]^4)/(162*mkin^3) - 
    (apinl1mus^3*mkin*TAGnc*Log[2]^4)/81 + 
    (11*apinl1mus^2*mcMSmusmc^24*TAGnc*Log[mcMSmusmc/mkin])/(28980*mkin^23) + 
    (194881611047644800407*apinl1mus^3*mcMSmusmc^24*TAGnc*
      Log[mcMSmusmc/mkin])/(20069484527233911369600*mkin^23) + 
    (20*apinl1mus^2*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mkin])/(39501*mkin^21) + 
    (84250035134168729123*apinl1mus^3*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mkin])/
     (7242811051244408640000*mkin^21) + 
    (9*apinl1mus^2*mcMSmusmc^20*TAGnc*Log[mcMSmusmc/mkin])/(12920*mkin^19) + 
    (2488663198975906531*apinl1mus^3*mcMSmusmc^20*TAGnc*Log[mcMSmusmc/mkin])/
     (175583298211985664000*mkin^19) + (16*apinl1mus^2*mcMSmusmc^18*TAGnc*
      Log[mcMSmusmc/mkin])/(16065*mkin^17) + 
    (287916747607903391*apinl1mus^3*mcMSmusmc^18*TAGnc*Log[mcMSmusmc/mkin])/
     (16342379002556006400*mkin^17) + (7*apinl1mus^2*mcMSmusmc^16*TAGnc*
      Log[mcMSmusmc/mkin])/(4680*mkin^15) + 
    (26396860982549*apinl1mus^3*mcMSmusmc^16*TAGnc*Log[mcMSmusmc/mkin])/
     (1178083838131200*mkin^15) + (12*apinl1mus^2*mcMSmusmc^14*TAGnc*
      Log[mcMSmusmc/mkin])/(5005*mkin^13) + 
    (1194668026861*apinl1mus^3*mcMSmusmc^14*TAGnc*Log[mcMSmusmc/mkin])/
     (40905688824000*mkin^13) + (5*apinl1mus^2*mcMSmusmc^12*TAGnc*
      Log[mcMSmusmc/mkin])/(1188*mkin^11) + 
    (241232375387*apinl1mus^3*mcMSmusmc^12*TAGnc*Log[mcMSmusmc/mkin])/
     (6224027040000*mkin^11) + (8*apinl1mus^2*mcMSmusmc^10*TAGnc*
      Log[mcMSmusmc/mkin])/(945*mkin^9) + 
    (2042394191*apinl1mus^3*mcMSmusmc^10*TAGnc*Log[mcMSmusmc/mkin])/
     (41150592000*mkin^9) + (3*apinl1mus^2*mcMSmusmc^8*TAGnc*
      Log[mcMSmusmc/mkin])/(140*mkin^7) + 
    (413297*apinl1mus^3*mcMSmusmc^8*TAGnc*Log[mcMSmusmc/mkin])/
     (12700800*mkin^7) + (4*apinl1mus^2*mcMSmusmc^6*TAGnc*
      Log[mcMSmusmc/mkin])/(45*mkin^5) - 
    (721229*apinl1mus^3*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mkin])/
     (1166400*mkin^5) - (13*apinl1mus^2*mcMSmusmc^4*TAGnc*
      Log[mcMSmusmc/mkin])/(18*mkin^3) - 
    (4321*apinl1mus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin])/(648*mkin^3) + 
    (4*apinl1mus^3*mcMSmusmc^2*TAGnc*Log[mcMSmusmc/mkin])/mkin - 
    (44*apinl1mus^3*mcMSmusmc^24*mufac*TAGnc*Log[mcMSmusmc/mkin])/
     (2835*mkin^24) - (320*apinl1mus^3*mcMSmusmc^22*mufac*TAGnc*
      Log[mcMSmusmc/mkin])/(16929*mkin^22) - 
    (2*apinl1mus^3*mcMSmusmc^20*mufac*TAGnc*Log[mcMSmusmc/mkin])/
     (85*mkin^20) - (256*apinl1mus^3*mcMSmusmc^18*mufac*TAGnc*
      Log[mcMSmusmc/mkin])/(8505*mkin^18) - 
    (14*apinl1mus^3*mcMSmusmc^16*mufac*TAGnc*Log[mcMSmusmc/mkin])/
     (351*mkin^16) - (64*apinl1mus^3*mcMSmusmc^14*mufac*TAGnc*
      Log[mcMSmusmc/mkin])/(1155*mkin^14) - 
    (20*apinl1mus^3*mcMSmusmc^12*mufac*TAGnc*Log[mcMSmusmc/mkin])/
     (243*mkin^12) - (128*apinl1mus^3*mcMSmusmc^10*mufac*TAGnc*
      Log[mcMSmusmc/mkin])/(945*mkin^10) - 
    (4*apinl1mus^3*mcMSmusmc^8*mufac*TAGnc*Log[mcMSmusmc/mkin])/(15*mkin^8) - 
    (64*apinl1mus^3*mcMSmusmc^6*mufac*TAGnc*Log[mcMSmusmc/mkin])/
     (81*mkin^6) + (8*apinl1mus^3*mcMSmusmc^4*mufac*TAGnc*
      Log[mcMSmusmc/mkin])/(3*mkin^4) - (11*apinl1mus^3*mcMSmusmc^24*mufac^2*
      TAGnc*Log[mcMSmusmc/mkin])/(1890*mkin^25) - 
    (40*apinl1mus^3*mcMSmusmc^22*mufac^2*TAGnc*Log[mcMSmusmc/mkin])/
     (5643*mkin^23) - (3*apinl1mus^3*mcMSmusmc^20*mufac^2*TAGnc*
      Log[mcMSmusmc/mkin])/(340*mkin^21) - 
    (32*apinl1mus^3*mcMSmusmc^18*mufac^2*TAGnc*Log[mcMSmusmc/mkin])/
     (2835*mkin^19) - (7*apinl1mus^3*mcMSmusmc^16*mufac^2*TAGnc*
      Log[mcMSmusmc/mkin])/(468*mkin^17) - 
    (8*apinl1mus^3*mcMSmusmc^14*mufac^2*TAGnc*Log[mcMSmusmc/mkin])/
     (385*mkin^15) - (5*apinl1mus^3*mcMSmusmc^12*mufac^2*TAGnc*
      Log[mcMSmusmc/mkin])/(162*mkin^13) - 
    (16*apinl1mus^3*mcMSmusmc^10*mufac^2*TAGnc*Log[mcMSmusmc/mkin])/
     (315*mkin^11) - (apinl1mus^3*mcMSmusmc^8*mufac^2*TAGnc*
      Log[mcMSmusmc/mkin])/(10*mkin^9) - 
    (8*apinl1mus^3*mcMSmusmc^6*mufac^2*TAGnc*Log[mcMSmusmc/mkin])/
     (27*mkin^7) + (apinl1mus^3*mcMSmusmc^4*mufac^2*TAGnc*
      Log[mcMSmusmc/mkin])/mkin^5 - (73801799*apinl1mus^3*mcMSmusmc^24*
      NLxMSOS*TAGnc*Log[mcMSmusmc/mkin])/(4231787807520*mkin^23) - 
    (14290513*apinl1mus^3*mcMSmusmc^22*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin])/
     (689665418442*mkin^21) - (1211963*apinl1mus^3*mcMSmusmc^20*NLxMSOS*TAGnc*
      Log[mcMSmusmc/mkin])/(50127997920*mkin^19) - 
    (197062*apinl1mus^3*mcMSmusmc^18*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin])/
     (7381208835*mkin^17) - (355*apinl1mus^3*mcMSmusmc^16*NLxMSOS*TAGnc*
      Log[mcMSmusmc/mkin])/(14455584*mkin^15) - 
    (1153*apinl1mus^3*mcMSmusmc^14*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin])/
     (225450225*mkin^13) + (391*apinl1mus^3*mcMSmusmc^12*NLxMSOS*TAGnc*
      Log[mcMSmusmc/mkin])/(4939704*mkin^11) + 
    (398*apinl1mus^3*mcMSmusmc^10*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin])/
     (893025*mkin^9) + (37*apinl1mus^3*mcMSmusmc^8*NLxMSOS*TAGnc*
      Log[mcMSmusmc/mkin])/(14700*mkin^7) + 
    (2*apinl1mus^3*mcMSmusmc^6*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin])/
     (75*mkin^5) + (apinl1mus^3*mcMSmusmc^4*NLxMSOS*TAGnc*
      Log[mcMSmusmc/mkin])/(12*mkin^3) + 
    (11*apinl1mus^3*mcMSmusmc*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/6 + 
    (533*apinl1mus^3*mcMSmusmc^25*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (277725*mkin^24) + (1620816161*apinl1mus^3*mcMSmusmc^24*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(822167470080*mkin^23) + 
    (445*apinl1mus^3*mcMSmusmc^23*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (192717*mkin^22) + (13926181*apinl1mus^3*mcMSmusmc^22*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(5839257600*mkin^21) + 
    (365*apinl1mus^3*mcMSmusmc^21*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (128877*mkin^20) + (156353*apinl1mus^3*mcMSmusmc^20*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(53084160*mkin^19) + 
    (293*apinl1mus^3*mcMSmusmc^19*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (82365*mkin^18) + (52013*apinl1mus^3*mcMSmusmc^18*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(13934592*mkin^17) + 
    (229*apinl1mus^3*mcMSmusmc^17*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (49725*mkin^16) + (565351*apinl1mus^3*mcMSmusmc^16*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(115605504*mkin^15) + 
    (173*apinl1mus^3*mcMSmusmc^15*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (27885*mkin^14) + (2161*apinl1mus^3*mcMSmusmc^14*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(322560*mkin^13) + 
    (125*apinl1mus^3*mcMSmusmc^13*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (14157*mkin^12) + (40553*apinl1mus^3*mcMSmusmc^12*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(4147200*mkin^11) + 
    (85*apinl1mus^3*mcMSmusmc^11*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (6237*mkin^10) + (17*apinl1mus^3*mcMSmusmc^10*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(1080*mkin^9) + 
    (53*apinl1mus^3*mcMSmusmc^9*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (2205*mkin^8) + (139*apinl1mus^3*mcMSmusmc^8*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(4608*mkin^7) + 
    (29*apinl1mus^3*mcMSmusmc^7*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (525*mkin^6) + (113*apinl1mus^3*mcMSmusmc^6*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(1296*mkin^5) + 
    (13*apinl1mus^3*mcMSmusmc^5*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/(45*mkin^4) + 
    (217*apinl1mus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (432*mkin^3) + (289*apinl1mus^3*mcMSmusmc^3*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(162*mkin^2) + 
    (3*apinl1mus^3*mcMSmusmc^2*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/(2*mkin) - 
    (apinl1mus^3*mcMSmusmc*NLxMSOS*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/9 - 
    (apinl1mus^3*mcMSmusmc^4*NLxMSOS*Pi^2*TAGnc*Log[mcMSmusmc/mkin])/
     (27*mkin^3) - (apinl1mus^3*mcMSmusmc^3*NLxMSOS*Pi^2*TAGnc*
      Log[mcMSmusmc/mkin])/(9*mkin^2) - 
    (apinl1mus^3*mcMSmusmc^2*Pi^4*TAGnc*Log[mcMSmusmc/mkin])/(12*mkin) - 
    (281941*apinl1mus^3*mcMSmusmc^24*TAGnc^2*Log[mcMSmusmc/mkin])/
     (629880300*mkin^23) - (318845*apinl1mus^3*mcMSmusmc^22*TAGnc^2*
      Log[mcMSmusmc/mkin])/(520109667*mkin^21) - 
    (1017529*apinl1mus^3*mcMSmusmc^20*TAGnc^2*Log[mcMSmusmc/mkin])/
     (1168484800*mkin^19) - (334304*apinl1mus^3*mcMSmusmc^18*TAGnc^2*
      Log[mcMSmusmc/mkin])/(258084225*mkin^17) - 
    (7469*apinl1mus^3*mcMSmusmc^16*TAGnc^2*Log[mcMSmusmc/mkin])/
     (3650400*mkin^15) - (263546*apinl1mus^3*mcMSmusmc^14*TAGnc^2*
      Log[mcMSmusmc/mkin])/(75150075*mkin^13) - 
    (9545*apinl1mus^3*mcMSmusmc^12*TAGnc^2*Log[mcMSmusmc/mkin])/
     (1411344*mkin^11) - (14048*apinl1mus^3*mcMSmusmc^10*TAGnc^2*
      Log[mcMSmusmc/mkin])/(893025*mkin^9) - 
    (193*apinl1mus^3*mcMSmusmc^8*TAGnc^2*Log[mcMSmusmc/mkin])/(3675*mkin^7) - 
    (8*apinl1mus^3*mcMSmusmc^6*TAGnc^2*Log[mcMSmusmc/mkin])/(75*mkin^5) + 
    (5*apinl1mus^3*mcMSmusmc^4*TAGnc^2*Log[mcMSmusmc/mkin])/(12*mkin^3) - 
    (apinl1mus^3*mcMSmusmc*Pi^2*TAGnc^2*Log[mcMSmusmc/mkin])/9 - 
    (apinl1mus^3*mcMSmusmc^4*Pi^2*TAGnc^2*Log[mcMSmusmc/mkin])/(27*mkin^3) - 
    (apinl1mus^3*mcMSmusmc^3*Pi^2*TAGnc^2*Log[mcMSmusmc/mkin])/(9*mkin^2) + 
    (apinl1mus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[2]*Log[mcMSmusmc/mkin])/
     (9*mkin^3) + (941278536953*apinl1mus^3*mcMSmusmc^24*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(46646042002560*mkin^23) + 
    (3036958280711*apinl1mus^3*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (124450902576000*mkin^21) + (45436526101*apinl1mus^3*mcMSmusmc^20*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(1508495788800*mkin^19) + 
    (3391310411*apinl1mus^3*mcMSmusmc^18*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (88921857024*mkin^17) + (1141880431*apinl1mus^3*mcMSmusmc^16*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(22884301440*mkin^15) + 
    (4128997*apinl1mus^3*mcMSmusmc^14*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (60540480*mkin^13) + (14841023*apinl1mus^3*mcMSmusmc^12*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(149688000*mkin^11) + 
    (2585963*apinl1mus^3*mcMSmusmc^10*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (16329600*mkin^9) + (9019*apinl1mus^3*mcMSmusmc^8*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(30240*mkin^7) + 
    (5329*apinl1mus^3*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (6480*mkin^5) + (apinl1mus^2*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (3*mkin^3) + (95*apinl1mus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (18*mkin^3) - (16*apinl1mus^3*mcMSmusmc^4*mufac*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(9*mkin^4) - 
    (2*apinl1mus^3*mcMSmusmc^4*mufac^2*TAGnc*Log[mcMSmusmc/mkin]^2)/
     (3*mkin^5) - (13*apinl1mus^3*mcMSmusmc^4*NLxMSOS*TAGnc*
      Log[mcMSmusmc/mkin]^2)/(54*mkin^3) - 
    (apinl1mus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[mcMSmusmc/mkin]^2)/(4*mkin^3) + 
    (11*apinl1mus^3*mcMSmusmc^24*TAGnc^2*Log[mcMSmusmc/mkin]^2)/
     (43470*mkin^23) + (40*apinl1mus^3*mcMSmusmc^22*TAGnc^2*
      Log[mcMSmusmc/mkin]^2)/(118503*mkin^21) + 
    (3*apinl1mus^3*mcMSmusmc^20*TAGnc^2*Log[mcMSmusmc/mkin]^2)/
     (6460*mkin^19) + (32*apinl1mus^3*mcMSmusmc^18*TAGnc^2*
      Log[mcMSmusmc/mkin]^2)/(48195*mkin^17) + 
    (7*apinl1mus^3*mcMSmusmc^16*TAGnc^2*Log[mcMSmusmc/mkin]^2)/
     (7020*mkin^15) + (8*apinl1mus^3*mcMSmusmc^14*TAGnc^2*
      Log[mcMSmusmc/mkin]^2)/(5005*mkin^13) + 
    (5*apinl1mus^3*mcMSmusmc^12*TAGnc^2*Log[mcMSmusmc/mkin]^2)/
     (1782*mkin^11) + (16*apinl1mus^3*mcMSmusmc^10*TAGnc^2*
      Log[mcMSmusmc/mkin]^2)/(2835*mkin^9) + 
    (apinl1mus^3*mcMSmusmc^8*TAGnc^2*Log[mcMSmusmc/mkin]^2)/(70*mkin^7) - 
    (13*apinl1mus^3*mcMSmusmc^4*TAGnc^2*Log[mcMSmusmc/mkin]^2)/(54*mkin^3) - 
    (apinl1mus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]^3)/(3*mkin^3) + 
    (2*apinl1mus^3*mcMSmusmc^4*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]^3)/
     (27*mkin^3) + (2*apinl1mus^3*mcMSmusmc^4*TAGnc^2*Log[mcMSmusmc/mkin]^3)/
     (27*mkin^3) - (88*apinl1mus^2*mufac*Log[(2*mufac)/mus])/9 - 
    (3416*apinl1mus^3*mufac*Log[(2*mufac)/mus])/9 - 
    (11*apinl1mus^2*mufac^2*Log[(2*mufac)/mus])/(3*mkin) - 
    (733*apinl1mus^3*mufac^2*Log[(2*mufac)/mus])/(6*mkin) + 
    (16*apinl1mus^2*mufac*NLxOSKIN*Log[(2*mufac)/mus])/27 + 
    (3388*apinl1mus^3*mufac*NLxOSKIN*Log[(2*mufac)/mus])/81 + 
    (2*apinl1mus^2*mufac^2*NLxOSKIN*Log[(2*mufac)/mus])/(9*mkin) + 
    (715*apinl1mus^3*mufac^2*NLxOSKIN*Log[(2*mufac)/mus])/(54*mkin) - 
    (256*apinl1mus^3*mufac*NLxOSKIN^2*Log[(2*mufac)/mus])/243 - 
    (26*apinl1mus^3*mufac^2*NLxOSKIN^2*Log[(2*mufac)/mus])/(81*mkin) + 
    (88*apinl1mus^3*mufac*Pi^2*Log[(2*mufac)/mus])/9 + 
    (11*apinl1mus^3*mufac^2*Pi^2*Log[(2*mufac)/mus])/(3*mkin) - 
    (16*apinl1mus^3*mufac*NLxOSKIN*Pi^2*Log[(2*mufac)/mus])/27 - 
    (2*apinl1mus^3*mufac^2*NLxOSKIN*Pi^2*Log[(2*mufac)/mus])/(9*mkin) + 
    (484*apinl1mus^3*mufac*Log[(2*mufac)/mus]^2)/9 + 
    (121*apinl1mus^3*mufac^2*Log[(2*mufac)/mus]^2)/(6*mkin) - 
    (176*apinl1mus^3*mufac*NLxOSKIN*Log[(2*mufac)/mus]^2)/27 - 
    (22*apinl1mus^3*mufac^2*NLxOSKIN*Log[(2*mufac)/mus]^2)/(9*mkin) + 
    (16*apinl1mus^3*mufac*NLxOSKIN^2*Log[(2*mufac)/mus]^2)/81 + 
    (2*apinl1mus^3*mufac^2*NLxOSKIN^2*Log[(2*mufac)/mus]^2)/(27*mkin) - 
    (4*apinl1mus^3*mkin*TAGkinmc*Log[mus^2/mcMSmusmc^2])/9 - 
    (8*apinl1mus^2*mufac*TAGkinmc*Log[mus^2/mcMSmusmc^2])/27 - 
    (314*apinl1mus^3*mufac*TAGkinmc*Log[mus^2/mcMSmusmc^2])/27 - 
    (apinl1mus^2*mufac^2*TAGkinmc*Log[mus^2/mcMSmusmc^2])/(9*mkin) - 
    (15*apinl1mus^3*mufac^2*TAGkinmc*Log[mus^2/mcMSmusmc^2])/(4*mkin) + 
    (128*apinl1mus^3*mufac*NLxOSKIN*TAGkinmc*Log[mus^2/mcMSmusmc^2])/243 + 
    (13*apinl1mus^3*mufac^2*NLxOSKIN*TAGkinmc*Log[mus^2/mcMSmusmc^2])/
     (81*mkin) + (8*apinl1mus^3*mufac*Pi^2*TAGkinmc*Log[mus^2/mcMSmusmc^2])/
     27 + (apinl1mus^3*mufac^2*Pi^2*TAGkinmc*Log[mus^2/mcMSmusmc^2])/
     (9*mkin) + (88*apinl1mus^3*mufac*TAGkinmc*Log[(2*mufac)/mus]*
      Log[mus^2/mcMSmusmc^2])/27 + (11*apinl1mus^3*mufac^2*TAGkinmc*
      Log[(2*mufac)/mus]*Log[mus^2/mcMSmusmc^2])/(9*mkin) - 
    (16*apinl1mus^3*mufac*NLxOSKIN*TAGkinmc*Log[(2*mufac)/mus]*
      Log[mus^2/mcMSmusmc^2])/81 - (2*apinl1mus^3*mufac^2*NLxOSKIN*TAGkinmc*
      Log[(2*mufac)/mus]*Log[mus^2/mcMSmusmc^2])/(27*mkin) - 
    (2*apinl1mus^3*mkin*TAGkinmc*Log[mus^2/mcMSmusmc^2]^2)/27 + 
    (4*apinl1mus^3*mufac*TAGkinmc*Log[mus^2/mcMSmusmc^2]^2)/81 + 
    (apinl1mus^3*mufac^2*TAGkinmc*Log[mus^2/mcMSmusmc^2]^2)/(54*mkin) + 
    (2*apinl1mus^3*mkin*TAGkinmc^2*Log[mus^2/mcMSmusmc^2]^2)/27 - 
    apinl1mus*mkin*Log[mus^2/mkin^2] - 
    (461*apinl1mus^2*mkin*Log[mus^2/mkin^2])/72 - 
    (22273*apinl1mus^3*mkin*Log[mus^2/mkin^2])/324 - 
    (16*apinl1mus^2*mufac*Log[mus^2/mkin^2])/9 - 
    (2950*apinl1mus^3*mufac*Log[mus^2/mkin^2])/81 - 
    (2*apinl1mus^2*mufac^2*Log[mus^2/mkin^2])/(3*mkin) - 
    (1277*apinl1mus^3*mufac^2*Log[mus^2/mkin^2])/(108*mkin) + 
    (13*apinl1mus^2*mkin*NLxMSOS*Log[mus^2/mkin^2])/36 + 
    (1283*apinl1mus^3*mkin*NLxMSOS*Log[mus^2/mkin^2])/162 + 
    (4*apinl1mus^3*mufac*NLxMSOS*Log[mus^2/mkin^2])/81 + 
    (apinl1mus^3*mufac^2*NLxMSOS*Log[mus^2/mkin^2])/(54*mkin) - 
    (89*apinl1mus^3*mkin*NLxMSOS^2*Log[mus^2/mkin^2])/648 + 
    (128*apinl1mus^3*mufac*NLxOSKIN*Log[mus^2/mkin^2])/81 + 
    (13*apinl1mus^3*mufac^2*NLxOSKIN*Log[mus^2/mkin^2])/(27*mkin) - 
    (3*apinl1mus^3*mkin*Pi^2*Log[mus^2/mkin^2])/2 + 
    (8*apinl1mus^3*mufac*Pi^2*Log[mus^2/mkin^2])/9 + 
    (apinl1mus^3*mufac^2*Pi^2*Log[mus^2/mkin^2])/(3*mkin) + 
    (13*apinl1mus^3*mkin*NLxMSOS*Pi^2*Log[mus^2/mkin^2])/36 - 
    (apinl1mus^3*mkin*NLxMSOS^2*Pi^2*Log[mus^2/mkin^2])/54 - 
    (9727*apinl1mus^3*mcMSmusmc^24*TAGnc*Log[mus^2/mkin^2])/
     (41473600*mkin^23) - (39833*apinl1mus^3*mcMSmusmc^22*TAGnc*
      Log[mus^2/mkin^2])/(115579926*mkin^21) - 
    (352467*apinl1mus^3*mcMSmusmc^20*TAGnc*Log[mus^2/mkin^2])/
     (667705600*mkin^19) - (16277*apinl1mus^3*mcMSmusmc^18*TAGnc*
      Log[mus^2/mkin^2])/(19117350*mkin^17) - 
    (1587*apinl1mus^3*mcMSmusmc^16*TAGnc*Log[mus^2/mkin^2])/
     (1081600*mkin^15) - (138339*apinl1mus^3*mcMSmusmc^14*TAGnc*
      Log[mus^2/mkin^2])/(50100050*mkin^13) - 
    (1229*apinl1mus^3*mcMSmusmc^12*TAGnc*Log[mus^2/mkin^2])/
     (209088*mkin^11) - (997*apinl1mus^3*mcMSmusmc^10*TAGnc*
      Log[mus^2/mkin^2])/(66150*mkin^9) - 
    (4167*apinl1mus^3*mcMSmusmc^8*TAGnc*Log[mus^2/mkin^2])/(78400*mkin^7) - 
    (19*apinl1mus^3*mcMSmusmc^6*TAGnc*Log[mus^2/mkin^2])/(50*mkin^5) + 
    (151*apinl1mus^3*mcMSmusmc^4*TAGnc*Log[mus^2/mkin^2])/(48*mkin^3) + 
    (9*apinl1mus^3*mcMSmusmc^2*TAGnc*Log[mus^2/mkin^2])/(2*mkin) + 
    (13*apinl1mus^2*mkin*TAGnc*Log[mus^2/mkin^2])/36 + 
    (1283*apinl1mus^3*mkin*TAGnc*Log[mus^2/mkin^2])/162 + 
    (4*apinl1mus^3*mufac*TAGnc*Log[mus^2/mkin^2])/81 + 
    (apinl1mus^3*mufac^2*TAGnc*Log[mus^2/mkin^2])/(54*mkin) + 
    (9727*apinl1mus^3*mcMSmusmc^24*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/
     (559893600*mkin^23) + (39833*apinl1mus^3*mcMSmusmc^22*NLxMSOS*TAGnc*
      Log[mus^2/mkin^2])/(1560329001*mkin^21) + 
    (39163*apinl1mus^3*mcMSmusmc^20*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/
     (1001558400*mkin^19) + (16277*apinl1mus^3*mcMSmusmc^18*NLxMSOS*TAGnc*
      Log[mus^2/mkin^2])/(258084225*mkin^17) + 
    (529*apinl1mus^3*mcMSmusmc^16*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/
     (4867200*mkin^15) + (15371*apinl1mus^3*mcMSmusmc^14*NLxMSOS*TAGnc*
      Log[mus^2/mkin^2])/(75150075*mkin^13) + 
    (1229*apinl1mus^3*mcMSmusmc^12*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/
     (2822688*mkin^11) + (997*apinl1mus^3*mcMSmusmc^10*NLxMSOS*TAGnc*
      Log[mus^2/mkin^2])/(893025*mkin^9) + 
    (463*apinl1mus^3*mcMSmusmc^8*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/
     (117600*mkin^7) + (19*apinl1mus^3*mcMSmusmc^6*NLxMSOS*TAGnc*
      Log[mus^2/mkin^2])/(675*mkin^5) - (151*apinl1mus^3*mcMSmusmc^4*NLxMSOS*
      TAGnc*Log[mus^2/mkin^2])/(648*mkin^3) - 
    (apinl1mus^3*mcMSmusmc^2*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/(3*mkin) - 
    (89*apinl1mus^3*mkin*NLxMSOS*TAGnc*Log[mus^2/mkin^2])/324 - 
    (3*apinl1mus^3*mcMSmusmc*Pi^2*TAGnc*Log[mus^2/mkin^2])/4 + 
    (apinl1mus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[mus^2/mkin^2])/(4*mkin^3) - 
    (3*apinl1mus^3*mcMSmusmc^3*Pi^2*TAGnc*Log[mus^2/mkin^2])/(4*mkin^2) + 
    (13*apinl1mus^3*mkin*Pi^2*TAGnc*Log[mus^2/mkin^2])/36 + 
    (apinl1mus^3*mcMSmusmc*NLxMSOS*Pi^2*TAGnc*Log[mus^2/mkin^2])/18 - 
    (apinl1mus^3*mcMSmusmc^4*NLxMSOS*Pi^2*TAGnc*Log[mus^2/mkin^2])/
     (54*mkin^3) + (apinl1mus^3*mcMSmusmc^3*NLxMSOS*Pi^2*TAGnc*
      Log[mus^2/mkin^2])/(18*mkin^2) - (apinl1mus^3*mkin*NLxMSOS*Pi^2*TAGnc*
      Log[mus^2/mkin^2])/27 + (9727*apinl1mus^3*mcMSmusmc^24*TAGnc^2*
      Log[mus^2/mkin^2])/(559893600*mkin^23) + 
    (39833*apinl1mus^3*mcMSmusmc^22*TAGnc^2*Log[mus^2/mkin^2])/
     (1560329001*mkin^21) + (39163*apinl1mus^3*mcMSmusmc^20*TAGnc^2*
      Log[mus^2/mkin^2])/(1001558400*mkin^19) + 
    (16277*apinl1mus^3*mcMSmusmc^18*TAGnc^2*Log[mus^2/mkin^2])/
     (258084225*mkin^17) + (529*apinl1mus^3*mcMSmusmc^16*TAGnc^2*
      Log[mus^2/mkin^2])/(4867200*mkin^15) + 
    (15371*apinl1mus^3*mcMSmusmc^14*TAGnc^2*Log[mus^2/mkin^2])/
     (75150075*mkin^13) + (1229*apinl1mus^3*mcMSmusmc^12*TAGnc^2*
      Log[mus^2/mkin^2])/(2822688*mkin^11) + 
    (997*apinl1mus^3*mcMSmusmc^10*TAGnc^2*Log[mus^2/mkin^2])/
     (893025*mkin^9) + (463*apinl1mus^3*mcMSmusmc^8*TAGnc^2*
      Log[mus^2/mkin^2])/(117600*mkin^7) + 
    (19*apinl1mus^3*mcMSmusmc^6*TAGnc^2*Log[mus^2/mkin^2])/(675*mkin^5) - 
    (151*apinl1mus^3*mcMSmusmc^4*TAGnc^2*Log[mus^2/mkin^2])/(648*mkin^3) - 
    (apinl1mus^3*mcMSmusmc^2*TAGnc^2*Log[mus^2/mkin^2])/(3*mkin) - 
    (89*apinl1mus^3*mkin*TAGnc^2*Log[mus^2/mkin^2])/648 + 
    (apinl1mus^3*mcMSmusmc*Pi^2*TAGnc^2*Log[mus^2/mkin^2])/18 - 
    (apinl1mus^3*mcMSmusmc^4*Pi^2*TAGnc^2*Log[mus^2/mkin^2])/(54*mkin^3) + 
    (apinl1mus^3*mcMSmusmc^3*Pi^2*TAGnc^2*Log[mus^2/mkin^2])/(18*mkin^2) - 
    (apinl1mus^3*mkin*Pi^2*TAGnc^2*Log[mus^2/mkin^2])/54 - 
    (apinl1mus^3*mkin*Pi^2*Log[2]*Log[mus^2/mkin^2])/2 + 
    (apinl1mus^3*mkin*NLxMSOS*Pi^2*Log[2]*Log[mus^2/mkin^2])/27 + 
    (apinl1mus^3*mkin*Pi^2*TAGnc*Log[2]*Log[mus^2/mkin^2])/27 + 
    (11*apinl1mus^3*mcMSmusmc^24*TAGnc*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (6440*mkin^23) + (10*apinl1mus^3*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(4389*mkin^21) + 
    (81*apinl1mus^3*mcMSmusmc^20*TAGnc*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (25840*mkin^19) + (8*apinl1mus^3*mcMSmusmc^18*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(1785*mkin^17) + 
    (7*apinl1mus^3*mcMSmusmc^16*TAGnc*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (1040*mkin^15) + (54*apinl1mus^3*mcMSmusmc^14*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(5005*mkin^13) + 
    (5*apinl1mus^3*mcMSmusmc^12*TAGnc*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (264*mkin^11) + (4*apinl1mus^3*mcMSmusmc^10*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(105*mkin^9) + 
    (27*apinl1mus^3*mcMSmusmc^8*TAGnc*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (280*mkin^7) + (2*apinl1mus^3*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(5*mkin^5) - (13*apinl1mus^3*mcMSmusmc^4*TAGnc*
      Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/(4*mkin^3) - 
    (11*apinl1mus^3*mcMSmusmc^24*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(86940*mkin^23) - 
    (20*apinl1mus^3*mcMSmusmc^22*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(118503*mkin^21) - 
    (3*apinl1mus^3*mcMSmusmc^20*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(12920*mkin^19) - 
    (16*apinl1mus^3*mcMSmusmc^18*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(48195*mkin^17) - 
    (7*apinl1mus^3*mcMSmusmc^16*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(14040*mkin^15) - 
    (4*apinl1mus^3*mcMSmusmc^14*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(5005*mkin^13) - 
    (5*apinl1mus^3*mcMSmusmc^12*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(3564*mkin^11) - 
    (8*apinl1mus^3*mcMSmusmc^10*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(2835*mkin^9) - 
    (apinl1mus^3*mcMSmusmc^8*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(140*mkin^7) - 
    (4*apinl1mus^3*mcMSmusmc^6*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(135*mkin^5) + 
    (13*apinl1mus^3*mcMSmusmc^4*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(54*mkin^3) - (11*apinl1mus^3*mcMSmusmc^24*TAGnc^2*
      Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/(86940*mkin^23) - 
    (20*apinl1mus^3*mcMSmusmc^22*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(118503*mkin^21) - 
    (3*apinl1mus^3*mcMSmusmc^20*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(12920*mkin^19) - 
    (16*apinl1mus^3*mcMSmusmc^18*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(48195*mkin^17) - 
    (7*apinl1mus^3*mcMSmusmc^16*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(14040*mkin^15) - 
    (4*apinl1mus^3*mcMSmusmc^14*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(5005*mkin^13) - 
    (5*apinl1mus^3*mcMSmusmc^12*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(3564*mkin^11) - 
    (8*apinl1mus^3*mcMSmusmc^10*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(2835*mkin^9) - 
    (apinl1mus^3*mcMSmusmc^8*TAGnc^2*Log[mcMSmusmc/mkin]*Log[mus^2/mkin^2])/
     (140*mkin^7) - (4*apinl1mus^3*mcMSmusmc^6*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(135*mkin^5) + 
    (13*apinl1mus^3*mcMSmusmc^4*TAGnc^2*Log[mcMSmusmc/mkin]*
      Log[mus^2/mkin^2])/(54*mkin^3) + (3*apinl1mus^3*mcMSmusmc^4*TAGnc*
      Log[mcMSmusmc/mkin]^2*Log[mus^2/mkin^2])/(2*mkin^3) - 
    (apinl1mus^3*mcMSmusmc^4*NLxMSOS*TAGnc*Log[mcMSmusmc/mkin]^2*
      Log[mus^2/mkin^2])/(9*mkin^3) - (apinl1mus^3*mcMSmusmc^4*TAGnc^2*
      Log[mcMSmusmc/mkin]^2*Log[mus^2/mkin^2])/(9*mkin^3) + 
    (88*apinl1mus^3*mufac*Log[(2*mufac)/mus]*Log[mus^2/mkin^2])/9 + 
    (11*apinl1mus^3*mufac^2*Log[(2*mufac)/mus]*Log[mus^2/mkin^2])/(3*mkin) - 
    (16*apinl1mus^3*mufac*NLxOSKIN*Log[(2*mufac)/mus]*Log[mus^2/mkin^2])/27 - 
    (2*apinl1mus^3*mufac^2*NLxOSKIN*Log[(2*mufac)/mus]*Log[mus^2/mkin^2])/
     (9*mkin) - (apinl1mus^3*mkin*TAGkinmc*Log[mus^2/mcMSmusmc^2]*
      Log[mus^2/mkin^2])/3 + (8*apinl1mus^3*mufac*TAGkinmc*
      Log[mus^2/mcMSmusmc^2]*Log[mus^2/mkin^2])/27 + 
    (apinl1mus^3*mufac^2*TAGkinmc*Log[mus^2/mcMSmusmc^2]*Log[mus^2/mkin^2])/
     (9*mkin) - (apinl1mus^3*mkin*TAGkinmc*Log[mus^2/mcMSmusmc^2]^2*
      Log[mus^2/mkin^2])/18 + (apinl1mus^3*mkin*TAGkinmc^2*
      Log[mus^2/mcMSmusmc^2]^2*Log[mus^2/mkin^2])/18 - 
    (23*apinl1mus^2*mkin*Log[mus^2/mkin^2]^2)/24 - 
    (14275*apinl1mus^3*mkin*Log[mus^2/mkin^2]^2)/864 - 
    (46*apinl1mus^3*mufac*Log[mus^2/mkin^2]^2)/27 - 
    (23*apinl1mus^3*mufac^2*Log[mus^2/mkin^2]^2)/(36*mkin) + 
    (apinl1mus^2*mkin*NLxMSOS*Log[mus^2/mkin^2]^2)/12 + 
    (107*apinl1mus^3*mkin*NLxMSOS*Log[mus^2/mkin^2]^2)/48 + 
    (4*apinl1mus^3*mufac*NLxMSOS*Log[mus^2/mkin^2]^2)/27 + 
    (apinl1mus^3*mufac^2*NLxMSOS*Log[mus^2/mkin^2]^2)/(18*mkin) - 
    (13*apinl1mus^3*mkin*NLxMSOS^2*Log[mus^2/mkin^2]^2)/216 + 
    (apinl1mus^2*mkin*TAGnc*Log[mus^2/mkin^2]^2)/12 + 
    (107*apinl1mus^3*mkin*TAGnc*Log[mus^2/mkin^2]^2)/48 + 
    (4*apinl1mus^3*mufac*TAGnc*Log[mus^2/mkin^2]^2)/27 + 
    (apinl1mus^3*mufac^2*TAGnc*Log[mus^2/mkin^2]^2)/(18*mkin) - 
    (13*apinl1mus^3*mkin*NLxMSOS*TAGnc*Log[mus^2/mkin^2]^2)/108 - 
    (13*apinl1mus^3*mkin*TAGnc^2*Log[mus^2/mkin^2]^2)/216 - 
    (601*apinl1mus^3*mkin*Log[mus^2/mkin^2]^3)/432 + 
    (25*apinl1mus^3*mkin*NLxMSOS*Log[mus^2/mkin^2]^3)/108 - 
    (apinl1mus^3*mkin*NLxMSOS^2*Log[mus^2/mkin^2]^3)/108 + 
    (25*apinl1mus^3*mkin*TAGnc*Log[mus^2/mkin^2]^3)/108 - 
    (apinl1mus^3*mkin*NLxMSOS*TAGnc*Log[mus^2/mkin^2]^3)/54 - 
    (apinl1mus^3*mkin*TAGnc^2*Log[mus^2/mkin^2]^3)/108 + 
    (4*apinl1mus^3*mkin*TAGkinmc*Log[mus^2/musmc^2])/9 - 
    (16*apinl1mus^3*mufac*TAGkinmc*Log[mus^2/musmc^2])/27 - 
    (2*apinl1mus^3*mufac^2*TAGkinmc*Log[mus^2/musmc^2])/(9*mkin) + 
    (apinl1mus^3*mkin*TAGkinmc*Log[mus^2/mkin^2]*Log[mus^2/musmc^2])/3 + 
    (4*apinl1mus^3*mkin*TAGkinmc*Log[musmc^2/mcMSmusmc^2])/9 - 
    (10163*apinl1mus^3*mcMSmusmc^24*TAGnc*Log[musmc^2/mcMSmusmc^2])/
     (11664450*mkin^23) - (5066*apinl1mus^3*mcMSmusmc^22*TAGnc*
      Log[musmc^2/mcMSmusmc^2])/(4298427*mkin^21) - 
    (5507*apinl1mus^3*mcMSmusmc^20*TAGnc*Log[musmc^2/mcMSmusmc^2])/
     (3338528*mkin^19) - (7678*apinl1mus^3*mcMSmusmc^18*TAGnc*
      Log[musmc^2/mcMSmusmc^2])/(3186225*mkin^17) - 
    (283*apinl1mus^3*mcMSmusmc^16*TAGnc*Log[musmc^2/mcMSmusmc^2])/
     (76050*mkin^15) - (3166*apinl1mus^3*mcMSmusmc^14*TAGnc*
      Log[musmc^2/mcMSmusmc^2])/(511225*mkin^13) - 
    (899*apinl1mus^3*mcMSmusmc^12*TAGnc*Log[musmc^2/mcMSmusmc^2])/
     (78408*mkin^11) - (298*apinl1mus^3*mcMSmusmc^10*TAGnc*
      Log[musmc^2/mcMSmusmc^2])/(11907*mkin^9) - 
    (179*apinl1mus^3*mcMSmusmc^8*TAGnc*Log[musmc^2/mcMSmusmc^2])/
     (2450*mkin^7) - (94*apinl1mus^3*mcMSmusmc^6*TAGnc*
      Log[musmc^2/mcMSmusmc^2])/(225*mkin^5) + 
    (56*apinl1mus^3*mcMSmusmc^4*TAGnc*Log[musmc^2/mcMSmusmc^2])/(27*mkin^3) + 
    (2*apinl1mus^3*mcMSmusmc^2*TAGnc*Log[musmc^2/mcMSmusmc^2])/mkin - 
    (apinl1mus^3*mcMSmusmc*Pi^2*TAGnc*Log[musmc^2/mcMSmusmc^2])/6 + 
    (2*apinl1mus^3*mcMSmusmc^4*Pi^2*TAGnc*Log[musmc^2/mcMSmusmc^2])/
     (9*mkin^3) - (apinl1mus^3*mcMSmusmc^3*Pi^2*TAGnc*
      Log[musmc^2/mcMSmusmc^2])/(2*mkin^2) + 
    (22*apinl1mus^3*mcMSmusmc^24*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(2415*mkin^23) + 
    (40*apinl1mus^3*mcMSmusmc^22*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(3591*mkin^21) + 
    (9*apinl1mus^3*mcMSmusmc^20*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(646*mkin^19) + 
    (32*apinl1mus^3*mcMSmusmc^18*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(1785*mkin^17) + 
    (14*apinl1mus^3*mcMSmusmc^16*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(585*mkin^15) + 
    (24*apinl1mus^3*mcMSmusmc^14*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(715*mkin^13) + 
    (5*apinl1mus^3*mcMSmusmc^12*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(99*mkin^11) + 
    (16*apinl1mus^3*mcMSmusmc^10*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(189*mkin^9) + 
    (6*apinl1mus^3*mcMSmusmc^8*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(35*mkin^7) + 
    (8*apinl1mus^3*mcMSmusmc^6*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(15*mkin^5) - 
    (20*apinl1mus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]*
      Log[musmc^2/mcMSmusmc^2])/(9*mkin^3) + 
    (4*apinl1mus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]^2*
      Log[musmc^2/mcMSmusmc^2])/(3*mkin^3) + 
    (apinl1mus^3*mkin*TAGkinmc*Log[mus^2/mkin^2]*Log[musmc^2/mcMSmusmc^2])/
     3 + (220*apinl1mus^3*mkin*PolyLog[4, 1/2])/27 - 
    (8*apinl1mus^3*mkin*NLxMSOS*PolyLog[4, 1/2])/27 - 
    (20*apinl1mus^3*mcMSmusmc^4*TAGnc*PolyLog[4, 1/2])/(27*mkin^3) - 
    (8*apinl1mus^3*mkin*TAGnc*PolyLog[4, 1/2])/27 + 
    (apinl1mus^2*mkin*Zeta[3])/6 - (61*apinl1mus^3*mkin*Zeta[3])/27 - 
    (3070*apinl1mus^3*mufac*Zeta[3])/27 - (1535*apinl1mus^3*mufac^2*Zeta[3])/
     (36*mkin) + (241*apinl1mus^3*mkin*NLxMSOS*Zeta[3])/72 - 
    (7*apinl1mus^3*mkin*NLxMSOS^2*Zeta[3])/54 + 
    (140*apinl1mus^3*mufac*NLxOSKIN*Zeta[3])/27 + 
    (35*apinl1mus^3*mufac^2*NLxOSKIN*Zeta[3])/(18*mkin) + 
    (1439*apinl1mus^3*mkin*Pi^2*Zeta[3])/432 + 
    (2710689767*apinl1mus^3*mcMSmusmc^24*TAGnc*Zeta[3])/
     (1644334940160*mkin^23) + (23017987*apinl1mus^3*mcMSmusmc^22*TAGnc*
      Zeta[3])/(11678515200*mkin^21) + (254791*apinl1mus^3*mcMSmusmc^20*TAGnc*
      Zeta[3])/(106168320*mkin^19) + (83291*apinl1mus^3*mcMSmusmc^18*TAGnc*
      Zeta[3])/(27869184*mkin^17) + (885457*apinl1mus^3*mcMSmusmc^16*TAGnc*
      Zeta[3])/(231211008*mkin^15) + (3287*apinl1mus^3*mcMSmusmc^14*TAGnc*
      Zeta[3])/(645120*mkin^13) + (59231*apinl1mus^3*mcMSmusmc^12*TAGnc*
      Zeta[3])/(8294400*mkin^11) + (187*apinl1mus^3*mcMSmusmc^10*TAGnc*
      Zeta[3])/(17280*mkin^9) + (173*apinl1mus^3*mcMSmusmc^8*TAGnc*Zeta[3])/
     (9216*mkin^7) + (29*apinl1mus^3*mcMSmusmc^6*TAGnc*Zeta[3])/
     (648*mkin^5) + (2309*apinl1mus^3*mcMSmusmc^4*TAGnc*Zeta[3])/
     (864*mkin^3) + (11*apinl1mus^3*mcMSmusmc^2*TAGnc*Zeta[3])/(2*mkin) + 
    (241*apinl1mus^3*mkin*TAGnc*Zeta[3])/72 - 
    (apinl1mus^3*mcMSmusmc^4*NLxMSOS*TAGnc*Zeta[3])/(3*mkin^3) - 
    (7*apinl1mus^3*mkin*NLxMSOS*TAGnc*Zeta[3])/27 - 
    (3*apinl1mus^3*mcMSmusmc^2*Pi^2*TAGnc*Zeta[3])/(4*mkin) - 
    (7*apinl1mus^3*mkin*TAGnc^2*Zeta[3])/54 - 
    (14*apinl1mus^3*mcMSmusmc^4*TAGnc*Log[mcMSmusmc/mkin]*Zeta[3])/
     (9*mkin^3) + (19*apinl1mus^3*mkin*Log[mus^2/mkin^2]*Zeta[3])/12 + 
    (7*apinl1mus^3*mkin*NLxMSOS*Log[mus^2/mkin^2]*Zeta[3])/9 + 
    (7*apinl1mus^3*mkin*TAGnc*Log[mus^2/mkin^2]*Zeta[3])/9 - 
    (1975*apinl1mus^3*mkin*Zeta[5])/216 - 
    (5*apinl1mus^3*mcMSmusmc^2*TAGnc*Zeta[5])/(2*mkin)


(* ************************************************************ *)

(* Compute \\Lambda^(nf) using 3 methods:
   1. Explicit formulae for \\Lambda: lamexpl[alphas_,mu_,nn_,loops_]
   2. Implicit formulae for \\Lambda: lamimpl[alphas_,mu_,nn_,loops_]
   3. Solving the integral exactly:   lamexact[alphas_,mu_,nn_]
*)

(* Determine \Lamba for given \alpha_s				*)

LamExpl[alphas_,mu_,nn_,loops_] := Module[
    {logmu2lam2,nnn},
    If[(loops<1)||(loops>5), Print["LamExpl: Invalid # loops."];
       Return[];];
    (* \log (\mu^2/\Lambda^2) =	
     *)
    logmu2lam2 = 1/(as*b0) + (as*(-b1^2 + b0*b2))/b0^3 +
	(as^2*(b1^3 - 2*b0*b1*b2 + b0^2*b3))/(2*b0^4) + (b1*Log[as])/b0^2 +
	    (b1*Log[b0])/b0^2
       + (as^3*(-b1^4 + 3*b0*b1^2*b2 
                - b0^2*b2^2 - 2*b0^2*b1*b3 +  b0^3*b4))/  (3*b0^5)
       ;
    (* \Lambda/\mu = Exp[ lamomu[nn_] ]; nn: # active flavours
     *)
    lamomu[nnn_]:=Collect[Expand[ 
	-1/2*logmu2lam2/.setbeta/.nf->nnn
	],as]/.as:>ToExpression["api"<>ToString[nnn]];
    Return[N[mu*Exp[
	Expand[
	    (Expand[ lamomu[nn] * ToExpression["api"<>ToString[nn]]^2 ]
	     /.cut[ToExpression["api"<>ToString[nn]],loops]
	     )
		*1/ToExpression["api"<>ToString[nn]]^2
		]
	/.ToExpression["api"<>ToString[nn]]->alphas/Pi] 
	     ,numprec]];
];

(* ************************************************************ *)

LamImpl[alphas_,mu_,nn_,loops_] := Module[
    {lam,astmp},
    If[(loops<1)||(loops>5), Print["LamImpl: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==1, astmp=as0];
    If[loops==2, astmp=as1];
    If[loops==3, astmp=as2];
    If[loops==4, astmp=as3];
    If[loops==5, astmp=as4];
    Return[lam/.FindRoot[
	(Expand[
	    astmp/.setasL
		/.setbeta/.nf->nn]/.L->Log[mu^2/lam^2])==alphas/Pi,
	    {lam,LamExpl[alphas,mu,nn,loops]} 
    ]
       ]
];

(* ************************************************************ *)

(* Compute \alpha_s(\mu);
   Input: \Lambda, \mu, # active flavours and # loops
*)

AlphasLam[lambda_,mu_,nn_,loops_] := Module[
    {},
    If[(loops<1)||(loops>5), Print["AlphasLam: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[ N[mu/lambda] < rmulam ,
	Print["WARNING: the ratio \\mu/\\Lambda = ",
	      N[mu/lambda]," is very small!"];
    ];
    If[loops==1, Return[ N[Pi*as0/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
    If[loops==2, Return[ N[Pi*as1/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
    If[loops==3, Return[ N[Pi*as2/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
    If[loops==4, Return[ N[Pi*as3/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
    If[loops==5, Return[ N[Pi*as4/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
];

(* ************************************************************ *)

(* Compute \alpha_s(\mu);
   Solve differential equation numerically
   Input: \alpha_s(\mu0), \mu0, \mu, # active flavours and # loops
*)

AlphasExact[alphasmu0_,mu0_,muend_,nnf_,loops_] := Module[
    {rhs,solx,a,x,lambda,alphasmu0tmp},
    alphasmu0tmp = Rationalize[alphasmu0,10^-numprec];
    If[(loops<1)||(loops>5), Print["AlphasExact: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    lambda = LamExpl[alphasmu0,mu0,nnf,4];
    If[ N[muend/lambda] < rmulam ,
	Print["WARNING: the ratio \\mu/\\Lambda = ",
	      N[muend/lambda]," is very small!"];
    ];
    If[ mu0 == muend , Return[alphasmu0] ];
    If[ loops==1, rhs = -a[x]^2*(b0) ];
    If[ loops==2, rhs = -a[x]^2*(b0 + b1*a[x]) ];
    If[ loops==3, rhs = -a[x]^2*(b0 + b1*a[x] + b2*a[x]^2) ];
    If[ loops==4, rhs = -a[x]^2*(b0 + b1*a[x] + b2*a[x]^2 + b3*a[x]^3) ];
    If[ loops==5, rhs = -a[x]^2*(b0 + b1*a[x] + b2*a[x]^2 + b3*a[x]^3 
				 + b4*a[x]^4);
	(* Print["b4 is set to ",Global`b4num /. nf->nnf //N]; *) ];
    sol = N[(Pi*a[x]/.
	     NDSolve[{a[mu0^2]==alphasmu0tmp/Pi,
		      x*a'[x] == rhs
         /.setbeta(*/.{b4->Global`b4num}*)/.nf->nnf}, 
         a[x], 
         {x,muend^2,mu0^2},        WorkingPrecision->26,
                                   AccuracyGoal->16,
                                   PrecisionGoal->16,
                                   MaxSteps->1500]
/.x->muend^2)[[1]],numprec];
Return[sol];
];

(* ************************************************************ *)

(*
 Light-mass corrections at order \alpha_s^2 and \alpha_s^3

 deltamOS2mMS[] is the difference to the massless case, i.e. it vanishes for mq=0.
 Examples for format of argument mq:
 - {mb}                 -> light mass: mb(mb)
 - {mb,mc}                             mb(mb) and mc(mc) order is important
 - {{mb,mub},{mc,muc}}                 mb(mub) and mc(muc) order is important
 - {{mb,mub},mc}                       mb(mub) and mc(mc) order is important
 - {mb,{mc,muc}}                       mb(mb) and mc(muc) order is important
*)

deltamOS2mMS[mOS_,mq_,asmu_,mu_,nnf_,loops_] := Module[
    {len,ii,delout,deltazmxfSUB,api,xf,nl,mfMS,muf,Mq,nlini},
    delout = 0;
    len = Length[mq];
    Mq = mOS;
    nlini = nnf;
    For[ii=1, ii<=len, ii++,
      If[Length[ mq[[ii]] ]==0,
        mfMS = mq[[ii]]; muf  = mq[[ii]];
        ,
        mfMS = mq[[ii,1]]; muf  = mq[[ii,2]];
      ];
      xf = mfMS/Mq;
      nl = nlini-1; nlini--;
      deltazmxfSUB = (
        + api^2*(-(Pi^2*xf)/6 + xf^2/2 - (Pi^2*xf^3)/6 + 
     (Pi^2*xf^4)/18 + (xf^4*Log[xf]^2)/3 + 
     Log[xf]*(xf^2/3 + (-1/3 + xf/3 + xf^3/3 - xf^4/3)*Log[1 - xf] + 
       (-1/3 - xf/3 - xf^3/3 - xf^4/3)*Log[1 + xf]) + 
     (-1/3 - xf/3 - xf^3/3 - xf^4/3)*PolyLog[2, -xf] + 
     (-1/3 + xf/3 + xf^3/3 - xf^4/3)*PolyLog[2, xf])
        + api^3*(-24.05855452178852*xf + 
        0.9826666666666666*nl*xf - 1.9355435378308894*xf^2 + 
        0.30033333333333334*nl*xf^2 - 3.918074523431132*xf^3 + 
        2.2579995190968614*xf^4 - 0.7860673441568309*xf^5 + 
        (-6.610555555555555*xf + 0.5346666666666666*nl*xf + 
          2.4651111111111113*xf^2 - 0.22*nl*xf^2 - 0.7243333333333333*xf^3 + 
          0.067*nl*xf^3)*Log[mu^2/Mq^2] + (-1.6403325580080583*xf + 
          1.8093423466268328*xf^2 - 2.1719725592400154*xf^3 + 
          1.730499639322646*xf^4 - 0.5895505081176231*xf^5)*Log[muf^2/mfMS^2] + 
        (16.947666666666667*xf - 1.1013333333333333*nl*xf + 
          2.7875555555555556*xf^2 - 0.03433333333333333*xf^3)*Log[xf])
        /. cut[api,loops]);
      delout = delout + deltazmxfSUB;
(*
      Print["mfMS = ",mfMS," muf = ",muf," nl = ",nl];
      Print[deltazmxfSUB];
*)
    ];
    Return[delout /. api->asmu/Pi];
];

deltamMS2mOS[mMS_,mqMS_,asmu_,mu_,nnf_,loops_] := Module[
    {len,ii,delout,ideltazmxfSUB,api,xfq,nl,mfMS,muf,mq,nlini},
    delout = 0;
    len = Length[mqMS];
    mq = mMS;
    nlini = nnf;
    For[ii=1, ii<=len, ii++,
      If[Length[ mqMS[[ii]] ]==0,
        mfMS = mqMS[[ii]]; muf  = mqMS[[ii]];
        ,
        mfMS = mqMS[[ii,1]]; muf  = mqMS[[ii,2]];
      ];
      xfq = mfMS/mq;
      nl = nlini-1; nlini--;
      ideltazmxfSUB = (
      + (-1)*(api^2*(-(Pi^2*xfq)/6 + xfq^2/2 - (Pi^2*xfq^3)/6 + 
      (Pi^2*xfq^4)/18 + (xfq^4*Log[xfq]^2)/3 + 
      Log[xfq]*(xfq^2/3 + (-1/3 + xfq/3 + xfq^3/3 - xfq^4/3)*Log[1 - xfq] + 
        (-1/3 - xfq/3 - xfq^3/3 - xfq^4/3)*Log[1 + xfq]) + 
      (-1/3 - xfq/3 - xfq^3/3 - xfq^4/3)*PolyLog[2, -xfq] + 
      (-1/3 + xfq/3 + xfq^3/3 - xfq^4/3)*PolyLog[2, xfq]))
      + api^3*(22.968067155676593*xfq - 
      0.9826666666666666*nl*xfq + 4.014666666666667*xfq^2 - 
      0.30033333333333334*nl*xfq^2 + 1.0221111111111112*xfq^3 + 
      (Pi^2*xfq^3)/9 + 0.049333333333333326*xfq^4 - (Pi^2*xfq^4)/27 - 
      (2*xfq^4*Log[xfq]^2)/9 + Log[xfq]*(-16.947666666666667*xfq + 
        1.1013333333333333*nl*xfq - 3.009777777777778*xfq^2 + 
        0.03433333333333333*xfq^3 + (2/9 - (2*xfq)/9 - (2*xfq^3)/9 + 
          (2*xfq^4)/9)*Log[1 - xfq] + (2/9 + (2*xfq)/9 + (2*xfq^3)/9 + 
          (2*xfq^4)/9)*Log[1 + xfq]) + (2*PolyLog[2, -xfq])/9 + 
      (2*xfq*PolyLog[2, -xfq])/9 + (2*xfq^3*PolyLog[2, -xfq])/9 + 
      (2*xfq^4*PolyLog[2, -xfq])/9 + (2*PolyLog[2, xfq])/9 - 
      (2*xfq*PolyLog[2, xfq])/9 - (2*xfq^3*PolyLog[2, xfq])/9 + 
      (2*xfq^4*PolyLog[2, xfq])/9 + Log[mu^2/mq^2]*(6.610555555555555*xfq - 
        0.5346666666666666*nl*xfq + (Pi^2*xfq)/6 - 2.1317777777777778*xfq^2 + 
        0.22*nl*xfq^2 + 0.7243333333333333*xfq^3 - 0.067*nl*xfq^3 - 
        (Pi^2*xfq^3)/6 + (Pi^2*xfq^4)/9 + (2*xfq^4*Log[xfq]^2)/3 + 
        Log[xfq]*((2/3 - xfq/3 + xfq^3/3 - (2*xfq^4)/3)*Log[1 - xfq] + 
          (2/3 + xfq/3 - xfq^3/3 - (2*xfq^4)/3)*Log[1 + xfq]) + 
        (2*PolyLog[2, -xfq])/3 + (xfq*PolyLog[2, -xfq])/3 - 
        (xfq^3*PolyLog[2, -xfq])/3 - (2*xfq^4*PolyLog[2, -xfq])/3 + 
        (2*PolyLog[2, xfq])/3 - (xfq*PolyLog[2, xfq])/3 + 
        (xfq^3*PolyLog[2, xfq])/3 - (2*xfq^4*PolyLog[2, xfq])/3) + 
      Log[muf^2/mfMS^2]*((Pi^2*xfq)/6 - (4*xfq^2)/3 + (Pi^2*xfq^3)/2 - 
        (2*Pi^2*xfq^4)/9 - (4*xfq^4*Log[xfq]^2)/3 + 
        Log[xfq]*((-2*xfq^2)/3 + (-xfq/3 - xfq^3 + (4*xfq^4)/3)*
           Log[1 - xfq] + (xfq/3 + xfq^3 + (4*xfq^4)/3)*Log[1 + xfq]) + 
        (xfq*PolyLog[2, -xfq])/3 + xfq^3*PolyLog[2, -xfq] + 
        (4*xfq^4*PolyLog[2, -xfq])/3 - (xfq*PolyLog[2, xfq])/3 - 
        xfq^3*PolyLog[2, xfq] + (4*xfq^4*PolyLog[2, xfq])/3))
        /. cut[api,loops]);
      delout = delout + ideltazmxfSUB;
(*
      Print["mfMS = ",mfMS," muf = ",muf," nl = ",nl];
      Print[ideltazmxfSUB];
*)
    ];
    Return[delout /. api->asmu/Pi];
];


(*
 OLD: Light-mass corrections at order \alpha_s^2
*)

delta[mOS_,mq_] := Module[
    {i,del,delta,r},
    del=0;
    delta[r_] := If[(r<=1) && (r>=0),
		    Pi^2/8*r-0.597*r^2+0.230*r^3,
		    Print["\\Delta(",N[r],
                          ") IS CALLED; THE FUNCTION IS NOT IMPLEMENTED 
FOR ARGUMENTS OUTSIDE THE INTERVAL [0,1]."];
		    Abort[]
		];
    For[i=1,i<=Length[mq],i++,
	del=del+delta[mq[[i]]/mOS];
    ];
    Return[del];
];

(* ************************************************************ *)

(*
 m_MS = M_OS * ( 1 + ... )
*)

mOS2mMS[mOS_,mq_,asmu_,mu_,nnf_,loops_] := 
    mOS2mMS[mOS,mq,asmu,mu,nnf,loops,1];

mOS2mMS[mOS_,mq_,asmu_,mu_,nnf_,loops_,fdelm_] := Module[
    {extmp,ex4ltmp,Mu,deltalight=0},
    If[(loops<0)||(loops>4), Print["mOS2mMS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[mq=!={}, deltalight = deltamOS2mMS[mOS,mq,asmu,mu,nnf,loops]];
(*
    Print[deltalight];
    Print[api^2*(-4/3 * delta[mOS,mq])/.api->asmu/Pi];
    Abort[];
*)
    If[loops==0,extmp=1,
       extmp = 
	   ( msfromos 
	    + deltalight
	    + api^3*msfromos3 /.setMS2OS
	    )/.cut[api,loops] /.{lmM->Log[mu^2/mOS^2],nl->nnf-1
				 } /.num1 /.cf2num;
   ];

    If[ loops == 4, 
	(ex4ltmp = 
 - 3654.15040757339*fdelm - 1524.2292266911543*Log[Mu^2/mOS^2] - 
 288.778291935394*Log[Mu^2/mOS^2]^2 - 32.54735725308642*Log[Mu^2/mOS^2]^3 - 
 1.85546875*Log[Mu^2/mOS^2]^4 + (-1 + nnf)^3*(0. + 0.678141025604516*fdelm + 
   0.3205521521864135*Log[Mu^2/mOS^2] + 0.0800290327210927*
    Log[Mu^2/mOS^2]^2 + 0.010030864197530864*Log[Mu^2/mOS^2]^3) + 
 (-1 + nnf)^2*(0. - 43.48241924867489*fdelm - 
   19.82672048099557*Log[Mu^2/mOS^2] - 4.482957520194182*Log[Mu^2/mOS^2]^2 - 
   0.5270061728395061*Log[Mu^2/mOS^2]^3 - 0.04108796296296297*
    Log[Mu^2/mOS^2]^4) + (-1 + nnf)*(0. + 756.9421565599532*fdelm + 
   330.1770776731065*Log[Mu^2/mOS^2] + 67.99849534415492*Log[Mu^2/mOS^2]^2 + 
   7.595293209876542*Log[Mu^2/mOS^2]^3 + 0.48119212962962954*
    Log[Mu^2/mOS^2]^4)
	 );
	extmp = extmp + api^4 * ex4ltmp /. {nl->nnf-1, Mu->mu} // Expand;
      ];

    Return[ N[ mOS * (extmp/.api->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 m_OS = M_MS * ( 1 + ... )
*)

mMS2mOS[mMS_,mq_,asmu_,mu_,nnf_,loops_] := mMS2mOS[mMS,mq,asmu,mu,nnf,loops,1];
mMS2mOS[mMS_,mq_,asmu_,mu_,nnf_,loops_,fdelm_] := Module[
    {extmp,ex4ltmp,Mu,deltalight=0},
    If[(loops<0)||(loops>4), Print["mMS2mOS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[mq=!={}, deltalight = deltamMS2mOS[mMS,mq,asmu,mu,nnf,loops]];
(*
    Print[deltalight];
    Print[api^2*(+4/3 * delta[mMS,mq])/.api->asmu/Pi];
    Abort[];
*)
    If[loops==0,extmp=1,
       extmp = 
	   ( osfromms 
	    + deltalight
	    + api^3*osfromms3 /.setMS2OS
          )/.cut[api,loops] /.{lmm->Log[mu^2/mMS^2],nl->nnf-1} /.num1 /.cf2num;
   ];
    If[ loops == 4, 
	(ex4ltmp = 
3567.6071717381847*fdelm + 1727.2260148986106*Log[(1.*mu^2)/mMS^2] + 
 409.2429990574718*Log[(1.*mu^2)/mMS^2]^2 + 
 66.93663194444443*Log[(1.*mu^2)/mMS^2]^3 + 
 8.056278935185185*Log[(1.*mu^2)/mMS^2]^4 + 
 (-1 + nnf)^3*(-0.678141025604516*fdelm - 0.3205521521864134*
    Log[(1.*mu^2)/mMS^2] - 0.0800290327210927*Log[(1.*mu^2)/mMS^2]^2 - 
   0.010030864197530864*Log[(1.*mu^2)/mMS^2]^3) + 
 (-1 + nnf)*(-745.7206623740028*fdelm - 358.29765085086774*
    Log[(1.*mu^2)/mMS^2] - 87.39262571554698*Log[(1.*mu^2)/mMS^2]^2 - 
   11.883873456790122*Log[(1.*mu^2)/mMS^2]^3 - 
   1.2705439814814814*Log[(1.*mu^2)/mMS^2]^4) + 
 (-1 + nnf)^2*(43.396250117985666*fdelm + 20.528466368867228*
    Log[(1.*mu^2)/mMS^2] + 4.971905254812516*Log[(1.*mu^2)/mMS^2]^2 + 
   0.6304012345679011*Log[(1.*mu^2)/mMS^2]^3 + 0.06655092592592593*
    Log[(1.*mu^2)/mMS^2]^4)
	 );
	extmp = extmp + api^4 * ex4ltmp /. {nl->nnf-1, Mu->mu} // Expand;
      ];
    Return[ N[ mMS * (extmp/.api->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 Compute in a first step m_MS(m_MS) and afterwards m_MS(mu).
*)

mOS2mMSrun[mOS_,mq_,asmu_,mu_,nnf_,loops_] := Module[
    {mum,asM,asmum},
    asM = AlphasExact[asmu,mu,mOS,nnf,loops];
    mum = mOS2mSI[mOS,mq,asM,nnf,loops];
    asmum = AlphasExact[asmu,mu,mum,nnf,loops];
    Return[ N[ mMS2mMS[mum,asmum,asmu,nnf,loops] , numprec ] ];
];

(* ************************************************************ *)

(*
 Compute in a first step m_MS(m_MS). Then use m_OS = m_MS * ( 1 + ... )
 for mu=m_MS.
*)

mMS2mOSrun[mMS_,mq_,asmu_,mu_,nnf_,loops_] := Module[
    {asmMS,mMS2mSI2,xx,res1},
    asmMS = AlphasLam[LamImpl[asmu,mu,nnf,loops],xx,nnf,loops];
    mMS2mSI2 = mMS2mMS[mMS,asmu,asmMS,nnf,loops];
    res1 = FindRoot[ xx == mMS2mSI2 , {xx,mMS} ];
    (* print mMS2mSI2 and \alpha_s(mMS2mSI2)
     Print[ mMS2mSI2/.res1 ];
     Print[ asmMS/.res1 ];
     *)
    Return[ N[ mMS2mOS[mMS2mSI2/.res1,mq,asmMS/.res1,mMS2mSI2/.res1,nnf,loops] 
	       , numprec ] ];
];

(* ************************************************************ *)

mOS2mMSit[mOS_,mq_,asmu_,mu_,nnf_,loops_] := Module[
    {extmp},
    If[(loops<0)||(loops>5), Print["mOS2mMSit: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    extmp = mMS2mOS[mMS,mq,asmu,mu,nnf,loops];
    Return[ N[ mMS/.FindRoot[ mOS == extmp , {mMS,mOS}]
	       , numprec ] ];
];

(* ************************************************************ *)

(*
 \mu_m = M_OS * ( 1 + ... )
*)

mOS2mSI[mOS_,mq_,asM_,nnf_,loops_] := mOS2mSI[mOS,mq,asM,nnf,loops,1];
mOS2mSI[mOS_,mq_,asM_,nnf_,loops_,fdelm_] := Module[
    {extmp,deltalight=0},
    If[(loops<0)||(loops>4), Print["mOS2mSI: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[mq=!={}, deltalight = deltamOS2mMS[mOS,mq,asM,mOS,nnf,loops]];
    If[loops==0,extmp=1,
       extmp = 
	   ( mumfromos 
	    + deltalight
	    + apiM^3*mumfromos3 /.setMS2OS
          )/.cut[apiM,loops] /.{nl->nnf-1} /.num1 /.cf2num;
   ];
    If[ loops == 4, 
        (ex4ltmp =
fdelm*(-3214.227044839041 + 692.4809215366435*(-1 + nnf) - 
  41.95978562498058*(-1 + nnf)^2 + 0.678141025604516*(-1 + nnf)^3)
	 );
	extmp = extmp + apiM^4 * ex4ltmp /. {nl->nnf-1} // Expand;
      ];
    Return[ N[ mOS * (extmp/.apiM->asM/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 Running MS-bar mass
*)

mMS2mMS[mmu0_,asmu0_,asmu_,nnf_,loops_] := Module[
    {ccc},
    If[(loops<0)||(loops>5), Print["mMS2mMS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,
       Return[ N[mmu0, numprec] ]
   ];
    If[loops==1, 
       ccc = x^(g0/b0) /.setbeta /.setgamma /.nf->nnf;
       ,
       ccc = (x^(g0/b0)*
	      (Expand[ ms2msc /.setMS2OS ]/.cut[x,loops-1])
	      ) /.setbeta /.setgamma /.nf->nnf;
   ];
    Return[ N[ mmu0 * (ccc/.x->asmu/Pi)/(ccc/.x->asmu0/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 Compute in a first step m_MS(m_MS). Then use m_OS = m_MS * ( 1 + ... )
 for mu=m_MS.
*)

mMS2mSIold[mMS_,asmu_,mu_,nnf_,loops_] := Module[
    {asmMS,mMS2mSI2,xx,res1},
    asmMS = AlphasLam[LamImpl[asmu,mu,nnf,loops],xx,nnf,loops];
    mMS2mSI2 = mMS2mMS[mMS,asmu,asmMS,nnf,loops];
    res1 = FindRoot[ xx == mMS2mSI2 , {xx,mMS} ];

    Return[ N[ mMS2mSI2/.res1 , numprec ] ];
];


(mMS2mSI[mqmu_,asmu_,mu_,nfnum_,nolo_] := Module[
    {mcmc,i,imax,asmc},
    imax = 30;
    mcmc = mqmu;
    For[i=1,i<=imax,i++,
	asmc = AlphasExact[asmu ,mu ,mcmc ,nfnum ,nolo];
        mcmc = mMS2mMS[mqmu, asmu, asmc, nfnum, nolo];
(*         Print[mcmc];						  *)
    ];
    Return[mcmc];
]);

(* ************************************************************ *)

(*
 m_RI = m_MS * ( 1 + ... )
*)

mMS2mRI[mMS_,asmu_,nnf_,loops_] := Module[
    {extmp},
    If[(loops<0)||(loops>4), Print["mMS2mRI: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( rifromms 
	    )/.setMS2OS/.cut[api,loops] /.{nf->nnf} /.num1 /.cf2num;
   ];
    Return[ N[ mMS * (extmp/.api->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 m_MS = m_RI * ( 1 + ... )
*)

mRI2mMS[mRI_,asmu_,nnf_,loops_] := Module[
    {extmp},
    If[(loops<0)||(loops>4), Print["mRI2mMS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( msfromri 
	    )/.setMS2OS/.cut[api,loops] /.{nf->nnf} /.num1 /.cf2num;
   ];
    Return[ N[ mRI * (extmp/.api->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)
(*
 m_RGI = m_MS * ( 1 + ... )
*)
(*
mMS2mRGI[1,0.1,10,5,1]
*)
mMS2mRGI[mMS_,asmu_,nnf_,loops_] := Module[
    {ccc},
    If[(loops<0)||(loops>5), Print["mMS2mRGI: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,
       Return[ N[mMS, numprec] ]
   ];
    If[loops==1, 
       ccc = x^(g0/b0) /.setbeta /.setgamma /.nf->nnf;
       ,
       ccc = (x^(g0/b0)*
	      (Expand[ ms2msc /.setMS2OS ]/.cut[x,loops-1])
	      ) /.setbeta /.setgamma /.nf->nnf;
   ];
    Return[ N[ mMS / (ccc/.x->asmu/Pi) , numprec ] ];
];

mMS2mRGImod[mMS_,asmu_,nnf_,loops_] := Module[
    {ccc,bet0},
    bet0 = 11/4 - nnf/6;
    If[(loops<0)||(loops>5), Print["mMS2mRGImod: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,
       Return[ N[mmu0, numprec] ]
   ];
    If[loops==1, 
       ccc = x^(g0/b0) /.setbeta /.setgamma /.nf->nnf;
       ,
       ccc = (x^(g0/b0)*
	      (Expand[ ms2msc /.setMS2OS ]/.cut[x,loops-1])
	      ) /.setbeta /.setgamma /.nf->nnf;
   ];
    Return[ N[ mMS / (ccc/.x-> 2*bet0*asmu/Pi ) , numprec ] ];
(*                             ^^^^^^ *)
];


(* ************************************************************ *)

(*
 m_MS = m_RGI * ( 1 + ... )
*)

mRGI2mMS[mRGI_,asmu_,nnf_,loops_] := Module[
    {ccc},
    If[(loops<0)||(loops>4), Print["mRGI2mMS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,
       Return[ N[mmu0, numprec] ]
   ];
    If[loops==1, 
       ccc = x^(g0/b0) /.setbeta /.setgamma /.nf->nnf;
       ,
       ccc = (x^(g0/b0)*
	      (Expand[ ms2msc /.setMS2OS ]/.cut[x,loops-1])
	      ) /.setbeta /.setgamma /.nf->nnf;
   ];
    Return[ N[ mRGI * (ccc/.x->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 decoupling of \alpha_s: on-shell scheme (for heavy mass)
 compute: \alpha_s^(nf+1)
 input: \alpha_s^(nf), ...
 the decoupling is performed at "loops-1" order
*)

DecAsUpOS[alsp_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["DecAsUpOS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as6to5os/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as6to5os/.setdec/.api->0 ];
    Return[N[
	alsp*(dectmp/.num1/.nl->nnl/.api->alsp/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(* 
 decoupling of \alpha_s: MS-bar scheme (for heavy mass)
*)

DecAsUpMS[alsp_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["DecAsUpMS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as6to5ms/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as6to5ms/.setdec/.api->0 ];
    Return[N[
	alsp*(dectmp/.num1/.nl->nnl/.api->alsp/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(* 
 decoupling of \alpha_s: scale invariant mass (for heavy mass)
*)

DecAsUpSI[alsp_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["DecAsUpSI: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as6to5si/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as6to5si/.setdec/.api->0 ];
    Return[N[
	alsp*(dectmp/.num1/.nl->nnl/.api->alsp/Pi/.lmmu->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of \alpha_s: on-shell scheme (for heavy mass)
 compute: \alpha_s^(nf-1)
 input: \alpha_s^(nf), ...
 the decoupling is performed at "loops-1" order
*)

DecAsDownOS[als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["DecAsDownOS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as5to6os/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as5to6os/.setdec/.api->0 ];
    Return[N[
	als*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(* 
 decoupling of \alpha_s: MS-bar scheme (for heavy mass)
*)

DecAsDownMS[als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["DecAsDownMS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as5to6ms/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as5to6ms/.setdec/.api->0 ];
    Return[N[
	als*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(* 
 decoupling of \alpha_s: scale invariant mass (for heavy mass)
*)

DecAsDownSI[als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["DecAsDownSI: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as5to6si/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as5to6si/.setdec/.api->0 ];
    Return[N[
	als*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmmu->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of quark mass: on-shell scheme (for heavy mass)
 compute: m_q^(nf)
 input: m_q^(nf-1), ...
 the decoupling is performed at "loops-1" order
*)

DecMqUpOS[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["DecMqUpOS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq6to5os/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq6to5os/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of quark mass: MS-bar scheme (for heavy mass)
 compute: m_q^(nf)
 input: m_q^(nf-1), ...
*)

DecMqUpMS[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["DecMqUpMS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq6to5ms/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq6to5ms/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of quark mass: scale invariant mass (for heavy mass)
 compute: m_q^(nf)
 input: m_q^(nf-1), ...
*)

DecMqUpSI[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["DecMqUpSI: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq6to5si/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq6to5si/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmmu->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decouplingof quark mass: on-shell scheme (for heavy mass)
 compute: m_q^(nf-1)
 input: m_q^(nf), ...
*)

DecMqDownOS[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["DecMqDownOS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq5to6os/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq5to6os/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];


(* ************************************************************ *)

(*
 decoupling of quark mass: MS-bar scheme (for heavy mass)
 compute: m_q^(nf-1)
 input: m_q^(nf), ...
*)

DecMqDownMS[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["DecMqDownMS: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq5to6ms/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq5to6ms/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of quark mass: scale invariant mass (for heavy mass)
 compute: m_q^(nf-1)
 input: m_q^(nf), ...
*)

DecMqDownSI[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>5), Print["DecMqDownSI: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq5to6si/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq5to6si/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmmu->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of \Lambda
 compute: \Lambda^(nl+1)
 input: \Lambda^(nl), ...
*)

DecLambdaUp[lam_,massth_,nnl_,loops_] := Module[
    {exra,exra2,iLb0m,ex1,ex2,ex3,ex4},
    If[(loops<1)||(loops>5), Print["DecLambdaUp: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    exra = Expand[
(L*(1/2 - b0[-1 + nf]/b0times2[nf]) + 
 (b1p[-1 + nf]^4/(6*b0[-1 + nf]^3) + b1p[nf]^4/(3*b0[-1 + nf]^3) + 
   (2*b2p[-1 + nf]^2)/(3*b0[-1 + nf]^3) - (b2p[-1 + nf]*b2p[nf])/
    b0[-1 + nf]^3 + b2p[nf]^2/(3*b0[-1 + nf]^3) - 
   (b1p[-1 + nf]*b3p[-1 + nf])/(6*b0[-1 + nf]^3) + 
   b4p[-1 + nf]/(3*b0[-1 + nf]^3) - b4p[nf]/(3*b0[-1 + nf]^3) + 
   (b2p[-1 + nf]*c2[nf, -1 + nf])/b0[-1 + nf]^3 - 
   (b2p[nf]*c2[nf, -1 + nf])/b0[-1 + nf]^3 - c2[nf, -1 + nf]^2/
    b0[-1 + nf]^3 + b1p[-1 + nf]^2*(-(b2p[-1 + nf]/b0[-1 + nf]^3) + 
     b2p[nf]/b0[-1 + nf]^3 - c2[nf, -1 + nf]/b0[-1 + nf]^3) + 
   b1p[nf]^2*(-(b1p[-1 + nf]^2/b0[-1 + nf]^3) + b2p[-1 + nf]/b0[-1 + nf]^3 - 
     b2p[nf]/b0[-1 + nf]^3 + c2[nf, -1 + nf]/b0[-1 + nf]^3) + 
   b1p[nf]*(b1p[-1 + nf]^3/(2*b0[-1 + nf]^3) - 
     b3p[-1 + nf]/(2*b0[-1 + nf]^3) + (2*b3p[nf])/(3*b0[-1 + nf]^3) - 
     c3[nf, -1 + nf]/b0[-1 + nf]^3) + c4[nf, -1 + nf]/b0[-1 + nf]^3 + 
   (b1p[-1 + nf]^4/b0[-1 + nf]^3 - (b1p[-1 + nf]^2*b1p[nf]^2)/b0[-1 + nf]^3 + 
     (b1p[-1 + nf]*b1p[nf]^3)/b0[-1 + nf]^3 + b1p[-1 + nf]^2*
      (-(b2p[-1 + nf]/b0[-1 + nf]^3) + b2p[nf]/b0[-1 + nf]^3 - 
       c2[nf, -1 + nf]/b0[-1 + nf]^3) + 
     b1p[nf]*(-(b1p[-1 + nf]^3/b0[-1 + nf]^3) + b1p[-1 + nf]*
        ((2*b2p[-1 + nf])/b0[-1 + nf]^3 - (2*b2p[nf])/b0[-1 + nf]^3 + 
         (2*c2[nf, -1 + nf])/b0[-1 + nf]^3)) + 
     b1p[-1 + nf]*(-(b3p[-1 + nf]/b0[-1 + nf]^3) + b3p[nf]/b0[-1 + nf]^3 - 
       (2*c3[nf, -1 + nf])/b0[-1 + nf]^3))*Log[L] + 
   (b1p[-1 + nf]^4/(2*b0[-1 + nf]^3) - (3*b1p[-1 + nf]^3*b1p[nf])/
      (2*b0[-1 + nf]^3) + (b1p[-1 + nf]^2*b1p[nf]^2)/b0[-1 + nf]^3 + 
     b1p[-1 + nf]^2*(b2p[-1 + nf]/b0[-1 + nf]^3 - b2p[nf]/b0[-1 + nf]^3 + 
       c2[nf, -1 + nf]/b0[-1 + nf]^3))*Log[L]^2 + 
   (-b1p[-1 + nf]^4/(3*b0[-1 + nf]^3) + (b1p[-1 + nf]^3*b1p[nf])/
      (3*b0[-1 + nf]^3))*Log[L]^3)/(L^3*b0times2[nf]) + 
 ((-b1p[-1 + nf] + b1p[nf])*Log[L] + 
   (-b1p[-1 + nf]^2 + b1p[nf]^2 + b2p[-1 + nf] - b2p[nf] + c2[nf, -1 + nf] + 
     (-b1p[-1 + nf]^2 + b1p[-1 + nf]*b1p[nf])*Log[L])/Lb0m + 
   (-b1p[-1 + nf]^3/2 - b1p[nf]^3/2 + b3p[-1 + nf]/2 - b3p[nf]/2 + 
     b1p[nf]*(b1p[-1 + nf]^2 - b2p[-1 + nf] + b2p[nf] - c2[nf, -1 + nf]) + 
     c3[nf, -1 + nf] + (b1p[-1 + nf]^2*b1p[nf] - b1p[-1 + nf]*b1p[nf]^2 + 
       b1p[-1 + nf]*(-b2p[-1 + nf] + b2p[nf] - c2[nf, -1 + nf]))*Log[L] + 
     (b1p[-1 + nf]^3/2 - (b1p[-1 + nf]^2*b1p[nf])/2)*Log[L]^2)/Lb0m^2 + 
   b1p[nf]*Log[b0[-1 + nf]/b0[nf]])/b0times2[nf]
)/.{b0times2[nn_]->2*b0[nn]}
		  /.{Lb0m->(L*b0[nf-1])}
	      ];
    ex1 = Coefficient[exra,L,1];
    ex2 = Coefficient[exra,L,0];
    ex3 = Coefficient[exra,L,-1];
    ex4 = Coefficient[exra,L,-2];
    ex5 = Coefficient[exra,L,-3];

    If[ loops==5 , exra = ex1*L+ex2+ex3/L+ex4/L^2+ex5/L^3 ];
    If[ loops==4 , exra = ex1*L+ex2+ex3/L+ex4/L^2 ];
    If[ loops==3 , exra = ex1*L+ex2+ex3/L ];
    If[ loops==2 , exra = ex1*L+ex2 ];
    If[ loops==1 , exra = ex1*L ];

    exra2 = Expand[ exra  
		    /.{nf->nnl+1}
		    /.{b0[nn_] -> (b0/.{setbeta/.nf->nn}),
		       b1p[nn_]-> ((b1/b0)/.{setbeta/.nf->nn}),
		       b2p[nn_]-> ((b2/b0)/.{setbeta/.nf->nn}),
		       b3p[nn_]-> ((b3/b0)/.{setbeta/.nf->nn}),
		       b4p[nn_]-> ((b4/b0)/.{setbeta/.nf->nn})}
		    /.{c2[nn__]   -> -Coefficient[as5to6si /. setdec /. lmmu->0,api,2]/.num1 }
		    /.{c3[mm_,nn_]-> -Coefficient[as5to6si /. setdec /. lmmu->0/. nl->nn,api,3]/.num1 }
		    /.{c4[mm_,nn_]-> -Coefficient[as5to6si /. setdec /. lmmu->0/. nl->nn,api,4]/.num1 }
		][[1]];
    Return[ lam*Exp[exra2 /. L->Log[massth^2/lam^2]] ] ;
];

(* ************************************************************ *)

(*
 decoupling of \Lambda
 compute: \Lambda^(nl)
 input: \Lambda^(nl+1), ...
*)

DecLambdaDown[lam_,massth_,nnl_,loops_] := Module[
    {exra,exra2,iLb0},
    If[(loops<1)||(loops>5), Print["DecLambdaDown: PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    exra = Expand[
(
L*(1/2 - b0[nf]/b0times2[-1 + nf]) + 
 (b1p[-1 + nf]^4/(3*b0[nf]^3) + b1p[nf]^4/(6*b0[nf]^3) + 
   b2p[-1 + nf]^2/(3*b0[nf]^3) - (b2p[-1 + nf]*b2p[nf])/b0[nf]^3 + 
   (2*b2p[nf]^2)/(3*b0[nf]^3) - (b1p[nf]*b3p[nf])/(6*b0[nf]^3) - 
   b4p[-1 + nf]/(3*b0[nf]^3) + b4p[nf]/(3*b0[nf]^3) - 
   (b2p[-1 + nf]*c2[-1 + nf, nf])/b0[nf]^3 + (b2p[nf]*c2[-1 + nf, nf])/
    b0[nf]^3 - c2[-1 + nf, nf]^2/b0[nf]^3 + 
   b1p[nf]^2*(b2p[-1 + nf]/b0[nf]^3 - b2p[nf]/b0[nf]^3 - 
     c2[-1 + nf, nf]/b0[nf]^3) + b1p[-1 + nf]^2*(-(b1p[nf]^2/b0[nf]^3) - 
     b2p[-1 + nf]/b0[nf]^3 + b2p[nf]/b0[nf]^3 + c2[-1 + nf, nf]/b0[nf]^3) + 
   b1p[-1 + nf]*(b1p[nf]^3/(2*b0[nf]^3) + (2*b3p[-1 + nf])/(3*b0[nf]^3) - 
     b3p[nf]/(2*b0[nf]^3) - c3[-1 + nf, nf]/b0[nf]^3) + 
   c4[-1 + nf, nf]/b0[nf]^3 + ((b1p[-1 + nf]^3*b1p[nf])/b0[nf]^3 - 
     (b1p[-1 + nf]^2*b1p[nf]^2)/b0[nf]^3 + b1p[nf]^4/b0[nf]^3 + 
     b1p[nf]^2*(b2p[-1 + nf]/b0[nf]^3 - b2p[nf]/b0[nf]^3 - 
       c2[-1 + nf, nf]/b0[nf]^3) + b1p[-1 + nf]*(-(b1p[nf]^3/b0[nf]^3) + 
       b1p[nf]*((-2*b2p[-1 + nf])/b0[nf]^3 + (2*b2p[nf])/b0[nf]^3 + 
         (2*c2[-1 + nf, nf])/b0[nf]^3)) + 
     b1p[nf]*(b3p[-1 + nf]/b0[nf]^3 - b3p[nf]/b0[nf]^3 - 
       (2*c3[-1 + nf, nf])/b0[nf]^3))*Log[L] + 
   ((b1p[-1 + nf]^2*b1p[nf]^2)/b0[nf]^3 - (3*b1p[-1 + nf]*b1p[nf]^3)/
      (2*b0[nf]^3) + b1p[nf]^4/(2*b0[nf]^3) + 
     b1p[nf]^2*(-(b2p[-1 + nf]/b0[nf]^3) + b2p[nf]/b0[nf]^3 + 
       c2[-1 + nf, nf]/b0[nf]^3))*Log[L]^2 + 
   ((b1p[-1 + nf]*b1p[nf]^3)/(3*b0[nf]^3) - b1p[nf]^4/(3*b0[nf]^3))*Log[L]^3)/
  (L^3*b0times2[-1 + nf]) + ((b1p[-1 + nf] - b1p[nf])*Log[L] + 
   (b1p[-1 + nf]^2 - b1p[nf]^2 - b2p[-1 + nf] + b2p[nf] + c2[-1 + nf, nf] + 
     (b1p[-1 + nf]*b1p[nf] - b1p[nf]^2)*Log[L])/Lb0 + 
   (-b1p[-1 + nf]^3/2 - b1p[nf]^3/2 - b3p[-1 + nf]/2 + b3p[nf]/2 + 
     b1p[-1 + nf]*(b1p[nf]^2 + b2p[-1 + nf] - b2p[nf] - c2[-1 + nf, nf]) + 
     c3[-1 + nf, nf] + (-(b1p[-1 + nf]^2*b1p[nf]) + b1p[-1 + nf]*b1p[nf]^2 + 
       b1p[nf]*(b2p[-1 + nf] - b2p[nf] - c2[-1 + nf, nf]))*Log[L] + 
     (-(b1p[-1 + nf]*b1p[nf]^2)/2 + b1p[nf]^3/2)*Log[L]^2)/Lb0^2 + 
   b1p[-1 + nf]*Log[b0[nf]/b0[-1 + nf]])/b0times2[-1 + nf]
)/.{b0times2[nn_]->2*b0[nn]}
		  /.{Lb0->(L*b0[nf])}
	      ];
    
    ex1 = Coefficient[exra,L,1];
    ex2 = Coefficient[exra,L,0];
    ex3 = Coefficient[exra,L,-1];
    ex4 = Coefficient[exra,L,-2];
    ex5 = Coefficient[exra,L,-3];

    If[ loops==5 , exra = ex1*L+ex2+ex3/L+ex4/L^2+ex5/L^3 ];
    If[ loops==4 , exra = ex1*L+ex2+ex3/L+ex4/L^2 ];
    If[ loops==3 , exra = ex1*L+ex2+ex3/L ];
    If[ loops==2 , exra = ex1*L+ex2 ];
    If[ loops==1 , exra = ex1*L ];
    
    exra2 = Expand[ exra  /.{iLb0->1/(L*b0[nf])}
		    /.{nf->nnl+1}
		    /.{b0[nn_] -> (b0/.{setbeta/.nf->nn}),
		       b1p[nn_]-> ((b1/b0)/.{setbeta/.nf->nn}),
		       b2p[nn_]-> ((b2/b0)/.{setbeta/.nf->nn}),
		       b3p[nn_]-> ((b3/b0)/.{setbeta/.nf->nn}),
		       b4p[nn_]-> ((b4/b0)/.{setbeta/.nf->nn})}
		    /.{c2[nn__]   -> Coefficient[as5to6si /. setdec /. lmmu->0,api,2]/.num1 }
		    /.{c3[mm_,nn_]-> Coefficient[as5to6si /. setdec /. lmmu->0/. nl->mm,api,3]/.num1 }
		    /.{c4[mm_,nn_]-> Coefficient[as5to6si /. setdec /. lmmu->0/. nl->mm,api,4]/.num1 }
		][[1]];
    
    Return[lam*Exp[exra2 /. L->Log[massth^2/lam^2]] ] ;
];

(* ************************************************************ *)

(* 
 running-decoupling-running-decoupling-running-...
 input: \alpha_s^(l)(mu1), ...
 output: \alpha_s^(h)(mu2)
*)

AlL2AlH[asl_,mu1_,decpar_,mu2_,loops_] := Module[
    {i,asini,muini,res1,res2,res3,decpar2},
    asini=asl;
    muini=mu1;
    decpar2 = Sort[decpar];
    For[i=2,i<=Length[decpar2],i++,
	If[(decpar2[[i]][[1]]-decpar2[[i-1]][[1]]=!=1) ||
	   (decpar2[[i]][[1]]-decpar2[[i-1]][[1]]==0) ,
	   Print["WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."];
	   Abort[];
       ];
    ];
    For[i=1,i<=Length[decpar2],i++,
	res1 = AlphasExact[asini,muini,decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	res2 = DecAsUpOS[res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			 decpar2[[i]][[1]]-1,loops];
	asini = res2;
	muini = decpar2[[i]][[3]];
    ];
    res3 = AlphasExact[asini,muini,mu2,decpar2[[i-1]][[1]],loops];
    Return[res3];
];

(* ************************************************************ *)

(* 
 running-decoupling-running-decoupling-running-...
 input: \alpha_s^(h)(mu1), ...
 output: \alpha_s^(l)(mu2)
*)

AlH2AlL[ash_,mu1_,decpar_,mu2_,loops_] := Module[
    {i,asini,muini,res1,res2,res3,decpar2},
    asini=ash;
    muini=mu1;
    decpar2 = Reverse[Sort[decpar]];
    For[i=2,i<=Length[decpar2],i++,
	If[(decpar2[[i]][[1]]-decpar2[[i-1]][[1]]=!=-1) ||
	   (decpar2[[i]][[1]]-decpar2[[i-1]][[1]]==0) ,
	   Print["WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."];
	   Abort[];
       ];
    ];
    For[i=1,i<=Length[decpar2],i++,
	res1 = AlphasExact[asini,muini,decpar2[[i]][[3]],
			   decpar2[[i]][[1]],loops];
	res2 = DecAsDownOS[res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	asini = res2;
	muini = decpar2[[i]][[3]];
    ];
    res3 = AlphasExact[asini,muini,mu2,decpar2[[i-1]][[1]]-1,loops];
    Return[res3];
];

(* ************************************************************ *)

(* 
 example: running-decoupling-running-decoupling-running-...
 input: m_q^(l)(mu1), ...
 output: m_q^(h)(mu2)
*)

mL2mH[mql_,asl_,mu1_,decpar_,mu2_,loops_] := Module[
    {i,asini,muini,res1,res2,res3,res4,res5,res6,decpar2},
    asini=asl;
    muini=mu1;
    mqini=mql;
    decpar2 = Sort[decpar];
    For[i=2,i<=Length[decpar2],i++,
	If[(decpar2[[i]][[1]]-decpar2[[i-1]][[1]]=!=1) ||
	   (decpar2[[i]][[1]]-decpar2[[i-1]][[1]]==0) ,
	   Print["WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."];
	   Abort[];
       ];
    ];
    For[i=1,i<=Length[decpar2],i++,
	(* res1 and res2: alpha_s *)
	(* res3 and res4: masses  *)
	res1 = AlphasExact[asini,muini,decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	res3 = mMS2mMS[mqini,asini,res1,decpar2[[i]][[1]]-1,loops];
	res2 = DecAsUpOS[res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			 decpar2[[i]][[1]]-1,loops];
	res4 = DecMqUpOS[res3,res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			 decpar2[[i]][[1]]-1,loops];
	asini = res2;
	mqini = res4;
	muini = decpar2[[i]][[3]];
    ];
    res5 = AlphasExact[asini,muini,mu2,decpar2[[i-1]][[1]],loops];
    res6 = mMS2mMS[mqini,asini,res5,decpar2[[i-1]][[1]],loops];
    Return[res6];
];

(* ************************************************************ *)

(* 
 example: running-decoupling-running-decoupling-running-...
 input: m_q^(h)(mu1), ...
 output: m_q^(l)(mu2)
*)

mH2mL[mqh_,ash_,mu1_,decpar_,mu2_,loops_] := Module[
    {i,asini,muini,res1,res2,res3,res4,res5,res6,decpar2},
    asini=ash;
    muini=mu1;
    mqini=mqh;
    decpar2 = Reverse[Sort[decpar]];
    For[i=2,i<=Length[decpar2],i++,
	If[(decpar2[[i]][[1]]-decpar2[[i-1]][[1]]=!=-1) ||
	   (decpar2[[i]][[1]]-decpar2[[i-1]][[1]]==0) ,
	   Print["WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."];
	   Abort[];
       ];
    ];
    For[i=1,i<=Length[decpar2],i++,
	(* res1 and res2: alpha_s *)
	(* res3 and res4: masses  *)
	res1 = AlphasExact[asini,muini,decpar2[[i]][[3]],
			   decpar2[[i]][[1]],loops];
	res3 = mMS2mMS[mqini,asini,res1,decpar2[[i]][[1]],loops];
	res2 = DecAsDownOS[res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	res4 = DecMqDownOS[res3,res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	asini = res2;
	mqini = res4;
	muini = decpar2[[i]][[3]];
    ];
    res5 = AlphasExact[asini,muini,mu2,decpar2[[i-1]][[1]]-1,loops];
    res6 = mMS2mMS[mqini,asini,res5,decpar2[[i-1]][[1]]-1,loops];
    Return[res6];
];

(* ************************************************************ *)

(*
 AsRunDec[als, mu0, mu, l] computes \alpha_s^(m)(mu) from the knowledge
 of \alpha_s^(n)(mu0)  
*)

AsRunDec[als_,mu0_,mu_,loops_] := Module[
    {m,n,kk,decpar},
    If[mu>=Global`Mt/.Global`NumDef,m=6,
       If[mu>=Global`Mb/.Global`NumDef,m=5,
	  If[mu>=Global`Mc/.Global`NumDef,m=4,
	     If[mu<Global`Mc/.Global`NumDef,m=3
		]]]];
    If[mu0>=Global`Mt/.Global`NumDef,n=6,
       If[mu0>=Global`Mb/.Global`NumDef,n=5,
	  If[mu0>=Global`Mc/.Global`NumDef,n=4,
	     If[mu0<Global`Mc/.Global`NumDef,n=3
		]]]];
    If[m==n,
       (* Return[{AlphasExact[als,mu0,mu,n,loops],m,n}]; *)
       Return[AlphasExact[als,mu0,mu,n,loops]];
   ];
    decpar = {};
    If[m>n,
       For[kk=n+1,kk<=m,kk++,
	   If[kk==4,decpar=Join[decpar,{{4,Global`Mc/.Global`NumDef,
					 Global`Mc/.Global`NumDef}}]];
	   If[kk==5,decpar=Join[decpar,{{5,Global`Mb/.Global`NumDef,
					 Global`Mb/.Global`NumDef}}]];
	   If[kk==6,decpar=Join[decpar,{{6,Global`Mt/.Global`NumDef,
					 Global`Mt/.Global`NumDef}}]];
       ];
       (* Return[{AlL2AlH[als,mu0,decpar,mu,loops],m,n}]; *)
       Return[AlL2AlH[als,mu0,decpar,mu,loops]];
   ];
    If[m<n,
       For[kk=n,kk>=m+1,kk--,
	   If[kk==4,decpar=Join[decpar,{{4,Global`Mc/.Global`NumDef,
					 Global`Mc/.Global`NumDef}}]];
	   If[kk==5,decpar=Join[decpar,{{5,Global`Mb/.Global`NumDef,
					 Global`Mb/.Global`NumDef}}]];
	   If[kk==6,decpar=Join[decpar,{{6,Global`Mt/.Global`NumDef,
					 Global`Mt/.Global`NumDef}}]];
       ];
       (* Print[decpar]; *)
       (* Return[{AlH2AlL[als,mu0,decpar,mu,loops],m,n}]; *)
       Return[AlH2AlL[als,mu0,decpar,mu,loops]];
   ];
];

(* ************************************************************ *)

(* 
AsmMSrunexact: solve simultaneously RGE for \alpha_s and m_q.
*)

AsmMSrunexact[mmu0_, asmu0_, mu0_, mu_, nnf_, l_] := Module[
    {res,beta,gammam,x,xl,api},
    If[l>5,Print["AsmMSrunexact: l>5 not implemented."];Abort[];];

    beta   = Expand[-api[x]^2*(xl*b0+xl^2*api[x]*b1+
			       xl^3*api[x]^2*b2+xl^4*api[x]^3*b3+xl^5*api[x]^4*b4)
		    ]/.cut[xl,l]/.{xl->1};
    gammam = Expand[-api[x]  *(xl*g0+xl^2*api[x]*g1+
			       xl^3*api[x]^2*g2+xl^4*api[x]^3*g3+xl^5*api[x]^4*g4)
		    ]/.cut[xl,l]/.{xl->1};

    res = NDSolve[ {x*api'[x]      == beta /. setbeta  (*/.{b4->Global`b4num}*) /. nf->nnf ,
		    x/mq[x]*mq'[x] == gammam /. setgamma /. nf->nnf,
		    api[mu0^2]     == asmu0/Pi, 
		    mq[mu0^2]      == mmu0},    {api[x],mq[x]}, {x,mu^2,mu0^2}
(*			,
		   WorkingPrecision->26,
		   AccuracyGoal->16,
		   PrecisionGoal->16,
		   MaxSteps->5000 *)
		   ];
    Return[{mq[x],Pi*api[x]}/.res[[1]]/.{x->mu^2}];
];

(* ************************************************************ *)

(*
 mOS2mMS[mOS, nf, l] computes MS-bar mass corresponding to mOS
*)

mOS2mMS[mOS_,nnf_,loops_] := Module[
    {as},
    as=AsRunDec[Global`asMz/.Global`NumDef,Global`Mz/.Global`NumDef,mOS,loops];
    Return[mOS2mMS[mOS,{},as,mOS,nnf,loops]];
];

(* ************************************************************ *)

(*
 mMS2mOS[mOS, nf, l] computes the on-shell mass corresponding to mMS
*)

mMS2mOS[mMS_,nnf_,loops_] := Module[
    {as},
    as=AsRunDec[Global`asMz/.Global`NumDef,Global`Mz/.Global`NumDef,mMS,loops];
    Return[mMS2mOS[mMS,{},as,mMS,nnf,loops]];
];

(* ************************************************************ *)

(* 
 input: \mu_c^(4)
 output: m_c(MZ)^(5)
*)

Mc5Mzfrommuc4[asMz_,muc4_,Mb_,mub_,Mz_,loops_] := Module[
    {alsmuth,alsmuthp,alsmuc,mcthp,mcth,mcMZ},

    alsmuth = AlphasExact[asMz,Mz,mub,5,loops];
    alsmuthp = DecAsDownOS[alsmuth,Mb,mub,4,loops];
    alsmuc = AlphasExact[alsmuthp,mub,muc4,4,loops];

    mcthp = mMS2mMS[muc4,alsmuc,alsmuthp,4,loops];
    mcth = DecMqUpOS[mcthp,alsmuthp,Mb,mub,4,loops];
    mcMZ = mMS2mMS[mcth,alsmuth,asMz,5,loops];

(*    Print["alsmuc, alsmuthp, alsmuth: ", alsmuc," ",alsmuthp," ",alsmuth]; *)
(*    Print["mcth, mcthp, mcMZ: ", mcth," ",mcthp," ",mcMZ]; *)
    
    Return[mcMZ];
];

(* ************************************************************ *)

(* 
 compare running of "exact" vs "\lambda" 
 input: \alpha_s^(5)(Mz)
 output: \alpha_s^(5)(Mb)
*)

AsMbfromAsMz[asMz_,Mb_,loops_] := Module[
    {res1,res2},
    res1 = AlphasExact[asMz,Global`Mz/.Global`NumDef,Mb,5,loops];
    res2 = AlphasLam[LamImpl[asMz,Global`Mz/.Global`NumDef,5,loops],Mb,5,loops];
    Return[ {res1, res2}];
];

AsMbfromAsMz[asMz_,loops_] := AsMbfromAsMz[asMz,Global`Mb/.Global`NumDef,
					   loops];

(* ************************************************************ *)

(* 
 running-decoupling-running
 input: \alpha_s^(5)(Mz)
 output: \alpha_s^(4)(Mc)
*)

As4McfromAs5Mz[asMz_,Mb_,mub_,Mc_,loops_] := Module[
    {res1,res2,res3},
    res1 = AlphasExact[asMz,Global`Mz/.Global`NumDef,mub,5,loops];
    res2 = DecAsDownOS[res1,Mb,mub,4,loops];
    res3 = AlphasExact[res2,mub,Mc,4,loops];
    Return[res3];
];

(* ************************************************************ *)

(* Additions August 2015: *)

(* ************************************************************ *)

(*
 mOS2mPS[mOS_, mq_, asmu_, Mu_, Muf_, nl_, loops_] computes the PS mass from the OS mass
*)

mOS2mPS[mOS_, mq_, asmu_, Mu_, Muf_, nl_, loops_] := 
    Module[{extmp, api}, If[loops < 0 || loops > 4, 
       Print["mOS2mPS: PROCEDURE IS NOT IMPLEMENTED FOR ", loops, " LOOPS."]; 
        Return[]]; If[loops == 0, extmp = 0, 
       extmp = Expand[deltamf[api, Muf, Mu, nl]] /. cut[api, loops]; ]; 
      Return[N[mOS - (extmp /. api -> asmu/Pi), numprec]]; ];
  
deltamf[api_, Muf_, Mu_, nl_] := 
    Muf*((4*api)/3 + api^2*(97/9 + nl*(-22/27 - (2*Log[Mu^2/Muf^2])/9) + 
        (11*Log[Mu^2/Muf^2])/3) + api^3*(33623/216 + 3*Pi^2 - (3*Pi^4)/16 + 
        (610*Log[Mu^2/Muf^2])/9 + (121*Log[Mu^2/Muf^2]^2)/12 + 
        nl^2*(157/243 + (22*Log[Mu^2/Muf^2])/81 + Log[Mu^2/Muf^2]^2/27) + 
        nl*(-7145/324 - (493*Log[Mu^2/Muf^2])/54 - (11*Log[Mu^2/Muf^2]^2)/9 - 
          (13*Zeta[3])/9) + (11*Zeta[3])/2) + 
      api^4*(3567.723056629293 + (7271*Log[Mu^2/Muf^2]^2)/24 + 
        (1331*Log[Mu^2/Muf^2]^3)/48 + nl^3*(-2951/4374 - 
          (157*Log[Mu^2/Muf^2])/486 - (11*Log[Mu^2/Muf^2]^2)/162 - 
          Log[Mu^2/Muf^2]^3/162) + nl*(-701.2303148875468 - 
          (8485*Log[Mu^2/Muf^2]^2)/144 - (121*Log[Mu^2/Muf^2]^3)/24 + 
          Log[Mu^2/Muf^2]*(-253189/864 - (3*Pi^2)/2 + (3*Pi^4)/32 - 
            (44*Zeta[3])/3)) + nl^2*(1751971/46656 + Pi^4/135 + 
          (773*Log[Mu^2/Muf^2]^2)/216 + (11*Log[Mu^2/Muf^2]^3)/36 + 
          Log[Mu^2/Muf^2]*(15355/864 + (13*Zeta[3])/18) + 
          (259*Zeta[3])/108) + Log[Mu^2/Muf^2]*(26125/18 + (99*Pi^2)/4 - 
          (99*Pi^4)/64 + (363*Zeta[3])/8)));

(* ************************************************************ *)

(*
 mMS2mPS[mMS, mq, asmu, mu, muf, nl, loops] computes the PS mass from the MSbar mass
*)

mMS2mPS[mMS_, mq_, asmu_, Mu_, Muf_, nl_, loops_] := 
    mMS2mPS[mMS, mq, asmu, Mu, Muf, nl, loops, 1, "no"];

mMS2mPS[mMS_, mq_, asmu_, Mu_, Muf_, nl_, loops_, marker_String] := 
    mMS2mPS[mMS, mq, asmu, Mu, Muf, nl, loops, 1, marker];

mMS2mPS[mMS_, mq_, asmu_, Mu_, Muf_, nl_, loops_, fdelm_/;NumericQ[fdelm]] := 
    mMS2mPS[mMS, mq, asmu, Mu, Muf, nl, loops, fdelm, "no"];

mMS2mPS[mMS_, mq_, asmu_, Mu_, Muf_, nl_, loops_, fdelm_, marker_String] := 
    Module[{extmp, api, delmuf, exmOS, exmOSprime, apiprime, exals, mark}, 
           If[ marker == "no", mark = 1, mark = "RunDecXMS2OS" ];
	   If[loops < 0 || loops > 4, 
	      Print["mMS2mPS: PROCEDURE IS NOT IMPLEMENTED FOR ",  loops, " LOOPS."]; 
	      Return[]]; 
	   If[loops == 0, extmp = mMS, 
	      delmuf = deltamf[api, Muf, Mu, nl]; 
	      exmOS = mMS2mOS[mMS, {}, api*Pi, Mu, nl + 1, loops, fdelm];
	      exals = (1*DecAsUpMS[apiprime*Pi, mMS, Mu, nl, 4])/(apiprime*Pi); 
	      exmOSprime = Normal[Series[
		  exmOS /. {api -> apiprime*exals}, {apiprime, 0, loops}]]; 
	      extmp = Expand[mark * (exmOSprime /. apiprime -> api) - delmuf] /. 
		  cut[api, loops]; ]; Return[N[extmp /. api -> asmu/Pi, numprec]]; 
       ];

(* ************************************************************ *)

(*
 mPS2mMS[] computes the MSbar mass from the PS mass
*)

mPS2mMS[mPS_, mq_, asnlmu_, Mu_, Muf_, nl_, loops_] := 
    mPS2mMS[mPS, mq, asnlmu, Mu, Muf, nl, loops, 1];
 
mPS2mMS[mPS_, mq_, asnlmu_, Mu_, Muf_, nl_, loops_, fdelm_] :=
    Module[{mtmp, mMS}, mtmp = mMS2mPS[mMS, mq, asnlmu, Mu, Muf, nl, loops, fdelm];
    Return[mMS /. FindRoot[mPS == mtmp, {mMS, mPS}]]; ];

(* ************************************************************ *)

(*
 mPS2mSI[] computes the SI mass from the PS mass
*)

mPS2mSI[mPS_, mq_, asfct_, Muf_, nl_, loops_] := 
    mPS2mSI[mPS, mq, asfct, Muf, nl, loops, 1];
 
mPS2mSI[mPS_, mq_, asfct_, Muf_, nl_, loops_, fdelm_] := 
    Module[{mMS1, mMS, acc}, mMS1 = 0; mMS = mPS; acc = 1.*^-6; 
      While[Abs[mMS1 - mMS] > acc, mMS1 = mMS; 
        mMS = mPS2mMS[mPS, mq, asfct[mMS1], mMS1, Muf, nl, loops, fdelm]; ]; 
      Return[mMS]; ];

(* ************************************************************ *)

(*
 mOS2mPS[mOS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, PRIME] computes the RS or RS' mass from the OS mass
*)

mOS2mRS[mOS_, mq_, asmu_, Mu_, nuf_, nl_, loops_] := 
    mOS2mRS[mOS, mq, asmu, Mu, nuf, nl, loops, "no"]

mOS2mRSp[mOS_, mq_, asmu_, Mu_, nuf_, nl_, loops_] := 
    mOS2mRS[mOS, mq, asmu, Mu, nuf, nl, loops, "yes"]
 
mOS2mRS[mOS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, PRIME_] := 
    Module[{extmp, api}, If[loops < 0 || loops > 4, 
       Print["mOS2mRS: PROCEDURE IS NOT IMPLEMENTED FOR ", loops, " LOOPS."]; 
        Return[]]; If[loops == 0, extmp = 0, 
       If[PRIME === "yes", extmp = exOS2RSp[api, Mu, nuf, nl]; , 
         extmp = exOS2RS[api, Mu, nuf, nl]; ]; 
        extmp = extmp /. cut[api, loops]; ]; 
      Return[N[mOS - (extmp /. api -> asmu/Pi), numprec]]; ];
 
exOS2RSp[aapi_, mmu_, nnuf_, nnl_] = 
    (aapi^2*(11 - (2*nnl)/3)*nnuf*Pi*((ctil[3, nnl]*Gamma[-1 + nu[nnl]])/
         Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[nu[nnl]])/
         Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[1 + nu[nnl]])/
         Gamma[nu[nnl]] + Gamma[2 + nu[nnl]]/Gamma[1 + nu[nnl]])*Nm[nnl])/2 + 
     aapi^3*(((11 - (2*nnl)/3)^2*nnuf*Pi*((ctil[3, nnl]*Gamma[nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[2 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[3 + nu[nnl]]/Gamma[1 + nu[nnl]])*Nm[nnl])/
        4 + ((11 - (2*nnl)/3)*nnuf*Pi*((ctil[3, nnl]*Gamma[-1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[2 + nu[nnl]]/Gamma[1 + nu[nnl]])*
         ((11*Log[mmu^2/nnuf^2])/2 - (nnl*Log[mmu^2/nnuf^2])/3)*Nm[nnl])/2) + 
     aapi^4*(((11 - (2*nnl)/3)^3*nnuf*Pi*((ctil[3, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[2 + nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[3 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[4 + nu[nnl]]/Gamma[1 + nu[nnl]])*Nm[nnl])/
        8 + ((11 - (2*nnl)/3)^2*nnuf*Pi*((ctil[3, nnl]*Gamma[nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[2 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[3 + nu[nnl]]/Gamma[1 + nu[nnl]])*
         ((33*Log[mmu^2/nnuf^2])/4 - (nnl*Log[mmu^2/nnuf^2])/2)*Nm[nnl])/4 + 
       ((11 - (2*nnl)/3)*nnuf*Pi*((ctil[3, nnl]*Gamma[-1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[2 + nu[nnl]]/Gamma[1 + nu[nnl]])*
         (612*Log[mmu^2/nnuf^2] - 76*nnl*Log[mmu^2/nnuf^2] + 
          1089*Log[mmu^2/nnuf^2]^2 - 132*nnl*Log[mmu^2/nnuf^2]^2 + 
          4*nnl^2*Log[mmu^2/nnuf^2]^2)*Nm[nnl])/96);

ctil[0,nf_] = 1;
(* 22jul20: ctil[3,{3,4,5}] updated *)
ctil[1,3] = -0.1638; ctil[2,3] = 0.2372; ctil[3,3] = 0.0217337; (* -0.1205;*) nu[3] = 0.3951;
ctil[1,4] = -0.1054; ctil[2,4] = 0.2736; ctil[3,4] = 0.0423924; (* -0.1610;*) nu[4] = 0.3696;
ctil[1,5] =  0.0238; ctil[2,5] = 0.3265; ctil[3,5] = 0.0595309; (* -0.2681;*) nu[5] = 0.3289;
Nm[3] = 0.563; (*\pm 26*)
Nm[4] = 0.547; (*\pm 33*)
Nm[5] = 0.527; (*\pm 51*)
 
exOS2RS[aapi_, mmu_, nnuf_, nnl_] = 
    aapi*nnuf*Pi*(1 + ctil[1, nnl] + ctil[2, nnl] + ctil[3, nnl])*Nm[nnl] + 
     aapi^3*(((11 - (2*nnl)/3)^2*nnuf*Pi*((ctil[3, nnl]*Gamma[nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[2 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[3 + nu[nnl]]/Gamma[1 + nu[nnl]])*Nm[nnl])/
        4 + nnuf*Pi*(1 + ctil[1, nnl] + ctil[2, nnl] + ctil[3, nnl])*
        Log[mmu^2/nnuf^2]*((102 - (38*nnl)/3)/16 + 
         ((11 - (2*nnl)/3)^2*Log[mmu^2/nnuf^2])/16)*Nm[nnl] + 
       ((11 - (2*nnl)/3)*nnuf*Pi*((ctil[3, nnl]*Gamma[-1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[2 + nu[nnl]]/Gamma[1 + nu[nnl]])*
         ((11*Log[mmu^2/nnuf^2])/2 - (nnl*Log[mmu^2/nnuf^2])/3)*Nm[nnl])/2) + 
     aapi^2*(((11 - (2*nnl)/3)*nnuf*Pi*((ctil[3, nnl]*Gamma[-1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[2 + nu[nnl]]/Gamma[1 + nu[nnl]])*Nm[nnl])/
        2 - (nnuf*Pi*(1 + ctil[1, nnl] + ctil[2, nnl] + ctil[3, nnl])*
         (-33*Log[mmu^2/nnuf^2] + 2*nnl*Log[mmu^2/nnuf^2])*Nm[nnl])/12) + 
     aapi^4*(((11 - (2*nnl)/3)^3*nnuf*Pi*((ctil[3, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[2 + nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[3 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[4 + nu[nnl]]/Gamma[1 + nu[nnl]])*Nm[nnl])/
        8 + ((11 - (2*nnl)/3)^2*nnuf*Pi*((ctil[3, nnl]*Gamma[nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[2 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[3 + nu[nnl]]/Gamma[1 + nu[nnl]])*
         ((33*Log[mmu^2/nnuf^2])/4 - (nnl*Log[mmu^2/nnuf^2])/2)*Nm[nnl])/4 + 
       nnuf*Pi*(1 + ctil[1, nnl] + ctil[2, nnl] + ctil[3, nnl])*
        Log[mmu^2/nnuf^2]*((2857/2 - (5033*nnl)/18 + (325*nnl^2)/54)/64 + 
         (5*(102 - (38*nnl)/3)*(11 - (2*nnl)/3)*Log[mmu^2/nnuf^2])/128 + 
         ((11 - (2*nnl)/3)^3*Log[mmu^2/nnuf^2]^2)/64)*Nm[nnl] + 
       ((11 - (2*nnl)/3)*nnuf*Pi*((ctil[3, nnl]*Gamma[-1 + nu[nnl]])/
           Gamma[-2 + nu[nnl]] + (ctil[2, nnl]*Gamma[nu[nnl]])/
           Gamma[-1 + nu[nnl]] + (ctil[1, nnl]*Gamma[1 + nu[nnl]])/
           Gamma[nu[nnl]] + Gamma[2 + nu[nnl]]/Gamma[1 + nu[nnl]])*
         (612*Log[mmu^2/nnuf^2] - 76*nnl*Log[mmu^2/nnuf^2] + 
          1089*Log[mmu^2/nnuf^2]^2 - 132*nnl*Log[mmu^2/nnuf^2]^2 + 
          4*nnl^2*Log[mmu^2/nnuf^2]^2)*Nm[nnl])/96);

(* ************************************************************ *)
(*
mMS2mRS[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, PRIME_, fdelm_]
   computes the RS mass form the MS mass
   PRIME: "yes": consider RS' mass else consider RS mass
*)

mMS2mRS[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, "no", 1];

mMS2mRS[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, "no", 1, "no"];

mMS2mRS[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, fdelm_/;NumericQ[fdelm]==True] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, "no", fdelm, "no"];

mMS2mRS[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, PRIME_String, fdelm_/;NumericQ[fdelm]==True] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, PRIME, fdelm, "no"];

mMS2mRSp[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, "yes", 1];

mMS2mRSp[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, "yes", 1, "no"];

mMS2mRSp[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, fdelm_/;NumericQ[fdelm]==True] := 
    mMS2mRS[mMS, mq, asmu, Mu, nuf, nl, loops, "yes", fdelm, "no"];

mMS2mRS[mMS_, mq_, asmu_, Mu_, nuf_, nl_, loops_, PRIME_String, fdelm_, marker_String] := 
    Module[{extmp, api, delmuf, exmOS, exmOSprime, apiprime, exals}, 
     If[ marker == "no", mark = 1, mark = "RunDecXMS2OS" ];
     If[loops < 0 || loops > 4, Print["mMS2mRS: PROCEDURE IS NOT IMPLEMENTED FOR ", 
         loops, " LOOPS."]; Return[]]; If[loops == 0, extmp = mMS, 
       If[PRIME === "yes", delmuf = exOS2RSp[api, Mu, nuf, nl]; , 
         delmuf = exOS2RS[api, Mu, nuf, nl]; ]; 
        exmOS = mMS2mOS[mMS, {}, api*Pi, Mu, nl + 1, loops, fdelm]; 
	exals = (DecAsUpMS[apiprime*Pi, mMS, Mu, nl, 4])/(apiprime*Pi); 
        exmOSprime = Normal[Series[exmOS /. {api -> apiprime*exals}, 
           {apiprime, 0, loops}]]; extmp = 
         Expand[mark * (exmOSprime /. apiprime -> api) - delmuf] /. 
          cut[api, loops]; ]; Return[N[extmp /. api -> asmu/Pi, numprec]]; ];


(* ************************************************************ *)
(*
mRS2mMS[mRS_, mq_, asnlmu_, Mu_, nuf_, nl_, loops_, PRIME_, fdelm_]
   computes the MS mass form the RS mass
   PRIME: "yes": consider RS' mass else consider RS mass
*)

mRS2mMS[mRS_, mq_, asnlmu_, Mu_, nuf_, nl_, loops_] := 
    mRS2mMS[mRS, mq, asnlmu, Mu, nuf, nl, loops, "no", 1];

mRS2mMS[mRS_, mq_, asnlmu_, Mu_, nuf_, nl_, loops_, fdelm_/;NumericQ[fdelm]==True] := 
    mRS2mMS[mRS, mq, asnlmu, Mu, nuf, nl, loops, "no", fdelm];

mRSp2mMS[mRS_, mq_, asnlmu_, Mu_, nuf_, nl_, loops_] := 
    mRS2mMS[mRS, mq, asnlmu, Mu, nuf, nl, loops, "yes", 1];

mRSp2mMS[mRS_, mq_, asnlmu_, Mu_, nuf_, nl_, loops_, fdelm_/;NumericQ[fdelm]==True] := 
    mRS2mMS[mRS, mq, asnlmu, Mu, nuf, nl, loops, "yes", fdelm];
 
mRS2mMS[mRS_, mq_, asnlmu_, Mu_, nuf_, nl_, loops_, PRIME_, fdelm_] :=  Module[
    {mtmp, mMS}, 
    mtmp = mMS2mRS[mMS, mq, asnlmu, Mu, nuf, nl, loops, PRIME, fdelm]; 
    Return[mMS /. FindRoot[mRS == mtmp, {mMS, mRS}]]; ];


(* ************************************************************ *)
(*
mRS2mSI[mRS_, mq_, asfct_, nuf_, nl_, loops_, PRIME_, fdelm_]
   computes the SI mass form the RS mass
   PRIME: "yes": consider RS' mass else consider RS mass
*)

mRS2mSI[mRS_, mq_, asfct_, nuf_, nl_, loops_] := 
    mRS2mSI[mRS, mq, asfct, nuf, nl, loops, "no", 1];

mRS2mSI[mRS_, mq_, asfct_, nuf_, nl_, loops_ , fdelm_/;NumericQ[fdelm]==True] := 
    mRS2mSI[mRS, mq, asfct, nuf, nl, loops, "no", fdelm];

mRSp2mSI[mRS_, mq_, asfct_, nuf_, nl_, loops_] := 
    mRS2mSI[mRS, mq, asfct, nuf, nl, loops, "yes", 1];

mRSp2mSI[mRS_, mq_, asfct_, nuf_, nl_, loops_ , fdelm_/;NumericQ[fdelm]==True] := 
    mRS2mSI[mRS, mq, asfct, nuf, nl, loops, "yes", fdelm];
 
mRS2mSI[mRS_, mq_, asfct_, nuf_, nl_, loops_, PRIME_, fdelm_] := 
    Module[{mMS1, mMS, acc}, 
	   mMS1 = 0; mMS = mRS; 
	   acc = 1.*^-6; 
	   While[Abs[mMS1 - mMS] > acc, mMS1 = mMS; 
		 mMS = mRS2mMS[mRS, mq, asfct[mMS1], 
			       mMS1, nuf, nl, loops, PRIME, fdelm]; ]; 
	   Return[mMS]; ];


(* ************************************************************ *)
(*
 mOS2m1S[mOS_, mq_, asmu_, mu_, nl_, loops_] computes the 1S mass from the OS mass
*)

mOS2m1S[mOS_, mq_, asmu_, mu_, nl_, loops_] := 
    Module[{xx, EnC, E1p, nn, eps}, If[loops < 0 || loops > 4, 
       Print["mOS2m1S: PROCEDURE IS NOT IMPLEMENTED FOR ", loops, " LOOPS."]; 
        Return[]]; E1p = 
       EnC*(1 + (asmu*xx*(97/6 + nl*(-11/9 - (2*Log[(3*mu)/(4*asmu*mOS)])/
                 3) + 11*Log[(3*mu)/(4*asmu*mOS)]))/Pi + 
           (asmu^2*xx^2*(1793/12 + (2917/216 - (11*nl)/18 + nl^2/54)*Pi^2 - 
              (9*Pi^4)/32 + (927*Log[(3*mu)/(4*asmu*mOS)])/4 + 
              (363*Log[(3*mu)/(4*asmu*mOS)]^2)/4 + nl*(-1693/72 - 
                (193*Log[(3*mu)/(4*asmu*mOS)])/6 - 11*Log[(3*mu)/(4*asmu*
                     mOS)]^2 - (19*Zeta[3])/2) + nl^2*(77/108 + 
                Log[(3*mu)/(4*asmu*mOS)] + Log[(3*mu)/(4*asmu*mOS)]^2/3 + 
                (2*Zeta[3])/9) + (275*Zeta[3])/4))/Pi^2 + 
           (asmu^3*xx^3*(1267919/1728 + Pi^4*(-723119/51840 + (11*nl^2)/
                 1080 - nl^3/4860 + nl*(59677/77760 + (3*Log[(3*mu)/(4*asmu*
                       mOS)])/8) - (99*Log[(3*mu)/(4*asmu*mOS)])/16) + 
              (4521*Log[(3*mu)/(4*asmu*mOS)]^2)/2 + (1331*
                Log[(3*mu)/(4*asmu*mOS)]^3)/2 + (114917*Zeta[3])/48 + 
              Log[(3*mu)/(4*asmu*mOS)]*(247675/96 + (3025*Zeta[3])/2) + 
              Pi^2*(265.389067842508 + (865*Log[mu/mOS])/18 + 
                (26897*Log[(3*mu)/(4*asmu*mOS)])/108 + nl^2*(905/432 + 
                  (11*Log[(3*mu)/(4*asmu*mOS)])/9 - (11*Zeta[3])/6) + 
                nl^3*(-19/486 - (2*Log[(3*mu)/(4*asmu*mOS)])/81 + 
                  Zeta[3]/27) + nl*(-397591/7776 - (5095*Log[(3*mu)/(4*asmu*
                       mOS)])/162 + (121*Zeta[3])/4)) + (13432.614375 - 
                3289.906669391583*nl - (1000*nl^3)/729 + nl^2*
                 ((14002/81 - (416*Zeta[3])/3)/3 + (3*(12541/243 + (64*Pi^4)/
                      135 + (368*Zeta[3])/3))/4))/32 + nl*(-52033/288 - 
                (10955*Log[(3*mu)/(4*asmu*mOS)]^2)/24 - 
                121*Log[(3*mu)/(4*asmu*mOS)]^3 + Log[(3*mu)/(4*asmu*mOS)]*
                 (-166309/288 - (902*Zeta[3])/3) - (8797*Zeta[3])/18 - 
                363*Zeta[5]) + nl^3*(-98/729 - (5*Log[(3*mu)/(4*asmu*mOS)]^2)/
                 9 - (4*Log[(3*mu)/(4*asmu*mOS)]^3)/27 + 
                Log[(3*mu)/(4*asmu*mOS)]*(-50/81 - (8*Zeta[3])/27) - 
                (44*Zeta[3])/81 - (4*Zeta[5])/9) + (3993*Zeta[5])/2 + 
              nl^2*(3073/288 + (1027*Log[(3*mu)/(4*asmu*mOS)]^2)/36 + 
                (22*Log[(3*mu)/(4*asmu*mOS)]^3)/3 + (3239*Zeta[3])/108 + 
                Log[(3*mu)/(4*asmu*mOS)]*(10351/288 + (158*Zeta[3])/9) + 
                22*Zeta[5])))/Pi^3) /. {xx -> 1, 
          EnC -> -((asmu^2*cf^2*mOS)/(4*1^2))} /. {cf -> 4/3}; 
      E1p = Expand[(1*(E1p /. asmu -> eps*asmu))/eps] /. 
        {eps^(nn_.) /; nn > loops :> 0}; Return[mOS + (1*E1p)/2 /. 
        eps -> 1]; ];
 
(* ************************************************************ *)
(*
 mMS2m1S[mMS_, mq_, asmu_, mu_, nl_, loops_, fdelm_] 
 computes the 1S mass from the MS mass
*)

mMS2m1S[mMS_, mq_, asmu_, mu_, nl_, loops_] := 
    mMS2m1S[mMS, mq, asmu, mu, nl, loops, 1, "no"];

mMS2m1S[mMS_, mq_, asmu_, mu_, nl_, loops_, marker_String] := 
    mMS2m1S[mMS, mq, asmu, mu, nl, loops, 1, marker];

mMS2m1S[mMS_, mq_, asmu_, mu_, nl_, loops_, fdelm_/;NumericQ[fdelm]] := 
    mMS2m1S[mMS, mq, asmu, mu, nl, loops, fdelm, "no"];
 
mMS2m1S[mMS_, mq_, asmu_, mu_, nl_, loops_, fdelm_, marker_String] := 
    Module[{mOS, m1S, res, epstmp, asmutmp, asmutmp2, keeplab, mark, markeps}, 
	   keeplab = "no";
     If[ marker == "no", mark = 1; markeps = 1, mark = "RunDecXMS2OS"; markeps = "RunDecXeps" ];
     If[loops == 0, Return[mMS]]; mOS = mMS2mOS[mMS, {}, asmutmp, mu, nl + 1, 
        loops, fdelm]; mOS = mOS /. {asmutmp -> DecAsUpMS[asmutmp2, mMS, mu, nl, 
           loops]}; mOS = Normal[Series[mOS, {asmutmp2, 0, loops}]] /. 
        {asmutmp2 -> asmu*eps*xmMS2mOS}; 
      m1S = Normal[Series[mOS2m1S[mOS, {}, asmu*epstmp, mu, nl, loops], 
         {epstmp, 0, loops + 3}]]; 
      m1S = Normal[Series[m1S /. {Log[epstmp] -> 0} /. 
          {epstmp^(nn_) :> eps^(nn - 1)*xmOS2m1S}, {eps, 0, loops}]]; 
      If[keeplab === "yes", 
       m1S = Expand[PowerExpand[m1S] /. Log[eps] -> 0 /. Log[xmOS2m1S] -> 
                0 /. Log[xmMS2mOS] -> 0 /. Log[xx_] :> Log[xx /. 
                {eps -> 1, xmMS2mOS -> 1}]] //. {xmOS2m1S^(nn_) :> xmOS2m1S, 
            xmMS2mOS^(nn_) :> xmMS2mOS} /. {xmMS2mOS*xmOS2m1S -> xmOS2m1S}; 
        res = Collect[{mOS, m1S}, {eps}, Expand]; Print["mOS = ", 
         InputForm[res[[1]]]]; m1S = res[[2]]; Print["m1S = ", 
         InputForm[m1S]]; 
         , 
         m1S = m1S /. {eps -> markeps, xmOS2m1S -> 1, xmMS2mOS -> mark}; ]; 
         Return[m1S]; ];


(* ************************************************************ *)
(*
 m1S2mMS[m1S_, mq_, asnlmu_, Mu_, nl_, loops_, fdelm_]
 computes the MS mass from the 1S mass
*)

m1S2mMS[m1S_, mq_, asnlmu_, Mu_, nl_, loops_] := m1S2mMS[m1S, mq, asnlmu, Mu, 
     nl, loops, 1];
 
m1S2mMS[m1S_, mq_, asnlmu_, Mu_, nl_, loops_, fdelm_] := 
    Module[{mtmp, mMS}, mtmp = mMS2m1S[mMS, mq, asnlmu, Mu, nl, loops, fdelm];
      Return[mMS /. FindRoot[m1S == mtmp, {mMS, m1S}]]; ];


(* ************************************************************ *)
(*
 m1S2mSI[m1S_, mq_, asfct_, nl_, loops_,fdelm_]
 computes the SI mass from the 1S mass
*)

m1S2mSI[m1S_, mq_, asfct_, nl_, loops_] := m1S2mSI[m1S, mq, asfct, nl, loops, 1];
 
m1S2mSI[m1S_, mq_, asfct_, nl_, loops_, fdelm_] := 
    Module[{mMS1, mMS, acc}, 
	   mMS1 = 0; mMS = m1S; 
	   acc = 1.*^-6; 
	   While[Abs[mMS1 - mMS] > acc, 
		 mMS1 = mMS; mMS = m1S2mMS[m1S, mq, asfct[mMS1], mMS1, nl, loops, 
					   fdelm]; ]; 
	   Return[mMS]; ];

(* ************************************************************ *)
(*
 mKIN2mMS[zmbKIN_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, znlMSOS_, znlOSKIN_, zloops_, CASE_, KINverbose_]
 computes the MSbar mass from the kinetic mass
*)

tagto1 = Map[ Rule[#,1]&, {RunDec`Modules`TAGdecmcMSOS, RunDec`Modules`TAGnc, RunDec`Modules`TAGxmc, RunDec`Modules`TAGkinmc}];

mKIN2mMS[zmbKIN_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, zloops_, "A"] := 
 mKIN2mMS[zmbKIN, {zmcMS,zmusmc}, zasmus, zmus, zmuf, 3, 3, zloops, "deccharm"] /. tagto1;
mKIN2mMS[zmbKIN_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, zloops_, "B"] := 
 mKIN2mMS[zmbKIN, {zmcMS,zmusmc}, zasmus, zmus, zmuf, 3, 3, zloops, ""]  /. tagto1;
mKIN2mMS[zmbKIN_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, zloops_, "C"] := 
 mKIN2mMS[zmbKIN, {zmcMS,zmusmc}, zasmus, zmus, zmuf, 3, 4, zloops, ""] /. TAGkinmc->0  /. tagto1;
mKIN2mMS[zmbKIN_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, zloops_, "D"] := 
 mKIN2mMS[zmbKIN, {zmcMS,zmusmc}, zasmus, zmus, zmuf, 3, 3, zloops, "deccharm","yes"] /. {TAGkinmc->0,TAGxmc->0,TAGdecmcMSOS->0,TAGnc->0};

 mKIN2mMS[zmbKIN_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, znlMSOS_, znlOSKIN_, zloops_, CASE_] := 
 mKIN2mMS[zmbKIN, {zmcMS,zmusmc}, zasmus, zmus, zmuf, znlMSOS, znlOSKIN, zloops, CASE, "yes"];

 mKIN2mMS[zmbKIN_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, znlMSOS_, znlOSKIN_, zloops_, CASE_, KINverbose_] := Module[
    {kk,xmc,numlist,numTAG,xmcMS,xmusmc},
    If[ (CASE =!= "deccharm") && (CASE =!= ""), Print["mKIN2mMS: illigal option CASE = ",CASE]; Abort[] ];
    {xmcMS,xmusmc} = {zmcMS,zmusmc};
    If[zmcMS==0,xmcMS=10^-20]; (* note: set mc=0; massless part in *)
    If[zmusmc==0,xmusmc=1];    (* nm term is needed since NLxMSOS->NL-1 is chosen *)
    (numlist = {mkin      -> zmbKIN,
		mus       -> zmus,
		mcMSmusmc -> xmcMS * xmc,
		musmc     -> xmusmc,
		mufac     -> zmuf,
		NLxMSOS   -> znlMSOS,
                NLxOSKIN  -> znlOSKIN});
    numTAG = {DUMMY -> DUMMY};
    If[KINverbose =!= "yes", numTAG = {TAGdecmcMSOS->1, TAGxmc->1, TAGkinmc->1, TAGnc->1} ];
    If[zmusmc==0,numTAG = {TAGdecmcMSOS->0, TAGxmc->0, TAGkinmc->0, TAGnc->0}];
    If[CASE === "deccharm",
       Return[ Expand[PowerExpand[RunDecmkin2mMSnl/.numlist]] 
                   /. {(apinlmus^kk_/;kk>zloops) :> 0}
		   //. {Log[xmc]->0,xmc^nn_:>xmc} /. xmc->TAGxmc
		   /.numTAG
                   //. {TAGnc^kk_ -> TAGnc}
		   /. {apinlmus -> zasmus/Pi} 
                   ];
    ,
       Return[ Expand[PowerExpand[RunDecmkin2mMSmc/.numlist]] 
                   /. {(apinl1mus^kk_/;kk>zloops) :> 0}
		   //. {Log[xmc]->0,xmc^nn_:>xmc} /. xmc->TAGxmc
		   /.numTAG
                   //. {TAGkinmc^kk_ -> TAGkinmc}
                   //. {TAGnc^kk_ -> TAGnc}
		   /. {apinl1mus -> zasmus/Pi} ];
   ];
];

(* ************************************************************ *)
(*
 mMS2mKIN[zmbMS_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, znlMSOS_, znlOSKIN_, zloops_, CASE_, KINverbose_]
 computes the kinetic mass from the MSbar mass
*)
mMS2mKIN[zmbMS_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, zloops_, "A"] := 
 mMS2mKIN[zmbMS, {zmcMS,zmusmc}, zasmus, zmus, zmuf, 3, 3, zloops, "deccharm"] /. tagto1;
mMS2mKIN[zmbMS_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, zloops_, "B"] := 
 mMS2mKIN[zmbMS, {zmcMS,zmusmc}, zasmus, zmus, zmuf, 3, 3, zloops, ""] /. tagto1;
mMS2mKIN[zmbMS_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, zloops_, "C"] := 
 mMS2mKIN[zmbMS, {zmcMS,zmusmc}, zasmus, zmus, zmuf, 3(**LABnlMSOS*), 4(**LABnlOSKIN*), zloops, ""] /. TAGkinmc->0 /. tagto1;
mMS2mKIN[zmbMS_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, zloops_, "D"] := 
 mMS2mKIN[zmbMS, {zmcMS,zmusmc}, zasmus, zmus, zmuf, 3, 3, zloops, "deccharm","yes"]  /. {TAGkinmc->0,TAGxmc->0,TAGdecmcMSOS->0,TAGnc->0};
 
 mMS2mKIN[zmbMS_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, znlMSOS_, znlOSKIN_, zloops_, CASE_] := 
 mMS2mKIN[zmbMS, {zmcMS,zmusmc}, zasmus, zmus, zmuf, znlMSOS, znlOSKIN, zloops, CASE, "yes"];

 mMS2mKIN[zmbMS_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, znlMSOS_, znlOSKIN_, zloops_, CASE_, KINverbose_] := Module[
    (*{numlist,kk, mbMSmus,mus,TAGdecmcMSOS,apinlmus,apinl1mus,mufac,NL},*)
    {kk,xmc,numlist,numTAG,xmcMS,xmusmc},
    If[ (CASE =!= "deccharm") && (CASE =!= ""), Print["mMS2mKIN: illigal option CASE = ",CASE]; Abort[] ];
    {xmcMS,xmusmc} = {zmcMS,zmusmc};
    If[zmcMS==0,xmcMS=10^-20]; (* note: set mc=0; massless part in *)
    If[zmusmc==0,xmusmc=1];    (* nm term is needed since NLxMSOS->NL-1 is chosen *)
    (numlist = {mbMSmus   -> zmbMS,
		mus       -> zmus,
		mcMSmusmc -> xmcMS * xmc,
		musmc     -> xmusmc,
		mufac     -> zmuf,
		NLxMSOS   -> znlMSOS,
                NLxOSKIN  -> znlOSKIN});
    numTAG = {DUMMY -> DUMMY};
    If[KINverbose =!= "yes", numTAG = {TAGdecmcMSOS->1, TAGxmc->1, TAGkinmc->1, TAGnc->1} ];
    If[zmusmc==0,numTAG = {TAGdecmcMSOS->0, TAGxmc->0, TAGkinmc->0, TAGnc->0}];
    If[CASE === "deccharm",
       Return[ Expand[PowerExpand[RunDecmMS2mkinnl/.numlist]] 
                   /. {(apinlmus^kk_/;kk>zloops) :> 0}
		   //. {Log[xmc]->0,xmc^nn_:>xmc} /. xmc->TAGxmc
		   /.numTAG
                   //. {TAGnc^kk_ -> TAGnc}
		   /. {apinlmus -> zasmus/Pi} 
                   ];
    ,
       Return[ Expand[PowerExpand[RunDecmMS2mkinmc/.numlist]] 
                   /. {(apinl1mus^kk_/;kk>zloops) :> 0}
		   //. {Log[xmc]->0,xmc^nn_:>xmc} /. xmc->TAGxmc
		   /.numTAG
                   //. {TAGkinmc^kk_ -> TAGkinmc}
                   //. {TAGnc^kk_ -> TAGnc}
		   /. {apinl1mus -> zasmus/Pi} ];
   ];
];

(* ************************************************************ *)

End[];

(* ************************************************************ *)

EndPackage[];

(* ************************************************************ *)
(* ************************************************************ *)

