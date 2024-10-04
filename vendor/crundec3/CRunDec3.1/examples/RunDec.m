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
    "mKIN2mMS[mbKIN, {mcMS,musmc}, asmus, mus, muf, nl, nloops, CASE]  TODO UPDATE
computes mMS(mus) from mKIN(muf). The routines are developed for bottom and charm quarks.
However, they can also be applied to other qurk flavours.
A theory with nf active flavours is assumed.
CASE defines the scheme which has to be adapted for the lighter quark (charm) mass.
CASE = \"A\" and \"D\": nl=nf-2, asmus = \\alpha_s^{(nl)}(mus)
CASE = \"B\" and \"C\": nl=nf-1, asmus = \\alpha_s^{(nl)}(mus)
nloops=1,2 or 3 is the number of loops.
{mcMS,musmc} specifies the lighter quark mass in the MSbar scheme: mcMS(musmc).
{zmcMS,zmusmc} = {0,0}: mass relations fro the light quark are computed; CASE is ignored.";

 mMS2mKIN::usage =
    "mMS2mKIN[mbMS, {mcMS,musmc}, asmus, mus, muf, nl, nloops, CASE]
computes mKIN(muf) from mMS(mus). The routines are developed for bottom and charm quarks.
However, they can also be applied to other qurk flavours.
A theory with nf active flavours is assumed.
CASE defines the scheme which has to be adapted for the lighter quark (charm) mass.
CASE = \"A\" and \"D\": nl=nf-2, asmus = \\alpha_s^{(nl)}(mus)
CASE = \"B\" and \"C\": nl=nf-1, asmus = \\alpha_s^{(nl)}(mus)
nloops=1,2 or 3 is the number of loops.
{mcMS,musmc} specifies the lighter quark mass in the MSbar scheme: mcMS(musmc).
{zmcMS,zmusmc} = {0,0}: mass relations fro the light quark are computed; CASE is ignored."

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

RunDecmMS2mkinnl = 1.*mbMSmus + apinlmus*(-1.7777777777777777*mufac - 
       (0.6666666666666666*mufac^2)/mbMSmus + 
       mbMSmus*(4/3 + Log[mus^2/mbMSmus^2])) + 
     apinlmus^2*(1.596*mcMSmusmc + (0.17773333333333333*mcMSmusmc^3)/
        mbMSmus^2 + mufac*(-23.078870161994644 + 
         NL*(1.5802469135802468 - 0.5925925925925926*Log[(2*mufac)/mus]) + 
         9.777777777777779*Log[(2*mufac)/mus]) + 
       (-0.6285333333333333*mcMSmusmc^2 + mufac^2*(-5.932354088525768 + 
           NL*(0.48148148148148145 - 0.2222222222222222*Log[(2*mufac)/mus]) + 
           3.6666666666666665*Log[(2*mufac)/mus] + 0.6666666666666666*
            Log[mus^2/mbMSmus^2]))/mbMSmus + 
       mbMSmus*(13.443429834770056 + 1.9583333333333333*Log[mus^2/mbMSmus^2]^
           2 + NL*(-1.0413669111716308 - 0.36111111111111105*
            Log[mus^2/mbMSmus^2] - 0.08333333333333333*Log[mus^2/mbMSmus^2]^
             2) + (2*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/9 + 
         Log[mus^2/mbMSmus^2]*(7.069444444444444 + 
           (TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])/6))) + 
     apinlmus^3*(23.91811111111111*mcMSmusmc + 
       (0.11848888888888888*mcMSmusmc^3*mufac^2)/mbMSmus^4 + 
       (0.04933333333333332*mcMSmusmc^4 - 0.4190222222222222*mcMSmusmc^2*
          mufac^2)/mbMSmus^3 - 18.049*mcMSmusmc*Log[mcMSmusmc/mbMSmus] + 
       mufac*(-362.1430910903762 + 276.53423844860777*Log[(2*mufac)/mus] - 
         53.77777777777778*Log[(2*mufac)/mus]^2 + 
         NL^2*(-1.4473655615553573 + 1.0534979423868311*Log[(2*mufac)/mus] - 
           0.19753086419753085*Log[(2*mufac)/mus]^2) + 
         NL*(50.92088541456441 - 35.583444305527294*Log[(2*mufac)/mus] + 
           6.518518518518517*Log[(2*mufac)/mus]^2)) + 
       9.273222222222222*mcMSmusmc*Log[mus^2/mbMSmus^2] + 
       NL*(-0.9826666666666665*mcMSmusmc + 1.101333333333333*mcMSmusmc*
          Log[mcMSmusmc/mbMSmus] - 0.5346666666666666*mcMSmusmc*
          Log[mus^2/mbMSmus^2]) + 0.532*mcMSmusmc*TAGdecmcMSOS*
        Log[mus^2/mcMSmusmc^2] + 1.596*mcMSmusmc*Log[musmc^2/mcMSmusmc^2] + 
       (4.229311111111111*mcMSmusmc^2 - 2.7875555555555556*mcMSmusmc^2*
          Log[mcMSmusmc/mbMSmus] - 2.8946222222222215*mcMSmusmc^2*
          Log[mus^2/mbMSmus^2] + NL*(-0.3003333333333333*mcMSmusmc^2 + 
           0.21999999999999997*mcMSmusmc^2*Log[mus^2/mbMSmus^2]) - 
         0.2095111111111111*mcMSmusmc^2*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2] + 
         mufac^2*(-56.99806410856273 - 20.166666666666664*Log[(2*mufac)/mus]^
             2 + NL^2*(-0.30819418434869106 + 0.32098765432098764*
              Log[(2*mufac)/mus] - 0.07407407407407407*Log[(2*mufac)/mus]^
               2) + Log[(2*mufac)/mus]*(78.64478386267234 - 
             3.6666666666666665*Log[mus^2/mbMSmus^2]) + 0.6388888888888888*
            Log[mus^2/mbMSmus^2]^2 + NL*(9.864994084863337 + 
             2.4444444444444446*Log[(2*mufac)/mus]^2 + Log[(2*mufac)/mus]*
              (-10.603050873831995 + 0.2222222222222222*Log[mus^2/
                  mbMSmus^2]) - 0.7222222222222222*Log[mus^2/mbMSmus^2] - 
             0.05555555555555555*Log[mus^2/mbMSmus^2]^2) + 
           0.14814814814814814*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2] + 
           Log[mus^2/mbMSmus^2]*(9.756428162599843 + 0.1111111111111111*
              TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2])) - 1.2570666666666666*
          mcMSmusmc^2*Log[musmc^2/mcMSmusmc^2])/mbMSmus + 
       (1.1406*mcMSmusmc^3 + 1.064*mcMSmusmc*mufac^2 + 0.03433333333333333*
          mcMSmusmc^3*Log[mcMSmusmc/mbMSmus] + 0.6728444444444446*mcMSmusmc^3*
          Log[mus^2/mbMSmus^2] - 0.06699999999999999*mcMSmusmc^3*NL*
          Log[mus^2/mbMSmus^2] + 0.05924444444444444*mcMSmusmc^3*TAGdecmcMSOS*
          Log[mus^2/mcMSmusmc^2] + 0.5332*mcMSmusmc^3*
          Log[musmc^2/mcMSmusmc^2])/mbMSmus^2 + 
       mbMSmus*(190.39103694423557 - (11*TAGdecmcMSOS)/54 + 
         4.30787037037037*Log[mus^2/mbMSmus^2]^3 + 
         NL^2*(0.6526907490815436 + 0.3201161308843708*Log[mus^2/mbMSmus^2] + 
           0.060185185185185175*Log[mus^2/mbMSmus^2]^2 + 0.009259259259259259*
            Log[mus^2/mbMSmus^2]^3) + 5.536698833812239*TAGdecmcMSOS*
          Log[mus^2/mcMSmusmc^2] + (TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]^2)/
          27 + Log[mus^2/mbMSmus^2]^2*(25.13310185185185 + 
           0.6527777777777777*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]) + 
         NL*(-26.655179856935014 - 0.39814814814814803*Log[mus^2/mbMSmus^2]^
             3 - 0.3471223037238769*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2] + 
           Log[mus^2/mbMSmus^2]*(-13.828790516488414 - 0.12037037037037035*
              TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]) + Log[mus^2/mbMSmus^2]^2*
            (-2.673611111111111 - 0.027777777777777776*TAGdecmcMSOS*
              Log[mus^2/mcMSmusmc^2])) - (4*TAGdecmcMSOS*
           Log[musmc^2/mcMSmusmc^2])/9 + Log[mus^2/mbMSmus^2]*
          (96.08143257463242 - (11*TAGdecmcMSOS)/72 + 3.1481481481481475*
            TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2] + 
           (TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]^2)/36 - 
           (TAGdecmcMSOS*Log[musmc^2/mcMSmusmc^2])/3)))
 
RunDecmMS2mkinmc = 1.596*apinl1mus^2*mcMSmusmc + 23.91811111111111*
      apinl1mus^3*mcMSmusmc + (0.11848888888888888*apinl1mus^3*mcMSmusmc^3*
       mufac^2)/mbMSmus^4 + (0.04933333333333332*apinl1mus^3*mcMSmusmc^4 - 
       0.4190222222222222*apinl1mus^3*mcMSmusmc^2*mufac^2)/mbMSmus^3 - 
     18.049*apinl1mus^3*mcMSmusmc*Log[mcMSmusmc/mbMSmus] + 
     9.273222222222222*apinl1mus^3*mcMSmusmc*Log[mus^2/mbMSmus^2] + 
     NL*(-0.9826666666666665*apinl1mus^3*mcMSmusmc + 1.101333333333333*
        apinl1mus^3*mcMSmusmc*Log[mcMSmusmc/mbMSmus] - 
       0.5346666666666666*apinl1mus^3*mcMSmusmc*Log[mus^2/mbMSmus^2]) + 
     mufac*(-1.7777777777777777*apinl1mus - 23.078870161994644*apinl1mus^2 - 
       362.1430910903762*apinl1mus^3 - 0.2716049382716049*apinl1mus^3*
        TAGkinmc - 53.77777777777778*apinl1mus^3*Log[(2*mufac)/mus]^2 + 
       NL^2*(-1.4473655615553573*apinl1mus^3 + 1.0534979423868311*apinl1mus^3*
          Log[(2*mufac)/mus] - 0.19753086419753085*apinl1mus^3*
          Log[(2*mufac)/mus]^2) + (0.2962962962962963*apinl1mus^2*TAGkinmc + 
         8.507771535479696*apinl1mus^3*TAGkinmc)*Log[mus^2/mcMSmusmc^2] - 
       0.04938271604938271*apinl1mus^3*TAGkinmc*Log[mus^2/mcMSmusmc^2]^2 + 
       Log[(2*mufac)/mus]*(9.777777777777779*apinl1mus^2 + 
         276.53423844860777*apinl1mus^3 - 3.2592592592592595*apinl1mus^3*
          TAGkinmc*Log[mus^2/mcMSmusmc^2]) + 
       NL*(1.5802469135802468*apinl1mus^2 + 50.92088541456441*apinl1mus^3 + 
         6.518518518518517*apinl1mus^3*Log[(2*mufac)/mus]^2 - 
         0.5267489711934156*apinl1mus^3*TAGkinmc*Log[mus^2/mcMSmusmc^2] + 
         Log[(2*mufac)/mus]*(-0.5925925925925926*apinl1mus^2 - 
           35.583444305527294*apinl1mus^3 + 0.19753086419753085*apinl1mus^3*
            TAGkinmc*Log[mus^2/mcMSmusmc^2]))) + 1.596*apinl1mus^3*mcMSmusmc*
      Log[musmc^2/mcMSmusmc^2] + (-0.6285333333333333*apinl1mus^2*
        mcMSmusmc^2 + 4.229311111111111*apinl1mus^3*mcMSmusmc^2 - 
       2.7875555555555556*apinl1mus^3*mcMSmusmc^2*Log[mcMSmusmc/mbMSmus] - 
       2.8946222222222215*apinl1mus^3*mcMSmusmc^2*Log[mus^2/mbMSmus^2] + 
       NL*(-0.3003333333333333*apinl1mus^3*mcMSmusmc^2 + 
         0.21999999999999997*apinl1mus^3*mcMSmusmc^2*Log[mus^2/mbMSmus^2]) + 
       mufac^2*(-0.6666666666666666*apinl1mus - 5.932354088525768*
          apinl1mus^2 - 56.99806410856273*apinl1mus^3 - 0.10185185185185186*
          apinl1mus^3*TAGkinmc - 20.166666666666664*apinl1mus^3*
          Log[(2*mufac)/mus]^2 + NL^2*(-0.30819418434869106*apinl1mus^3 + 
           0.32098765432098764*apinl1mus^3*Log[(2*mufac)/mus] - 
           0.07407407407407407*apinl1mus^3*Log[(2*mufac)/mus]^2) + 
         0.6388888888888888*apinl1mus^3*Log[mus^2/mbMSmus^2]^2 + 
         (0.1111111111111111*apinl1mus^2*TAGkinmc + 2.4311550665456267*
            apinl1mus^3*TAGkinmc)*Log[mus^2/mcMSmusmc^2] - 
         0.018518518518518517*apinl1mus^3*TAGkinmc*Log[mus^2/mcMSmusmc^2]^2 + 
         Log[(2*mufac)/mus]*(3.6666666666666665*apinl1mus^2 + 
           78.64478386267234*apinl1mus^3 - 3.6666666666666665*apinl1mus^3*
            Log[mus^2/mbMSmus^2] - 1.2222222222222219*apinl1mus^3*TAGkinmc*
            Log[mus^2/mcMSmusmc^2]) + Log[mus^2/mbMSmus^2]*
          (0.6666666666666666*apinl1mus^2 + 9.756428162599843*apinl1mus^3 - 
           0.1111111111111111*apinl1mus^3*TAGkinmc*Log[mus^2/mcMSmusmc^2]) + 
         NL*(0.48148148148148145*apinl1mus^2 + 9.864994084863337*
            apinl1mus^3 + 2.4444444444444446*apinl1mus^3*Log[(2*mufac)/mus]^
             2 - 0.7222222222222222*apinl1mus^3*Log[mus^2/mbMSmus^2] - 
           0.05555555555555555*apinl1mus^3*Log[mus^2/mbMSmus^2]^2 - 
           0.16049382716049382*apinl1mus^3*TAGkinmc*Log[mus^2/mcMSmusmc^2] + 
           Log[(2*mufac)/mus]*(-0.2222222222222222*apinl1mus^2 - 
             10.603050873831995*apinl1mus^3 + 0.2222222222222222*apinl1mus^3*
              Log[mus^2/mbMSmus^2] + 0.07407407407407407*apinl1mus^3*TAGkinmc*
              Log[mus^2/mcMSmusmc^2]))) - 1.2570666666666666*apinl1mus^3*
        mcMSmusmc^2*Log[musmc^2/mcMSmusmc^2])/mbMSmus + 
     (0.17773333333333333*apinl1mus^2*mcMSmusmc^3 + 1.1406*apinl1mus^3*
        mcMSmusmc^3 + 1.064*apinl1mus^3*mcMSmusmc*mufac^2 + 
       0.03433333333333333*apinl1mus^3*mcMSmusmc^3*Log[mcMSmusmc/mbMSmus] + 
       0.6728444444444446*apinl1mus^3*mcMSmusmc^3*Log[mus^2/mbMSmus^2] - 
       0.06699999999999999*apinl1mus^3*mcMSmusmc^3*NL*Log[mus^2/mbMSmus^2] + 
       0.5332*apinl1mus^3*mcMSmusmc^3*Log[musmc^2/mcMSmusmc^2])/mbMSmus^2 + 
     mbMSmus*(1. + (4*apinl1mus)/3 + 13.443429834770056*apinl1mus^2 + 
       190.39103694423557*apinl1mus^3 + (1.9583333333333333*apinl1mus^2 + 
         25.13310185185185*apinl1mus^3)*Log[mus^2/mbMSmus^2]^2 + 
       4.30787037037037*apinl1mus^3*Log[mus^2/mbMSmus^2]^3 + 
       NL*(-1.0413669111716308*apinl1mus^2 - 26.655179856935014*apinl1mus^3 + 
         (-0.36111111111111105*apinl1mus^2 - 13.828790516488414*apinl1mus^3)*
          Log[mus^2/mbMSmus^2] + (-0.08333333333333333*apinl1mus^2 - 
           2.673611111111111*apinl1mus^3)*Log[mus^2/mbMSmus^2]^2 - 
         0.39814814814814803*apinl1mus^3*Log[mus^2/mbMSmus^2]^3) + 
       NL^2*(0.6526907490815436*apinl1mus^3 + 0.3201161308843708*apinl1mus^3*
          Log[mus^2/mbMSmus^2] + 0.060185185185185175*apinl1mus^3*
          Log[mus^2/mbMSmus^2]^2 + 0.009259259259259259*apinl1mus^3*
          Log[mus^2/mbMSmus^2]^3) + 0.4444444444444444*apinl1mus^3*TAGkinmc*
        Log[mus^2/mcMSmusmc^2] + (0.07407407407407407*apinl1mus^3*TAGkinmc - 
         (2*apinl1mus^3*TAGkinmc^2)/27)*Log[mus^2/mcMSmusmc^2]^2 - 
       (4*apinl1mus^3*TAGkinmc*Log[musmc^2/mcMSmusmc^2])/9 + 
       Log[mus^2/mbMSmus^2]*(apinl1mus + 7.069444444444444*apinl1mus^2 + 
         96.08143257463242*apinl1mus^3 + 0.3333333333333333*apinl1mus^3*
          TAGkinmc*Log[mus^2/mcMSmusmc^2] + 
         (0.05555555555555555*apinl1mus^3*TAGkinmc - (apinl1mus^3*TAGkinmc^2)/
            18)*Log[mus^2/mcMSmusmc^2]^2 - (apinl1mus^3*TAGkinmc*
           Log[musmc^2/mcMSmusmc^2])/3))
 
RunDecmkin2mMSnl = 1.*mkin + apinlmus*(-1.3333333333333333*mkin + 
       1.7777777777777777*mufac + (0.6666666666666666*mufac^2)/mkin - 
       1.*mkin*Log[(1.*mus^2)/mkin^2]) + apinlmus^2*(0. - 1.596*mcMSmusmc - 
       (0.17773333333333333*mcMSmusmc^3)/mkin^2 + 
       (0.6285333333333333*mcMSmusmc^2)/mkin - 14.332318723658945*mkin + 
       24.264055347179827*mufac + (7.265687421859101*mufac^2)/mkin + 
       1.0413669111716308*mkin*NL - 1.5802469135802468*mufac*NL - 
       (0.48148148148148145*mufac^2*NL)/mkin - 9.777777777777779*mufac*
        Log[(2*mufac)/mus] - (3.6666666666666665*mufac^2*Log[(2*mufac)/mus])/
        mkin + 0.5925925925925926*mufac*NL*Log[(2*mufac)/mus] + 
       (0.2222222222222222*mufac^2*NL*Log[(2*mufac)/mus])/mkin - 
       0.2222222222222222*mkin*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2] - 
       6.402777777777777*mkin*Log[(1.*mus^2)/mkin^2] - 
       1.7777777777777777*mufac*Log[(1.*mus^2)/mkin^2] - 
       (0.6666666666666667*mufac^2*Log[(1.*mus^2)/mkin^2])/mkin + 
       0.36111111111111105*mkin*NL*Log[(1.*mus^2)/mkin^2] - 
       0.16666666666666666*mkin*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]*
        Log[(1.*mus^2)/mkin^2] - 0.9583333333333333*mkin*
        Log[(1.*mus^2)/mkin^2]^2 + 0.08333333333333333*mkin*NL*
        Log[(1.*mus^2)/mkin^2]^2) + 
     (apinlmus^3*(0. - 0.04933333333333332*mcMSmusmc^4*mkin^2 - 
        1.7330444444444446*mcMSmusmc^3*mkin^3 - 2.972244444444443*mcMSmusmc^2*
         mkin^4 - 24.98211111111111*mcMSmusmc*mkin^5 - 
        199.0954170543889*mkin^6 + 0.6319407407407407*mcMSmusmc^3*mkin^2*
         mufac - 1.1173925925925925*mcMSmusmc^2*mkin^3*mufac + 
        374.8147588995222*mkin^5*mufac + 0.2369777777777778*mcMSmusmc^3*mkin*
         mufac^2 - 0.4190222222222222*mcMSmusmc^2*mkin^2*mufac^2 + 
        80.56030314981179*mkin^4*mufac^2 + 0.3003333333333333*mcMSmusmc^2*
         mkin^4*NL + 0.9826666666666665*mcMSmusmc*mkin^5*NL + 
        26.923898212450222*mkin^6*NL - 51.40701502104119*mkin^5*mufac*NL - 
        11.309438529307783*mkin^4*mufac^2*NL - 0.6526907490815436*mkin^6*
         NL^2 + 1.4473655615553573*mkin^5*mufac*NL^2 + 0.30819418434869106*
         mkin^4*mufac^2*NL^2 + 0.2037037037037037*mkin^6*TAGdecmcMSOS - 
        0.03433333333333333*mcMSmusmc^3*mkin^3*Log[(1.*mcMSmusmc)/mkin] + 
        2.7875555555555556*mcMSmusmc^2*mkin^4*Log[(1.*mcMSmusmc)/mkin] + 
        18.049*mcMSmusmc*mkin^5*Log[(1.*mcMSmusmc)/mkin] - 
        1.101333333333333*mcMSmusmc*mkin^5*NL*Log[(1.*mcMSmusmc)/mkin] - 
        283.05275696712624*mkin^5*mufac*Log[(2*mufac)/mus] - 
        85.97811719600566*mkin^4*mufac^2*Log[(2*mufac)/mus] + 
        35.97850603392235*mkin^5*mufac*NL*Log[(2*mufac)/mus] + 
        11.047495318276438*mkin^4*mufac^2*NL*Log[(2*mufac)/mus] - 
        1.0534979423868311*mkin^5*mufac*NL^2*Log[(2*mufac)/mus] - 
        0.32098765432098764*mkin^4*mufac^2*NL^2*Log[(2*mufac)/mus] + 
        53.77777777777778*mkin^5*mufac*Log[(2*mufac)/mus]^2 + 
        20.166666666666664*mkin^4*mufac^2*Log[(2*mufac)/mus]^2 - 
        6.518518518518517*mkin^5*mufac*NL*Log[(2*mufac)/mus]^2 - 
        2.4444444444444446*mkin^4*mufac^2*NL*Log[(2*mufac)/mus]^2 + 
        0.19753086419753085*mkin^5*mufac*NL^2*Log[(2*mufac)/mus]^2 + 
        0.07407407407407407*mkin^4*mufac^2*NL^2*Log[(2*mufac)/mus]^2 - 
        0.05924444444444444*mcMSmusmc^3*mkin^3*TAGdecmcMSOS*
         Log[mus^2/mcMSmusmc^2] + 0.2095111111111111*mcMSmusmc^2*mkin^4*
         TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2] - 0.532*mcMSmusmc*mkin^5*
         TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2] - 5.832995130108536*mkin^6*
         TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2] + 0.19753086419753085*mkin^5*
         mufac*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2] + 0.07407407407407407*
         mkin^4*mufac^2*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2] + 
        0.3471223037238769*mkin^6*NL*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2] - 
        0.037037037037037035*mkin^6*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]^2 - 
        0.850577777777778*mcMSmusmc^3*mkin^3*Log[(1.*mus^2)/mkin^2] + 
        2.894622222222221*mcMSmusmc^2*mkin^4*Log[(1.*mus^2)/mkin^2] - 
        7.6772222222222215*mcMSmusmc*mkin^5*Log[(1.*mus^2)/mkin^2] - 
        85.06494327546267*mkin^6*Log[(1.*mus^2)/mkin^2] - 
        27.64677139656254*mkin^5*mufac*Log[(1.*mus^2)/mkin^2] - 
        8.534205940377621*mkin^4*mufac^2*Log[(1.*mus^2)/mkin^2] + 
        0.06699999999999999*mcMSmusmc^3*mkin^3*NL*Log[(1.*mus^2)/mkin^2] - 
        0.21999999999999997*mcMSmusmc^2*mkin^4*NL*Log[(1.*mus^2)/mkin^2] + 
        0.5346666666666666*mcMSmusmc*mkin^5*NL*Log[(1.*mus^2)/mkin^2] + 
        12.67198262007108*mkin^6*NL*Log[(1.*mus^2)/mkin^2] + 
        1.6296296296296293*mkin^5*mufac*NL*Log[(1.*mus^2)/mkin^2] + 
        0.49999999999999956*mkin^4*mufac^2*NL*Log[(1.*mus^2)/mkin^2] - 
        0.3201161308843708*mkin^6*NL^2*Log[(1.*mus^2)/mkin^2] + 
        0.1527777777777778*mkin^6*TAGdecmcMSOS*Log[(1.*mus^2)/mkin^2] + 
        9.777777777777779*mkin^5*mufac*Log[(2*mufac)/mus]*
         Log[(1.*mus^2)/mkin^2] + 3.666666666666667*mkin^4*mufac^2*
         Log[(2*mufac)/mus]*Log[(1.*mus^2)/mkin^2] - 0.5925925925925926*
         mkin^5*mufac*NL*Log[(2*mufac)/mus]*Log[(1.*mus^2)/mkin^2] - 
        0.2222222222222222*mkin^4*mufac^2*NL*Log[(2*mufac)/mus]*
         Log[(1.*mus^2)/mkin^2] - 2.925925925925925*mkin^6*TAGdecmcMSOS*
         Log[mus^2/mcMSmusmc^2]*Log[(1.*mus^2)/mkin^2] - 
        0.2962962962962963*mkin^5*mufac*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]*
         Log[(1.*mus^2)/mkin^2] - 0.1111111111111111*mkin^4*mufac^2*
         TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]*Log[(1.*mus^2)/mkin^2] + 
        0.12037037037037035*mkin^6*NL*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]*
         Log[(1.*mus^2)/mkin^2] - 0.027777777777777776*mkin^6*TAGdecmcMSOS*
         Log[mus^2/mcMSmusmc^2]^2*Log[(1.*mus^2)/mkin^2] - 
        16.521990740740744*mkin^6*Log[(1.*mus^2)/mkin^2]^2 - 
        1.7037037037037035*mkin^5*mufac*Log[(1.*mus^2)/mkin^2]^2 - 
        0.6388888888888893*mkin^4*mufac^2*Log[(1.*mus^2)/mkin^2]^2 + 
        2.229166666666667*mkin^6*NL*Log[(1.*mus^2)/mkin^2]^2 + 
        0.14814814814814814*mkin^5*mufac*NL*Log[(1.*mus^2)/mkin^2]^2 + 
        0.05555555555555555*mkin^4*mufac^2*NL*Log[(1.*mus^2)/mkin^2]^2 - 
        0.060185185185185175*mkin^6*NL^2*Log[(1.*mus^2)/mkin^2]^2 - 
        0.3194444444444444*mkin^6*TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]*
         Log[(1.*mus^2)/mkin^2]^2 + 0.027777777777777776*mkin^6*NL*
         TAGdecmcMSOS*Log[mus^2/mcMSmusmc^2]*Log[(1.*mus^2)/mkin^2]^2 - 
        1.391203703703704*mkin^6*Log[(1.*mus^2)/mkin^2]^3 + 
        0.2314814814814814*mkin^6*NL*Log[(1.*mus^2)/mkin^2]^3 - 
        0.009259259259259259*mkin^6*NL^2*Log[(1.*mus^2)/mkin^2]^3 - 
        0.5332*mcMSmusmc^3*mkin^3*Log[musmc^2/mcMSmusmc^2] + 
        1.2570666666666666*mcMSmusmc^2*mkin^4*Log[musmc^2/mcMSmusmc^2] - 
        1.596*mcMSmusmc*mkin^5*Log[musmc^2/mcMSmusmc^2] + 
        0.4444444444444444*mkin^6*TAGdecmcMSOS*Log[musmc^2/mcMSmusmc^2] + 
        0.3333333333333333*mkin^6*TAGdecmcMSOS*Log[(1.*mus^2)/mkin^2]*
         Log[musmc^2/mcMSmusmc^2]))/mkin^5
 
RunDecmkin2mMSmc = 1.*mkin + apinl1mus*(-1.3333333333333333*mkin + 
       1.7777777777777777*mufac + (0.6666666666666666*mufac^2)/mkin - 
       1.*mkin*Log[(1.*mus^2)/mkin^2]) + apinl1mus^2*
      (0. - 1.596*mcMSmusmc - (0.17773333333333333*mcMSmusmc^3)/mkin^2 + 
       (0.6285333333333333*mcMSmusmc^2)/mkin - 14.332318723658943*mkin + 
       24.264055347179827*mufac + (7.265687421859101*mufac^2)/mkin + 
       1.0413669111716308*mkin*NL - 1.5802469135802468*mufac*NL - 
       (0.48148148148148145*mufac^2*NL)/mkin - 9.777777777777779*mufac*
        Log[(2*mufac)/mus] - (3.6666666666666665*mufac^2*Log[(2*mufac)/mus])/
        mkin + 0.5925925925925926*mufac*NL*Log[(2*mufac)/mus] + 
       (0.2222222222222222*mufac^2*NL*Log[(2*mufac)/mus])/mkin - 
       0.2962962962962963*mufac*TAGkinmc*Log[mus^2/mcMSmusmc^2] - 
       (0.1111111111111111*mufac^2*TAGkinmc*Log[mus^2/mcMSmusmc^2])/mkin - 
       6.402777777777778*mkin*Log[(1.*mus^2)/mkin^2] - 
       1.7777777777777777*mufac*Log[(1.*mus^2)/mkin^2] - 
       (0.6666666666666665*mufac^2*Log[(1.*mus^2)/mkin^2])/mkin + 
       0.36111111111111105*mkin*NL*Log[(1.*mus^2)/mkin^2] - 
       0.9583333333333333*mkin*Log[(1.*mus^2)/mkin^2]^2 + 
       0.08333333333333333*mkin*NL*Log[(1.*mus^2)/mkin^2]^2) + 
     (apinl1mus^3*(0. - 0.04933333333333332*mcMSmusmc^4*mkin^2 - 
        1.7330444444444444*mcMSmusmc^3*mkin^3 - 2.9722444444444434*
         mcMSmusmc^2*mkin^4 - 24.98211111111111*mcMSmusmc*mkin^5 - 
        199.0954170543889*mkin^6 + 0.6319407407407407*mcMSmusmc^3*mkin^2*
         mufac - 1.1173925925925925*mcMSmusmc^2*mkin^3*mufac + 
        374.8147588995222*mkin^5*mufac + 0.2369777777777778*mcMSmusmc^3*mkin*
         mufac^2 - 0.4190222222222222*mcMSmusmc^2*mkin^2*mufac^2 + 
        80.5603031498118*mkin^4*mufac^2 + 0.3003333333333333*mcMSmusmc^2*
         mkin^4*NL + 0.9826666666666665*mcMSmusmc*mkin^5*NL + 
        26.923898212450222*mkin^6*NL - 51.40701502104118*mkin^5*mufac*NL - 
        11.309438529307782*mkin^4*mufac^2*NL - 0.6526907490815436*mkin^6*
         NL^2 + 1.4473655615553573*mkin^5*mufac*NL^2 + 0.30819418434869106*
         mkin^4*mufac^2*NL^2 + 0.2716049382716049*mkin^5*mufac*TAGkinmc + 
        0.10185185185185186*mkin^4*mufac^2*TAGkinmc - 0.03433333333333333*
         mcMSmusmc^3*mkin^3*Log[(1.*mcMSmusmc)/mkin] + 
        2.7875555555555556*mcMSmusmc^2*mkin^4*Log[(1.*mcMSmusmc)/mkin] + 
        18.049*mcMSmusmc*mkin^5*Log[(1.*mcMSmusmc)/mkin] - 
        1.101333333333333*mcMSmusmc*mkin^5*NL*Log[(1.*mcMSmusmc)/mkin] - 
        283.05275696712624*mkin^5*mufac*Log[(2*mufac)/mus] - 
        85.97811719600567*mkin^4*mufac^2*Log[(2*mufac)/mus] + 
        35.97850603392235*mkin^5*mufac*NL*Log[(2*mufac)/mus] + 
        11.04749531827644*mkin^4*mufac^2*NL*Log[(2*mufac)/mus] - 
        1.0534979423868311*mkin^5*mufac*NL^2*Log[(2*mufac)/mus] - 
        0.32098765432098764*mkin^4*mufac^2*NL^2*Log[(2*mufac)/mus] + 
        53.77777777777778*mkin^5*mufac*Log[(2*mufac)/mus]^2 + 
        20.166666666666664*mkin^4*mufac^2*Log[(2*mufac)/mus]^2 - 
        6.518518518518517*mkin^5*mufac*NL*Log[(2*mufac)/mus]^2 - 
        2.4444444444444446*mkin^4*mufac^2*NL*Log[(2*mufac)/mus]^2 + 
        0.19753086419753085*mkin^5*mufac*NL^2*Log[(2*mufac)/mus]^2 + 
        0.07407407407407407*mkin^4*mufac^2*NL^2*Log[(2*mufac)/mus]^2 - 
        0.4444444444444444*mkin^6*TAGkinmc*Log[mus^2/mcMSmusmc^2] - 
        8.705302399677226*mkin^5*mufac*TAGkinmc*Log[mus^2/mcMSmusmc^2] - 
        2.653377288767849*mkin^4*mufac^2*TAGkinmc*Log[mus^2/mcMSmusmc^2] + 
        0.5267489711934156*mkin^5*mufac*NL*TAGkinmc*Log[mus^2/mcMSmusmc^2] + 
        0.16049382716049382*mkin^4*mufac^2*NL*TAGkinmc*
         Log[mus^2/mcMSmusmc^2] + 3.2592592592592595*mkin^5*mufac*TAGkinmc*
         Log[(2*mufac)/mus]*Log[mus^2/mcMSmusmc^2] + 1.2222222222222219*
         mkin^4*mufac^2*TAGkinmc*Log[(2*mufac)/mus]*Log[mus^2/mcMSmusmc^2] - 
        0.19753086419753085*mkin^5*mufac*NL*TAGkinmc*Log[(2*mufac)/mus]*
         Log[mus^2/mcMSmusmc^2] - 0.07407407407407407*mkin^4*mufac^2*NL*
         TAGkinmc*Log[(2*mufac)/mus]*Log[mus^2/mcMSmusmc^2] - 
        0.07407407407407407*mkin^6*TAGkinmc*Log[mus^2/mcMSmusmc^2]^2 + 
        0.04938271604938271*mkin^5*mufac*TAGkinmc*Log[mus^2/mcMSmusmc^2]^2 + 
        0.018518518518518517*mkin^4*mufac^2*TAGkinmc*Log[mus^2/mcMSmusmc^2]^
          2 + 0.07407407407407407*mkin^6*TAGkinmc^2*Log[mus^2/mcMSmusmc^2]^
          2 - 0.8505777777777779*mcMSmusmc^3*mkin^3*Log[(1.*mus^2)/mkin^2] + 
        2.894622222222222*mcMSmusmc^2*mkin^4*Log[(1.*mus^2)/mkin^2] - 
        7.6772222222222215*mcMSmusmc*mkin^5*Log[(1.*mus^2)/mkin^2] - 
        85.06494327546267*mkin^6*Log[(1.*mus^2)/mkin^2] - 
        27.64677139656255*mkin^5*mufac*Log[(1.*mus^2)/mkin^2] - 
        8.534205940377625*mkin^4*mufac^2*Log[(1.*mus^2)/mkin^2] + 
        0.06699999999999999*mcMSmusmc^3*mkin^3*NL*Log[(1.*mus^2)/mkin^2] - 
        0.21999999999999997*mcMSmusmc^2*mkin^4*NL*Log[(1.*mus^2)/mkin^2] + 
        0.5346666666666666*mcMSmusmc*mkin^5*NL*Log[(1.*mus^2)/mkin^2] + 
        12.67198262007108*mkin^6*NL*Log[(1.*mus^2)/mkin^2] + 
        1.6296296296296293*mkin^5*mufac*NL*Log[(1.*mus^2)/mkin^2] + 
        0.5*mkin^4*mufac^2*NL*Log[(1.*mus^2)/mkin^2] - 
        0.3201161308843708*mkin^6*NL^2*Log[(1.*mus^2)/mkin^2] + 
        9.777777777777779*mkin^5*mufac*Log[(2*mufac)/mus]*
         Log[(1.*mus^2)/mkin^2] + 3.666666666666666*mkin^4*mufac^2*
         Log[(2*mufac)/mus]*Log[(1.*mus^2)/mkin^2] - 0.5925925925925926*
         mkin^5*mufac*NL*Log[(2*mufac)/mus]*Log[(1.*mus^2)/mkin^2] - 
        0.2222222222222222*mkin^4*mufac^2*NL*Log[(2*mufac)/mus]*
         Log[(1.*mus^2)/mkin^2] - 0.3333333333333333*mkin^6*TAGkinmc*
         Log[mus^2/mcMSmusmc^2]*Log[(1.*mus^2)/mkin^2] + 
        0.2962962962962963*mkin^5*mufac*TAGkinmc*Log[mus^2/mcMSmusmc^2]*
         Log[(1.*mus^2)/mkin^2] + 0.1111111111111111*mkin^4*mufac^2*TAGkinmc*
         Log[mus^2/mcMSmusmc^2]*Log[(1.*mus^2)/mkin^2] - 
        0.05555555555555555*mkin^6*TAGkinmc*Log[mus^2/mcMSmusmc^2]^2*
         Log[(1.*mus^2)/mkin^2] + 0.05555555555555555*mkin^6*TAGkinmc^2*
         Log[mus^2/mcMSmusmc^2]^2*Log[(1.*mus^2)/mkin^2] - 
        16.521990740740737*mkin^6*Log[(1.*mus^2)/mkin^2]^2 - 
        1.7037037037037033*mkin^5*mufac*Log[(1.*mus^2)/mkin^2]^2 - 
        0.6388888888888885*mkin^4*mufac^2*Log[(1.*mus^2)/mkin^2]^2 + 
        2.229166666666667*mkin^6*NL*Log[(1.*mus^2)/mkin^2]^2 + 
        0.14814814814814814*mkin^5*mufac*NL*Log[(1.*mus^2)/mkin^2]^2 + 
        0.05555555555555555*mkin^4*mufac^2*NL*Log[(1.*mus^2)/mkin^2]^2 - 
        0.060185185185185175*mkin^6*NL^2*Log[(1.*mus^2)/mkin^2]^2 - 
        1.391203703703704*mkin^6*Log[(1.*mus^2)/mkin^2]^3 + 
        0.2314814814814814*mkin^6*NL*Log[(1.*mus^2)/mkin^2]^3 - 
        0.009259259259259259*mkin^6*NL^2*Log[(1.*mus^2)/mkin^2]^3 - 
        0.5332*mcMSmusmc^3*mkin^3*Log[musmc^2/mcMSmusmc^2] + 
        1.2570666666666666*mcMSmusmc^2*mkin^4*Log[musmc^2/mcMSmusmc^2] - 
        1.596*mcMSmusmc*mkin^5*Log[musmc^2/mcMSmusmc^2] + 
        0.4444444444444444*mkin^6*TAGkinmc*Log[musmc^2/mcMSmusmc^2] + 
        0.3333333333333333*mkin^6*TAGkinmc*Log[(1.*mus^2)/mkin^2]*
         Log[musmc^2/mcMSmusmc^2]))/mkin^5


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
 mKIN2mMS[mbKIN_, {mcMS_,musmc_}, asmus_, mus_, muf_, nl_, nloops_, CASE_]
 computes the MSbar mass from the kinetic mass
*)
mKIN2mMS[zmbKIN_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, znl_, zloops_, CASE_] := Module[
    {kk,xmc,numlist,numTAG,xmcMS,xmusmc},
    {xmcMS,xmusmc} = {zmcMS,zmusmc};
    If[zmcMS==0,xmcMS=10^-20]; (* note: set mc=0; massless part in *)
    If[zmusmc==0,xmusmc=1];    (* nm term is needed since NLxMSOS->NL-1 is chosen *)
    (numlist = {mkin      -> zmbKIN,
                mus       -> zmus,
                mcMSmusmc -> xmcMS * xmc,
                musmc     -> xmusmc,
                mufac     -> zmuf,
                NL        -> znl});
    numTAG = {DUMMY -> DUMMY};
    If[CASE == "A", numTAG = {TAGdecmcMSOS->1, TAGxmc->1, TAGkinmc->1}];
    If[CASE == "B", numTAG = {TAGdecmcMSOS->1, TAGxmc->1, TAGkinmc->1}];
    If[CASE == "C", numTAG = {TAGdecmcMSOS->1, TAGxmc->1, TAGkinmc->0}];
    If[CASE == "D", numTAG = {TAGdecmcMSOS->0, TAGxmc->0, TAGkinmc->0}];
    If[zmusmc==0,numTAG = {TAGdecmcMSOS->0, TAGxmc->0, TAGkinmc->0}];
    If[(CASE === "A") || (CASE === "D"),
       Return[ Expand[PowerExpand[RunDecmkin2mMSnl/.numlist]]
                   /. {(apinlmus^kk_/;kk>zloops) :> 0}
                   //. {Log[xmc]->0,xmc^nn_:>xmc} /. xmc->TAGxmc
                   /.numTAG
                   /. {apinlmus -> zasmus/Pi}
                   ];
    ,
       Return[ Expand[PowerExpand[RunDecmkin2mMSmc/.numlist]]
                   /. {(apinl1mus^kk_/;kk>zloops) :> 0}
                   //. {Log[xmc]->0,xmc^nn_:>xmc} /. xmc->TAGxmc
                   /.numTAG
                   //. {TAGkinmc^kk_ -> TAGkinmc}
                   /. {apinl1mus -> zasmus/Pi} ];
   ];
];

(* ************************************************************ *)
(*
 mMS2mKIN[mbMS_, {mcMS_,musmc_}, asmus_, mus_, muf_, nl_, loops_, CASE_]
 computes the kinetic mass from the MSbar mass
*)
mMS2mKIN[zmbMS_, {zmcMS_,zmusmc_}, zasmus_, zmus_, zmuf_, znl_, zloops_, CASE_] := Module[
    (*{numlist,kk, mbMSmus,mus,TAGdecmcMSOS,apinlmus,apinl1mus,mufac,NL},*)
    {kk,xmc,numlist,numTAG,xmcMS,xmusmc},
    {xmcMS,xmusmc} = {zmcMS,zmusmc};
    If[zmcMS==0,xmcMS=10^-20]; (* note: set mc=0; massless part in *)
    If[zmusmc==0,xmusmc=1];    (* nm term is needed since NLxMSOS->NL-1 is chosen *)
    (numlist = {mbMSmus   -> zmbMS,
                mus       -> zmus,
                mcMSmusmc -> xmcMS * xmc,
                musmc     -> xmusmc,
                mufac     -> zmuf,
                NL        -> znl});
    numTAG = {DUMMY -> DUMMY};
    If[CASE == "A", numTAG = {TAGdecmcMSOS->1, TAGxmc->1, TAGkinmc->1}];
    If[CASE == "B", numTAG = {TAGdecmcMSOS->1, TAGxmc->1, TAGkinmc->1}];
    If[CASE == "C", numTAG = {TAGdecmcMSOS->1, TAGxmc->1, TAGkinmc->0}];
    If[CASE == "D", numTAG = {TAGdecmcMSOS->0, TAGxmc->0, TAGkinmc->0}];
    If[zmusmc==0,numTAG = {TAGdecmcMSOS->0, TAGxmc->0, TAGkinmc->0}];
    If[(CASE === "A") || (CASE === "D"),
       Return[ Expand[PowerExpand[RunDecmMS2mkinnl/.numlist]]
                   /. {(apinlmus^kk_/;kk>zloops) :> 0}
                   //. {Log[xmc]->0,xmc^nn_:>xmc} /. xmc->TAGxmc
                   /.numTAG
                   /. {apinlmus -> zasmus/Pi}
                   ];
    ,
       Return[ Expand[PowerExpand[RunDecmMS2mkinmc/.numlist]]
                   /. {(apinl1mus^kk_/;kk>zloops) :> 0}
                   //. {Log[xmc]->0,xmc^nn_:>xmc} /. xmc->TAGxmc
                   /.numTAG
                   //. {TAGkinmc^kk_ -> TAGkinmc}
                   /. {apinl1mus -> zasmus/Pi} ];
   ];
];

(* ************************************************************ *)

End[];

(* ************************************************************ *)

EndPackage[];

(* ************************************************************ *)
(* ************************************************************ *)

