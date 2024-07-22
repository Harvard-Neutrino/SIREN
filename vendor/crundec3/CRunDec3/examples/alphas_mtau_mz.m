(*
  RunDec.m, version 3
  alphas_mtau_mz.m
*)

<< "RunDec.m"

(* Numerical input values from arXiv:0801.1821: *)
mudecc = 3.0;
mudecb = mub /. NumDef;
asmtauerror = 0.016;

(* vary decoupling scale between mu/f and f*mu: *)
f = 3;

(* Calculate alpha_s(M_Z) from alpha_s(mu) using n-loop running and (n-1)-loop decoupling. *)
As5MZ[asin_?NumberQ, mu_?NumberQ, thr1_?NumberQ, thr2_?NumberQ, n_?IntegerQ] := Module[
    {as3c, as4c, as4b, as5b, as5, as5mz},
    as3c    = AlphasExact[asin, mu, thr1, 3, n];
    as4c    = DecAsUpSI[as3c, muc /. NumDef, thr1, 3, n];
    as4b    = AlphasExact[as4c, thr1, thr2, 4, n];
    as5b    = DecAsUpSI[as4b, mub /. NumDef, thr2, 4, n];
    as5mz   = AlphasExact[as5b, thr2, Mz /. NumDef, 5, n];
    Return[as5mz];
];

Print["Calculate \\alpha_s(M_Z) for two- to five-loop running."];
Do[
    as5          = As5MZ[asMtau /. NumDef, Mtau /. NumDef, mudecc, mudecb, i];
    truncerr     = Abs[as5 - As5MZ[asMtau /. NumDef, Mtau /. NumDef, mudecc, mudecb, i-1]];
    step = 1;
    tmp = Map[ As5MZ[asMtau /. NumDef, Mtau /. NumDef, #, mudecb, i]&, Range[mudecc/f, f*mudecc, step]];
    scaleuncertc = Max[ tmp ] - Min[tmp];
    tmp = Map[ As5MZ[asMtau /. NumDef, Mtau /. NumDef, mudecc, #, i]&, Range[mudecb/f, f*mudecb, step]];
    scaleuncertb = Max[ tmp ] - Min[tmp];
    totaluncert  = Sqrt[truncerr^2 + scaleuncertc^2 + scaleuncertb^2];
    Print[""];
    Print["\\alpha_s(M_Z) with ", i, "-loop running and ", i-1, "-loop decoupling:  ", 
          As5MZ[asMtau /. NumDef, Mtau /. NumDef, mudecc, mudecb, i]];
    Print["take into account uncertainty on \\alpha_s(M_\\tau):       ",
      "+", Abs[as5 - As5MZ[(asMtau /. NumDef) + asmtauerror, Mtau /. NumDef, mudecc, mudecb, i]]];
    Print["                                                         -", 
      Abs[as5 - As5MZ[(asMtau /.NumDef) - asmtauerror, Mtau /. NumDef, mudecc, mudecb, i]]];
    Print["uncertainty due to truncation of the perturbation series:    ", truncerr];
    Print["decoupling scale uncertainty for \\mu_dec,c:                  ", scaleuncertc];
    Print["decoupling scale uncertainty for \\mu_dec,b:                  ", scaleuncertb];
    Print["total uncertainty (truncation and dec. scale):         +/-", totaluncert],
    {i,2,5}
]

