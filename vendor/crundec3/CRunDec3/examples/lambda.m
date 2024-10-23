(*
  RunDec.m, version 3
  lambda.m
*)

<< "RunDec.m"

mudecc = 3.0;
mudecb = mub /. NumDef;

(* vary decoupling scale between mu/f and f*mu: *)
f = 3.0;

(* Calculate alpha_s^{(3)}(mu) from Lambda^{(3)}: *)
As3Lambda[lamin_?NumberQ, scale1_?NumberQ, scale2_?NumberQ, n_?IntegerQ] := Module[
    {as31, as32},
    as31 = AlphasLam[lamin, scale1, 3, n];
    as32 = AlphasExact[as31, scale1, scale2, 3, n];
    Return[as32];
];

(* Calculate alpha_s^{(5)}(M_Z) from alpha_s^{(3)}(mu): *)
As5MZ[asin_?NumberQ, mu_?NumberQ, thr1_?NumberQ, thr2_?NumberQ, n_?IntegerQ] := Module[
    {as3c, as4c, as4b, as5, as5mz},
    as3c  = AlphasExact[asin, mu, thr1, 3, n];
    as4c  = DecAsUpSI[as3c, muc /. NumDef, thr1, 3, n];
    as4b  = AlphasExact[as4c, thr1, thr2, 4, n];
    as5b  = DecAsUpSI[as4b, mub /. NumDef, thr2, 4, n];
    as5mz = AlphasExact[as5b, thr2, Mz /. NumDef, 5, n];
    Return[as5mz];
];

(* Calculate alpha_s^{(5)}(M_Z) from Lambda^{(3)}: *)
As5MZfromLambda[lamin_?NumberQ, scale1_?NumberQ, thr1_?NumberQ, thr2_?NumberQ, n_?IntegerQ] := 
    Module[{as3, as5},
    as3 = As3Lambda[lamin, scale1, thr1, n];
    as5 = As5MZ[as3, thr1, thr1, thr2, n];
    Return[as5];
];

(* Calculate alpha_s^{(5)}(M_Z) from Lambda^{(3)} including uncertainties: *)
As5MZfromLambda3[lamin_?NumberQ, lamerrp_?NumberQ, lamerrm_?NumberQ, mudecc_?NumberQ, mudecb_?NumberQ, n_?IntegerQ] := 
    Module[{as5, tmp, step, truncuncert, expuncert, scaleuncertlam, scaleuncertc, scaleuncertb, totalscaleuncert},
    as5         = As5MZfromLambda[lamin, mudecc, mudecc, mudecb, n];

    truncuncert = Abs[as5 - As5MZfromLambda[lamin, mudecc, mudecc, mudecb, n-1]];
    expuncert   = Max[Abs[as5 - As5MZfromLambda[lamin+lamerrp, mudecc, mudecc, mudecb, n]],
            Abs[as5 - As5MZfromLambda[lamin-lamerrm, mudecc, mudecc, mudecb, n]]];
    step = 1;
    tmp = Map[ As5MZfromLambda[lamin, #, mudecc, mudecb, n]&, Range[mudecc/f, f*mudecc, step]];
    scaleuncertlam = Max[ tmp ] - Min[tmp];
    tmp = Map[ As5MZfromLambda[lamin, mudecc, #, mudecb, n]&, Range[mudecc/f, f*mudecc, step]];
    scaleuncertc = Max[ tmp ] - Min[tmp];
    tmp = Map[ As5MZfromLambda[lamin, mudecc, mudecc, #, n]&, Range[mudecb/f, f*mudecb, step]];
    scaleuncertb = Max[ tmp ] - Min[tmp];
    totalscaleuncert = Sqrt[scaleuncertc^2 + scaleuncertb^2 + scaleuncertlam^2];
    Return[{as5, truncuncert, expuncert, totalscaleuncert, scaleuncertlam, scaleuncertc, scaleuncertb}];
];

(* Calculate alpha_s^{(6)}(M_t) from alpha_s^{(5)}(M_Z): *)
As6Mt[asin_?NumberQ, thr_?NumberQ, n_?IntegerQ] := Module[{as5thr, as6mt, as6thr},
    as5thr = AlphasExact[asin, Mz /. NumDef, thr, 5, n];
    as6thr = DecAsUpOS[as5thr, Mt /. NumDef, thr, 5, n];
    as6mt  = AlphasExact[as6thr, thr, Mt /. NumDef, 6, n];
    Return[as6mt];
];

(* Calculate alpha_s^{(4)}(mu_b) from alpha_s^{(5)}(M_Z): *)
As4mub[asin_?NumberQ, thr_?NumberQ, n_?IntegerQ] := Module[{as5thr, as4thr, as4mb},
    as5thr = AlphasExact[asin, Mz /. NumDef, thr, 5, n];
    as4thr = DecAsDownSI[as5thr, mub /. NumDef, thr, 4, n];
    as4mb  = AlphasExact[as4thr, thr, mub /. NumDef, 4, n];
    Return[as4mb];
];

(* Calculate alpha_s^{(3)}(3 GeV) from alpha_s^{(5)}(M_Z): *)
As33gev[asin_?NumberQ, thr1_?NumberQ, thr2_?NumberQ, n_?IntegerQ] := Module[
    {as5thr1, as4thr1, as4thr2, as3thr2, as33gev},
    as5thr1 = AlphasExact[asin, Mz /. NumDef, thr1, 5, n];
    as4thr1 = DecAsDownSI[as5thr1, mub /. NumDef, thr1, 4, n];
    as4thr2 = AlphasExact[as4thr1, thr1, thr2, 4, n];
    as3thr2 = DecAsDownSI[as4thr2, muc /. NumDef, thr2, 3, n];
    as33gev = AlphasExact[as3thr2, thr2, 3.0, 3, n];
    Return[as33gev];
];

(* Calculate Lambda^{(nf)} from alpha_s^{(nf)}(mu) and its uncertainties. *)
CalcLam[as_, aserr_, mu_, nf_, nloops_] := Module[{expl4, impl4, expl, impl,
        explplus, implplus, explminus, implminus, mean, mean4, truncuncert, diffuncert, expuncert},
    expl4       = LamExpl[as, mu, nf, nloops-1];
    impl4       = LamImpl[as, mu, nf, nloops-1];
    expl        = LamExpl[as, mu, nf, nloops];
    impl        = LamImpl[as, mu, nf, nloops];
    explplus    = Abs[LamExpl[as + aserr, mu, nf, nloops] - expl];
    implplus    = Abs[LamImpl[as + aserr, mu, nf, nloops] - impl];
    explminus   = Abs[LamExpl[as - aserr, mu, nf, nloops] - expl];
    implminus   = Abs[LamImpl[as - aserr, mu, nf, nloops] - impl];
    mean        = (expl + impl)/2;
    mean4       = (expl4 + impl4)/2;
    truncuncert = Abs[mean - mean4];
    diffuncert  = Abs[expl - impl]/2;
    expuncert   = Max[explplus, implplus, explminus, implminus];
    Return[{mean, expuncert, truncuncert, diffuncert}];
];

Print["Comparison between LamImpl, LamExpl for nf = 5:"];
Print["#loops\t LamExpl\t LamImpl\t difference (MeV)"];
Do[
    expl[i] = LamExpl[asMz /. NumDef, Mz /. NumDef, 5, i];
    impl[i] = LamImpl[asMz /. NumDef, Mz /. NumDef, 5, i];
    diff[i] = Abs[expl[i] - impl[i]];
    Print[i,"\t",10^3*expl[i],"\t\t",10^3*impl[i],"\t\t",10^3*diff[i]],
    {i,1,5}
];

Print[""];
Print["Calculate values in table 3."];
Print["In the last column in addition the uncertainty from the variation"];
Print["of the matching scale is shown."];
Print["nf\t \\Lambda^(nf)\t exp. uncert.\t \\delta_trunc\t \\delta_diff\t \\delta_scale"];
lam5 = CalcLam[asMz /. NumDef, 0.0011, Mz /. NumDef, 5, 5];
Print[5,"\t",lam5[[1]],"\t",lam5[[2]],"\t",lam5[[3]],"\t", lam5[[4]],"\t",0];
as6      = As6Mt[asMz /. NumDef, Mt /. NumDef, 5];
as6l4    = As6Mt[asMz /. NumDef, Mt /. NumDef, 4];
as6plus  = Abs[As6Mt[(asMz /. NumDef) + 0.0011, Mt /. NumDef, 5] - as6];
as6minus = Abs[As6Mt[(asMz /. NumDef) - 0.0011, Mt /. NumDef, 5] - as6];
lam6     = CalcLam[as6, Max[as6plus, as6minus], Mt /. NumDef, 6, 5];
lam6l4   = CalcLam[as6l4, 0, Mt /. NumDef, 6, 4];
step = 10;
tmp = Map[ As6Mt[N[asMz /. NumDef], #, 5]&, Range[N[Mt/f /. NumDef], N[f*Mt /. NumDef], step]];
max = Max[ tmp ]; min = Min[tmp];
Print[6,"\t",lam6[[1]],"\t",lam6[[2]],"\t",Abs[lam6[[1]] - lam6l4[[1]]],"\t",lam6[[4]], max - min];
as4      = As4mub[asMz /. NumDef, mub /. NumDef, 5];
as4l4    = As4mub[asMz /. NumDef,mub /. NumDef, 4];
as4plus  = Abs[As4mub[(asMz /. NumDef) + 0.0011, mub /. NumDef, 5] - as4];
as4minus = Abs[As4mub[(asMz /. NumDef) - 0.0011, mub /. NumDef, 5] - as4];
lam4     = CalcLam[as4, Max[as4plus, as4minus], mub /. NumDef, 4, 5];
lam4l4   = CalcLam[as4l4, 0, mub /. NumDef, 4, 4];
step = 1;
tmp = Map[ As4mub[N[asMz /. NumDef], #, 5]&, Range[N[mub/f /. NumDef], N[f*mub /. NumDef], step]];
max = Max[ tmp ]; min = Min[tmp];
Print[4,"\t",lam4[[1]],"\t",lam4[[2]],"\t",Abs[lam4[[1]] - lam4l4[[1]]],"\t",lam4[[4]],"\t",max - min];
as3      = As33gev[asMz /. NumDef, mub /. NumDef, 3.0, 5];
as3l4    = As33gev[asMz /. NumDef, mub /. NumDef, 3.0, 4];
as3plus  = Abs[As33gev[(asMz /. NumDef) + 0.0011, mub /. NumDef, 3.0, 5] - as3];
as3minus = Abs[As33gev[(asMz /. NumDef) - 0.0011, mub /. NumDef, 3.0, 5] - as3];
lam3     = CalcLam[as3, Max[as3plus, as3minus], 3.0, 3, 5];
lam3l4   = CalcLam[as3l4, 0, 3.0, 3, 4];
step = 1;
tmp = Map[ As33gev[N[asMz /. NumDef], mub /. NumDef, #, 5]&, Range[3./f, f*3., step]];
max = Max[ tmp ]; min = Min[tmp];
Print[3,"\t",lam3[[1]],"\t",lam3[[2]],"\t",Abs[lam3[[1]] - lam3l4[[1]]],"\t",lam3[[4]],"\t",max - min];


Print[""];
Print["Compute \\alpha_s similar to arXiv:1701.03075:"];
Lambda3      = 0.332;
Lambda3plus  = Lambda3 + 0.014;
Lambda3minus = Lambda3 - 0.014;

nloops = 5;

Lambda4      = DecLambdaUp[Lambda3, muc /. NumDef, 3, nloops];
Lambda5      = DecLambdaUp[Lambda4, mub /. NumDef, 4, nloops];
aslam        = AlphasLam[Lambda5, Mz /. NumDef, 5, nloops];

Lambda4l4    = DecLambdaUp[Lambda3, muc /. NumDef, 3, nloops-1];
Lambda5l4    = DecLambdaUp[Lambda4l4, mub /. NumDef, 4, nloops-1];
aslaml4      = AlphasLam[Lambda5l4, Mz /. NumDef, 5, nloops-1];

Lambda4l3    = DecLambdaUp[Lambda3, muc /. NumDef, 3, nloops-2];
Lambda5l3    = DecLambdaUp[Lambda4l3, mub /. NumDef, 4, nloops-2];
aslaml3      = AlphasLam[Lambda5l3, Mz /. NumDef, 5, nloops-2];

Lambda4plus  = DecLambdaUp[Lambda3plus, muc /. NumDef, 3, nloops];
Lambda5plus  = DecLambdaUp[Lambda4plus, mub /. NumDef, 4, nloops];
aslamplus    = AlphasLam[Lambda5plus, Mz /. NumDef, 5, nloops];

Lambda4minus = DecLambdaUp[Lambda3minus, muc /. NumDef, 3, nloops];
Lambda5minus = DecLambdaUp[Lambda4minus, mub /. NumDef, 4, nloops];
aslamminus   = AlphasLam[Lambda5minus, Mz /. NumDef, 5, nloops];

Print["\\Lambda^{(4)} = ", Lambda4, "+", Abs[Lambda4 - Lambda4plus], "-", Abs[Lambda4 - Lambda4minus]];
Print["\\Lambda^{(5)} = ", Lambda5, "+", Abs[Lambda5 - Lambda5plus], "-", Abs[Lambda5 - Lambda5minus]];
Print["\\alpha_s^{(5)}(M_Z) = ", aslam];
Print[" + ", Abs[aslam - aslamplus], " - ",Abs[aslam - aslamminus]," (\\delta\\Lambda^(3))"];
Print[" +/- ", Abs[aslam - aslaml4] + Abs[aslaml3 - aslaml4]," (truncation)"];
(* last uncertainty: see comment after eq.(5.4) of arXiv:1701.03075 *)

Print[""];
Print["Compute \\alpha_s as proposed in the paper:"];
as5mz = As5MZfromLambda3[Lambda3, 0.014, 0.014, 3.0, mub /. NumDef, 5];

Print["\\alpha_s(M_Z) with ", 5, "-loop running and ", 4, "-loop decoupling:  ", as5mz[[1]]];
Print["uncertainty due to \\delta\\Lambda:                +/-", as5mz[[2]]];
Print["uncertainty due to truncation of the perturbation series: +/-", as5mz[[3]]];
Print["scale uncertainty for \\mu_lam,c:                              ", as5mz[[5]]];
Print["scale uncertainty for \\mu_dec,c:                              ", as5mz[[6]]];
Print["scale uncertainty for \\mu_dec,b:                              ", as5mz[[7]]];
Print["uncertainty introduced through evolution:                 +/-", as5mz[[4]]];

Print[""];
Print["Compute \\alpha_s^{(3)}(M_\\tau) using \\Lambda^{(3)} from arXiv:1701.03075:"];
Lambda3      = 0.332;
Lambda3plus  = Lambda3 + 0.014;
Lambda3minus = Lambda3 - 0.014;

astau       = AlphasLam[Lambda3, Mtau /. NumDef, 3, 5];
astau4      = AlphasLam[Lambda3, Mtau /. NumDef, 3, 4];
astauplus   = Abs[AlphasLam[Lambda3plus, Mtau /. NumDef, 3, 5] - astau];
astauminus  = Abs[AlphasLam[Lambda3minus, Mtau /. NumDef, 3, 5] - astau];
truncuncert = Abs[astau - astau4];
expuncert   = Max[astauplus, astauminus];
Print["\\Lambda^{(3)} -> \\alpha_s^{(3)}(M_\\tau)"];
Print["\\alpha_s^{(3)}(M_\\tau) =                ", astau];
Print["Uncertainty due to \\delta\\Lambda^{(3)}: +/-",expuncert];
Print["Uncertainty due to truncation:          +/-", truncuncert];

Print[""];
as3Gev      = AlphasLam[Lambda3, 3.0, 3, 5];
as3Gev4     = AlphasLam[Lambda3, 3.0, 3, 4];
as3Gevplus  = AlphasLam[Lambda3plus, 3.0, 3, 5];
as3Gevminus = AlphasLam[Lambda3minus, 3.0, 3, 5];
astau       = AlphasExact[as3Gev, 3.0, Mtau /. NumDef, 3, 5];
astau4      = AlphasExact[as3Gev4, 3.0, Mtau /. NumDef, 3, 4];
astauplus   = Abs[AlphasExact[as3Gevplus, 3.0, Mtau /. NumDef, 3, 5] - astau];
astauminus  = Abs[AlphasExact[as3Gevplus, 3.0, Mtau /. NumDef, 3, 5] - astau];
truncuncert = Abs[astau - astau4];
expuncert    = Max[astauplus, astauminus];
Print["\\Lambda^{(3)} -> \\alpha_s^{(3)}(3 geV) -> \\alpha_s^{(3)}(M_\\tau)"];
Print["\\alpha_s^{(3)}(M_\\tau) =                ", astau];
Print["uncertainty due to \\delta\\Lambda^{(3)}: +/-", expuncert];
Print["uncertainty due to truncation:          +/-", truncuncert];
