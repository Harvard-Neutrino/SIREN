(*
  RunDec.m, version 3
  alphas_cms.m
*)

<< "RunDec.m"

(* Numerical input for alpha_s^(6)(mu) and mu from arXiv:1304.7498: *)
muin  = 896;
as6mu = 0.0889;

Print[""];
as5Mz = AlH2AlL[as6mu,muin,{{6,Mt/.NumDef,2*Mt/.NumDef}}, Mz/.NumDef, 5];
Print["\\alpha_s^(5)(M_Z) interpreting \\alpha_s(896 GeV) as \\alpha_s^{(6)}: ", as5Mz];
as5Mz = AlphasExact[as6mu, muin, Mz/.NumDef, 5, 5];
Print["\\alpha_s^(5)(M_Z) interpreting \\alpha_s(896 GeV) as \\alpha_s^{(5)}: ", as5Mz];


(* Numerical input values for alpha_s^(5)(Mz) and 'Q' from arXiv:1609.05331: *)
asmz      = 0.1162;
asmzplus  = asmz + 0.0070;
asmzminus = asmz - 0.0062;
Q         = 1508.04;

(* vary decoupling scale between mu/f and f*mu: *)
f = 4;

(* Compute alpha_s^(5)(Q): *)

(* 'nloops': number of loops used for running: *)
nloops = 2;

asQ            = AlphasExact[asmz, Mz /. NumDef, Q, 5, nloops];
delasQplus     = Abs[asQ - AlphasExact[asmzplus, Mz /. NumDef, Q, 5, nloops]];
delasQminus    = Abs[asQ - AlphasExact[asmzminus, Mz /. NumDef, Q, 5, nloops]];
truncuncert    = Abs[asQ - AlphasExact[asmz, Mz /. NumDef, Q, 5, nloops-1]];
totuncertplus  = Sqrt[truncuncert^2 + delasQplus^2];
totuncertminus = Sqrt[truncuncert^2 + delasQminus^2];

Print[""];
Print["\\alpha_s^(5)(1508.04 GeV) using 2-loop evolution: ", asQ];
Print["uncertainty from \\delta\\alpha_s:                 +",delasQplus, " - ", delasQminus];
Print["uncertainty due to truncation of pert. series: +/-", truncuncert];
Print["total uncertainty:                               +", totuncertplus, " -", totuncertminus];

nloops = 5;

asQ            = AlphasExact[asmz, Mz /. NumDef, Q, 5, nloops];
delasQplus     = Abs[asQ - AlphasExact[asmzplus, Mz /. NumDef, Q, 5, nloops]];
delasQminus    = Abs[asQ - AlphasExact[asmzminus, Mz /. NumDef, Q, 5, nloops]];
truncuncert    = Abs[asQ - AlphasExact[asmz, Mz /. NumDef, Q, 5, nloops-1]];
totuncertplus  = Sqrt[truncuncert^2 + delasQplus^2];
totuncertminus = Sqrt[truncuncert^2 + delasQminus^2];

Print[""];
Print["\\alpha_s^(5)(1508.04 GeV) using 5-loop evolution: ",asQ];
Print["uncertainty from \\delta\\alpha_s:                 +",delasQplus, " - ", delasQminus];
Print["uncertainty due to truncation of pert. series: +/-", truncuncert];
Print["total uncertainty:                               +", totuncertplus, " -", totuncertminus];

(* Compute alpha_s^(6)(Q): *)

nloops = 5;
(* nloops = 2; *)

dec            = {{6, Mt /. NumDef, 2*Mt /. NumDef}};
asQ            = AlL2AlH[asmz, Mz /. NumDef, dec, Q, nloops];
delasQplus     = Abs[asQ - AlL2AlH[asmzplus, Mz /. NumDef, dec, Q, nloops]];
delasQminus    = Abs[asQ - AlL2AlH[asmzminus, Mz /. NumDef, dec, Q, nloops]];
truncuncert    = Abs[asQ - AlL2AlH[asmz, Mz /. NumDef, dec, Q, nloops-1]];
step = 5;
tmp = Map[ AlL2AlH[asmz, N[Mz /. NumDef], {{6, N[Mt /. NumDef], #}}, Q, nloops]&, 
   Range[N[Mt /. NumDef]/f, f*N[Mt /. NumDef], step]];
max = Max[ tmp ]; min = Min[tmp];
scaleuncert    = max - min;
totuncertplus  = Sqrt[truncuncert^2 + delasQplus^2 + scaleuncert^2];
totuncertminus = Sqrt[truncuncert^2 + delasQminus^2 + scaleuncert^2];

Print[""];
Print["\\alpha_s^(6)(1508.04 GeV) using ",nloops,"-loop evolution: ", asQ];
Print["uncertainty from \\delta\\alpha_s:                 +", delasQplus, " - ", delasQminus];
Print["uncertainty due to truncation of pert. series:   +/-", truncuncert];
Print["minimal and maximal values of alpha_s in interval"];
Print["       ",N[Mt /. NumDef]/f," <= mu_dec <= ",f*N[Mt /. NumDef],"."];
Print["max: ", max, ", min: ", min];
Print["uncertainty due to variation of decoupling scale: +/-", scaleuncert];
Print["total uncertainty:                                  +", totuncertplus, " -", totuncertminus];

Print[""];
asQ = AlL2AlH[asmz, Mz /. NumDef, {{6, Mt /. NumDef, Mz /. NumDef}}, Q, nloops];
Print["alpha_s^(6)(Q) for decoupling scale M_Z:   ", asQ];
asQ = AlL2AlH[asmz, Mz /. NumDef, {{6, Mt /. NumDef, Mt /. NumDef}}, Q, nloops];
Print["alpha_s^(6)(Q) for decoupling scale M_t:   ", asQ];
asQ = AlL2AlH[asmz, Mz /. NumDef, {{6, Mt /. NumDef, 2*Mt /. NumDef}}, Q, nloops];
Print["alpha_s^(6)(Q) for decoupling scale 2*M_t: ", asQ];
asQ = AlL2AlH[asmz, Mz /. NumDef, {{6, Mt /. NumDef, 5*Mt /. NumDef}}, Q, nloops];
Print["alpha_s^(6)(Q) for decoupling scale 5*M_t: ", asQ];
asQ = AlL2AlH[asmz, Mz /. NumDef, {{6, Mt /. NumDef, Q}}, Q, nloops];
Print["alpha_s^(6)(Q) for decoupling scale Q:     ", asQ];

