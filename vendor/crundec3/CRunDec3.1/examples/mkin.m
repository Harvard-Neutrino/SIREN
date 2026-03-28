(*
Jul2020
Oct2020
Further examples for additional routines implemented in RunDec_v3.1.m
*)

<<RunDec.m

??mKIN2mMS
??mMS2mKIN

(* definitions: *)

(Mznum = 91.1876;
asMzin = 0.1179;
mbMSin = 4.163;
mc3MSin  = 0.993;
mc2MSin  = 1.0986725567495788;
mcmcMSin = 1.2785614182835914;
mucdec = 3;)

(
as5[mu_] := AlphasExact[asMzin, Mznum, mu, 5, 4];
as4[mu_] := AlphasExact[DecAsDownSI[AlphasExact[asMzin, Mznum, 2*mbMSin, 5, 4], mbMSin, 2*mbMSin, 4, 4], 2*mbMSin, mu, 4, 4];
as3[mu_] := AlphasExact[DecAsDownMS[as4[mucdec], mc3MSin, mucdec, 3, 4], mucdec, mu, 3, 4];
mbMS[muact_] := mMS2mMS[mbMSin, as5[mbMSin], as5[muact], 5, 5];
mcMS[muact_] := mMS2mMS[mcMSmbMSin, as4[mbMSin], as4[muact], 4, 5];
)

(* charm: *)
Print["Charm quark"];

mufin    = 0.5;

(
murenin = 3;        Print[ tmp=mMS2mKIN[mc3MSin, {0,0}, AA*as3[murenin], murenin, mufin, 3, 3, 3, ""]," = ",tmp/.AA->1 ];
murenin = 2;        Print[ tmp=mMS2mKIN[mc2MSin, {0,0}, AA*as3[murenin], murenin, mufin, 3, 3, 3, ""]," = ",tmp/.AA->1 ];
murenin = mcmcMSin; Print[ tmp=mMS2mKIN[mcmcMSin, {0,0}, AA*as3[murenin], murenin, mufin, 3, 3, 3, ""]," = ",tmp/.AA->1 ];
)

(* bottom: *)
Print["Bottom quark"];

(* mMS(mMS) -> mKIN(mufin) *)
Print["mMS -> mKIN"];

(musmcact = 3;
musact = mbMSin;   mbMSact = mbMS[musact];
mufin = 1;)

mcMSmbMSin = mMS2mMS[mc3MSin,as4[3],as4[mbMSin],4, 5]

Print["mb(mb) = ",mbMSin];

(* (A) *)
tmpA = Collect[mMS2mKIN[mbMSact, {mcMS[musmcact],musmcact}, AA*as3[musact], musact, mufin, 3, "A"],AA,Expand]
(* (B) *)
tmpB = Collect[mMS2mKIN[mbMSact, {mcMS[musmcact],musmcact}, AA*as4[musact], musact, mufin, 3, "B"],AA,Expand]
(* (C) *)
tmpC = Collect[mMS2mKIN[mbMSact, {mcMS[musmcact],musmcact}, AA*as4[musact], musact, mufin, 3, "C"],AA,Expand]
(* (D) *)
tmpD = Collect[mMS2mKIN[mbMSact, {mcMS[musmcact],musmcact}, AA*as3[musact], musact, mufin, 3, "D"],AA,Expand]

Print[{tmpA,tmpB,tmpC,tmpD} // MatrixForm];
Print[{tmpA,tmpB,tmpC,tmpD}/. AA->1 // MatrixForm];

Print["mKIN -> mMS"];

(* mKIN(mufin) -> mMS(mbkinin) *)

(mbkinin   = 4.550;
musmcact  = 2.;
mufin     = 1.;)

Print["mbkin = ",mbkinin];
Print["mb(",mbkinin,") ="];

(* (A) *)
tmpA = Collect[mKIN2mMS[mbkinin, {mcMS[musmcact],musmcact}, AA*as3[mbkinin], mbkinin, mufin, 3, "A"],AA,Expand]
(* (B) *)
tmpB = Collect[mKIN2mMS[mbkinin, {mcMS[musmcact],musmcact}, AA*as4[mbkinin], mbkinin, mufin, 3, "B"],AA,Expand]
(* (C) *)
tmpC = Collect[mKIN2mMS[mbkinin, {mcMS[musmcact],musmcact}, AA*as4[mbkinin], mbkinin, mufin, 3, "C"],AA,Expand]
(* (D) *)
tmpD = Collect[mKIN2mMS[mbkinin, {mcMS[musmcact],musmcact}, AA*as3[mbkinin], mbkinin, mufin, 3, "D"],AA,Expand]

Print[{tmpA,tmpB,tmpC,tmpD} // MatrixForm];
Print[{tmpA,tmpB,tmpC,tmpD}/. AA->1 // MatrixForm];
