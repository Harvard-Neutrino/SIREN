(*
  RunDec.m, version 3
  thr_masses.m
*)

<< "RunDec.m"


xerr = 0.002;
mbSI = mub /. NumDef;
mtOS = Mt /. NumDef;
mcSI = muc /. NumDef;

(* Module to calculate the scale-invariant mass for given
alpha_s(M_Z) and 1S, PS, RS or RS' mass including their uncertainties. *)
mThr2mSI[mThr_,merr_,asmz_,aserr_,nl_,nloop_,scheme_] := Module[
    {as, asp, asm, muf, mSI, mcentral, msoserr, expuncert, muncert, asuncert},
    as[mu_] := Switch[nl,
                5, AlphasExact[asmz, Mz /. NumDef, mu, 5, 4],
                4, AlphasExact[DecAsDownSI[AlphasExact[asmz, Mz /. NumDef, 2*mbSI, 5, 4],
                   mbSI, 2*mbSI, 4, 4], 2*mbSI, mu, 4, 4],
                3, AlphasExact[DecAsDownSI[AlphasExact[DecAsDownSI[AlphasExact[
                    asmz, Mz /. NumDef, 2*mbSI, 5, 4], mbSI, 2*mbSI, 4, 4],
                   2*mbSI, 3.0, 4, 4], mcSI, 3.0, 3, 4], 3.0, mu, 3, 4]
                ];
    asp[mu_] := Switch[nl,
                5, AlphasExact[asmz + aserr, Mz /. NumDef, mu, 5, 4],
                4, AlphasExact[DecAsDownSI[AlphasExact[asmz + aserr, Mz /. NumDef,
                   2*mbSI, 5, 4], mbSI, 2*mbSI, 4, 4], 2*mbSI, mu, 4, 4],
                3, AlphasExact[DecAsDownSI[AlphasExact[DecAsDownSI[AlphasExact[
                    asmz + aserr, Mz /. NumDef, 2*mbSI, 5, 4], mbSI, 2*mbSI, 4, 4],
                    2*mbSI, 3.0, 4, 4], mcSI, 3.0, 3, 4], 3.0, mu, 3, 4]
                ];
    asm[mu_] := Switch[nl,
                5, AlphasExact[asmz - aserr, Mz /. NumDef, mu, 5, 4],
                4, AlphasExact[DecAsDownSI[AlphasExact[asmz - aserr, Mz /. NumDef,
                   2*mbSI, 5, 4], mbSI, 2*mbSI, 4, 4], 2*mbSI, mu, 4, 4],
                3, AlphasExact[DecAsDownSI[AlphasExact[DecAsDownSI[AlphasExact[
                    asmz - aserr, Mz /. NumDef, 2*mbSI, 5, 4], mbSI, 2*mbSI, 4, 4],
                    2*mbSI, 3.0, 4, 4], mcSI, 3.0, 3, 4], 3.0, mu, 3, 4]
                ];
    muf := Switch[nl, 5, 80.0, 4, 2.0, 3, 2.0];
    
    mSI[m_,asnl_,err_] := Switch[scheme,
                "1S", m1S2mSI[m, {}, asnl, nl, nloop, 1+err],
                "PS", mPS2mSI[m, {}, asnl, muf, nl, nloop, 1+err],
                "RS", mRS2mSI[m, {}, asnl, muf, nl, nloop, 1+err],
                "RSp", mRSp2mSI[m, {}, asnl, muf, nl, nloop, 1+err]];


    mcentral := mSI[mThr, as, 0];
    msoserr := Abs[mSI[mThr, as, xerr] - mcentral];
    If[merr =!= 0,
        muncertplus := Abs[mcentral - mSI[mThr + merr[[1]], as, 0]],
        muncertplus = 0;
    ];
    If[Length[merr] > 1,
        muncertminus := Abs[mcentral - mSI[mThr - merr[[2]], as, 0]],
        muncertminus = muncertplus;
    ];
    If[aserr =!= 0,
        asuncert := Max[Abs[mcentral - mSI[mThr, asp, 0]], Abs[mcentral - mSI[mThr, asm, 0]]],
        asuncert = 0;
    ];
    Return[{mcentral, muncertplus, muncertminus, asuncert, msoserr}];
];

(* Module to calculate the MS-bar mass at a scale mu for given
alpha_s(M_Z) and 1S, PS, RS or RS' mass including their uncertainties. *)
mThr2mMS[mThr_,merr_,asmz_,aserr_,scale_,nl_,nloop_,scheme_] := Module[
    {as, asp, asm, muf, mMS, mcentral, msoserr, expuncert, muncert, asuncert},
    as := Switch[nl,
                5, AlphasExact[asmz, Mz /. NumDef, scale, 5, 4],
                4, AlphasExact[DecAsDownSI[AlphasExact[asmz, Mz /. NumDef, 2*mbSI, 5, 4],
                   mbSI, 2*mbSI, 4, 4], 2*mbSI, scale, 4, 4],
                3, AlphasExact[DecAsDownSI[AlphasExact[DecAsDownSI[AlphasExact[
                    asmz, Mz /. NumDef, 2*mbSI, 5, 4], mbSI, 2*mbSI, 4, 4],
                    2*mbSI, 3.0, 4, 4], mcSI, 3.0, 3, 4], 3.0, scale, 3, 4]
                ];
    asp := Switch[nl,
                5, AlphasExact[asmz + aserr, Mz /. NumDef, scale, 5, 4],
                4, AlphasExact[DecAsDownSI[AlphasExact[asmz + aserr, Mz /. NumDef,
                   2*mbSI, 5, nloop], mbSI, 2*mbSI, 4, 4], 2*mbSI, scale, 4, 4],
                3, AlphasExact[DecAsDownSI[AlphasExact[DecAsDownSI[AlphasExact[
                    asmz + aserr, Mz /. NumDef, 2*mbSI, 5, 4], mbSI, 2*mbSI, 4, 4],
                    2*mbSI, 3.0, 4, 4], mcSI, 3.0, 3, 4], 3.0, scale, 3, 4]
                ];
    asm := Switch[nl,
                5, AlphasExact[asmz - aserr, Mz /. NumDef, scale, 5, 4],
                4, AlphasExact[DecAsDownSI[AlphasExact[asmz - aserr, Mz /. NumDef, 2*mbSI, 5, 4],
                   mbSI, 2*mbSI, 4, 4], 2*mbSI, scale, 4, 4],
                3, AlphasExact[DecAsDownSI[AlphasExact[DecAsDownSI[AlphasExact[
                   asmz - aserr, Mz /. NumDef, 2*mbSI, 5, 4], mbSI, 2*mbSI, 4, 4],
                   2*mbSI, 3.0, 4, 4], mcSI, 3.0, 3, 4], 3.0, scale, 3, 4]
                ];
    muf := Switch[nl, 5, 80.0, 4, 2.0, 3, 2.0];
    
    mMS[m_,asnl_,err_] := Switch[scheme,
                "1S", m1S2mMS[m, {}, asnl, scale, nl, nloop, 1+err],
                "PS", mPS2mMS[m, {}, asnl, scale, muf, nl, nloop, 1+err],
                "RS", mRS2mMS[m, {}, asnl, scale, muf, nl, nloop, 1+err],
                "RSp", mRSp2mMS[m, {}, asnl, scale, muf, nl, nloop, 1+err]];


    mcentral := mMS[mThr, as, 0];
    msoserr := Abs[mMS[mThr, as, xerr] - mcentral];
    If[merr != 0,
        muncertplus := Abs[mcentral - mMS[mThr + merr[[1]], as, 0]],
        muncertplus = 0;
    ];
    If[Length[merr] > 1,
        muncertminus := Abs[mcentral - mMS[mThr - merr[[2]], as, 0]],
        muncertminus = muncertplus;
    ];
    If[aserr != 0,
        asuncert := Max[Abs[mcentral - mMS[mThr, asp, 0]], Abs[mcentral - mMS[mThr, asm, 0]]],
        asuncert = 0;
    ];
    Return[{mcentral, muncertplus, muncertminus, asuncert, msoserr}];
];

Print[""];
Print["Reproduce table 3 of arXiv:1606.06754:"];
Print["(Deviations in the last digit might occur due to rounding"];
Print["effects in the input.)"];
Print[""];
Print["mt(mt)"];
Print["input    mtPS =  mt1S =  mtRS =  mtRSp ="];
Print["#loops   168.049 172.060 166.290 171.785"];
Do[
    mfromPS = mThr2mSI[168.049, 0, asMz /. NumDef, 0, 5, i, "PS"];
    mfrom1S = mThr2mSI[172.060, 0, asMz /. NumDef, 0, 5, i, "1S"];
    mfromRS = mThr2mSI[166.290, 0, asMz /. NumDef, 0, 5, i, "RS"];
    mfromRSp = mThr2mSI[171.785, 0, asMz /. NumDef, 0, 5, i, "RSp"];
    Print[i, "       ", mfromPS[[1]], "  ", mfrom1S[[1]], "  ", mfromRS[[1]], "  ", mfromRSp[[1]]];
    ,{i,1,4}];
Print["4 (x1.002) ", mfromPS[[1]] - mfromPS[[5]], "  ", mfrom1S[[1]] - mfrom1S[[5]],
      "  ", mfromRS[[1]] - mfromRS[[5]], "  ", mfromRSp[[1]] - mfromRSp[[5]]];

Print["mb(mb)"];
Print["input    mbPS =  mb1S =  mbRS =  mbRSp ="];
Print["#loops   4.481   4.668   4.364   4.692"];
Do[
    mfromPS = mThr2mSI[4.481, 0, asMz /. NumDef, 0, 4, i, "PS"];
    mfrom1S = mThr2mSI[4.668, 0, asMz /. NumDef, 0, 4, i, "1S"];
    mfromRS = mThr2mSI[4.364, 0, asMz /. NumDef, 0, 4, i, "RS"];
    mfromRSp = mThr2mSI[4.692, 0, asMz /. NumDef, 0, 4, i, "RSp"];
    Print[i, "       ", mfromPS[[1]], "  ", mfrom1S[[1]], "  ", mfromRS[[1]], "  ", mfromRSp[[1]]];
    ,{i,1,4}];
Print["4 (x1.002) ", mfromPS[[1]] - mfromPS[[5]], "  ", mfrom1S[[1]] - mfrom1S[[5]],
      "  ", mfromRS[[1]] - mfromRS[[5]], "  ", mfromRSp[[1]] - mfromRSp[[5]]];

Print["mc(mc) table commented"];
(* ***
Print["mc(mc)"];
Print["input    mcPS =  mc1S =  mcRS =  mcRSp ="];
Print["#loops   1.130   1.513   1.035   1.351"];
Do[
    mfromPS = mThr2mSI[1.130, 0, asMz /. NumDef, 0, 3, i, "PS"];
    mfrom1S = mThr2mSI[1.513, 0, asMz /. NumDef, 0, 3, i, "1S"];
    mfromRS = mThr2mSI[1.035, 0, asMz /. NumDef, 0, 3, i, "RS"];
    mfromRSp = mThr2mSI[1.351, 0, asMz /. NumDef, 0, 3, i, "RSp"];
    Print[i, "       ", mfromPS[[1]], "  ", mfrom1S[[1]], "  ", mfromRS[[1]], "  ", mfromRSp[[1]]];
    ,{i,1,4}];
Print["4 (x1.002) ", mfromPS[[1]] - mfromPS[[5]], "  ", mfrom1S[[1]] - mfrom1S[[5]],
      "  ", mfromRS[[1]] - mfromRS[[5]], "  ", mfromRSp[[1]] - mfromRSp[[5]]];
*** *)

Print["mc(3 GeV)"];
Print["input    mcPS =  mc1S =  mcRS =  mcRSp ="];
Print["#loops   1.153   1.5145  1.043   1.357"];
Do[
    mfromPS = mThr2mMS[1.153, 0, asMz /. NumDef, 0, 3.0, 3, i, "PS"]; 
    mfrom1S = mThr2mMS[1.5145, 0, asMz /. NumDef, 0, 3.0,3, i, "1S"];
    mfromRS = mThr2mMS[1.043, 0, asMz /. NumDef, 0, 3.0,3, i, "RS"];
    mfromRSp = mThr2mMS[1.357, 0, asMz /. NumDef, 0, 3.0,3, i, "RSp"];
    Print[i, "       ", mfromPS[[1]], "  ", mfrom1S[[1]], "  ", mfromRS[[1]], "  ", mfromRSp[[1]]];
    ,{i,1,4}];
Print["4 (x1.002) ", mfromPS[[1]] - mfromPS[[5]], "  ", mfrom1S[[1]] - mfrom1S[[5]],
      "  ", mfromRS[[1]] - mfromRS[[5]], "  ", mfromRSp[[1]] - mfromRSp[[5]]];

Print[""];
Print["Calculate m_t(m_t) for a given PS mass:"];
mtMS = mThr2mSI[168.049, {0.1}, asMz /. NumDef, 0.0011, 5, 4, "PS"];
Print["mt_MS = ", mtMS[[1]]];
Print["\\delta_PS:       +/-", mtMS[[2]]];
Print["\\delta_\\alpha_s: +/-", mtMS[[4]]];


Print[""];
Print["Calculate m_b(m_b) from PS mass provided in arXiv:1411.3132:"];
mbPS      = 4.532;
mbPSplus  = 0.013;
mbPSminus = 0.039;
alsuncert = 0.0011;
resmMS = mThr2mSI[mbPS, {mbPSplus, mbPSminus}, asMz /. NumDef, alsuncert, 4, 4, "PS"];
Print[resmMS[[1]]," + ",resmMS[[2]]," - ",resmMS[[3]],"(\\delta m_PS)"];
Print[" +/- ",resmMS[[4]],"(\\delta\\alpha_s) +/- ",resmMS[[5]],"(4l-MS-OS)"];
Print[" = ",resmMS[[1]],
      " + ",Sqrt[resmMS[[2]]^2+resmMS[[4]]^2+resmMS[[5]]^2],
      " - ",Sqrt[resmMS[[3]]^2+resmMS[[4]]^2+resmMS[[5]]^2]];
