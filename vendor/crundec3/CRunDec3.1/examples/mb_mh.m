(*
  RunDec.m, version 3
  mb_mh.m
*)

<< "RunDec.m"

(* Module to calculate mb^{(5)}(M_h) using alpha_s^{(5)}(M_Z) and its uncertainty as input. *)
mbMh[as_,aserr_,nloop_] := Module[{mb10, mb10mberrplus, mb10mberrminus, mb10aserrplus, mb10aserrminus,
        as10, asMh, as10plus, asMhplus, as10minus, asMhminus, mbmh,
	    mbmherr, mbmhplus, mbmhminus, mbmhasplus, mbmhasminus},
    (* Calculate mb(10 GeV) and its uncertainties. *)
    mb10           = 3.610 - 12/1000*(as - 0.1189)/0.002;
    mb10mberrplus  = mb10 + 11/1000;
    mb10mberrminus = mb10 - 11/1000;
    mb10aserrplus  = mb10 + 12/1000*aserr/0.002;
    mb10aserrminus = mb10 - 12/1000*aserr/0.002;

    (* Calculate alpha_s at 10 GeV and Mh with its uncertainties. *)
    as10      = AlphasExact[as, Mz /. NumDef, 10.0, 5, nloop];
    asMh      = AlphasExact[as, Mz /. NumDef, Mh /. NumDef, 5, nloop];
    as10plus  = AlphasExact[as + aserr, Mz /. NumDef, 10.0, 5, nloop];
    asMhplus  = AlphasExact[as + aserr, Mz /. NumDef, Mh /. NumDef, 5, nloop];
    as10minus = AlphasExact[as - aserr, Mz /. NumDef, 10.0, 5, nloop];
    asMhminus = AlphasExact[as - aserr, Mz /. NumDef, Mh /. NumDef, 5, nloop];

    (* Calculate mb(Mh) and its uncertainty. *)
    mbmh        = mMS2mMS[mb10, as10, asMh, 5, nloop];
    mbmhplus    = Abs[mMS2mMS[mb10mberrplus, as10, asMh, 5, nloop] - mbmh];
    mbmhminus   = Abs[mMS2mMS[mb10mberrminus, as10, asMh, 5, nloop] - mbmh];
    mbmhasplus  = Abs[mMS2mMS[mb10aserrplus, as10plus, asMhplus, 5, nloop] - mbmh];
    mbmhasminus = Abs[mMS2mMS[mb10aserrminus, as10minus, asMhminus, 5, nloop] - mbmh];
    mbmherr     = Sqrt[Max[mbmhplus, mbmhminus]^2 + Max[mbmhasplus, mbmhasminus]^2];

    Return[{mbmh, mbmherr}];
];


Print["Reproduce mb(Mh) from arXiv:1502.00509:"];
asmzuncert = 0.0011;
asold      = 0.1189;
asolderr   = 0.002;

mbmh = mbMh[asold, asolderr, 5];
Print["mb(Mh) = ", mbmh[[1]], " \\pm ", mbmh[[2]]];

Print["Use input values from RunDec.m:"];
mbmh = mbMh[asMz /. NumDef, asmzuncert, 5];
Print["mb(Mh) = ", mbmh[[1]], " \\pm ", mbmh[[2]]];
Print[""];

Print["Compute mb^(6)(mt) from mb^(5)(10 GeV):"];
mb10   = 3.610 - 12/1000*((asMz /. NumDef) - 0.1189)/0.002;
mu10   = 10;
nloops = 5;
as10   = AlphasExact[asMz /. NumDef, Mz /. NumDef, mu10, 5, nloops];
mbmt   = mL2mH[mb10, as10, mu10, {{6, Mt/. NumDef, 2*Mt/. NumDef}}, Mt /. NumDef, nloops];
Print["mb(Mt) = ", mbmt];
