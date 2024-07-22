(*
  RunDec.m, version 3
  msos_lightq.m
*)

<< "RunDec.m"

as6[mu_] := AlL2AlH[asMz /. NumDef,Mz /. NumDef, {{6, Mt /. NumDef, 2*Mt /. NumDef}}, mu, 5];
as5[mu_] := AlphasExact[asMz /. NumDef, Mz /. NumDef, mu, 5, 5];

Print[""];
Print["Calculate on-shell top quark mass from MS-bar top quark mass.
   Separate effects of finite bottom and charm masses. 
   alpha_s^k corrections are marked by XXX^k."];
mut  = 163; (* mt(mt) in GeV *)
Mt4  = mMS2mOS[mut, {}, as6[mut]*XXX, mut, 6, 4];
Mtb  = mMS2mOS[mut, {mub /. NumDef}, as6[mut]*XXX, mut, 6, 4];
Mtbc = mMS2mOS[mut, {{mub /. NumDef, mub /. NumDef}, {mc3 /. NumDef, 3.0}}, as6[mut]*XXX, mut, 6, 4];
delb = Expand[(Mtb - Mt4)];
delc = Expand[(Mtbc - Mtb)];

Print["m_t = ", mut];
Print["M_t = ", tmp=Expand[Mtbc]//N, " = ", tmp/.XXX->1];
Print["without light quark mass effects: ", tmp=Expand[Mt4]//N, " = ", tmp/.XXX->1];
Print["influence of bottom quark: ", delb, " = ", delb/.XXX->1];
Print["influence of charm quark: ", delc, " = ", delc/.XXX->1];

Print["Calculate on-shell bottom quark mass from MS-bar bottom quark mass.
   Separate finite charm mass effects."];
Mb4 = mMS2mOS[mub /. NumDef, {}, as5[mub /. NumDef]*XXX, mub /. NumDef, 5, 4];
Mbc = mMS2mOS[mub /. NumDef, {{mc3 /. NumDef, 3.0}}, as5[mub /. NumDef]*XXX, mub /. NumDef, 5, 4];
delc = Expand[(Mbc - Mb4)];

Print["M_b = ", tmp=Expand[Mbc], " = ", tmp/.XXX->1];
Print["without charm quark mass effects: ", tmp=Expand[Mb4], " = ", tmp/.XXX->1];
Print["influence of charm quark: ", delc, " = ", delc/.XXX->1];
