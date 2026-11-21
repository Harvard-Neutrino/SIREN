"""Reproduce the charm dimuon kinematics from the IceCube tau/charm update slides
with the full SIREN+Pythia chain (PythiaDISCrossSection -> CharmMesonDecay).

For E_nu = 100 GeV nu_mu charged-current charm DIS: PythiaDISCrossSection samples
the primary muon + D meson, then the D decays semileptonically to a muon
(CharmMesonDecay, muon channel). This is the 'pythiaSIREN' curve in the
Pythia-vs-pythiaSIREN comparison (Diffuse WG). Saves the five slide observables
to an .npz for plotting against the slide histograms:

  th_mumu : opening angle (primary mu, decay mu)        [deg]
  zmu     : semileptonic D decay energy fraction E_mu/E_D
  y       : Bjorken-y = 1 - E_mu,prim / E_nu
  th_Dmu  : opening angle (D meson, primary mu)         [deg]
  Q2      : Q^2 = 2 E_nu E_mu,prim (1 - cos th_nu,mu)   [GeV^2]

Requires a Pythia8 build + LHAPDF + a total charm spline (SIREN_PYTHIA_WIDE_SIGMA;
see generate_charm_pythia_splines.py and README_charm.md). The means reproduce
the slide panels (theta_mumu ~9 deg core, <z> ~0.25, <y> ~0.56, theta_Dmu ~7 deg,
<Q2> ~15 GeV^2).
"""
import os
import math
import sys

import numpy as np

import siren
import siren.interactions
import siren.dataclasses
import siren.utilities

PT = siren.dataclasses.Particle.ParticleType
M_N = 0.9314943
E_NU = 100.0


def _vmag(p):
    return math.sqrt(p[1] ** 2 + p[2] ** 2 + p[3] ** 2)


def _angle_deg(a, b):
    na, nb = _vmag(a), _vmag(b)
    if na == 0 or nb == 0:
        return float("nan")
    c = (a[1] * b[1] + a[2] * b[2] + a[3] * b[3]) / (na * nb)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "siren_dimuon.npz"
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 2500
    sigma_spline = os.environ["SIREN_PYTHIA_WIDE_SIGMA"]
    pdf = os.environ.get("SIREN_PYTHIA_PDF", "LHAPDF6:CT18NLO")
    pythia_data = os.environ.get("PYTHIA8DATA", "")

    rng = siren.utilities.SIREN_random(20260506)
    xs = siren.interactions.PythiaDISCrossSection(
        "", sigma_spline, 1, M_N, 1.0, [PT.NuMu], [PT.PPlus], pythia_data, pdf, "cm")
    sigs = list(xs.GetPossibleSignaturesFromParents(PT.NuMu, PT.PPlus))

    th_mumu, zmu, yv, th_Dmu, Q2v = [], [], [], [], []
    decays = {}
    n_done = attempts = 0
    while n_done < N and attempts < N * 30:
        attempts += 1
        ir = siren.dataclasses.InteractionRecord()
        ir.signature = sigs[n_done % len(sigs)]
        ir.primary_momentum = [E_NU, 0.0, 0.0, E_NU]
        ir.primary_mass = 0.0
        ir.target_mass = M_N
        cdr = siren.dataclasses.CrossSectionDistributionRecord(ir)
        try:
            xs.SampleFinalState(cdr, rng)
        except RuntimeError:
            continue
        ir_out = siren.dataclasses.InteractionRecord()
        ir_out.signature = ir.signature
        ir_out.primary_momentum = [E_NU, 0.0, 0.0, E_NU]
        ir_out.primary_mass = 0.0
        cdr.finalize(ir_out)
        secs = list(ir_out.secondary_momenta)
        smass = list(ir_out.secondary_masses)
        p_lep, p_D = secs[0], secs[2]          # primary muon, D meson
        D_type = ir_out.signature.secondary_types[2]
        E_lep, E_D = p_lep[0], p_D[0]
        if E_D <= 0:
            continue

        if D_type not in decays:
            try:
                decays[D_type] = siren.interactions.CharmMesonDecay(primary_type=D_type)
            except Exception:
                decays[D_type] = None
        dec = decays[D_type]
        if dec is None:
            continue
        mu_sig = None
        for ds in dec.GetPossibleSignaturesFromParent(D_type):
            st = list(ds.secondary_types)
            if len(st) == 3 and st[1] in (PT.MuPlus, PT.MuMinus):
                mu_sig = ds
                break
        if mu_sig is None:
            continue
        drec = siren.dataclasses.InteractionRecord()
        drec.signature = mu_sig
        drec.primary_momentum = [p_D[0], p_D[1], p_D[2], p_D[3]]
        drec.primary_mass = smass[2]
        dcdr = siren.dataclasses.CrossSectionDistributionRecord(drec)
        try:
            dec.SampleFinalState(dcdr, rng)
        except RuntimeError:
            continue
        dout = siren.dataclasses.InteractionRecord()
        dout.signature = mu_sig
        dout.primary_momentum = [p_D[0], p_D[1], p_D[2], p_D[3]]
        dout.primary_mass = smass[2]
        dcdr.finalize(dout)
        p_mu_dec = list(dout.secondary_momenta)[1]    # decay muon

        th_mumu.append(_angle_deg(p_lep, p_mu_dec))
        zmu.append(p_mu_dec[0] / E_D)
        yv.append(1.0 - E_lep / E_NU)
        th_Dmu.append(_angle_deg(p_D, p_lep))
        Q2v.append(2.0 * E_NU * E_lep * (1.0 - p_lep[3] / _vmag(p_lep)))
        n_done += 1

    np.savez(out, th_mumu=np.array(th_mumu), zmu=np.array(zmu), y=np.array(yv),
             th_Dmu=np.array(th_Dmu), Q2=np.array(Q2v))
    print(f"wrote {out} with {n_done} events")


if __name__ == "__main__":
    main()
