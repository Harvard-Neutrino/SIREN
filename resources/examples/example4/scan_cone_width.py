"""
Scan pi->V1 cone half-angle and compute weighted event counts.

If biasing is correct, sum(weights) should be independent of cone width.

Requires a dk2nu beam-simulation ROOT file (not shipped with SIREN) at
sources/G4BNB/nubeam_reference.dk2nu.root relative to the run directory.
"""

import os, sys, math, numpy as np, time
os.environ["SIREN_QUIET"] = "1"

import siren
from siren._util import load_detector, GenerateEvents
from siren import _util as _su
from siren.geometry import Box

PT = siren.dataclasses.Particle.ParticleType
_dt = os.path.join(_su.resource_package_dir(), "processes", "DarkNewsTables")
_mp = _su.load_module("mp", os.path.join(_dt, "MesonProduction.py"))
_vp = _su.load_module("vp", os.path.join(_dt, "VectorPortal.py"))
from siren import dk2nu as _dk

# Load data and detector
dk = _dk.read_dk2nu("sources/G4BNB/nubeam_reference.dk2nu.root", parent_pdg=[211])
model = load_detector("SBN", detector="SBND")
fid = Box(widths=(4.0, 4.0, 5.0))
dist = _dk.dk2nu_to_primary_distribution(dk, model, parent_pdg=[211])

pion_t = PT(211); v1p = PT(5922); v1s = PT(5923); chi_t = PT(5917)

# Shared processes
offshell = _vp.VectorPortalOffShellXS(
    m_chi=8e-3, m_chi_prime=50e-3, m_V1=17e-3, m_V2=200e-3,
    g_D=1.0, epsilon=1e-4, pdgid_V1=5923)
v1ee = _vp.DarkPhotonDecay(0.017, 7e-5, pdgid_V1=5923)

# Physical (unbiased) processes for weighting
decay_phys = _mp.MesonThreeBodySIRENDecay(
    m_mediator=0.017, g_mu=1e-3, pdgid_mediator=5922)
v1chi_phys = _vp.DarkPhotonToChiDecay(
    0.017, 8e-3, 1.0, pdgid_V1=5922, pdgid_chi=5917)

def stop(tree, d, i):
    st = d.record.signature.secondary_types[i]
    pt = d.record.signature.primary_type
    if st == v1p: return False
    if st == chi_t and pt != chi_t: return False
    if st == v1s: return False
    return True

N = 500

print(f"{'cone':>8s} {'n_sig':>6s} {'eff':>6s} {'n_w':>5s}"
      f" {'sum_w':>13s} {'mean_w':>13s} {'std/m':>7s} {'time':>6s}")
print("-" * 75)

for cone_deg in [180, 45, None]:
    t0 = time.time()

    kw_pi = {"detector_position": (0,0,0), "detector_radius": 2.5}
    if cone_deg is not None:
        kw_pi["cone_half_angle"] = math.radians(cone_deg)

    decay_inj = _mp.BiasedMesonThreeBodyDecay(
        m_mediator=0.017, g_mu=1e-3, pdgid_mediator=5922, **kw_pi)
    v1chi_inj = _vp.BiasedDarkPhotonToChiDecay(
        0.017, 8e-3, 1.0, detector_position=(0,0,0), detector_radius=2.5,
        pdgid_V1=5922, pdgid_chi=5917)

    inj = siren.injection.Injector()
    inj.number_of_events = N
    inj.detector_model = model
    inj.primary_type = pion_t
    inj.primary_interactions = [decay_inj]
    inj.primary_injection_distributions = [dist]
    inj.secondary_interactions = {
        v1p: [v1chi_inj], v1s: [v1ee], chi_t: [offshell]}
    inj.secondary_injection_distributions = {
        v1p: [siren.distributions.SecondaryBoundedVertexDistribution()],
        v1s: [siren.distributions.SecondaryBoundedVertexDistribution()],
        chi_t: [siren.distributions.SecondaryBoundedVertexDistribution(fid)]}
    inj.stopping_condition = stop

    events, _ = GenerateEvents(inj)
    n_sig = sum(1 for ev in events
                for d in ev.tree if int(d.record.signature.primary_type) == 5917)

    weights = []
    n_failed = 0
    first_failure_type = None
    if n_sig > 0:
        try:
            w = siren.injection.Weighter(
                injectors=[inj],
                detector_model=model,
                primary_type=pion_t,
                primary_interactions=[decay_phys],
                primary_physical_distributions=[dist],
                secondary_interactions={
                    v1p: [v1chi_phys], v1s: [v1ee], chi_t: [offshell]},
                secondary_physical_distributions={})
        except Exception as e:
            print(f"  Weighter error: {e}")
        else:
            for ev in events:
                try:
                    wt = w.event_weight(ev)
                    if np.isfinite(wt) and wt > 0:
                        weights.append(wt)
                except Exception as e:
                    if first_failure_type is None:
                        first_failure_type = type(e)
                    elif type(e) is not first_failure_type:
                        # A different exception type than the first weighting
                        # failure indicates a bug, not a per-event edge case.
                        raise
                    n_failed += 1

    dt = time.time() - t0
    wa = np.array(weights) if weights else np.array([0.0])
    sw = np.sum(wa)
    mw = np.mean(wa) if len(wa) > 0 else 0.0
    rel = np.std(wa) / mw if mw > 0 else float("inf")
    eff = n_sig / N

    label = f"{cone_deg}" if cone_deg is not None else "auto"
    line = (f"{label:>8s} {n_sig:>6d} {eff:>6.3f} {len(weights):>5d}"
            f" {sw:>13.5e} {mw:>13.5e} {rel:>7.2f} {dt:>6.1f}s")
    print(line, flush=True)
    if n_failed > 0:
        print(f"  weighting failures: {n_failed}/{len(events)} "
              f"(exception type: {first_failure_type.__name__})", flush=True)
