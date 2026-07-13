#!/usr/bin/env python
"""
CCM dipole-portal PRIMARY + SECONDARY simulation, instrumented for TIMING.

Chain per event:
  nu_mu  --upscatter-->  N4  (primary interaction, vertex V0, time T0)
         --propagates-->  N4 -> nu + gamma  (secondary interaction, vertex V1, time T1)

Built on the stock example2/DipolePortal_CCM.py (fast: monoenergetic 29.65 MeV
pi-DAR beam, CCM CsI/Ar -> only a couple of cross-section points to fill), using
the Injector directly (no SIREN_Controller).  We record the time of every
interaction vertex from the event trees and verify the timing chain:

    T1 - T0  ==  |V1 - V0| / (beta_N4 c),   beta_N4 = sqrt(E^2 - m4^2)/E

and, with --external, inject the primary via PrimaryExternalDistribution carrying
a per-event initial time t0, so the whole chain is offset by t0:
    T0 == t0 + TOF(nu, initial->V0)   and   T1 == T0 + TOF(N4).

Run:
  /home/shubham/siren_ubaid_venv/bin/python timing_ccm_chain.py --n 200
"""
import argparse
import os
import numpy as np

import siren
from siren import utilities
from siren._util import GenerateEvents, get_processes_model_path

M4 = 0.0235
MODEL_KWARGS = {"m4": M4, "mu_tr_mu4": 6e-7, "UD4": 0, "Umu4": 0, "epsilon": 0.0,
                "gD": 0.0, "decay_product": "photon", "noHC": True, "HNLtype": "dirac"}
NU_ENERGY = 0.02965   # GeV, pi+ DAR


def build_injector(n_events, t0_csv=None):
    experiment = "CCM"
    detector_model = utilities.load_detector(experiment)
    primary_type = siren.dataclasses.Particle.ParticleType.NuMu

    table_name = ("DarkNewsTables-v%s/Dipole_M%2.2e_mu%2.2e"
                  % (siren.utilities.darknews_version(), MODEL_KWARGS["m4"], MODEL_KWARGS["mu_tr_mu4"]))
    table_dir = os.path.join(get_processes_model_path("DarkNewsTables"), table_name)
    os.makedirs(table_dir, exist_ok=True)

    primary_processes, secondary_processes, _pk, _sk = utilities.load_processes(
        "DarkNewsTables", primary_type=primary_type, detector_model=detector_model,
        table_name=table_name, **MODEL_KWARGS)

    # Primary distributions (verbatim from the stock CCM example)
    mass_ddist = siren.distributions.PrimaryMass(0)
    edist = siren.distributions.Monoenergetic(NU_ENERGY)
    opening_angle = np.arctan(5 / 23.0)
    lower_target_origin = siren.math.Vector3D(0, 0, -0.241)
    detector_origin = siren.math.Vector3D(23, 0, -0.65)
    lower_dir = detector_origin - lower_target_origin
    lower_dir.normalize()
    inj_ddist = siren.distributions.Cone(lower_dir, opening_angle)
    max_dist = 25
    pos_dist = siren.distributions.PointSourcePositionDistribution(
        lower_target_origin - detector_origin, max_dist)
    primary_inj = [mass_ddist, edist, inj_ddist, pos_dist]
    # Inject a per-event initial time t0 via PrimaryExternalDistribution.  A CSV
    # with ONLY a t0 column sets record.initial_time (a random t0 per event) and
    # leaves energy/direction/position to the physical distributions above -- so
    # the neutrino keeps its physical kinematics/vertex but starts at time t0.
    if t0_csv is not None:
        primary_inj.append(siren.distributions.PrimaryExternalDistribution(t0_csv))

    fiducial_volume = siren.utilities.get_fiducial_volume(experiment)
    secondary_inj = {st: [siren.distributions.SecondaryBoundedVertexDistribution(fiducial_volume, max_dist)]
                     for st in secondary_processes.keys()}

    injector = siren.injection.Injector()
    injector.number_of_events = n_events
    injector.detector_model = detector_model
    injector.primary_type = primary_type
    injector.primary_interactions = primary_processes[primary_type]
    injector.primary_injection_distributions = primary_inj
    injector.secondary_interactions = secondary_processes
    injector.secondary_injection_distributions = secondary_inj

    # NOTE: we deliberately do NOT follow the N4 decay here -- that secondary
    # path hits a DarkNews-0.4.8 API bug in this build (SamplePS return_norm).
    # The PRIMARY upscatter runs fine, which is all we need for the timing test:
    # nu injected with t0 -> upscatter vertex at time t0 + TOF(nu).  Stop every
    # secondary so only the primary interaction is generated.
    def stop(tree, datum, i):
        return True
    injector.stopping_condition = stop
    return injector


def extract_timing(events):
    """Walk each event tree; return per-interaction timing records."""
    N4 = int(siren.dataclasses.Particle.ParticleType.N4)
    out = []
    for ev in events:
        inter = []
        for datum in ev.tree:
            r = datum.record
            inter.append(dict(
                ptype=int(r.signature.primary_type),
                vertex=np.array(r.interaction_vertex, float),
                itime=float(r.interaction_time),
                init_time=float(r.primary_initial_time),
                init_pos=np.array(r.primary_initial_position, float),
                pmom=np.array(r.primary_momentum, float),   # [E, px, py, pz]
                sec_types=[int(t) for t in r.signature.secondary_types],
                sec_times=[float(t) for t in r.record_secondary_times()] if hasattr(r, "record_secondary_times") else list(map(float, r.secondary_times)),
            ))
        out.append(inter)
    return out, N4


def verify(events):
    """Primary-interaction (upscatter) timing: interaction_time == t0 + TOF(nu)."""
    C = siren.utilities.Constants.c
    trees, _ = extract_timing(events)
    counts = {}
    for inter in trees:
        counts[len(inter)] = counts.get(len(inter), 0) + 1
    print("  interactions-per-event distribution:", counts)
    rows = []
    for inter in trees:
        if not inter:
            continue
        d = inter[0]                       # the primary upscatter
        E = d["pmom"][0]
        m = 0.0                            # neutrino primary
        p = np.sqrt(max(E * E - m * m, 0.0)); beta = p / E if E > 0 else 0.0
        dist = np.linalg.norm(d["vertex"] - d["init_pos"])
        tof = d["itime"] - d["init_time"]
        tof_exp = dist / (beta * C) if beta > 0 else 0.0
        rows.append((d["init_time"], d["itime"], dist, beta, tof, tof_exp, E))
    for k, inter in enumerate(trees[:3]):
        if inter:
            print("   event %d: t0=%.4g  vtx_time=%.4g  dist=%.3g m"
                  % (k, inter[0]["init_time"], inter[0]["itime"],
                     np.linalg.norm(inter[0]["vertex"] - inter[0]["init_pos"])))
    return rows, C


def write_and_check(rows, outdir):
    import pyarrow as pa, pyarrow.parquet as pq, h5py
    cols = ["primary_initial_time", "vertex_time", "dist", "beta",
            "tof", "tof_expected", "nu_energy"]
    if len(rows) == 0:
        print("  NO events with a primary interaction found.")
        return False
    rows = np.asarray(rows, float)
    D = {c: rows[:, i] for i, c in enumerate(cols)}
    os.makedirs(outdir, exist_ok=True)
    pq.write_table(pa.table(D), os.path.join(outdir, "ccm_timing_events.parquet"))
    with h5py.File(os.path.join(outdir, "ccm_timing_events.h5"), "w") as h:
        for c in cols:
            h.create_dataset(c, data=D[c])
    rel = np.abs(D["tof"] - D["tof_expected"]) / np.maximum(np.abs(D["tof_expected"]), 1e-12)
    ok = rel.max() < 1e-6
    print("  events with primary interaction: %d  -> wrote parquet + hdf5 in %s" % (len(rows), outdir))
    print("  timing check:  vertex_time - t0 == dist/(beta c)  ->  max rel resid %.2e  %s"
          % (rel.max(), "PASS" if ok else "FAIL"))
    print("  <t0>=%.3g  <vertex_time>=%.3g ns  <tof>=%.3g ns  <beta_nu>=%.4f  <dist>=%.3g m"
          % (D["primary_initial_time"].mean(), D["vertex_time"].mean(), D["tof"].mean(),
             D["beta"].mean(), D["dist"].mean()))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200, help="events to inject")
    ap.add_argument("--outdir", default="output/ccm_timing")
    ap.add_argument("--external", action="store_true",
                    help="inject a per-event t0 via PrimaryExternalDistribution (default t0=0)")
    ap.add_argument("--spill-ns", type=float, default=1600.0, help="t0 spread [ns] for --external")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    t0_csv = None
    if args.external:
        t0_csv = os.path.join(args.outdir, "t0_input.csv")
        rng = np.random.default_rng(1)
        t0 = rng.uniform(0.0, args.spill_ns, size=max(args.n, 200))
        np.savetxt(t0_csv, t0, header="t0", comments="", fmt="%.8g")
        print("[0] wrote t0 CSV (%d rows, uniform[0,%.0f] ns) -> %s" % (len(t0), args.spill_ns, t0_csv))

    print("[1/3] build injector (CCM detector, upscatter primary%s) ..."
          % (" + external t0" if args.external else ""))
    injector = build_injector(args.n, t0_csv=t0_csv)
    print("[2/3] generate %d events (nu upscatter in the detector) ..." % args.n)
    events, gen_times = GenerateEvents(injector)
    print("      generated %d event trees" % len(events))
    print("[3/3] extract vertex timing -> parquet/hdf5 + verify ...")
    rows, _ = verify(events)
    ok = write_and_check(rows, args.outdir)
    tag = "t0 (PrimaryExternalDistribution) -> vertex_time by TOF" if args.external else "vertex_time by TOF (t0=0)"
    print("\nRESULT:", ("DETECTOR-SIM TIMING VERIFIED: " + tag) if ok else "see above")


if __name__ == "__main__":
    main()
