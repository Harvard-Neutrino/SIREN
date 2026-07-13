#!/usr/bin/env python
"""
example1 DIS timing test: a DarkNews-free DETECTOR simulation with the time factor.

Primary nu_mu CC DIS (CSMSDISSplines -- pre-built cross-section splines, no DarkNews)
in a real detector (IceCube by default).  A per-event initial time t0 is injected via
PrimaryExternalDistribution (a t0-only CSV added to the physical distributions), and we
verify SIREN propagates it to the interaction vertex by time of flight:

    vertex_time == t0 + |vertex - initial_position| / (beta c),   beta = |p|/E (nu -> 1)

This mirrors the CCM timing sim but uses example1's spline-based DIS -- so it runs with
NO DarkNews (no table generation, none of the DarkNews decay bugs).

Run:
  /home/shubham/siren_ubaid_venv/bin/python timing_dis_chain.py --n 200 --external
"""
import argparse
import os
import numpy as np

import siren
from siren._util import GenerateEvents


def build_injector(n_events, experiment="IceCube", t0_csv=None):
    detector_model = siren.utilities.load_detector(experiment)
    primary_type = siren.dataclasses.Particle.ParticleType.NuMu
    primary_processes, _ = siren.utilities.load_processes(
        "CSMSDISSplines",
        primary_types=[primary_type],
        target_types=[siren.dataclasses.Particle.ParticleType.Nucleon],
        isoscalar=True, process_types=["CC"])

    dists = [
        siren.distributions.PrimaryMass(0),
        siren.distributions.PowerLaw(2, 1e3, 1e6),               # energy [GeV]
        siren.distributions.IsotropicDirection(),
        siren.distributions.ColumnDepthPositionDistribution(
            600, 600.0, siren.distributions.LeptonDepthFunction()),
    ]
    # per-event t0 via a t0-only external CSV -- sets record.initial_time and leaves
    # energy/direction/position to the physical distributions above.
    if t0_csv is not None:
        dists.append(siren.distributions.PrimaryExternalDistribution(t0_csv))

    injector = siren.injection.Injector()
    injector.number_of_events = n_events
    injector.detector_model = detector_model
    injector.primary_type = primary_type
    injector.primary_interactions = primary_processes[primary_type]
    injector.primary_injection_distributions = dists
    return injector


def verify_and_write(events, outdir):
    import pyarrow as pa, pyarrow.parquet as pq, h5py
    C = siren.utilities.Constants.c
    rows = []
    for ev in events:
        tree = list(ev.tree)
        if not tree:
            continue
        r = tree[0].record                      # primary DIS interaction
        E = float(np.array(r.primary_momentum)[0]); beta = 1.0 if E > 0 else 0.0  # nu massless
        ip = np.array(r.primary_initial_position, float)
        vx = np.array(r.interaction_vertex, float)
        dist = float(np.linalg.norm(vx - ip))
        t0 = float(r.primary_initial_time); vt = float(r.interaction_time)
        rows.append((t0, vt, dist, beta, vt - t0, dist / (beta * C) if beta > 0 else 0.0, E))
    if not rows:
        print("  NO events with a primary interaction."); return False
    a = np.asarray(rows, float)
    cols = ["primary_initial_time", "vertex_time", "dist", "beta", "tof", "tof_expected", "nu_energy"]
    D = {c: a[:, i] for i, c in enumerate(cols)}
    os.makedirs(outdir, exist_ok=True)
    pq.write_table(pa.table(D), os.path.join(outdir, "dis_timing_events.parquet"))
    with h5py.File(os.path.join(outdir, "dis_timing_events.h5"), "w") as h:
        for c in cols:
            h.create_dataset(c, data=D[c])
    rel = np.abs(D["tof"] - D["tof_expected"]) / np.maximum(np.abs(D["tof_expected"]), 1e-12)
    ok = rel.max() < 1e-6
    print("  events: %d -> parquet + hdf5 in %s" % (len(rows), outdir))
    print("  timing check: vertex_time - t0 == dist/(beta c) -> max rel resid %.2e  %s"
          % (rel.max(), "PASS" if ok else "FAIL"))
    print("  <t0>=%.4g  <vertex_time>=%.4g ns  <tof>=%.4g ns  <dist>=%.4g m"
          % (D["primary_initial_time"].mean(), D["vertex_time"].mean(), D["tof"].mean(), D["dist"].mean()))
    for i in range(min(3, len(rows))):
        print("   event %d: t0=%.4g  vtx_time=%.4g  dist=%.4g m" % (i, D["primary_initial_time"][i], D["vertex_time"][i], D["dist"][i]))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--experiment", default="IceCube")
    ap.add_argument("--outdir", default="output/dis_timing")
    ap.add_argument("--external", action="store_true", help="inject per-event t0 via PrimaryExternalDistribution")
    ap.add_argument("--spill-ns", type=float, default=1600.0)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    t0_csv = None
    if args.external:
        t0_csv = os.path.join(args.outdir, "t0_input.csv")
        t0 = np.random.default_rng(1).uniform(0.0, args.spill_ns, size=max(args.n, 200))
        np.savetxt(t0_csv, t0, header="t0", comments="", fmt="%.8g")
        print("[0] t0 CSV (%d rows) -> %s" % (len(t0), t0_csv))

    print("[1/3] build DIS injector (%s, CSMSDISSplines%s) ..." % (args.experiment, " + external t0" if args.external else ""))
    injector = build_injector(args.n, args.experiment, t0_csv=t0_csv)
    print("[2/3] generate %d DIS events ..." % args.n)
    events, _ = GenerateEvents(injector)
    print("      generated %d event trees" % len(events))
    print("[3/3] extract vertex timing -> verify ...")
    ok = verify_and_write(events, args.outdir)
    tag = "t0 (PrimaryExternalDistribution) -> vertex_time by TOF" if args.external else "vertex_time by TOF (t0=0)"
    print("\nRESULT:", ("DIS DETECTOR-SIM TIMING VERIFIED: " + tag) if ok else "see above")


if __name__ == "__main__":
    main()
