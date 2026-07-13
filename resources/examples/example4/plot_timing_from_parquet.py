#!/usr/bin/env python
"""
Plot the timing observables from a SIMULATION OUTPUT PARQUET file.

Reads the per-event timing columns written by the timing sims (e.g.
timing_sim_to_parquet.py, timing_dis_chain.py, timing_ccm_chain.py) and produces:
  (1) arrival-time spectrum  (vertex_time), grouped by particle mass if present
  (2) time delay vs energy   (delay = (vertex_time - t0) - dist/c), with the
      analytic  dt = (L/c)(1/beta - 1)  overlay (L inferred from <dist>)

Recognized columns (with fallbacks):
  mass                                  -> group populations (optional)
  energy | nu_energy                    -> x-axis of plot 2
  primary_initial_time | t0             -> beam start time
  vertex_time                           -> arrival time
  dist | length                         -> baseline (for L and the light-time)
  delay                                 -> used directly if present, else computed

Usage:
  python plot_timing_from_parquet.py <events.parquet> [-o out.png]
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyarrow.parquet as pq

C = 0.299792458   # m/ns


def col(t, *names):
    for n in names:
        if n in t.column_names:
            return t[n].to_numpy()
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parquet", help="simulation output parquet with timing columns")
    ap.add_argument("-o", "--out", default=None, help="output PNG (default: alongside the parquet)")
    args = ap.parse_args()

    t = pq.read_table(args.parquet)
    vertex_time = col(t, "vertex_time")
    t0 = col(t, "primary_initial_time", "t0")
    dist = col(t, "dist", "length")
    energy = col(t, "energy", "nu_energy")
    delay = col(t, "delay")
    mass = col(t, "mass")

    if vertex_time is None or t0 is None or dist is None:
        raise SystemExit("parquet lacks required timing columns (need vertex_time, "
                         "primary_initial_time/t0, dist/length). Found: %s" % t.column_names)
    if delay is None:
        delay = (vertex_time - t0) - dist / C           # delay vs a light-speed particle

    # group by mass if available, else a single population
    if mass is not None:
        groups = sorted(set(np.round(mass, 6)))
    else:
        groups = [None]
    cmap = plt.get_cmap("viridis")
    color = {g: ("k" if (g == 0.0) else cmap(0.15 + 0.7 * i / max(len(groups) - 1, 1)))
             for i, g in enumerate(groups)}

    L = float(np.median(dist))
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))

    # ---- Plot 1: arrival-time spectrum ----
    lo, hi = np.percentile(vertex_time, [0.5, 99.5])
    tbins = np.linspace(lo, hi, 90)
    for g in groups:
        sel = np.ones(len(vertex_time), bool) if g is None else (np.round(mass, 6) == g)
        lab = ("all events" if g is None else
               (r"$\nu$ ($m=0$)" if g == 0.0 else r"$m=%g$ MeV" % (g * 1e3)))
        ax[0].hist(vertex_time[sel], bins=tbins, histtype="step", lw=2, color=color[g], label=lab)
    ax[0].axvspan(L / C, L / C + (t0.max() - t0.min()), color="0.85", zorder=0, label="prompt window")
    ax[0].set_xlabel("arrival (vertex) time [ns]"); ax[0].set_ylabel("events / bin")
    ax[0].set_title("Arrival-time spectrum  (L$\\approx$%.0f m)" % L)
    ax[0].legend(fontsize=9)

    # ---- Plot 2: adaptive ----
    beta = col(t, "beta")
    massive = ((mass is not None and len(groups) > 1 and any(g and g > 0 for g in groups))
               or (beta is not None and np.nanmin(beta) < 0.999))
    if massive and energy is not None:
        # delay vs energy: the m^2/E^2 signature of a slow massive particle
        for g in groups:
            if g == 0.0 or g is None:
                continue
            sel = (np.round(mass, 6) == g)
            ax[1].scatter(energy[sel], np.maximum(delay[sel], 1e-3), s=6, alpha=0.25, color=color[g])
            Eg = np.linspace(g + 1e-3, energy.max(), 200)
            b = np.sqrt(Eg**2 - g**2) / Eg
            ax[1].plot(Eg, (L / C) * (1.0 / b - 1.0), color=color[g], lw=2, label=r"$m=%g$ MeV" % (g * 1e3))
        ax[1].set_xlabel("energy [GeV]"); ax[1].set_yscale("log")
        ax[1].set_ylabel(r"delay $\Delta t=(L/c)(1/\beta-1)$ [ns]")
        ax[1].set_title(r"Delay vs energy  ($\Delta t\!\approx\!\frac{L}{c}\frac{m^2}{2E^2}$)")
        if ax[1].get_legend_handles_labels()[0]:
            ax[1].legend(fontsize=9)
    else:
        # beta ~ 1 (neutrino) data: show the time-of-flight relation TOF vs distance.
        tof = col(t, "tof")
        if tof is None:
            tof = vertex_time - t0
        ax[1].scatter(dist, tof, s=8, alpha=0.35, color="tab:blue", label="events")
        xs = np.array([dist.min(), dist.max()])
        ax[1].plot(xs, xs / C, "r-", lw=2, label=r"$t=L/c$  (slope $1/c$, $\beta=1$)")
        ax[1].set_xlabel("distance  $|$vertex $-$ initial$|$  [m]")
        ax[1].set_ylabel("time of flight  $t_{vtx}-t_0$  [ns]")
        ax[1].set_title(r"Time of flight vs distance  ($\beta\!\approx\!%.4f$)"
                        % (float(np.nanmean(beta)) if beta is not None else 1.0))
        ax[1].legend(fontsize=9)

    fig.suptitle("Timing from simulation parquet: %s" % os.path.basename(args.parquet), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = args.out or os.path.splitext(args.parquet)[0] + "_timing.png"
    fig.savefig(out, dpi=140); plt.close(fig)
    print("read %d events from %s" % (t.num_rows, args.parquet))
    print("saved plot -> %s" % out)


if __name__ == "__main__":
    main()
