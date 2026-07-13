#!/usr/bin/env python
"""
Six plots you can make from ONLY the 4 timing columns of a timing parquet/h5:
  primary_initial_time (t0), vertex_time, tof, tof_expected.

Usage: python plot_timing_4param.py <events.parquet> [-o out.png]
"""
import argparse, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyarrow.parquet as pq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parquet")
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args()
    t = pq.read_table(a.parquet)
    g = lambda n: t[n].to_numpy()
    t0 = g("primary_initial_time"); vt = g("vertex_time")
    tof = g("tof"); tofx = g("tof_expected")
    resid = tof - tofx

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))

    # 1) VALIDATION: measured TOF vs analytic TOF  (should be y=x)
    ax[0,0].scatter(tofx, tof, s=14, color="tab:blue")
    lo, hi = min(tofx.min(), tof.min()), max(tofx.max(), tof.max())
    ax[0,0].plot([lo,hi],[lo,hi],"r-",lw=1.5,label="y = x")
    ax[0,0].set_xlabel("tof_expected  dist/(βc) [ns]"); ax[0,0].set_ylabel("tof  (vertex_time−t0) [ns]")
    ax[0,0].set_title("(1) TOF closure: measured vs analytic"); ax[0,0].legend(fontsize=8)

    # 2) CLOSURE residual histogram  (accuracy of the timing feature)
    ax[0,1].hist(resid, bins=30, color="tab:purple")
    ax[0,1].set_xlabel("tof − tof_expected [ns]"); ax[0,1].set_ylabel("events")
    ax[0,1].set_title("(2) Timing residual (max |Δ|=%.1e ns)" % np.abs(resid).max())

    # 3) ARRIVAL-TIME spectrum  (what the detector sees)
    ax[0,2].hist(vt, bins=30, color="tab:green", histtype="stepfilled", alpha=0.8)
    ax[0,2].set_xlabel("vertex_time (arrival) [ns]"); ax[0,2].set_ylabel("events")
    ax[0,2].set_title("(3) Arrival-time spectrum")

    # 4) BEAM SPILL: t0 distribution
    ax[1,0].hist(t0, bins=30, color="tab:orange", histtype="stepfilled", alpha=0.8)
    ax[1,0].set_xlabel("primary_initial_time  t0 [ns]"); ax[1,0].set_ylabel("events")
    ax[1,0].set_title("(4) Injected beam-spill (t0)")

    # 5) TIME-OF-FLIGHT distribution
    ax[1,1].hist(tof, bins=30, color="tab:red", histtype="stepfilled", alpha=0.8)
    ax[1,1].set_xlabel("tof [ns]"); ax[1,1].set_ylabel("events")
    ax[1,1].set_title("(5) Time-of-flight distribution")

    # 6) CUMULATIVE arrival time  (timing-cut efficiency)
    s = np.sort(vt); cum = np.arange(1, len(s)+1)/len(s)
    ax[1,2].step(s, cum, where="post", color="k")
    ax[1,2].set_xlabel("vertex_time [ns]"); ax[1,2].set_ylabel("cumulative fraction")
    ax[1,2].set_title("(6) Arrival-time CDF (timing-cut efficiency)"); ax[1,2].set_ylim(0,1.02)

    fig.suptitle("Timing plots from 4 columns  --  %s  (N=%d)" % (os.path.basename(a.parquet), t.num_rows), fontsize=13)
    fig.tight_layout(rect=[0,0,1,0.96])
    out = a.out or os.path.splitext(a.parquet)[0] + "_4param.png"
    fig.savefig(out, dpi=130); plt.close(fig)
    print("saved ->", out)


if __name__ == "__main__":
    main()
