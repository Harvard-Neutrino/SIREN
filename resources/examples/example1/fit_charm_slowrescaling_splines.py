"""Fit the slow-rescaling charm grid (from generate_charm_slowrescaling_splines)
into the FITS splines QuarkDISFromSpline consumes.

Reads xidiff.txt (3D log10 d2sigma/dxi dy on a log10(E), log10(xi), log10(y) grid)
and xisigma.txt (1D total sigma vs E), and writes

    dsdxidy_nu-N-cc-charm-CT14nlo_central.fits   (3D differential)
    sigma_nu-N-cc-charm-CT14nlo_central.fits     (1D total, cm^2)

Point SIREN_CHARM_SPLINE_DIR at the output directory to run
tests/python/test_quarkdis_slow_rescaling.py (kinematic bounds, differential
positivity, sample==density closure, and the charm-fraction normalization).
The differential stores log10(d2sigma/dxi dy) directly (QuarkDISFromSpline's
DifferentialCrossSection returns 10**value); the total stores log10(sigma_cm2).
"""
import os
import sys

import numpy as np
import photospline


def knots(c, order=2):
    c = np.asarray(c, float)
    d = c[1] - c[0]
    return np.concatenate([c[0] - d * np.arange(order, 0, -1), c,
                           c[-1] + d * np.arange(1, order + 1)])


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "."
    out = sys.argv[2] if len(sys.argv) > 2 else "."

    with open(os.path.join(src, "xidiff.txt")) as f:
        nE, nxi, ny = map(int, f.readline().split())
        Ev = [float(f.readline()) for _ in range(nE)]
        lxi = [float(f.readline()) for _ in range(nxi)]
        ly = [float(f.readline()) for _ in range(ny)]
        Z = np.full((nE, nxi, ny), -99.0)
        W = np.zeros((nE, nxi, ny))
        for i in range(nE):
            for j in range(nxi):
                for k in range(ny):
                    v, w = f.readline().split()
                    Z[i, j, k] = float(v)
                    W[i, j, k] = float(w)

    logE = np.log10(Ev)
    z, w = photospline.ndsparse.from_data(Z, W)
    sp = photospline.glam_fit(
        z, w, [logE, np.array(lxi), np.array(ly)],
        [knots(logE), knots(lxi), knots(ly)], [2, 2, 2], [1e-2, 1e-2, 1e-2], [2, 2, 2])
    diff_out = os.path.join(out, "dsdxidy_nu-N-cc-charm-CT14nlo_central.fits")
    sp.write(diff_out)

    E2, sig = np.loadtxt(os.path.join(src, "xisigma.txt"), unpack=True)
    zt, wt = photospline.ndsparse.from_data(np.log10(sig), np.ones_like(sig))
    spt = photospline.glam_fit(
        zt, wt, [np.log10(E2)], [knots(np.log10(E2))], [2], [1e-3], [2])
    tot_out = os.path.join(out, "sigma_nu-N-cc-charm-CT14nlo_central.fits")
    spt.write(tot_out)

    i100 = int(np.argmin(np.abs(E2 - 100.0)))
    print(f"wrote {diff_out}")
    print(f"wrote {tot_out}; sigma(100 GeV)={sig[i100]:.3e} cm^2  "
          f"charm_fraction={sig[i100] / (0.677e-38 * 100):.4f}")


if __name__ == "__main__":
    main()
