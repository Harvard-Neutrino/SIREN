"""Generate the Pythia charm-DIS splines used by PythiaDISCrossSection.

PythiaDISCrossSection samples the charm-DIS final state directly from Pythia and
uses photospline FITS tables only for the cross-section *rate* (a 1D total
sigma(E), always required) and, optionally, a 3D differential d2sigma/dx dy that
lets FinalStateProbability report a real density for reweighting. The splines
are built from the SAME Pythia configuration the sampler uses, so they close
against it. Pythia does NOT read these splines and the splines are NOT committed
to the repo (they depend on the PDF set and Pythia version) -- regenerate them
locally with this script.

Requirements (see README_charm.md):
  * SIREN built with -DSIREN_WITH_PYTHIA8=ON
  * LHAPDF_DATA_PATH set, with the requested PDF set installed
  * PYTHIA8DATA pointing at the Pythia xmldoc directory
  * DYLD_LIBRARY_PATH (macOS) including the prefix lib + the Pythia lib

Usage:
  python generate_charm_pythia_splines.py --out ./splines \
      --pdf LHAPDF6:CT18NLO --emin 100 --emax 1e6 --npoints 17

The default energy grid (100 GeV - 1 PeV) covers the IceCube diffuse
multi-cascade charm analysis band.
"""
import argparse
import os

import numpy as np

import siren.pythia_charm_splines as pcs

# nu / target PDG defaults (charged-current nu_mu on an isoscalar nucleon).
ISO_MASS = 0.9314943  # GeV


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=".", help="output directory")
    ap.add_argument("--pdf", default="LHAPDF6:CT18NLO", help="LHAPDF set (must be installed)")
    ap.add_argument("--interaction-type", type=int, default=1, help="1=CC (only CC is supported)")
    ap.add_argument("--primary-pdg", type=int, default=14, help="primary neutrino PDG (e.g. 14 = nu_mu)")
    ap.add_argument("--target-pdg", type=int, default=2212, help="target nucleon PDG (2212 p / 2112 n)")
    ap.add_argument("--emin", type=float, default=100.0, help="min neutrino energy [GeV]")
    ap.add_argument("--emax", type=float, default=1.0e6, help="max neutrino energy [GeV]")
    ap.add_argument("--npoints", type=int, default=17, help="log-spaced energy points")
    ap.add_argument("--minimum-q2", type=float, default=1.0, help="Pythia PhaseSpace:Q2Min [GeV^2]")
    ap.add_argument("--n-total", type=int, default=4000, help="events/energy for the total spline")
    ap.add_argument("--n-diff", type=int, default=20000, help="events/energy for the differential spline")
    ap.add_argument("--no-differential", action="store_true", help="only build the total spline")
    args = ap.parse_args()

    if not os.environ.get("LHAPDF_DATA_PATH"):
        raise SystemExit("LHAPDF_DATA_PATH is not set (parent of the PDF set folders).")
    pythia_data = os.environ.get("PYTHIA8DATA", "")
    os.makedirs(args.out, exist_ok=True)

    energies = np.logspace(np.log10(args.emin), np.log10(args.emax), args.npoints).tolist()
    print("energies (GeV):", [round(e, 1) for e in energies], flush=True)

    sigma_out = os.path.join(args.out, "pythia_charm_sigma.fits")
    pcs.generate_total_spline(
        sigma_out, interaction_type=args.interaction_type, primary_pdg=args.primary_pdg,
        target_pdg=args.target_pdg, target_mass=ISO_MASS, pdf_set=args.pdf,
        pythia_data_path=pythia_data, energies=energies, minimum_Q2=args.minimum_q2,
        n_events=args.n_total)
    print("wrote", sigma_out, flush=True)

    if not args.no_differential:
        diff_out = os.path.join(args.out, "pythia_charm_dsdxdy.fits")
        pcs.generate_differential_spline(
            diff_out, interaction_type=args.interaction_type, primary_pdg=args.primary_pdg,
            target_pdg=args.target_pdg, target_mass=ISO_MASS, pdf_set=args.pdf,
            pythia_data_path=pythia_data, energies=energies, minimum_Q2=args.minimum_q2,
            n_events=args.n_diff)
        print("wrote", diff_out, flush=True)


if __name__ == "__main__":
    main()
