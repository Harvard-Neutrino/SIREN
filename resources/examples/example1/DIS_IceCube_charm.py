#!/usr/bin/env python3
"""
Charm-DIS chain in IceCube on the Layer-2 Injector/Weighter interface.

Demonstrates the full charm chain as spec-form vertices:
  primary NuE CC+NC DIS (QuarkDISFromSpline, charm-target splines)
      emits {charged lepton, Hadrons (shower), D meson} directly -- no
      intermediate "Charm" quark, no separate hadronization step.
      For a neutrino primary QuarkDISFromSpline::DTypesForPrimary emits
      {D0, D+, Ds+} (the charge conjugates {D0bar, D-, Ds-} are emitted for
      an antineutrino primary).
  one secondary vertex per D species
      DMesonELoss  (propagation energy loss, re-emits the same D)
      CharmMesonDecay  (decay into leptons + K/pi)

10,000 events on IceCube by default (--events/--seed/--output), volume
injection inside the icecube sector, astrophysical power-law weighting.

The chain is built from ``siren.Vertex`` objects and driven directly through
``siren.injection.Injector``, ``siren.injection.Weighter`` and
``siren.generate``; ``siren.Simulation`` wraps the same three calls for the
common case (see DIS_IceCube.py).

Splines are read from the SIREN_CHARM_SPLINE_DIR environment variable; point it
at your own set of QuarkDIS charm-target spline files before running.

Usage:
    export SIREN_CHARM_SPLINE_DIR=/path/to/M_Muon_New
    python3 DIS_IceCube_charm.py
"""

import argparse
import os

import numpy as np

import siren


# ----------------------------------------------------------------------------
# Config (edit for your setup)
# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--events", type=int, default=10_000)
parser.add_argument("--seed", type=int, default=1)
parser.add_argument("--output", default="output/charm_example")
args = parser.parse_args()

# Spline directory: read from the SIREN_CHARM_SPLINE_DIR environment variable.
# These QuarkDIS charm-target .fits splines are large and machine-specific, so
# they are not bundled with SIREN. Point the variable at your own spline set,
# e.g.  export SIREN_CHARM_SPLINE_DIR=/path/to/M_Muon_New
SPLINES_DIR = os.environ.get("SIREN_CHARM_SPLINE_DIR")
if not SPLINES_DIR:
    raise RuntimeError(
        "SIREN_CHARM_SPLINE_DIR is not set. Set it to the directory containing "
        "the QuarkDIS charm-target spline files "
        "(dsdxidy_nu-N-{cc,nc}-charm-*.fits and sigma_nu-N-{cc,nc}-charm-*.fits) "
        "before running this example."
    )
EXPERIMENT        = "IceCube"
PRIMARY_TYPE      = siren.particles.NuE
GEN_EMIN, GEN_EMAX = 1e2, 1e6  # generation energy range [GeV]
OXYGEN_PDF        = "EPPS21nlo_CT18Anlo_O16_central"
HYDROGEN_PDF      = "HERAPDF20_NLO_EIG_central"

PT = siren.particles


# ----------------------------------------------------------------------------
# Primary interactions
# ----------------------------------------------------------------------------

def make_quark_dis_xs(pdf, target, current_type):
    """Build one QuarkDISFromSpline for a given PDF / nuclear target / CC or NC."""
    int_type = 1 if current_type == "cc" else 2
    isoscalar_mass = (0.938272 + 0.939565) / 2
    return siren.interactions.QuarkDISFromSpline(
        os.path.join(SPLINES_DIR, f"dsdxidy_nu-N-{current_type}-charm-{pdf}.fits"),
        os.path.join(SPLINES_DIR, f"sigma_nu-N-{current_type}-charm-{pdf}.fits"),
        int(int_type),         # interaction type: 1=CC, 2=NC
        isoscalar_mass,
        1,                     # min Q^2
        [PRIMARY_TYPE],
        [target],
        "m",                   # mass units
    )


detector_model = siren.load_detector(EXPERIMENT)

primary_interactions = [
    make_quark_dis_xs(OXYGEN_PDF,   PT.O16Nucleus, "cc"),
    make_quark_dis_xs(HYDROGEN_PDF, PT.HNucleus,   "cc"),
    make_quark_dis_xs(OXYGEN_PDF,   PT.O16Nucleus, "nc"),
    make_quark_dis_xs(HYDROGEN_PDF, PT.HNucleus,   "nc"),
]

# Generation energy spectrum: flat power law with index 1 over [emin, emax]
# Astrophysical reweight: E^-2.58 normalized to the HESE flux at 100 TeV
edist_gen = siren.dist.PowerLaw(1, GEN_EMIN, GEN_EMAX)
edist_phy = siren.dist.PowerLaw(2.58, 1e2, 1e6)
edist_phy.SetNormalizationAtEnergy(1.68e-18 * 1e4 * 4 * np.pi, 1e5)

direction = siren.dist.IsotropicDirection()
position = siren.get_volume_position_distribution_from_sector(
    detector_model, "icecube")


# ----------------------------------------------------------------------------
# The chain as spec-form vertices
#
# QuarkDISFromSpline emits a D meson directly as one of the primary-DIS
# secondaries (not a bare charm quark), so no CharmHadronization step is
# needed. Each D species gets its own vertex with energy loss and decay.
#
# IMPORTANT: every D type that QuarkDISFromSpline emits MUST have a vertex
# here and be named in the primary's expand rules. A secondary type without a
# registered process is dropped by the injector (Injector.cxx:
# `if(it == secondary_process_map.end()) continue;`), so an unregistered
# species would be generated by the primary DIS and then never decay or lose
# energy -- a dangling secondary that distorts the event sample.
#
# For the NuE (neutrino) primary used here, DTypesForPrimary emits {D0, D+, Ds+}.
# (An antineutrino primary would instead emit {D0bar, D-, Ds-}; this example
# would then need D_TYPES updated to the conjugates.) Both the vertices and
# the expand rules are built from D_TYPES so the registered set cannot
# silently desync from the emitter. DMesonELoss() handles all D species; the
# 2-body CharmMesonDecay supports D0/D+/Ds+ (and their conjugates).
# ----------------------------------------------------------------------------

# Must match QuarkDISFromSpline::DTypesForPrimary for a neutrino primary.
D_TYPES = [PT.D0, PT.DPlus, PT.DsPlus]

primary = siren.Vertex(
    PRIMARY_TYPE, primary_interactions,
    distributions=[
        siren.dist.PrimaryMass(0),
        siren.dist.PrimaryNeutrinoHelicityDistribution(),
        edist_gen,
        direction,
        position,
    ],
    physical=[
        siren.dist.PrimaryMass(0),
        siren.dist.PrimaryNeutrinoHelicityDistribution(),
        edist_phy,
        direction,
    ],
    # Recurse only into the D mesons; the charged lepton and the hadronic
    # shower are final.
    expand=tuple(siren.expand.child(d) for d in D_TYPES),
)

# One vertex per D species. Energy loss re-emits the same D, which recurses
# through this vertex again until it decays; the decay products are final.
# The vertices carry no physical distributions: the secondary vertex is
# placed at its physical position and the decay is weighted by the models.
secondaries = tuple(
    siren.Vertex(
        d,
        [
            siren.interactions.DMesonELoss(),
            siren.interactions.CharmMesonDecay(primary_type=d),
        ],
        position=siren.dist.SecondaryPhysicalVertexDistribution(),
        expand=(siren.expand.child(d),),
    )
    for d in D_TYPES
)


# ----------------------------------------------------------------------------
# Injector and Weighter -- Layer-2 spec form
# ----------------------------------------------------------------------------

injector = siren.injection.Injector(
    detector=detector_model,
    primary=primary,
    secondaries=secondaries,
    events=args.events,
    seed=args.seed,
)

# The spec-form Weighter inherits the detector, the interaction models, and
# each vertex's `physical` declarations from the injector.
weighter = siren.injection.Weighter(injector)

results = siren.generate(injector, weighter, events=args.events)
print(f"Generated {len(results)} events")
# Secondaries pruned by the expand rules (leptons, Hadrons, K/pi) appear in
# the ledger as unregistered secondary types; that is the declared pruning,
# not a misconfiguration.
print(injector.report())
results.summary()


# ----------------------------------------------------------------------------
# Save events (.siren_events + .hdf5 + .parquet)
# ----------------------------------------------------------------------------

fid_vol = siren.get_fiducial_volume(EXPERIMENT)
results.save(args.output, fid_vol=fid_vol)
print(f"Saved output to {args.output}.*")
