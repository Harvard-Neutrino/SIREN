#!/usr/bin/env python3
"""
Charm-DIS example using the upstream PR #71 Injector/Weighter idiom.
No SIREN_Controller.

Demonstrates the full modern charm chain:
  primary NuE CC+NC DIS (QuarkDISFromSpline, charm-target splines)
      emits {charged lepton, Hadrons (shower), D meson (D0/D+/D0Bar/D+Bar)}
      directly — no intermediate "Charm" quark, no separate hadronization step
  secondary on each D meson
      DMesonELoss  (propagation energy loss)
      CharmMesonDecay  (decay into leptons + K/pi)

10,000 events on IceCube, seed=1, volume injection inside the icecube sector,
astrophysical power-law weighting.

Spline paths below are cluster-specific; edit `SPLINES_DIR` to point at your
own set of QuarkDIS charm-target spline files.

Usage:
    python3 charm_example_new_idiom.py
"""

import os
import numpy as np

import siren
from siren._util import GenerateEvents, SaveEvents


# ----------------------------------------------------------------------------
# Config (edit for your setup)
# ----------------------------------------------------------------------------

SPLINES_DIR = (
    "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/miaochenjin/"
    "DBSearch/Simulation/Resources/Splines/M_Muon_New"
)
EXPERIMENT        = "IceCube"
PRIMARY_TYPE      = siren.dataclasses.ParticleType.NuE
NUMBER_OF_EVENTS  = 10_000
SEED              = 1
GEN_EMIN, GEN_EMAX = 1e2, 1e6  # generation energy range [GeV]
OXYGEN_PDF        = "EPPS21nlo_CT18Anlo_O16_central"
HYDROGEN_PDF      = "HERAPDF20_NLO_EIG_central"
OUTPUT_PREFIX     = "output/charm_example"


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

def cylinder_volume_position_distribution(detector_model, sector_name):
    """Build a CylinderVolumePositionDistribution for a named detector sector.

    Reproduces `SIREN_Controller.GetCylinderVolumePositionDistributionFromSector`
    without the Controller. Mirrors the inline pattern used in upstream
    `resources/examples/example1/DIS_ATLAS.py`.
    """
    geo = None
    for sector in detector_model.Sectors:
        if sector.name == sector_name:
            geo = sector.geo
            break
    if geo is None:
        raise ValueError(f"Detector sector {sector_name!r} not found")
    det_position = detector_model.GeoPositionToDetPosition(
        siren.detector.GeometryPosition(geo.placement.Position)
    )
    det_rotation  = geo.placement.Quaternion
    det_placement = siren.geometry.Placement(det_position.get(), det_rotation)
    cylinder = siren.geometry.Cylinder(
        det_placement, geo.Radius, geo.InnerRadius, geo.Z
    )
    return siren.distributions.CylinderVolumePositionDistribution(cylinder)


def make_quark_dis_xs(pdf, target, current_type):
    """Build one QuarkDISFromSpline for a given PDF / nuclear target / CC or NC."""
    int_type = 1 if current_type == "cc" else 2
    isoscalar_mass = (0.938272 + 0.939565) / 2
    return siren.interactions.QuarkDISFromSpline(
        os.path.join(SPLINES_DIR, f"dsdxdy_nu-N-{current_type}-charm-{pdf}.fits"),
        os.path.join(SPLINES_DIR, f"sigma_nu-N-{current_type}-charm-{pdf}.fits"),
        int(int_type),         # interaction type: 1=CC, 2=NC
        isoscalar_mass,
        1,                     # min Q^2
        [PRIMARY_TYPE],
        [target],
        "m",                   # mass units
    )


# ----------------------------------------------------------------------------
# Primary interactions
# ----------------------------------------------------------------------------

PT = siren.dataclasses.ParticleType

detector_model = siren.utilities.load_detector(EXPERIMENT)

primary_interactions = [
    make_quark_dis_xs(OXYGEN_PDF,   PT.O16Nucleus, "cc"),
    make_quark_dis_xs(HYDROGEN_PDF, PT.HNucleus,   "cc"),
    make_quark_dis_xs(OXYGEN_PDF,   PT.O16Nucleus, "nc"),
    make_quark_dis_xs(HYDROGEN_PDF, PT.HNucleus,   "nc"),
]

# Generation energy spectrum: flat power law with index 1 over [emin, emax]
# Astrophysical reweight: E^-2.58 normalized to the HESE flux at 100 TeV
edist_gen = siren.distributions.PowerLaw(1, GEN_EMIN, GEN_EMAX)
edist_phy = siren.distributions.PowerLaw(2.58, 1e2, 1e6)
edist_phy.SetNormalizationAtEnergy(1.68e-18 * 1e4 * 4 * np.pi, 1e5)

direction = siren.distributions.IsotropicDirection()
position  = cylinder_volume_position_distribution(detector_model, "icecube")

primary_injection_distributions = [
    siren.distributions.PrimaryMass(0),
    siren.distributions.PrimaryNeutrinoHelicityDistribution(),
    edist_gen,
    direction,
    position,
]
primary_physical_distributions = [
    siren.distributions.PrimaryMass(0),
    siren.distributions.PrimaryNeutrinoHelicityDistribution(),
    edist_phy,
    direction,
]


# ----------------------------------------------------------------------------
# Secondary interactions (D+ / D0)
#
# QuarkDISFromSpline emits a D meson directly as one of the primary-DIS
# secondaries (not a bare charm quark), so no CharmHadronization step is
# needed. We attach energy-loss and decay to each D species.
# ----------------------------------------------------------------------------

secondary_interactions = {
    PT.DPlus: [
        siren.interactions.DMesonELoss(),
        siren.interactions.CharmMesonDecay(primary_type=PT.DPlus),
    ],
    PT.D0: [
        siren.interactions.DMesonELoss(),
        siren.interactions.CharmMesonDecay(primary_type=PT.D0),
    ],
}

sec_vertex_dist = [siren.distributions.SecondaryPhysicalVertexDistribution()]
secondary_injection_distributions = {
    PT.DPlus: sec_vertex_dist,
    PT.D0:    sec_vertex_dist,
}
secondary_physical_distributions = {
    PT.DPlus: [],
    PT.D0:    [],
}


# ----------------------------------------------------------------------------
# Injector — new idiom, property-setter API
# ----------------------------------------------------------------------------

injector = siren.injection.Injector()
injector.seed                              = SEED
injector.number_of_events                  = NUMBER_OF_EVENTS
injector.detector_model                    = detector_model
injector.primary_type                      = PRIMARY_TYPE
injector.primary_interactions              = primary_interactions
injector.primary_injection_distributions   = primary_injection_distributions
injector.secondary_interactions            = secondary_interactions
injector.secondary_injection_distributions = secondary_injection_distributions
injector.stopping_condition                = lambda datum, i: False

events, gen_times = GenerateEvents(injector)
print(f"Generated {len(events)} events")


# ----------------------------------------------------------------------------
# Weighter — new idiom, built after the injector is materialised
# ----------------------------------------------------------------------------

weighter = siren.injection.Weighter()
weighter.injectors                       = [injector]
weighter.detector_model                  = detector_model
weighter.primary_type                    = PRIMARY_TYPE
weighter.primary_interactions            = list(primary_interactions)
weighter.primary_physical_distributions  = primary_physical_distributions
weighter.secondary_interactions          = secondary_interactions
weighter.secondary_physical_distributions = secondary_physical_distributions


# ----------------------------------------------------------------------------
# Save events (.siren_events + .hdf5 + .parquet)
# ----------------------------------------------------------------------------

os.makedirs(os.path.dirname(OUTPUT_PREFIX), exist_ok=True)
fid_vol = siren.utilities.get_fiducial_volume(EXPERIMENT)
SaveEvents(events, weighter, gen_times, fid_vol=fid_vol, output_filename=OUTPUT_PREFIX)
print(f"Saved output to {OUTPUT_PREFIX}.*")
