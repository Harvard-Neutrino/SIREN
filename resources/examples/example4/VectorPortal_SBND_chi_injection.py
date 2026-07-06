"""
Vector-portal dark matter at SBND: chi-first injection.

Injects chi directly into the detector with a pre-computed flux,
bypassing the geometric acceptance of the pion decay -> V1 -> chi chain.
The off-shell process chi + Ar -> chi + V1 + Ar produces the signal
V1 which decays to e+e-.

This is the simplest working example for computing signal rates.
The pion-first version requires additional biasing on the V1->chi chi
decay to achieve non-zero injection efficiency.

Reference: Dutta et al., PRL 129, 111803 (2022) [arXiv:2110.11944]

Requires dk2nu beam-simulation ROOT files (not shipped with SIREN); set
DK2NU_DIR to their directory.
"""

import os
import sys
import glob
import math
import numpy as np

import siren
from siren._util import load_detector, GenerateEvents
from siren import _util as _siren_util
from siren.math import Vector3D
from siren.geometry import Box

# ---------------------------------------------------------------------------
# Import model classes
# ---------------------------------------------------------------------------
_dt_base = os.path.join(
    _siren_util.resource_package_dir(), "processes", "DarkNewsTables",
)
_mod_vp = _siren_util.load_module(
    "siren.resources.processes.DarkNewsTables.VectorPortal",
    os.path.join(_dt_base, "VectorPortal.py"),
)
_mod_dk = _siren_util.load_module(
    "siren.resources.processes.DarkNewsTables.Dk2nuReader",
    os.path.join(_dt_base, "Dk2nuReader.py"),
)

VectorPortalOffShellXS = _mod_vp.VectorPortalOffShellXS
DarkPhotonDecay = _mod_vp.DarkPhotonDecay
compute_chi_flux_from_dk2nu = _mod_vp.compute_chi_flux_from_dk2nu
read_dk2nu = _mod_dk.read_dk2nu
print_summary = _mod_dk.print_summary

# ---------------------------------------------------------------------------
# Model parameters
# ---------------------------------------------------------------------------
M_CHI       = 8e-3
M_CHI_PRIME = 50e-3
M_V1        = 17e-3
M_V2        = 200e-3
G_D         = 1.0
EPSILON_1   = 7e-5
EPSILON_2   = 1e-4

M_PION = 0.13957039
M_MUON = 0.10565837

PDGID_V1  = 5922
PDGID_CHI = 5917

events_to_inject = 10000

# ---------------------------------------------------------------------------
# 1. Read dk2nu and compute chi flux
# ---------------------------------------------------------------------------
dk2nu_dir = os.environ.get(
    "DK2NU_DIR",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), r"../../../../G4BNB"),
)
dk2nu_files = sorted(glob.glob(os.path.join(dk2nu_dir, "*dk2nu*.root")))
if not dk2nu_files:
    dk2nu_files = sorted(glob.glob(os.path.join(dk2nu_dir, "nubeam*.root")))

print(f"Reading {len(dk2nu_files)} dk2nu file(s) ...")
dk2nu_data = read_dk2nu(dk2nu_files, parent_pdg=[211])
print_summary(dk2nu_data)

# ---------------------------------------------------------------------------
# 2. Load SBND detector
# ---------------------------------------------------------------------------
print("\nLoading SBND detector model ...")
detector_model = load_detector("SBN", detector="SBND")

chi_type = siren.dataclasses.Particle.ParticleType(PDGID_CHI)
v1_type  = siren.dataclasses.Particle.ParticleType(PDGID_V1)

# Fiducial volume
fiducial = Box(widths=(4.0, 4.0, 5.0))

# ---------------------------------------------------------------------------
# 3. Compute chi flux from dk2nu pion spectrum
# ---------------------------------------------------------------------------
print("Computing chi flux ...")

# Threshold for chi scattering
E_threshold = ((M_CHI_PRIME + 37.215)**2 - M_CHI**2 - 37.215**2) / (2.0 * 37.215)

chi_flux = compute_chi_flux_from_dk2nu(
    dk2nu_data=dk2nu_data,
    parent_pdg=211,
    m_meson=M_PION,
    m_lepton=M_MUON,
    m_V1=M_V1,
    m_chi=M_CHI,
    m_chi_prime=M_CHI_PRIME,
    g_D=G_D,
    epsilon=EPSILON_1,
    min_energy=E_threshold,
    max_energy=3.0,
    n_bins=100,
)

chi_flux_gen = compute_chi_flux_from_dk2nu(
    dk2nu_data=dk2nu_data,
    parent_pdg=211,
    m_meson=M_PION,
    m_lepton=M_MUON,
    m_V1=M_V1,
    m_chi=M_CHI,
    m_chi_prime=M_CHI_PRIME,
    g_D=G_D,
    epsilon=EPSILON_1,
    min_energy=E_threshold,
    max_energy=3.0,
    n_bins=100,
    physically_normalized=False,
)

# ---------------------------------------------------------------------------
# 4. Set up processes
# ---------------------------------------------------------------------------
print("Setting up processes ...")

offshell_xs = VectorPortalOffShellXS(
    m_chi=M_CHI,
    m_chi_prime=M_CHI_PRIME,
    m_V1=M_V1,
    m_V2=M_V2,
    g_D=G_D,
    epsilon=EPSILON_2,
    pdgid_V1=PDGID_V1,
    nuclear_pdgid=1000180400,
    nuclear_mass=37.215,
    nuclear_name="Ar40",
    A=40, Z=18,
)
print(f"  chi threshold: {E_threshold:.4f} GeV")

v1_decay = DarkPhotonDecay(M_V1, EPSILON_1, pdgid_V1=PDGID_V1)
print(f"  V1->ee width: {v1_decay._total_width:.4e} GeV")

primary_processes = {chi_type: [offshell_xs]}
secondary_processes = {
    v1_type: [v1_decay],
}

# ---------------------------------------------------------------------------
# 5. Distributions
# ---------------------------------------------------------------------------
mass_dist = siren.distributions.PrimaryMass(M_CHI)
direction_dist = siren.distributions.FixedDirection(Vector3D(0, 0, 1.0))

# Position: force scattering within fiducial volume
# DecayRangeFunction(mass, decay_width, n_decay_lengths, baseline_m)
# For chi (stable), use a very large decay length
decay_range = siren.distributions.DecayRangeFunction(
    M_CHI, 1e-30, 3, 110.0  # SBND baseline ~110m from BNB
)
position_dist = siren.distributions.RangePositionDistribution(
    2.0, 2.0, decay_range,
    set(detector_model.GetAvailableTargets(
        siren.detector.DetectorPosition(Vector3D(0, 0, 0))
    ))
)

primary_injection_distributions = [
    mass_dist,
    chi_flux_gen,
    direction_dist,
    position_dist,
]

primary_physical_distributions = [
    chi_flux,
    direction_dist,
]

secondary_injection_distributions = {
    v1_type: [siren.distributions.SecondaryBoundedVertexDistribution()],
}

def stop(tree, datum, i):
    st = datum.record.signature.secondary_types[i]
    if st == v1_type:
        return False
    return True

# ---------------------------------------------------------------------------
# 6. Run
# ---------------------------------------------------------------------------
print(f"\nInjecting {events_to_inject} events ...")

injector = siren.injection.Injector()
injector.number_of_events = events_to_inject
injector.detector_model = detector_model
injector.primary_type = chi_type
injector.primary_interactions = primary_processes[chi_type]
injector.primary_injection_distributions = primary_injection_distributions
injector.secondary_interactions = secondary_processes
injector.secondary_injection_distributions = secondary_injection_distributions
injector.stopping_condition = stop

events, gen_times = GenerateEvents(injector)

n_total = len(events)
total_time = sum(gen_times)
print(f"Generated {n_total} events in {total_time:.1f} s")

# Count signal events
n_scatter = sum(1 for ev in events for d in ev.tree
                if int(d.record.signature.primary_type) == PDGID_CHI and d.depth() == 0)
n_v1_decay = sum(1 for ev in events for d in ev.tree
                 if int(d.record.signature.primary_type) == PDGID_V1)

print(f"\n--- Results ---")
print(f"chi scattering events: {n_scatter}")
print(f"V1 -> e+e- events: {n_v1_decay}")
print(f"Injection efficiency: {n_scatter / events_to_inject:.2e}")

total_pot = dk2nu_data.get("pot", 0.0)
if total_pot > 0:
    n_pions = len(dk2nu_data["E"])
    pot_per_pion = total_pot / n_pions
    print(f"\nPOT per dk2nu pion: {pot_per_pion:.3e}")
