"""
Vector-portal dark matter at SBND, starting from pions.

Injection chain (pion-first, off-shell chi'):
    1. Read pi+ from BNB dk2nu ROOT files
    2. pi+ -> mu+ nu_mu V1            (three-body decay)
    3. V1 -> chi chi                   (dark photon -> DM pair)
    4. chi propagates to SBND
    5. chi + Ar -> chi + V1 + Ar       (off-shell chi', single vertex)
    6. V1 -> e+ e-                     (visible signal)

Biasing:
    - Pion decay vertex from dk2nu (PrimaryExternalDistribution)
    - chi scattering forced inside detector (SecondaryBoundedVertexDistribution)
    - Event weight accounts for interaction probability

Reference: Dutta et al., PRL 129, 111803 (2022) [arXiv:2110.11944]
"""

import os
import sys
import glob
import numpy as np

import siren
from siren import utilities
from siren._util import GenerateEvents
from siren.math import Vector3D
from siren.geometry import Box

# ---------------------------------------------------------------------------
# Import model classes via SIREN's module loader
# ---------------------------------------------------------------------------
from siren import _util as _siren_util

_dt_base = os.path.join(
    _siren_util.resource_package_dir(), "processes", "DarkNewsTables",
)

_mod_mp = _siren_util.load_module(
    "siren.resources.processes.DarkNewsTables.MesonProduction",
    os.path.join(_dt_base, "MesonProduction.py"),
)
_mod_vp = _siren_util.load_module(
    "siren.resources.processes.DarkNewsTables.VectorPortal",
    os.path.join(_dt_base, "VectorPortal.py"),
)
_mod_dk = _siren_util.load_module(
    "siren.resources.processes.DarkNewsTables.Dk2nuReader",
    os.path.join(_dt_base, "Dk2nuReader.py"),
)

MesonThreeBodySIRENDecay = _mod_mp.MesonThreeBodySIRENDecay
BiasedMesonThreeBodyDecay = _mod_mp.BiasedMesonThreeBodyDecay
VectorPortalOffShellXS = _mod_vp.VectorPortalOffShellXS
DarkPhotonDecay = _mod_vp.DarkPhotonDecay
BiasedDarkPhotonToChiDecay = _mod_vp.BiasedDarkPhotonToChiDecay
read_dk2nu = _mod_dk.read_dk2nu
dk2nu_to_primary_distribution = _mod_dk.dk2nu_to_primary_distribution
print_summary = _mod_dk.print_summary
PTYPE_PIPLUS = _mod_dk.PTYPE_PIPLUS

# ---------------------------------------------------------------------------
# Model parameters (Dutta et al. Table I, double-mediator)
# ---------------------------------------------------------------------------
M_CHI       = 8e-3    # GeV
M_CHI_PRIME = 50e-3   # GeV
M_V1        = 17e-3   # GeV  (light dark photon)
M_V2        = 200e-3  # GeV  (heavy upscattering mediator)
G_D         = 1.0
EPSILON_1   = 7e-5    # kinetic mixing for V1
EPSILON_2   = 1e-4    # kinetic mixing for V2
G_MU        = 1e-3    # Yukawa coupling of V1/phi to muon

M_PION = 0.13957039
M_MUON = 0.10565837

PDGID_PION      = 211
PDGID_MUPLUS    = -13
PDGID_NUMU      = 14
PDGID_V1_PROD   = 5922  # V1 from pion decay (decays to chi chi)
PDGID_V1_SIG    = 5923  # V1 from off-shell scattering (decays to e+e-)
PDGID_CHI       = 5917

events_to_inject = 10000

# ---------------------------------------------------------------------------
# 1. Read dk2nu files
# ---------------------------------------------------------------------------
dk2nu_dir = os.environ.get(
    "DK2NU_DIR",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), r"../../../../../../g4bnb/sources/G4BNB"),
)
dk2nu_files = sorted(glob.glob(os.path.join(dk2nu_dir, "*dk2nu*.root")))
if not dk2nu_files:
    dk2nu_files = sorted(glob.glob(os.path.join(dk2nu_dir, "nubeam*.root")))

if not dk2nu_files:
    print(f"No dk2nu files found in {dk2nu_dir}")
    sys.exit(1)

print(f"Reading {len(dk2nu_files)} dk2nu file(s) ...")
dk2nu_data = read_dk2nu(dk2nu_files, parent_pdg=[PTYPE_PIPLUS])
print_summary(dk2nu_data)

total_pot = dk2nu_data["pot"]
print(f"Total POT: {total_pot:.3e}")

# ---------------------------------------------------------------------------
# 2. Load SBND detector
# ---------------------------------------------------------------------------
print("\nLoading SBND detector model ...")
detector_model = utilities.load_detector("SBN", detector="SBND")

pion_type = siren.dataclasses.Particle.ParticleType(PDGID_PION)
v1_prod_type = siren.dataclasses.Particle.ParticleType(PDGID_V1_PROD)
v1_sig_type  = siren.dataclasses.Particle.ParticleType(PDGID_V1_SIG)
chi_type  = siren.dataclasses.Particle.ParticleType(PDGID_CHI)

# Fiducial volume: Box centered at detector origin
# SBND TPC active volume is roughly 4m x 4m x 5m
fiducial = Box(4.0, 4.0, 5.0)

# ---------------------------------------------------------------------------
# 3. Set up processes
# ---------------------------------------------------------------------------
print("Setting up processes ...")

# Physical pion decay (for weighting)
pion_decay_physical = MesonThreeBodySIRENDecay(
    m_meson=M_PION,
    m_lepton=M_MUON,
    m_mediator=M_V1,
    g_mu=G_MU,
    mediator_type="scalar",
    pdgid_meson=PDGID_PION,
    pdgid_lepton=PDGID_MUPLUS,
    pdgid_neutrino=PDGID_NUMU,
    pdgid_mediator=PDGID_V1_PROD,
)
print(f"  Pion 3-body width: {pion_decay_physical._total_width:.4e} GeV")

# Biased pion decay (for injection): V1 directed toward detector
# Detector origin is at (0,0,0) in detector coordinates
pion_decay_biased = BiasedMesonThreeBodyDecay(
    m_meson=M_PION,
    m_lepton=M_MUON,
    m_mediator=M_V1,
    m_chi=M_CHI,
    g_mu=G_MU,
    mediator_type="scalar",
    detector_position=(0.0, 0.0, 0.0),
    detector_radius=2.5,
    pdgid_meson=PDGID_PION,
    pdgid_lepton=PDGID_MUPLUS,
    pdgid_neutrino=PDGID_NUMU,
    pdgid_mediator=PDGID_V1_PROD,
)

# V1 from pion decay -> chi chi (biased toward detector)
v1_to_chi = BiasedDarkPhotonToChiDecay(
    M_V1, M_CHI, G_D,
    detector_position=(0.0, 0.0, 0.0),
    detector_radius=2.5,
    pdgid_V1=PDGID_V1_PROD,
    pdgid_chi=PDGID_CHI,
)
print(f"  V1->chi chi width: {v1_to_chi._total_width:.4e} GeV")

# V1 from off-shell scattering -> e+e- (signal)
v1_to_ee = DarkPhotonDecay(M_V1, EPSILON_1, pdgid_V1=PDGID_V1_SIG)
print(f"  V1->e+e- width: {v1_to_ee._total_width:.4e} GeV")

# chi + Ar -> chi + V1 + Ar  (off-shell chi', single vertex)
offshell_xs = VectorPortalOffShellXS(
    m_chi=M_CHI,
    m_chi_prime=M_CHI_PRIME,
    m_V1=M_V1,
    m_V2=M_V2,
    g_D=G_D,
    epsilon=EPSILON_2,
    pdgid_V1=PDGID_V1_SIG,
    nuclear_pdgid=1000180400,
    nuclear_mass=37.215,
    nuclear_name="Ar40",
    A=40, Z=18,
)
print(f"  chi scattering threshold: {offshell_xs._ups.Ethreshold:.4f} GeV")

primary_processes = {pion_type: [pion_decay_biased]}

secondary_processes = {
    v1_prod_type: [v1_to_chi],      # production V1 -> chi chi
    v1_sig_type: [v1_to_ee],         # signal V1 -> e+e-
    chi_type: [offshell_xs],          # chi + Ar -> chi + V1_sig + Ar
}

# ---------------------------------------------------------------------------
# 4. Distributions
# ---------------------------------------------------------------------------
primary_dist = dk2nu_to_primary_distribution(
    dk2nu_data, detector_model, parent_pdg=[PTYPE_PIPLUS]
)
print(f"  Loaded {primary_dist.GetPhysicalNumEvents()} pion events")

primary_injection_distributions = [primary_dist]

# Chi scattering biased into fiducial volume
chi_bounded = siren.distributions.SecondaryBoundedVertexDistribution(fiducial)

secondary_injection_distributions = {
    v1_prod_type: [siren.distributions.SecondaryBoundedVertexDistribution()],
    v1_sig_type: [siren.distributions.SecondaryBoundedVertexDistribution()],
    chi_type: [chi_bounded],
}

# ---------------------------------------------------------------------------
# 5. Stopping condition
# ---------------------------------------------------------------------------
def stop(tree, datum, i):
    secondary_type = datum.record.signature.secondary_types[i]
    parent_type = datum.record.signature.primary_type
    # Always stop: mu, nu, e (final state particles)
    if secondary_type in [
        siren.dataclasses.Particle.ParticleType.EMinus,
        siren.dataclasses.Particle.ParticleType.EPlus,
        siren.dataclasses.Particle.ParticleType.MuPlus,
        siren.dataclasses.Particle.ParticleType.MuMinus,
        siren.dataclasses.Particle.ParticleType(PDGID_NUMU),
    ]:
        return True
    # Stop the outgoing chi from the off-shell scattering
    # (parent is chi_type = the scattering interaction)
    if secondary_type == chi_type and parent_type == chi_type:
        return True
    # Stop the recoiling nucleus from off-shell scattering
    if secondary_type not in [chi_type, v1_prod_type, v1_sig_type]:
        if parent_type == chi_type:
            return True
    return False

# ---------------------------------------------------------------------------
# 6. Run injection
# ---------------------------------------------------------------------------
print(f"\nInjecting {events_to_inject} events ...")

injector = siren.injection.Injector()
injector.number_of_events = events_to_inject
injector.detector_model = detector_model
injector.primary_type = pion_type
injector.primary_interactions = primary_processes[pion_type]
injector.primary_injection_distributions = primary_injection_distributions
injector.secondary_interactions = secondary_processes
injector.secondary_injection_distributions = secondary_injection_distributions
injector.stopping_condition = stop

events, gen_times = GenerateEvents(injector)

n_total = len(events)
n_nonempty = sum(1 for e in events if len(e.tree) > 1)
total_time = sum(gen_times)

print(f"\nGenerated {n_total} events in {total_time:.1f} s")
print(f"  Non-trivial trees (>1 node): {n_nonempty}")
print(f"  Injection rate: {n_total / total_time:.0f} events/s")

# ---------------------------------------------------------------------------
# 7. Compute injection efficiency
# ---------------------------------------------------------------------------
print("\n--- Injection Efficiency ---")
print(f"Total pion decays sampled: {events_to_inject}")
print(f"Total events with interaction tree: {n_total}")

# Count how many events have a chi scattering vertex
n_chi_scatter = 0
n_v1_decay = 0
for event in events:
    for datum in event.tree:
        rec = datum.record
        primary_pdg = int(rec.signature.primary_type)
        if primary_pdg == PDGID_CHI:
            n_chi_scatter += 1
        if primary_pdg == PDGID_V1_SIG:
            n_v1_decay += 1

print(f"Events with chi scattering: {n_chi_scatter}")
print(f"Events with V1->ee signal: {n_v1_decay}")
print(f"Injection efficiency (chi scatter / injected): {n_chi_scatter / events_to_inject:.2e}")

if total_pot > 0:
    n_pions_in_file = len(dk2nu_data["E"])
    pot_per_event = total_pot / n_pions_in_file
    print(f"\n--- POT Scaling ---")
    print(f"POT per dk2nu pion: {pot_per_event:.3e}")
    print(f"POT sampled: {events_to_inject * pot_per_event:.3e}")
