"""
Vector-portal dark matter at MiniBooNE.

Physics:
    pi+ -> mu+ nu_mu V1  (three-body production)
    V1 -> chi chi'        (dark photon decay)
    chi N -> chi' N       (upscattering in detector)
    chi' -> chi V1        (de-excitation)
    V1 -> e+ e-           (visible signal)

Reference: Dutta et al., PRL 129, 111803 (2022) [arXiv:2110.11944]
Benchmark: double-mediator scenario (Table I).
"""

import os
import siren
from siren import utilities
from siren._util import get_processes_model_path
from siren import _util as _siren_util

_dt_base = os.path.join(
    _siren_util.resource_package_dir(), "processes", "DarkNewsTables",
)
_mod_processes = _siren_util.load_module(
    "siren.resources.processes.DarkNewsTables.processes",
    os.path.join(_dt_base, "processes.py"),
)

load_vector_portal = _mod_processes.load_vector_portal
PyDarkNewsCrossSection = _mod_processes.PyDarkNewsCrossSection
SaveDarkNewsProcesses = _mod_processes.SaveDarkNewsProcesses

# ---------------------------------------------------------------------------
# Model parameters (Dutta et al. Table I, double-mediator)
# ---------------------------------------------------------------------------
M_CHI       = 8e-3    # GeV  chi   (DM ground state)
M_CHI_PRIME = 50e-3   # GeV  chi'  (DM excited state)
M_V1        = 17e-3   # GeV  V1    (light dark photon)
M_V2        = 200e-3  # GeV  V2    (heavy upscattering mediator)
G_D         = 1.0
EPSILON_1   = 7e-5    # kinetic mixing for V1
EPSILON_2   = 1e-4    # kinetic mixing for V2

PDGID_CHI       = 5917
PDGID_CHI_PRIME = 5918
PDGID_V1        = 5922

# Production channel
M_MESON  = 0.13957   # GeV (pion)
M_LEPTON = 0.10566   # GeV (muon)
FLUX_TAG = "pion_numu"

events_to_inject = 10_000
experiment = "MiniBooNE"

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
detector_model = utilities.load_detector("SBN", detector=experiment)
chi_type       = siren.dataclasses.Particle.ParticleType(PDGID_CHI)
chi_prime_type = siren.dataclasses.Particle.ParticleType(PDGID_CHI_PRIME)
v1_type        = siren.dataclasses.Particle.ParticleType(PDGID_V1)

# ---------------------------------------------------------------------------
# Load processes via factory
# ---------------------------------------------------------------------------
primary_processes, secondary_processes, chi_flux = load_vector_portal(
    m_chi=M_CHI,
    m_chi_prime=M_CHI_PRIME,
    m_V1=M_V1,
    m_V2=M_V2,
    g_D=G_D,
    epsilon_1=EPSILON_1,
    epsilon_2=EPSILON_2,
    detector_model=detector_model,
    flux_tag=FLUX_TAG,
    m_meson=M_MESON,
    m_lepton=M_LEPTON,
    max_energy=3.0,
)

print(f"Primary cross sections: {len(primary_processes.get(chi_type, []))} targets")
print(f"Secondary decays: chi' ({len(secondary_processes.get(chi_prime_type, []))}), "
      f"V1 ({len(secondary_processes.get(v1_type, []))})")

# ---------------------------------------------------------------------------
# Injection distributions
# ---------------------------------------------------------------------------
# MiniBooNE fiducial volume: 5.0 m sphere about the tank center
# (detector-local origin). The SBN composite has no densities.dat, so
# the fiducial is built inline (cf. the SBND example).
fiducial_volume = siren.geometry.Sphere(5.0, 0.0)

primary_injection_distributions = [chi_flux]
primary_physical_distributions  = [chi_flux]

secondary_injection_distributions = {}
secondary_physical_distributions  = {}

def stop(datum, i):
    sec_type = cycler[i] if i < len(cycler) else None
    if sec_type == siren.dataclasses.Particle.ParticleType.EMinus:
        return True
    if sec_type == siren.dataclasses.Particle.ParticleType.EPlus:
        return True
    return i >= 3

# ---------------------------------------------------------------------------
# Injector
# ---------------------------------------------------------------------------
injector = siren.injection.Injector()
injector.number_of_events                  = events_to_inject
injector.detector_model                    = detector_model
injector.primary_type                      = chi_type
injector.primary_interactions              = primary_processes[chi_type]
injector.primary_injection_distributions   = primary_injection_distributions
injector.secondary_interactions            = secondary_processes
injector.secondary_injection_distributions = secondary_injection_distributions

print(f"Generating {events_to_inject} events ...")
events = injector.generate_events()
print(f"Done. Generated {len(events)} events.")
