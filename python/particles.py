"""
Curated particle type aliases for convenient access.

Usage::

    import siren
    siren.particles.NuMu
    siren.particles.Electron

The full enum remains available at
``siren.dataclasses.Particle.ParticleType`` for less common types.
"""

from . import dataclasses as _dc

_PT = _dc.ParticleType


def _safe(name):
    """Get a ParticleType by name, returning None if it does not exist."""
    return getattr(_PT, name, None)


# Neutrinos
NuE = _PT.NuE
NuMu = _PT.NuMu
NuTau = _PT.NuTau
NuEBar = _PT.NuEBar
NuMuBar = _PT.NuMuBar
NuTauBar = _PT.NuTauBar

# Charged leptons
Electron = _PT.EMinus
Positron = _PT.EPlus
EMinus = _PT.EMinus
EPlus = _PT.EPlus
Muon = _PT.MuMinus
MuMinus = _PT.MuMinus
MuPlus = _PT.MuPlus
TauMinus = _PT.TauMinus
TauPlus = _PT.TauPlus

# Nucleons
Nucleon = _PT.Nucleon
Neutron = _PT.Neutron
PPlus = _PT.PPlus

# Photon
Gamma = _PT.Gamma

# Mesons
Pi0 = _PT.Pi0
PiPlus = _PT.PiPlus
PiMinus = _PT.PiMinus
KPlus = _PT.KPlus
KMinus = _PT.KMinus
K0_Long = _PT.K0_Long
K0_Short = _PT.K0_Short
Eta = _safe("Eta")
EtaPrime = _safe("EtaPrime")

# BSM -- these may not exist in all builds
N4 = _safe("N4")
N5 = _safe("N5")
N6 = _safe("N6")
NuF4 = _safe("NuF4")
NuF4Bar = _safe("NuF4Bar")
ALP = _safe("ALP")

# Build reverse lookup: string name -> ParticleType
# Includes both our curated aliases and all enum members
_name_to_type = {}

# All enum members first (canonical names)
for _name in dir(_PT):
    if _name.startswith("_"):
        continue
    _val = getattr(_PT, _name)
    if isinstance(_val, type(_PT.NuMu)):
        _name_to_type[_name] = _val

# Add our convenience aliases
_name_to_type["Electron"] = _PT.EMinus
_name_to_type["Positron"] = _PT.EPlus
_name_to_type["Muon"] = _PT.MuMinus


def resolve(name_or_type):
    """Resolve a string name or ParticleType enum to a ParticleType.

    Parameters
    ----------
    name_or_type : str or ParticleType
        Either a string like ``"NuMu"`` or a ``ParticleType`` enum value.

    Returns
    -------
    ParticleType

    Raises
    ------
    ValueError
        If the string does not match any known particle type.
    TypeError
        If the argument is neither a string nor a ParticleType.
    """
    if isinstance(name_or_type, type(_PT.NuMu)):
        return name_or_type
    if isinstance(name_or_type, str):
        if name_or_type in _name_to_type:
            return _name_to_type[name_or_type]
        raise ValueError(
            f"Unknown particle type: {name_or_type!r}. "
            f"Available: {', '.join(sorted(_name_to_type.keys()))}"
        )
    raise TypeError(
        f"Expected a ParticleType enum or string name, got {type(name_or_type).__name__}"
    )
