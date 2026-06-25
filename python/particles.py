"""
Particle type aliases auto-generated from the C++ ``ParticleType`` enum.

Usage::

    import siren
    siren.particles.NuMu
    siren.particles.Electron   # convenience alias for EMinus

All members of the C++ enum are available as module-level attributes.
Human-friendly aliases (``Electron`` for ``EMinus``, etc.) are added
on top.  The full enum remains available at
``siren.dataclasses.ParticleType``.
"""

from . import dataclasses as _dc

_PT = _dc.ParticleType

# ------------------------------------------------------------------ #
#  Auto-export every enum member as a module-level attribute           #
# ------------------------------------------------------------------ #

_name_to_type = {}
_sentinel = _PT.NuMu  # known member to identify the enum type

for _name in dir(_PT):
    if _name.startswith("_"):
        continue
    _val = getattr(_PT, _name)
    if isinstance(_val, type(_sentinel)):
        globals()[_name] = _val
        _name_to_type[_name] = _val

# ------------------------------------------------------------------ #
#  Human-friendly aliases                                              #
# ------------------------------------------------------------------ #

_aliases = {
    "Electron": "EMinus",
    "Positron": "EPlus",
    "Muon": "MuMinus",
    "Proton": "PPlus",
}

for _alias, _canonical in _aliases.items():
    if _canonical in _name_to_type:
        globals()[_alias] = _name_to_type[_canonical]
        _name_to_type[_alias] = _name_to_type[_canonical]

# Clean up loop variables from module namespace
del _name, _val, _alias, _canonical, _sentinel


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
