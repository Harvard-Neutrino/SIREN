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

from .errors import ConfigurationError

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


# ------------------------------------------------------------------ #
#  User-defined (BSM) particle registration                            #
# ------------------------------------------------------------------ #

# name -> pdg code and name -> mass (GeV) for particles added via define().
_name_to_pdg = {}
_name_to_mass = {}


def define(name, pdg, mass):
    """Register a named particle so ``resolve(name)`` returns a ParticleType.

    The C++ ``ParticleType`` enum is PDG-coded (its members are the PDG
    integers themselves, e.g. ``N4 = 5914``), so a new BSM particle is
    represented by constructing ``ParticleType(pdg)`` from its PDG code --
    no new enum member is created.  ``mass`` (GeV) is recorded as metadata
    alongside the pdg code, for consumers that need it (e.g. phase-space
    templates), but the C++ enum itself carries no mass.

    Registering ``name``/``pdg`` with the values already on file (the exact
    same pdg and mass) is a no-op.  Registering ``name`` or ``pdg`` with a
    DIFFERENT pdg/mass than already registered raises ConfigurationError --
    this catches accidental redefinition rather than silently rebinding a
    name or code already in use.

    Parameters
    ----------
    name : str
        The identifier to register (becomes ``siren.particles.<name>`` and
        resolvable via ``resolve(name)``).
    pdg : int
        The PDG code identifying this particle's ParticleType.
    mass : float
        Mass in GeV, recorded in ``_name_to_mass``.

    Returns
    -------
    ParticleType
        The resolved (constructed-from-pdg) ParticleType, also stored at
        ``siren.particles.<name>``.
    """
    pdg = int(pdg)
    mass = float(mass)
    ptype = _PT(pdg)

    if name in _name_to_type:
        existing = _name_to_type[name]
        existing_pdg = _name_to_pdg.get(name, int(existing))
        if name in _name_to_mass:
            # A define()-registered name carries pdg and mass: both must match
            # for the re-registration to be an idempotent no-op.
            same = existing_pdg == pdg and _name_to_mass[name] == mass
        else:
            # A built-in enum member carries no mass; matching its pdg (the
            # bare enum value) is an idempotent no-op.
            same = existing_pdg == pdg
        if not same:
            raise ConfigurationError(
                f"particle name {name!r} is already registered with "
                f"pdg={existing_pdg}, mass={_name_to_mass.get(name)} "
                f"(requested pdg={pdg}, mass={mass})")
        return existing

    # A pdg code already in use -- under a different define()d name or a
    # built-in enum member -- is a conflict.
    for other_name, other_type in _name_to_type.items():
        other_pdg = _name_to_pdg.get(other_name, int(other_type))
        if other_pdg == pdg and other_name != name:
            raise ConfigurationError(
                f"pdg={pdg} is already registered under name {other_name!r} "
                f"(requested name {name!r})")

    _name_to_type[name] = ptype
    _name_to_pdg[name] = pdg
    _name_to_mass[name] = mass
    globals()[name] = ptype
    return ptype
