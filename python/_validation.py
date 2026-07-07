"""
Input validation helpers for SIREN's Python interface.

Provides early, human-readable error messages instead of letting
misconfiguration propagate into the C++ layer.
"""

from . import distributions as _d

# Categories of concrete distribution types
_ENERGY_TYPES = (
    _d.PowerLaw,
    _d.Monoenergetic,
    _d.TabulatedFluxDistribution,
)

_DIRECTION_TYPES = (
    _d.IsotropicDirection,
    _d.FixedDirection,
    _d.Cone,
)

_POSITION_TYPES = (
    _d.ColumnDepthPositionDistribution,
    _d.CylinderVolumePositionDistribution,
    _d.PointSourcePositionDistribution,
    _d.RangePositionDistribution,
    _d.VertexPositionDistribution,
)

_MASS_TYPES = (
    _d.PrimaryMass,
)


def classify_distribution(dist):
    """Return a string label for what role a distribution fills.

    Returns one of: ``"energy"``, ``"direction"``, ``"position"``,
    ``"mass"``, ``"helicity"``, or ``"unknown"``.
    """
    if isinstance(dist, _ENERGY_TYPES):
        return "energy"
    if isinstance(dist, _DIRECTION_TYPES):
        return "direction"
    if isinstance(dist, _POSITION_TYPES):
        return "position"
    if isinstance(dist, _MASS_TYPES):
        return "mass"
    if isinstance(dist, _d.PrimaryNeutrinoHelicityDistribution):
        return "helicity"
    if isinstance(dist, _d.NormalizationConstant):
        return "normalization"
    return "unknown"


def validate_injection_distributions(distributions):
    """Check that a list of injection distributions covers all required roles.

    Parameters
    ----------
    distributions : list
        List of PrimaryInjectionDistribution objects.

    Raises
    ------
    ValueError
        If required distribution types are missing, with a message listing
        what is needed and suggesting concrete classes.
    """
    roles = {classify_distribution(d) for d in distributions}

    missing = []
    if "energy" not in roles:
        missing.append(
            "energy (e.g. PowerLaw, Monoenergetic, TabulatedFluxDistribution)"
        )
    if "direction" not in roles:
        missing.append(
            "direction (e.g. IsotropicDirection, FixedDirection, Cone)"
        )
    if "position" not in roles:
        missing.append(
            "position (e.g. ColumnDepthPositionDistribution, "
            "CylinderVolumePositionDistribution, PointSourcePositionDistribution)"
        )

    if missing:
        raise ValueError(
            "Missing injection distributions:\n"
            + "\n".join(f"  - {m}" for m in missing)
        )


def validate_physical_distributions(distributions):
    """Check that a list of physical distributions covers required roles.

    Physical distributions need energy and direction but NOT position or mass.

    Parameters
    ----------
    distributions : list
        List of WeightableDistribution objects.

    Raises
    ------
    ValueError
        If required distribution types are missing.
    """
    roles = {classify_distribution(d) for d in distributions}

    missing = []
    if "energy" not in roles:
        missing.append(
            "energy (e.g. PowerLaw, Monoenergetic, TabulatedFluxDistribution)"
        )
    if "direction" not in roles:
        missing.append(
            "direction (e.g. IsotropicDirection, FixedDirection, Cone)"
        )

    if missing:
        raise ValueError(
            "Missing physical distributions:\n"
            + "\n".join(f"  - {m}" for m in missing)
        )
