"""
Input validation helpers for SIREN's Python interface.

Uses the C++ distribution base class hierarchy to classify
distributions by role.  Adding a new distribution in C++ that
inherits from the correct base class will be recognized here
automatically -- no manual plumbing required.
"""

from . import distributions as _d


def classify_distribution(dist):
    """Return a set of role labels for a distribution.

    Uses the C++ class hierarchy (``PrimaryEnergyDistribution``,
    ``PrimaryDirectionDistribution``, ``VertexPositionDistribution``,
    etc.) rather than enumerating concrete types.

    Returns a set that may contain multiple roles (e.g.
    ``{"energy", "direction"}`` for ``PrimaryEnergyDirectionDistribution``
    subclasses).
    """
    roles = set()
    if isinstance(dist, _d.PrimaryEnergyDistribution):
        roles.add("energy")
    if isinstance(dist, _d.PrimaryDirectionDistribution):
        roles.add("direction")
    if isinstance(dist, _d.PrimaryEnergyDirectionDistribution):
        roles.add("energy")
        roles.add("direction")
    if isinstance(dist, _d.VertexPositionDistribution):
        roles.add("position")
    if isinstance(dist, _d.PrimaryMass):
        roles.add("mass")
    if isinstance(dist, _d.PrimaryNeutrinoHelicityDistribution):
        roles.add("helicity")
    if isinstance(dist, _d.NormalizationConstant):
        roles.add("normalization")
    if isinstance(dist, _d.PrimaryAreaDistribution):
        roles.add("area")
    if not roles and isinstance(dist, _d.PrimaryInjectionDistribution):
        roles.add("unknown_primary")
    if not roles and isinstance(dist, _d.SecondaryInjectionDistribution):
        roles.add("unknown_secondary")
    return roles if roles else {"unknown"}


def collect_roles(distributions):
    """Return the union of all roles covered by a list of distributions."""
    roles = set()
    for d in distributions:
        roles |= classify_distribution(d)
    return roles


def validate_injection_distributions(distributions):
    """Check that a list of injection distributions covers all required roles.

    Raises
    ------
    ValueError
        If required distribution types are missing, with a message
        suggesting the relevant C++ base class.
    """
    roles = collect_roles(distributions)

    missing = []
    if "energy" not in roles:
        missing.append(
            "energy (any subclass of PrimaryEnergyDistribution "
            "or PrimaryEnergyDirectionDistribution)"
        )
    if "direction" not in roles:
        missing.append(
            "direction (any subclass of PrimaryDirectionDistribution "
            "or PrimaryEnergyDirectionDistribution)"
        )
    if "position" not in roles:
        missing.append(
            "position (any subclass of VertexPositionDistribution)"
        )

    if missing:
        raise ValueError(
            "Missing injection distributions:\n"
            + "\n".join(f"  - {m}" for m in missing)
        )


def validate_physical_distributions(distributions):
    """Check that a list of physical distributions covers required roles.

    Physical distributions need energy and direction but NOT position
    or mass.

    Raises
    ------
    ValueError
        If required distribution types are missing.
    """
    roles = collect_roles(distributions)

    missing = []
    if "energy" not in roles:
        missing.append(
            "energy (any subclass of PrimaryEnergyDistribution "
            "or PrimaryEnergyDirectionDistribution)"
        )
    if "direction" not in roles:
        missing.append(
            "direction (any subclass of PrimaryDirectionDistribution "
            "or PrimaryEnergyDirectionDistribution)"
        )

    if missing:
        raise ValueError(
            "Missing physical distributions:\n"
            + "\n".join(f"  - {m}" for m in missing)
        )
