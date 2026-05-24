"""
Input validation helpers for SIREN's Python interface.

Uses the C++ distribution methods ``SetVariables()``,
``RequiredVariables()``, and ``DensityVariables()`` to validate
distribution lists.  No hardcoded type lists -- adding a new
distribution in C++ with the correct overrides makes it work here
automatically.
"""

from . import distributions as _d

DV = _d.DistributionVariable

_REQUIRED_INJECTION_VARIABLES = {
    DV.PrimaryEnergy,
    DV.PrimaryDirection,
}

_REQUIRED_INJECTION_VERTEX = {
    DV.InitialPosition,
    DV.InteractionVertex,
}

_REQUIRED_PHYSICAL_VARIABLES = {
    DV.PrimaryEnergy,
    DV.PrimaryDirection,
}


def collect_set_variables(distributions):
    """Return the union of all SetVariables across a list of distributions."""
    result = set()
    for d in distributions:
        if hasattr(d, "SetVariables"):
            result |= d.SetVariables()
    return result


def collect_density_variables(distributions):
    """Return the union of all DensityVariables across a list of distributions."""
    result = set()
    for d in distributions:
        result.update(d.DensityVariables())
    return result


def validate_ordering(distributions):
    """Check that each distribution's RequiredVariables are satisfied
    by the SetVariables of preceding distributions.

    Parameters
    ----------
    distributions : list
        Ordered list of PrimaryInjectionDistribution objects.

    Raises
    ------
    ValueError
        If a distribution requires variables that haven't been set yet,
        with a message naming the distribution and the missing variables.
    """
    available = set()
    for i, dist in enumerate(distributions):
        if not hasattr(dist, "RequiredVariables"):
            continue
        required = dist.RequiredVariables()
        missing = required - available
        if missing:
            missing_names = sorted(v.name for v in missing)
            raise ValueError(
                f"Distribution {type(dist).__name__} (index {i}) requires "
                f"variables {missing_names} to be set first, but they are "
                f"not covered by preceding distributions. "
                f"Available so far: {sorted(v.name for v in available)}"
            )
        if hasattr(dist, "SetVariables"):
            available |= dist.SetVariables()


def validate_injection_distributions(distributions):
    """Check that injection distributions are complete and correctly ordered.

    Validates:
    1. All required variable roles are covered (energy, direction, vertex).
    2. Ordering constraints are satisfied (RequiredVariables).

    Parameters
    ----------
    distributions : list
        Ordered list of PrimaryInjectionDistribution objects.

    Raises
    ------
    ValueError
        If required distributions are missing or ordering is wrong.
    """
    covered = collect_set_variables(distributions)

    missing = []
    for var in _REQUIRED_INJECTION_VARIABLES:
        if var not in covered:
            missing.append(var.name)

    has_vertex = bool(covered & _REQUIRED_INJECTION_VERTEX)
    if not has_vertex:
        missing.append("InitialPosition/InteractionVertex "
                       "(any VertexPositionDistribution subclass)")

    if missing:
        raise ValueError(
            "Missing injection distributions for:\n"
            + "\n".join(f"  - {m}" for m in missing)
        )

    validate_ordering(distributions)


def validate_physical_distributions(distributions):
    """Check that physical distributions cover required roles.

    Physical distributions need energy and direction but NOT position
    or mass.

    Parameters
    ----------
    distributions : list
        List of WeightableDistribution objects.

    Raises
    ------
    ValueError
        If required distribution types are missing.
    """
    covered = collect_set_variables(distributions)

    missing = []
    for var in _REQUIRED_PHYSICAL_VARIABLES:
        if var not in covered:
            missing.append(var.name)

    if missing:
        raise ValueError(
            "Missing physical distributions for:\n"
            + "\n".join(f"  - {m}" for m in missing)
        )


def validate_reweighting_compatibility(injection_distributions, physical_distributions):
    """Check that physical distributions can reweight injection distributions.

    Physical DensityVariables must be a subset of injection DensityVariables.
    Variables that are delta functions in both injection and physical are OK.

    Parameters
    ----------
    injection_distributions : list
        Injection distributions.
    physical_distributions : list
        Physical distributions.

    Raises
    ------
    ValueError
        If physical distributions are differential in variables that the
        injection distributions don't cover.
    """
    inj_density = collect_density_variables(injection_distributions)
    phys_density = collect_density_variables(physical_distributions)

    extra = phys_density - inj_density
    if extra:
        raise ValueError(
            f"Physical distributions are differential in variables "
            f"{sorted(extra)} that are not covered by injection distributions. "
            f"Injection density variables: {sorted(inj_density)}, "
            f"Physical density variables: {sorted(phys_density)}"
        )
