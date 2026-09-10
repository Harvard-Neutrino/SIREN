"""
Input validation helpers for SIREN's Python interface.

Uses the C++ distribution methods ``SetVariables()``,
``RequiredVariables()``, and ``DensityVariables()`` to validate
distribution lists.  No hardcoded type lists -- adding a new
distribution in C++ with the correct overrides makes it work here
automatically.
"""

from collections import Counter

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

# Primary position measures need a little more structure than the opaque
# strings returned by DensityVariables().  A vertex is three-dimensional: two
# coordinates in the plane normal to the primary direction and one coordinate
# along that direction.  A primary-area density is expressed in the two
# point-of-closest-approach coordinates in that same transverse plane.
_PRIMARY_VERTEX_DENSITY = "InteractionVertexPosition"
_PRIMARY_AREA_DENSITY = "PointOfClosestApproach"
_PRIMARY_TRANSVERSE_DENSITIES = (
    "PrimaryPositionTransverse[0]",
    "PrimaryPositionTransverse[1]",
)
_PRIMARY_LONGITUDINAL_DENSITY = "PrimaryPositionLongitudinal"


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
            newly_set = dist.SetVariables()
            overlap = newly_set & available
            if overlap:
                raise ValueError(
                    f"Distribution {type(dist).__name__} (index {i}) sets "
                    f"variables {sorted(v.name for v in overlap)} that have "
                    f"already been set by preceding distributions. "
                    f"Available so far: {sorted(v.name for v in available)}"
                )
            available |= newly_set


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


def _effective_density_variables(distributions):
    """Return differential components, preserving repeated factors."""
    result = Counter()
    for distribution in distributions:
        for variable in distribution.DensityVariables():
            if variable == _PRIMARY_VERTEX_DENSITY:
                result.update(_PRIMARY_TRANSVERSE_DENSITIES)
                result[_PRIMARY_LONGITUDINAL_DENSITY] += 1
            elif variable == _PRIMARY_AREA_DENSITY:
                result.update(_PRIMARY_TRANSVERSE_DENSITIES)
            else:
                result[variable] += 1
    return result


def _format_density_variables(variables):
    return sorted(
        name if count == 1 else f"{name} (x{count})"
        for name, count in variables.items()
        if count > 0
    )


def validate_reweighting_compatibility(
    injection_distributions,
    physical_distributions,
    *,
    compute_interaction_probability=True,
    compute_position_probability=True,
):
    """Check that physical distributions can reweight injection distributions.

    Effective physical density variables must be a subset of the injection
    density variables.  Variables that are delta functions in both injection
    and physical are OK.

    Primary position measures are compared in transverse/longitudinal
    coordinates.  A ``VertexPositionDistribution`` supplies all three.  A
    ``PrimaryAreaDistribution`` supplies the two transverse coordinates, and
    the downstream normalized-position factor supplies the longitudinal
    coordinate when ``compute_position_probability`` is true.  Consequently,
    ordinary propagated weighting may leave two injection coordinates for an
    external flux-per-area factor, while a fixed-target area distribution plus
    normalized position density fully matches a fixed-target vertex
    distribution.  The interaction-probability factor is dimensionless, so
    ``compute_interaction_probability`` does not add a density variable.

    Parameters
    ----------
    injection_distributions : list
        Injection distributions.
    physical_distributions : list
        Physical distributions.
    compute_interaction_probability : bool, optional
        Whether downstream weighting includes the dimensionless interaction
        probability.  Accepted as part of the complete weighting-mode context;
        it does not change differential-variable compatibility.
    compute_position_probability : bool, optional
        Whether downstream weighting includes normalized position density.

    Raises
    ------
    ValueError
        If physical distributions are differential in variables that the
        injection distributions don't cover.
    """
    inj_density = _effective_density_variables(injection_distributions)
    phys_density = _effective_density_variables(physical_distributions)

    if compute_position_probability:
        phys_density[_PRIMARY_LONGITUDINAL_DENSITY] += 1

    extra = phys_density - inj_density
    if extra:
        raise ValueError(
            "Physical weighting is differential in variables "
            f"{_format_density_variables(extra)} that are not covered by "
            "injection distributions. "
            "Effective injection density variables: "
            f"{_format_density_variables(inj_density)}, effective physical "
            "density variables (including downstream weighting factors): "
            f"{_format_density_variables(phys_density)}"
        )
