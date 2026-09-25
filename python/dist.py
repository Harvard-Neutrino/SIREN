"""
Distribution shorthand aliases auto-generated from C++ bindings.

All concrete distribution classes from ``siren.distributions`` are
available here under their original names.  A small set of shorter
aliases is provided for the most common ones.

Usage::

    import siren
    siren.dist.PowerLaw(2, 1e3, 1e6)
    siren.dist.IsotropicDirection()
    siren.dist.ColumnDepth(600, 600.0)

Adding a new distribution in C++ and binding it to Python will make
it available here automatically.
"""

from . import distributions as _d

# ------------------------------------------------------------------ #
#  Auto-export every class from siren.distributions                    #
# ------------------------------------------------------------------ #

# Base classes that we expose but are not "concrete" distributions
_BASE_CLASSES = {
    "PrimaryInjectionDistribution",
    "SecondaryInjectionDistribution",
    "WeightableDistribution",
    "PrimaryEnergyDistribution",
    "PrimaryDirectionDistribution",
    "PrimaryEnergyDirectionDistribution",
    "PrimaryAreaDistribution",
    "VertexPositionDistribution",
    "SecondaryVertexPositionDistribution",
    "PhysicallyNormalizedDistribution",
    "DepthFunction",
    "RangeFunction",
}

for _name in dir(_d):
    if _name.startswith("_"):
        continue
    _obj = getattr(_d, _name)
    if isinstance(_obj, type):
        globals()[_name] = _obj

# Clean up loop variables
del _name, _obj

# ------------------------------------------------------------------ #
#  Short aliases for common distributions                              #
# ------------------------------------------------------------------ #

_short_aliases = {
    "ColumnDepth": "ColumnDepthPositionDistribution",
    "CylinderVolume": "CylinderVolumePositionDistribution",
    "SphereVolume": "SphereVolumePositionDistribution",
    "PointSource": "PointSourcePositionDistribution",
    "RangePosition": "RangePositionDistribution",
    "DecayRangePosition": "DecayRangePositionDistribution",
    "DecayRangeVertex": "SecondaryDecayRangePositionDistribution",
    "BoundedVertex": "SecondaryBoundedVertexDistribution",
    "SecondaryVertex": "SecondaryVertexPositionDistribution",
    "PhysicalVertex": "SecondaryPhysicalVertexDistribution",
    "Mass": "PrimaryMass",
    "Helicity": "PrimaryNeutrinoHelicityDistribution",
    "DecayRange": "DecayRangeFunction",
    "TabulatedFlux": "TabulatedFluxDistribution",
    "Tabulated2DFlux": "Tabulated2DFluxDistribution",
    "PiDARNuE": "PiDARNuEDistribution",
    "FixedTargetPosition": "FixedTargetPositionDistribution",
    "FixedTargetArea": "FixedTargetAreaDistribution",
    "BoundedPrimaryVertex": "PrimaryBoundedVertexDistribution",
    "PhysicalPrimaryVertex": "PrimaryPhysicalVertexDistribution",
    "External": "PrimaryExternalDistribution",
}

for _alias, _canonical in _short_aliases.items():
    if hasattr(_d, _canonical):
        globals()[_alias] = getattr(_d, _canonical)

del _alias, _canonical


# FixedDirection, Cone, and PointSourcePositionDistribution accept
# list/tuple natively via pybind11 overloads (std::array<double,3>
# constructors added in distributions.cxx).  No Python wrapper needed.


# Columns PrimaryExternalDistribution interprets itself (see
# PrimaryExternalDistribution.cxx); a metadata column must use another name.
_EXTERNAL_RESERVED_COLUMNS = frozenset(
    ("E", "m", "px", "py", "pz", "x", "y", "z", "x0", "y0", "z0", "t0", "weight"))


def on_shell_rows(keys, rows, *, tolerance=2e-5, keep_input_energy="E_table"):
    """Return external-table rows whose energy is rebuilt from mass and momentum.

    An opt-in policy for tabulated decaying parents. Rounding a row (for example
    to float32) moves E^2 - |p|^2 far from m^2 at high boost: a float32 pi0 row
    at gamma 1000 is about 6% off in invariant mass. SIREN's decay channels use
    the declared mass and three-momentum; canonicalizing the table makes every
    other consumer (the source's own density and ``emin`` cut, physical models
    that read the recorded energy) see that same parent.

    The ``m``, ``px``, ``py`` and ``pz`` columns are authoritative and ``E``
    becomes ``hypot(m, |p|)``. A row whose input energy differs from that by
    more than ``tolerance`` (relative), or whose rebuilt energy is not finite,
    raises ValueError instead of being moved, because that is not rounding.
    Each of these columns must appear exactly once. Rows keep their order, weights and all other
    columns. The input energy is kept in the column ``keep_input_energy``
    (``None`` drops it), which must be an ordinary metadata name: the
    distribution stores such a column in every record's interaction parameters,
    whereas a name it interprets (``weight``, ``t0``, a kinematic column or its
    own ``PrimaryExternalDistribution_*`` bookkeeping) would change the
    simulation. Apply this before constructing
    ``PrimaryExternalDistribution`` so that ``emin`` and the sampling weights
    use the canonical energies. Returns ``(keys, rows)``.
    """
    import math

    keys = list(keys)
    missing = [k for k in ("E", "m", "px", "py", "pz") if k not in keys]
    if missing:
        raise ValueError("on_shell_rows needs columns E, m, px, py, pz; missing %s" % missing)
    # The distribution reads the last of repeated columns; this helper would
    # rewrite the first. Refuse the ambiguity.
    repeated = [k for k in ("E", "m", "px", "py", "pz") if keys.count(k) > 1]
    if repeated:
        raise ValueError("columns %s appear more than once" % repeated)
    if keep_input_energy is not None:
        if not isinstance(keep_input_energy, str) or not keep_input_energy:
            raise ValueError("keep_input_energy must be a nonempty column name or None")
        if (keep_input_energy in _EXTERNAL_RESERVED_COLUMNS
                or keep_input_energy.startswith("PrimaryExternalDistribution_")):
            raise ValueError("%r is a column PrimaryExternalDistribution interprets; "
                             "use an ordinary metadata name" % keep_input_energy)
        if keep_input_energy in keys:
            raise ValueError("column %r already exists" % keep_input_energy)
    if not (tolerance >= 0 and math.isfinite(tolerance)):
        raise ValueError("tolerance must be finite and nonnegative")
    ie, im, ix, iy, iz = (keys.index(k) for k in ("E", "m", "px", "py", "pz"))
    out, bad = [], []
    for index, row in enumerate(rows):
        row = [float(value) for value in row]
        if len(row) != len(keys):
            # An appended input-energy column would land in the wrong place.
            raise ValueError("row %d has %d values for %d columns" % (index, len(row), len(keys)))
        energy, mass, px, py, pz = row[ie], row[im], row[ix], row[iy], row[iz]
        if not all(math.isfinite(v) for v in (energy, mass, px, py, pz)) or mass < 0:
            bad.append(index)
            continue
        # Same order of operations as SIREN's native |p| and hypot(m, |p|).
        on_shell = math.hypot(mass, math.sqrt(px * px + py * py + pz * pz))
        if not math.isfinite(on_shell) or abs(energy - on_shell) > tolerance * on_shell:
            bad.append(index)
            continue
        row[ie] = on_shell
        if keep_input_energy is not None:
            row.append(energy)
        out.append(row)
    if bad:
        raise ValueError(
            "%d row(s) are not within %g of their mass shell or have invalid "
            "values (first indices: %s)" % (len(bad), tolerance, bad[:10]))
    if keep_input_energy is not None:
        keys.append(keep_input_energy)
    return keys, out
