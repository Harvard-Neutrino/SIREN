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
from . import math as _math

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

# ------------------------------------------------------------------ #
#  FixedDirection wrapper: accept list/tuple as well as Vector3D       #
#                                                                      #
#  NOTE: This convenience wrapper belongs at the pybind level          #
#  (accepting py::list/py::tuple in the C++ constructor).  It lives    #
#  here as a stopgap until that change is made.                        #
# ------------------------------------------------------------------ #

_OrigFixedDirection = _d.FixedDirection


def FixedDirection(direction):
    """Create a fixed-direction distribution.

    Parameters
    ----------
    direction : list, tuple, or Vector3D
        The direction vector. Lists/tuples of length 3 are auto-converted.
    """
    if isinstance(direction, (list, tuple)):
        direction = _math.Vector3D(*direction)
    return _OrigFixedDirection(direction)


FixedDirection.__wrapped__ = _OrigFixedDirection
