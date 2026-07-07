"""
Distribution shorthand aliases.

Usage::

    import siren
    siren.dist.PowerLaw(2, 1e3, 1e6)
    siren.dist.IsotropicDirection()
    siren.dist.ColumnDepth(600, 600.0)

The full names remain available at ``siren.distributions.*``.
"""

from . import distributions as _d
from . import math as _math
from . import geometry as _geometry

# ---- Energy distributions ----
PowerLaw = _d.PowerLaw
Monoenergetic = _d.Monoenergetic
TabulatedFlux = _d.TabulatedFluxDistribution

# ---- Direction distributions ----
IsotropicDirection = _d.IsotropicDirection
Cone = _d.Cone
NormalizationConstant = _d.NormalizationConstant

# FixedDirection wrapper: accept list/tuple as well as Vector3D
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

# ---- Position distributions ----
ColumnDepth = _d.ColumnDepthPositionDistribution
CylinderVolume = _d.CylinderVolumePositionDistribution
PointSource = _d.PointSourcePositionDistribution
RangePosition = _d.RangePositionDistribution

# ---- Depth / range functions ----
LeptonDepthFunction = _d.LeptonDepthFunction
DecayRange = _d.DecayRangeFunction
DepthFunction = _d.DepthFunction
RangeFunction = _d.RangeFunction

# ---- Secondary distributions ----
BoundedVertex = _d.SecondaryBoundedVertexDistribution
SecondaryVertex = _d.SecondaryVertexPositionDistribution
PhysicalVertex = _d.SecondaryPhysicalVertexDistribution

# ---- Mass / helicity ----
Mass = _d.PrimaryMass
Helicity = _d.PrimaryNeutrinoHelicityDistribution

# ---- Base classes (for isinstance checks) ----
PrimaryInjectionDistribution = _d.PrimaryInjectionDistribution
SecondaryInjectionDistribution = _d.SecondaryInjectionDistribution
WeightableDistribution = _d.WeightableDistribution
PhysicallyNormalizedDistribution = _d.PhysicallyNormalizedDistribution
