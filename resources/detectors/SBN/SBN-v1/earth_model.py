"""
PREM Earth model for the SBN detector at Fermilab.

Adds concentric PREM spherical shells and an exponential atmosphere to a
DetectorModel that was loaded via LoadGDML. The GDML volumes always take
priority over these Earth sectors wherever both exist.

The PREM layers and the atmosphere model are the shared ones defined in
``siren.earth``; this module only supplies the SBN-specific pieces: the BNB
"up" direction, the ground level, and the local crustal-thickening correction.

The PREM layers follow Dziewonski & Anderson, PEPI 25 (1981) 297, with the
standard global-average Moho at 24.4 km depth. A local crustal thickening
correction overrides the upper mantle with crust (ROCK) in the Midwest region
from Fermilab to the DUNE far detector at Sanford Lab, reflecting the ~45 km Moho
depth measured by CRUST1.0 (Laske et al. 2013) for this stable continental
interior. The correction is an Earth-centered spherical shell spanning depths
24.4-45 km with a polar-angle cut (~17 deg around Fermilab), covering the entire
Fermilab-to-Sanford-Lab baseline (1285 km) with margin.

Earth center in BNB coordinates:
  The BNB frame has y = up, origin at the BNB target. Fermilab grade
  is _GRADE_Y_BNB = 7.62 m above the target. We equate grade with
  the PREM surface radius (6371 km), so the Earth center sits at
  (0, -(R_PREM - _GRADE_Y_BNB), 0) in BNB coordinates.

References:
  PREM: Dziewonski & Anderson, PEPI 25 (1981) 297
  Atmosphere: US Standard Atmosphere 1976
  CRUST1.0: Laske et al. 2013 (1x1 degree global crustal model)
  DUNE beam path: Roe 2016, arXiv:1608.03802
  Fermilab geology: FNAL ESH regional site characterization
"""

from __future__ import annotations

import math
from typing import Any

# Canonical shared Earth model (pure-data imports; the sector builder imports
# the compiled siren.detector module lazily, so importing these constants stays
# cheap).
from siren.earth import (
    R_PREM,
    CANONICAL_PREM_LAYERS as _PREM_LAYERS,
    USStandard1976,
    AtmosphereForm,
    DEFAULT_SHELL_ALTS as _ATMO_SHELL_ALTS,
    constant_shell_layers,
    build_earth_sectors,
    PositionedSector,
)

_GRADE_Y_BNB = 7.62
_PREM_MOHO_DEPTH = 24400.0
_LOCAL_MOHO_DEPTH = 45000.0

# Angular half-opening of the local crustal thickening correction, measured from
# the radial direction through Fermilab. 0.3 rad ~ 1900 km on the surface,
# comfortably covering the 1285 km Fermilab-to-Sanford-Lab baseline.
_LOCAL_MOHO_THETA = 0.3

# Shared atmosphere model (US Standard Atmosphere 1976, isothermal).
_ATMOSPHERE = USStandard1976()

# Constant-shell reference description of the atmosphere. The sectors are
# built as true exponential-density shells; these constants give the
# equivalent mass-conserving constant-shell approximation.
_RHO_0 = _ATMOSPHERE.rho0
_SCALE_HEIGHT = _ATMOSPHERE.scale_height


def _atmo_avg_density(h1: float, h2: float) -> float:
    return _ATMOSPHERE.shell_mean_density(h1, h2)


_ATMO_LAYERS = constant_shell_layers(_ATMOSPHERE, R_PREM, _ATMO_SHELL_ALTS)


def _all_layers():
    """All layers from innermost to outermost."""
    return list(_PREM_LAYERS) + list(_ATMO_LAYERS)


def add_earth_model(model: Any) -> None:
    """Add PREM + atmosphere sectors to a DetectorModel loaded from GDML.

    Must be called after model.LoadGDML() so that the GDML volumes already have
    positive levels. The PREM/atmosphere sectors are shifted below them (the GDML
    geometry always wins wherever it reaches).

    A local crustal thickening sector is added between the PREM upper mantle and
    crustal layers: an Earth-centered spherical shell (depths 24.4-45 km) with a
    polar-angle cut (~17 deg around Fermilab) that overrides the PREM mantle with
    ROCK in the Midwest, where the continental Moho is ~45 km deep (CRUST1.0).
    """
    from siren.math import Vector3D, Quaternion
    from siren.geometry import Sphere, Placement

    ground_level = _GRADE_Y_BNB

    # Earth center along the BNB +y (radially-outward through Fermilab) axis.
    earth_center = Vector3D(0.0, -(R_PREM - ground_level), 0.0)

    # Rotate so the Sphere's +z axis (theta=0) aligns with BNB +y.
    rot = Quaternion.rotation_between(Vector3D(0, 0, 1), Vector3D(0, 1, 0))
    rotated_placement = Placement(earth_center, rot)

    local_moho = PositionedSector(
        name="local_thick_crust",
        material="ROCK",
        density_const=2.9,
        geo=Sphere(
            rotated_placement,
            R_PREM - _PREM_MOHO_DEPTH,   # outer: PREM Moho depth
            R_PREM - _LOCAL_MOHO_DEPTH,  # inner: local Moho depth
            0.0, 2.0 * math.pi,          # full azimuth
            0.0, _LOCAL_MOHO_THETA,      # polar cap around Fermilab
        ),
        insert_before="moho_boundary",
    )

    build_earth_sectors(
        model,
        surface_offset=ground_level,
        prem_layers=_PREM_LAYERS,
        atmosphere=AtmosphereForm(model=_ATMOSPHERE, mode="exponential"),
        extra_sectors=[local_moho],
        up_axis=(0.0, 1.0, 0.0),
        r_prem=R_PREM,
    )
