"""
PREM Earth model for the SBN detector at Fermilab.

Adds concentric PREM spherical shells and a multi-shell atmosphere to
a DetectorModel that was loaded via LoadGDML. The GDML volumes (which
carry positive levels) always take priority over these Earth sectors
(which carry negative levels) wherever both exist.

The PREM layers follow Dziewonski & Anderson, PEPI 25 (1981) 297,
with an adjusted Moho depth of 50 km for the Illinois/Midwest crust
(vs. the PREM global average of 24.4 km).

The atmosphere approximates the US Standard Atmosphere 1976 as 6
constant-density spherical shells.

Earth center in BNB coordinates:
  The BNB frame has y = up, origin at the BNB target. Fermilab grade
  is _GRADE_Y_BNB = 7.62 m above the target. We equate grade with
  the PREM surface radius (6371 km), so the Earth center sits at
  (0, -(R_PREM - _GRADE_Y_BNB), 0) in BNB coordinates.

References:
  PREM: Dziewonski & Anderson, PEPI 25 (1981) 297
  Atmosphere: US Standard Atmosphere 1976
  Illinois Moho: ~50 km from receiver-function analysis (USArray)
  Fermilab geology: FNAL ESH regional site characterization
"""

from __future__ import annotations

import math
import os
from typing import Any

_THIS_DIR = os.path.dirname(os.path.realpath(__file__))

R_PREM = 6371000.0
_GRADE_Y_BNB = 7.62
_MOHO_DEPTH = 50000.0


# PREM Earth layers, ordered from innermost to outermost.
# Each entry: (name, outer_radius_m, material, density_type, density_params)
#   density_type: "constant" or "polynomial"
#   density_params: single float for constant, list of coefficients for polynomial
#     polynomial: rho(r) = p[0] + p[1]*r + p[2]*r^2 + ...  (r in meters, rho in g/cm3)
_PREM_LAYERS = [
    ("innercore_boundary",    1221500, "INNERCORE", "polynomial",
     [13.0885, 0.0, -2.17742748697875934e-13]),
    ("coremantle_boundary",   3480000, "OUTERCORE", "polynomial",
     [12.5815, -1.98367603202009108e-07, -8.97421093229181259e-14,
      -2.13773109929070169e-20]),
    ("lowermantle_boundary",  5701000, "MANTLE", "polynomial",
     [7.9565, -1.01649662533354259e-06, 1.36199775701391389e-13,
      -1.19131495406828110e-20]),
    ("lower_transition",      5771000, "MANTLE", "polynomial",
     [5.3197, -2.32867681682624407e-07]),
    ("middle_transition",     5971000, "MANTLE", "polynomial",
     [11.2494, -1.26036728927954783e-06]),
    ("upper_transition",      6151000, "MANTLE", "polynomial",
     [7.1089, -5.97159001726573544e-07]),
    ("moho_boundary",    int(R_PREM - _MOHO_DEPTH), "MANTLE", "polynomial",
     [2.691, 1.08679956050855438e-07]),
    ("inner_crust",           6356000, "ROCK", "constant", 2.900),
    ("upper_crust",     int(R_PREM),  "ROCK", "constant", 2.600),
]

# Atmosphere shells approximating US Standard Atmosphere 1976.
# Average density per shell from: rho_avg = rho_0*H*(exp(-h1/H) - exp(-h2/H))/(h2-h1)
# where rho_0 = 1.225e-3 g/cm3 and H = 8500 m (scale height).
_RHO_0 = 1.225e-3
_SCALE_HEIGHT = 8500.0

_ATMO_SHELL_ALTS = [
    (0,      4000),
    (4000,   10000),
    (10000,  20000),
    (20000,  40000),
    (40000,  80000),
    (80000,  180000),
]


def _atmo_avg_density(h1: float, h2: float) -> float:
    return (_RHO_0 * _SCALE_HEIGHT
            * (math.exp(-h1 / _SCALE_HEIGHT) - math.exp(-h2 / _SCALE_HEIGHT))
            / (h2 - h1))


_ATMO_LAYERS = []
for _i, (_h1, _h2) in enumerate(_ATMO_SHELL_ALTS):
    _ATMO_LAYERS.append((
        f"atmo_{_h1 // 1000}_{_h2 // 1000}km",
        int(R_PREM + _h2),
        "AIR",
        "constant",
        _atmo_avg_density(_h1, _h2),
    ))


def _all_layers():
    """All layers from innermost to outermost."""
    return list(_PREM_LAYERS) + list(_ATMO_LAYERS)


def add_earth_model(model: Any) -> None:
    """Add PREM + atmosphere sectors to a DetectorModel loaded from GDML.

    Must be called after model.LoadGDML() so that the GDML volumes already
    have positive levels. The PREM sectors get negative levels and serve
    as the background Earth wherever the GDML geometry does not reach.
    """
    from siren.math import Vector3D
    from siren.geometry import Sphere, Placement
    from siren.detector import (
        DetectorSector,
        RadialAxis1D,
        PolynomialDistribution1D,
        RadialAxisPolynomialDensityDistribution,
        ConstantDensityDistribution,
    )

    materials_path = os.path.join(_THIS_DIR, "materials.dat")
    model.LoadMaterialModel(materials_path)

    earth_center = Vector3D(0.0, -(R_PREM - _GRADE_Y_BNB), 0.0)
    placement = Placement(earth_center)
    radial_axis = RadialAxis1D(earth_center)

    layers = _all_layers()

    for i, (name, outer_radius, material, dtype, params) in enumerate(layers):
        sector = DetectorSector()
        sector.name = name
        sector.material_id = model.GetMaterials().GetMaterialId(material)
        # Innermost layer gets level -1 (highest priority among Earth sectors).
        sector.level = -(i + 1)
        sector.geo = Sphere(placement, float(outer_radius), 0.0)

        if dtype == "constant":
            sector.density = ConstantDensityDistribution(float(params))
        else:
            poly = PolynomialDistribution1D(list(params))
            sector.density = RadialAxisPolynomialDensityDistribution(
                radial_axis, poly)

        model.AddSector(sector)
