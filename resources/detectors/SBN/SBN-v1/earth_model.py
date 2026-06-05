"""
PREM Earth model for the SBN detector at Fermilab.

Adds concentric PREM spherical shells and a multi-shell atmosphere to
a DetectorModel that was loaded via LoadGDML. The GDML volumes always
take priority over these Earth sectors wherever both exist.

The PREM layers follow Dziewonski & Anderson, PEPI 25 (1981) 297,
with the standard global-average Moho at 24.4 km depth. A local
crustal thickening correction overrides the upper mantle with crust
(ROCK) in the Midwest region from Fermilab to the DUNE far detector
at Sanford Lab, reflecting the ~45 km Moho depth measured by
CRUST1.0 (Laske et al. 2013) for this stable continental interior.

The local correction is an Earth-centered spherical shell spanning
the depth range from 24.4 km (PREM Moho) to 45 km (local Moho),
with a polar-angle cut limiting it to ~17 degrees (~1900 km) around
the radial direction through Fermilab. This covers the entire
Fermilab-to-Sanford-Lab baseline (1285 km) with margin.

The atmosphere approximates the US Standard Atmosphere 1976 as 6
constant-density spherical shells.

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

R_PREM = 6371000.0
_PREM_MOHO_DEPTH = 24400.0
_LOCAL_MOHO_DEPTH = 45000.0

# Angular half-opening of the local crustal thickening correction,
# measured from the radial direction through Fermilab.  0.3 rad
# corresponds to ~1900 km on the surface, comfortably covering the
# 1285 km Fermilab-to-Sanford-Lab baseline with ~600 km margin.
_LOCAL_MOHO_THETA = 0.3


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
    ("moho_boundary",    int(R_PREM - _PREM_MOHO_DEPTH), "MANTLE", "polynomial",
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


def add_earth_model(model: Any, ground_level: float) -> None:
    """Add PREM + atmosphere sectors to a DetectorModel loaded from GDML.

    Must be called after model.LoadGDML() so that the GDML volumes already
    have positive levels. The PREM sectors get negative levels and serve
    as the background Earth wherever the GDML geometry does not reach.

    A local crustal thickening sector is added between the PREM upper
    mantle and crustal layers. It is an Earth-centered spherical shell
    (depths 24.4--45 km) with a polar-angle cut (~17 deg around
    Fermilab) that overrides the PREM mantle with ROCK in the Midwest,
    where the continental Moho is ~45 km deep (CRUST1.0).
    """
    from siren.math import Vector3D, Quaternion
    from siren.geometry import Sphere, Placement
    from siren.detector import (
        DetectorSector,
        RadialAxis1D,
        PolynomialDistribution1D,
        RadialAxisPolynomialDensityDistribution,
        ConstantDensityDistribution,
    )

    # PREM materials (ROCK, MANTLE, AIR, etc.) are defined in the
    # composite GDML alongside the detector materials, so they go
    # through the GDML MergeFrom pipeline and get proper collision
    # handling.  No separate LoadMaterialModel call is needed.

    earth_center = Vector3D(0.0, -(R_PREM - ground_level), 0.0)
    placement = Placement(earth_center)
    radial_axis = RadialAxis1D(earth_center)

    # --- Build the full layer list with local correction inserted ---
    # The local crustal thickening correction goes between
    # upper_transition and moho_boundary so that its level sits
    # between them in the priority ordering.
    prem = _all_layers()
    moho_idx = next(i for i, (n, *_) in enumerate(prem) if n == "moho_boundary")

    # Rotate so the Sphere's +z axis (theta=0 direction) aligns with
    # the BNB +y axis (radially outward through Fermilab).
    rot = Quaternion.rotation_between(Vector3D(0, 0, 1), Vector3D(0, 1, 0))
    rotated_placement = Placement(earth_center, rot)

    local_moho_entry = {
        "name": "local_thick_crust",
        "material": "ROCK",
        "density_type": "constant",
        "density_params": 2.9,
        "geo": Sphere(
            rotated_placement,
            R_PREM - _PREM_MOHO_DEPTH,   # outer: PREM Moho depth
            R_PREM - _LOCAL_MOHO_DEPTH,   # inner: local Moho depth
            0.0, 2.0 * math.pi,          # full azimuth
            0.0, _LOCAL_MOHO_THETA,       # polar cap around Fermilab
        ),
    }

    sectors = model.Sectors

    prem_names = {name for name, *_ in prem}
    prem_names.add("local_thick_crust")
    for s in sectors:
        if s.name in prem_names:
            raise RuntimeError(
                "add_earth_model has already been applied to this model "
                "(sector '{}' already exists)".format(s.name))

    min_level = min([s.level for s in sectors if s.name != "UNIVERSE"], default=0)

    num_new_sectors = len(prem) + 1  # +1 for local correction

    # Shift existing levels down to make room for the new sectors
    for s in sectors:
        if s.name != "UNIVERSE" and s.level >= min_level:
            s.level = s.level - min_level + num_new_sectors

    base_level = num_new_sectors - 1

    # Prepend the new sectors in order from innermost to outermost
    # Decrementing the level each time
    for i, (name, outer_radius, material, dtype, params) in enumerate(prem):
        # Insert local correction right before moho_boundary
        if i == moho_idx:
            lm = local_moho_entry
            sector = DetectorSector()
            sector.name = lm["name"]
            sector.material_id = model.GetMaterials().GetMaterialId(lm["material"])
            sector.level = base_level
            sector.geo = lm["geo"]
            sector.density = ConstantDensityDistribution(float(lm["density_params"]))
            sectors.append(sector)
            base_level -= 1

        sector = DetectorSector()
        sector.name = name
        sector.material_id = model.GetMaterials().GetMaterialId(material)
        sector.level = base_level
        sector.geo = Sphere(placement, float(outer_radius), 0.0)

        if dtype == "constant":
            sector.density = ConstantDensityDistribution(float(params))
        else:
            poly = PolynomialDistribution1D(list(params))
            sector.density = RadialAxisPolynomialDensityDistribution(
                radial_axis, poly)

        base_level -= 1
        sectors.append(sector)
    if base_level != -1:
        raise RuntimeError(
            "Earth model level assignment error: expected base_level == -1, "
            "got {}".format(base_level))

    # Sort the sectors by level before putting them back in the model
    sectors = sorted(sectors, key=lambda s: s.level)
    model.Sectors = sectors
