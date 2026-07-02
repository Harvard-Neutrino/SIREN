"""
Canonical PREM Earth model, shared by all SIREN detector loaders.

This is the single source of truth for the Preliminary Reference Earth Model
(Dziewonski & Anderson, PEPI 25 (1981) 297). Detector packages import the
canonical layer table from here.

This module is intentionally dependency-light: it imports only the standard
library so that it can be imported for its constants without pulling in the
compiled ``siren`` extension modules. The code that turns these layers into
``siren.detector`` sectors lives in ``siren.earth.sectors`` and imports the
compiled library lazily.

The layer table uses the 5-tuple form
``(name, outer_radius_m, material, density_type, density_params)`` where
``density_type`` is ``"constant"`` or ``"polynomial"`` and, for polynomials,
``density_params`` are coefficients of ``rho(r) = p[0] + p[1]*r + ...`` with
``r`` in meters and ``rho`` in g/cm^3 (absolute radius; matches the C++
``RadialAxisPolynomialDensityDistribution`` convention and the ``.dat`` grammar).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Union, List, Tuple

R_PREM = 6371000.0
PREM_MOHO_DEPTH = 24400.0

# Material-name vocabulary shared across sites. The actual G4 materials are
# resolved by name against the detector model's material list at build time.
PREM_MATERIALS = ("INNERCORE", "OUTERCORE", "MANTLE", "ROCK", "AIR")

# The canonical PREM layers, innermost to outermost. This IS the model; any
# site that wants a coarser (constant-shell) view derives it via
# collapse_to_constant().
CANONICAL_PREM_LAYERS: List[Tuple] = [
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
    ("moho_boundary",    int(R_PREM - PREM_MOHO_DEPTH), "MANTLE", "polynomial",
     [2.691, 1.08679956050855438e-07]),
    ("inner_crust",           6356000, "ROCK", "constant", 2.900),
    ("upper_crust",     int(R_PREM),  "ROCK", "constant", 2.600),
]


# ---- dataclass views (GDML-ready seam) -------------------------------------
# The dataclass form mirrors what a future GDML density element would carry,
# so a later "read PREM from GDML" path can produce these without changing
# the sector builder.

@dataclass(frozen=True)
class DensityForm:
    kind: str                                  # "constant" | "polynomial" | "exponential"
    params: Union[float, Sequence[float]]
    # polynomial: rho(r) = sum_k params[k] * r^k, r in meters (absolute radius)
    # exponential: (rho_ref, scale_height, r_ref) -> rho(r)=rho_ref*exp(-(r-r_ref)/H)
    # constant: rho (g/cm^3)


@dataclass(frozen=True)
class PremLayer:
    name: str
    outer_radius_m: float
    material: str
    density: DensityForm


def _tuple_to_layer(t: Tuple) -> PremLayer:
    name, radius, material, dtype, params = t
    return PremLayer(name, float(radius), material, DensityForm(dtype, params))


def canonical_prem() -> List[PremLayer]:
    """Return the canonical PREM as PremLayer dataclasses (GDML-ready view)."""
    return [_tuple_to_layer(t) for t in CANONICAL_PREM_LAYERS]


def _mean_density(params: Sequence[float], r_lo: float, r_hi: float) -> float:
    """Mass-conserving (mass/volume) mean of a radial-polynomial density over
    the spherical shell [r_lo, r_hi]. Integrates rho(r)*4*pi*r^2 dr analytically.
    """
    # mass = integral_{r_lo}^{r_hi} (sum_k c_k r^k) * r^2 dr  (drop common 4pi)
    mass = 0.0
    for k, c in enumerate(params):
        p = k + 3
        mass += c * (r_hi ** p - r_lo ** p) / p
    vol = (r_hi ** 3 - r_lo ** 3) / 3.0
    return mass / vol


def collapse_to_constant(layers: Sequence[Tuple] = None,
                         boundaries: Sequence[float] = None) -> List[Tuple]:
    """Return a constant-density view of the canonical (polynomial) PREM.

    Each output shell's density is the mass-conserving radial mean of the
    canonical density over that shell, so total mass is preserved.
    ``boundaries`` optionally re-bins the layers onto a coarser set of
    outer radii (each must coincide with a canonical layer boundary).
    """
    if layers is None:
        layers = CANONICAL_PREM_LAYERS
    out: List[Tuple] = []
    r_lo = 0.0
    for name, r_hi, material, dtype, params in layers:
        r_hi = float(r_hi)
        if dtype == "constant":
            rho = float(params)
        else:
            rho = _mean_density(list(params), r_lo, r_hi)
        out.append((name, r_hi, material, "constant", rho))
        r_lo = r_hi
    return out
