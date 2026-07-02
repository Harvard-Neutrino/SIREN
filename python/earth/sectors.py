"""
Shared sector-assembly for SIREN earth models.

Turns the canonical PREM layers (siren.earth.prem) and an atmosphere model
(siren.earth.atmosphere) into DetectorSectors on a GDML-loaded DetectorModel,
using the predict-count-and-shift-up level scheme: the number of new background
sectors is counted up front, existing GDML sectors are shifted up to make room,
and the new earth sectors occupy the freed low levels (innermost = highest
priority among the background, GDML volumes always above all of them).

The compiled siren.detector/geometry/math modules are imported lazily inside the
builder so this module (and the pure prem/atmosphere data modules) can be
imported cheaply for their constants.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, Tuple


@dataclass
class PositionedSector:
    """A site-specific background sector inserted into the earth stack.

    ``insert_before`` names a PREM layer; this sector is emitted just before it
    in the innermost->outermost walk, giving it higher priority than that layer
    (e.g. a local crustal thickening that overrides the PREM mantle).
    ``geo`` is a Geometry; supply either ``density`` (a DensityDistribution) or
    ``density_const`` (a scalar g/cm^3).
    """
    name: str
    material: str
    geo: Any
    density: Any = None
    density_const: Optional[float] = None
    insert_before: Optional[str] = None


def shift_levels_up(sectors, count, exempt=("UNIVERSE",)):
    """Shift existing sector levels up by ``count`` to open room below them.

    Normalizes the non-exempt sectors so their lowest level lands at ``count``
    (freeing levels 0..count-1 for new background sectors). Mutates in place and
    returns the min level that was used as the shift baseline.
    """
    movable = [s for s in sectors if s.name not in exempt]
    min_level = min((s.level for s in movable), default=0)
    for s in movable:
        if s.level >= min_level:
            s.level = s.level - min_level + count
    return min_level


def insert_sectors_above(model, anchor_name, new_sectors):
    """Insert already-built DetectorSectors just above ``anchor_name``.

    Uses the predict-count-and-shift-up scheme: every existing sector with a
    higher level than the anchor is shifted up by ``len(new_sectors)`` to open a
    gap, then the new sectors are assigned the freed levels (in ascending
    priority: the first entry lands just above the anchor). This places
    site-specific overburden between a background volume (e.g. "World") and the
    detector volumes without any absolute/negative level constants.
    """
    sectors = model.Sectors
    anchor = next((s for s in sectors if s.name == anchor_name), None)
    if anchor is None:
        raise RuntimeError(
            "insert_sectors_above: no sector named {!r}".format(anchor_name))
    anchor_level = anchor.level
    n = len(new_sectors)
    for s in sectors:
        if s.level > anchor_level:
            s.level = s.level + n
    level = anchor_level + 1
    for s in new_sectors:
        s.level = level
        level += 1
        sectors.append(s)
    model.Sectors = sorted(sectors, key=lambda s: s.level)


def build_earth_sectors(model, *, surface_offset, prem_layers=None,
                        atmosphere=None, extra_sectors=(),
                        up_axis=(0.0, 1.0, 0.0), r_prem=6371000.0):
    """Add PREM + atmosphere (+ site extras) to a GDML-loaded DetectorModel.

    Parameters
    ----------
    model : DetectorModel
        Must already have GDML volumes loaded (positive levels).
    surface_offset : float
        Detector-frame coordinate (along up_axis) of the local ground surface;
        the Earth center is placed at up_axis * -(r_prem - surface_offset).
    prem_layers : list of 5-tuples, optional
        (name, outer_radius_m, material, density_type, params). Defaults to the
        canonical PREM.
    atmosphere : AtmosphereForm, optional
        Atmosphere to append above r_prem (see siren.earth.atmosphere). If None,
        no atmosphere sectors are added.
    extra_sectors : sequence of PositionedSector
        Site-specific background sectors inserted into the PREM stack.
    up_axis : (x, y, z)
        Local "up" (radially-outward) direction in detector coordinates.
    """
    from siren.math import Vector3D
    from siren.geometry import Sphere, Placement
    from siren.detector import (
        DetectorSector,
        RadialAxis1D,
        PolynomialDistribution1D,
        RadialAxisPolynomialDensityDistribution,
        ConstantDensityDistribution,
        ExponentialDistribution1D,
        RadialAxisExponentialDensityDistribution,
    )

    from .prem import CANONICAL_PREM_LAYERS

    if prem_layers is None:
        prem_layers = CANONICAL_PREM_LAYERS

    ux, uy, uz = up_axis
    earth_center = Vector3D(ux * -(r_prem - surface_offset),
                            uy * -(r_prem - surface_offset),
                            uz * -(r_prem - surface_offset))
    placement = Placement(earth_center)
    radial_axis = RadialAxis1D(earth_center)
    materials = model.GetMaterials()

    def _make_sphere_sector(name, outer_radius, material, density_dist):
        s = DetectorSector()
        s.name = name
        s.material_id = materials.GetMaterialId(material)
        s.geo = Sphere(placement, float(outer_radius), 0.0)
        s.density = density_dist
        return s

    def _prem_density(dtype, params):
        if dtype == "constant":
            return ConstantDensityDistribution(float(params))
        poly = PolynomialDistribution1D(list(params))
        return RadialAxisPolynomialDensityDistribution(radial_axis, poly)

    def _extra_sector(ex):
        s = DetectorSector()
        s.name = ex.name
        s.material_id = materials.GetMaterialId(ex.material)
        s.geo = ex.geo
        if ex.density is not None:
            s.density = ex.density
        else:
            s.density = ConstantDensityDistribution(float(ex.density_const))
        return s

    def _atmo_sector(band):
        if band.mode == "exponential":
            # rho(r) = rho_base * exp(-(r - r_base)/H), anchored at the band base
            # so the exponent stays small (no overflow/underflow at earth radius).
            expo = ExponentialDistribution1D(
                -1.0 / band.scale_height, band.rho_base, band.r_base)
            dens = RadialAxisExponentialDensityDistribution(radial_axis, expo)
        else:
            dens = ConstantDensityDistribution(float(band.rho_base))
        return _make_sphere_sector(band.name, band.r_top, band.material, dens)

    # --- build the ordered new-sector list (innermost -> outermost) ---
    extras_by_before = {}
    for ex in extra_sectors:
        extras_by_before.setdefault(ex.insert_before, []).append(ex)

    ordered_new = []
    for (name, outer_radius, material, dtype, params) in prem_layers:
        for ex in extras_by_before.get(name, ()):
            ordered_new.append(_extra_sector(ex))
        ordered_new.append(
            _make_sphere_sector(name, outer_radius, material,
                                _prem_density(dtype, params)))
    # any extras keyed to a missing layer name -> append at outer edge of PREM
    for ex in extras_by_before.get(None, ()):
        ordered_new.append(_extra_sector(ex))

    if atmosphere is not None:
        for band in atmosphere.resolve(r_prem):
            ordered_new.append(_atmo_sector(band))

    # --- guard against double application ---
    sectors = model.Sectors
    new_names = {s.name for s in ordered_new}
    for s in sectors:
        if s.name in new_names:
            raise RuntimeError(
                "add_earth_model has already been applied to this model "
                "(sector '{}' already exists)".format(s.name))

    count = len(ordered_new)
    shift_levels_up(sectors, count)

    # innermost gets the highest background level (wins inside overlapping shells)
    base_level = count - 1
    for s in ordered_new:
        s.level = base_level
        base_level -= 1
        sectors.append(s)
    if base_level != -1:
        raise RuntimeError(
            "Earth model level assignment error: expected base_level == -1, "
            "got {}".format(base_level))

    model.Sectors = sorted(sectors, key=lambda s: s.level)
