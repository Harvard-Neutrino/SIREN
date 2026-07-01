"""
Shared atmosphere density model for SIREN detector loaders.

Promotes the previously per-detector (and, for DUNE, missing) atmosphere code
into a single shared module. Provides a pluggable barometric air-density model
(default US Standard Atmosphere 1976) and describes the atmosphere as a set of
altitude bands that ``siren.earth.sectors`` turns into detector sectors -- either
as true exponential-density shells (preferred) or, for backward compatibility,
as constant mass-conserving shells.

Dependency-light: standard library only. Sector construction (which needs the
compiled ``siren.detector`` module) lives in ``siren.earth.sectors``.

Altitudes here are measured above a reference surface (the loader anchors that
surface to a radius: R_PREM for SBN, the DEM site-surface radius for DUNE).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Tuple

# Default altitude bands (m above the reference surface), matching the six-shell
# US-Std-1976 approximation SBN has used. Band count sets the atmosphere sector
# count; each band becomes one sector.
DEFAULT_SHELL_ALTS: List[Tuple[float, float]] = [
    (0,      4000),
    (4000,   10000),
    (10000,  20000),
    (20000,  40000),
    (40000,  80000),
    (80000,  180000),
]


class Atmosphere:
    """Isothermal barometric air-density model: rho(h) = rho0 * exp(-h/H).

    Subclass and override ``density`` for a non-isothermal profile; the band
    machinery below only relies on ``density`` and ``scale_height``.
    """

    def __init__(self, rho0: float = 1.225e-3, scale_height: float = 8500.0):
        self.rho0 = rho0                # g/cm^3 at the reference surface
        self.scale_height = scale_height  # m

    def density(self, h: float) -> float:
        """Air density (g/cm^3) at altitude h (m) above the reference surface."""
        return self.rho0 * math.exp(-h / self.scale_height)

    def shell_mean_density(self, h1: float, h2: float) -> float:
        """Mass-conserving (column-preserving) mean density over band [h1, h2].

        rho_avg = rho0 * H * (exp(-h1/H) - exp(-h2/H)) / (h2 - h1)
        """
        H = self.scale_height
        return (self.rho0 * H
                * (math.exp(-h1 / H) - math.exp(-h2 / H))
                / (h2 - h1))


class USStandard1976(Atmosphere):
    """US Standard Atmosphere 1976 (isothermal approximation).

    rho0 = 1.225e-3 g/cm^3, scale height 8500 m -- the standard sea-level
    density and mean tropospheric scale height.
    """

    def __init__(self, rho0: float = 1.225e-3, scale_height: float = 8500.0):
        super().__init__(rho0=rho0, scale_height=scale_height)


def default_boundaries(top_m: float = 180000.0) -> List[Tuple[float, float]]:
    """Altitude bands (m) up to ``top_m``, clipped to the default set."""
    bands = [(h1, h2) for (h1, h2) in DEFAULT_SHELL_ALTS if h1 < top_m]
    if bands and bands[-1][1] > top_m:
        h1, _ = bands[-1]
        bands[-1] = (h1, top_m)
    return bands


@dataclass
class AtmosphereBand:
    """One atmosphere shell descriptor, resolved to absolute radii.

    ``mode`` selects how the sector builder realizes the density:
      - "exponential": rho(r) = rho_base * exp(-(r - r_base) / H) within the band
      - "constant":    rho(r) = rho_const (mass-conserving mean over the band)
    """
    name: str
    r_base: float          # inner radius (m): reference-surface radius + h1
    r_top: float           # outer radius (m): reference-surface radius + h2
    material: str
    mode: str
    rho_base: float        # g/cm^3 at r_base (exponential) or the constant value
    scale_height: float    # m (exponential only)


@dataclass
class AtmosphereForm:
    """How to build the atmosphere for a site.

    ``mode`` is "exponential" (true per-band exponential density, preferred) or
    "shells" (constant mass-conserving shells; legacy). ``surface_radius`` is the
    absolute radius (m) that altitude 0 maps to.
    """
    model: Atmosphere = field(default_factory=USStandard1976)
    mode: str = "exponential"
    top_m: float = 180000.0
    bands: List[Tuple[float, float]] = None

    def _bands(self) -> List[Tuple[float, float]]:
        return self.bands if self.bands is not None else default_boundaries(self.top_m)

    def resolve(self, surface_radius: float, name_prefix: str = "atmo") -> List[AtmosphereBand]:
        """Resolve to absolute-radius band descriptors for the sector builder."""
        out: List[AtmosphereBand] = []
        for (h1, h2) in self._bands():
            r_base = surface_radius + h1
            r_top = surface_radius + h2
            name = "{}_{}_{}km".format(name_prefix, int(h1) // 1000, int(h2) // 1000)
            if self.mode == "exponential":
                out.append(AtmosphereBand(
                    name=name, r_base=r_base, r_top=r_top, material="AIR",
                    mode="exponential", rho_base=self.model.density(h1),
                    scale_height=self.model.scale_height))
            else:
                out.append(AtmosphereBand(
                    name=name, r_base=r_base, r_top=r_top, material="AIR",
                    mode="constant", rho_base=self.model.shell_mean_density(h1, h2),
                    scale_height=self.model.scale_height))
        return out


def constant_shell_layers(model: Atmosphere = None,
                          surface_radius: float = 6371000.0,
                          bands: List[Tuple[float, float]] = None) -> List[Tuple]:
    """Legacy 5-tuple constant atmosphere layers, for compatibility re-exports.

    Returns ``(name, outer_radius_m, "AIR", "constant", rho_mean)`` entries.
    """
    if model is None:
        model = USStandard1976()
    if bands is None:
        bands = DEFAULT_SHELL_ALTS
    out = []
    for (h1, h2) in bands:
        out.append((
            "atmo_{}_{}km".format(int(h1) // 1000, int(h2) // 1000),
            int(surface_radius + h2),
            "AIR",
            "constant",
            model.shell_mean_density(h1, h2),
        ))
    return out
