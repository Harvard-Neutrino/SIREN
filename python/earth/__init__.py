"""
siren.earth -- canonical Earth (PREM) and atmosphere models shared by all
SIREN detector loaders.

The submodules ``prem`` and ``atmosphere`` are pure-data (standard library only);
``sectors`` builds DetectorSectors and imports the compiled siren.detector module
lazily inside its builder function. Import this package for the canonical
constants and the shared ``build_earth_sectors`` builder.
"""

from .prem import (
    R_PREM,
    PREM_MATERIALS,
    CANONICAL_PREM_LAYERS,
    PremLayer,
    DensityForm,
    canonical_prem,
    collapse_to_constant,
)
from .atmosphere import (
    Atmosphere,
    USStandard1976,
    AtmosphereForm,
    AtmosphereBand,
    DEFAULT_SHELL_ALTS,
    default_boundaries,
    constant_shell_layers,
)
from .sectors import (
    build_earth_sectors,
    PositionedSector,
    shift_levels_up,
    insert_sectors_above,
)

__all__ = [
    "R_PREM", "PREM_MATERIALS", "CANONICAL_PREM_LAYERS", "PremLayer",
    "DensityForm", "canonical_prem", "collapse_to_constant",
    "Atmosphere", "USStandard1976", "AtmosphereForm", "AtmosphereBand",
    "DEFAULT_SHELL_ALTS", "default_boundaries", "constant_shell_layers",
    "build_earth_sectors", "PositionedSector", "shift_levels_up",
    "insert_sectors_above",
]
