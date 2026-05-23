"""
Shared GDML composition engine for SBN detector loaders.

Generates a stub GDML that references source GDML files via the <file>
element. The SIREN GDML parser handles name collision resolution
automatically via GDMLData::MergeFrom.
"""

from __future__ import annotations

import os
from typing import Any


def _ensure_gdml_files(abs_dir: str, sources: list[dict[str, Any]]) -> None:
    """Download missing GDML files using siren.download."""
    from siren.download import ensure_files
    specs = []
    for src in sources:
        if not src.get("file"):
            continue
        specs.append({
            "path": os.path.join(abs_dir, src["file"]),
            "url": src.get("url", ""),
            "sha256": src.get("sha256", ""),
        })
    ensure_files(specs)


# Isotope table: (Z, A, molar_mass_g_per_mol).
# Each entry produces an <isotope> and a single-fraction <element>.
_ISOTOPES = [
    (1,   1,   1.008),
    (6,  12,  12.0),
    (7,  14,  14.0),
    (8,  16,  16.0),
    (12, 24,  24.0),
    (13, 27,  26.982),
    (14, 28,  28.086),
    (20, 40,  40.0),
    (26, 56,  55.845),
]

# Material table: (name, state, density_g_cm3, {element_key: mass_fraction}).
# element_key is (Z, A) referencing an entry in _ISOTOPES.
_MATERIALS = [
    # -- Site geology (used by GDML volumes) --
    ("env_GlacialTill", "solid", 1.9, {
        (14, 28): 0.303830, (8, 16): 0.534889, (13, 27): 0.079387,
        (26, 56): 0.034971, (20, 40): 0.035734, (1, 1): 0.011190,
    }),
    ("env_Dolomite", "solid", 2.6, {
        (20, 40): 0.2171, (12, 24): 0.1318,
        (6, 12): 0.1304, (8, 16): 0.5207,
    }),
    ("env_Air", "gas", 0.001225, {(7, 14): 0.7562, (8, 16): 0.2438}),
    # -- PREM Earth model (used by AddSector layers; density here is
    #    nominal -- actual density comes from each sector's
    #    DensityDistribution). --
    ("ROCK",      "solid",  2.6,   {(14, 28): 0.4674, (8, 16): 0.5326}),
    ("MANTLE",    "solid",  3.3,   {(14, 28): 0.4674, (8, 16): 0.5326}),
    ("OUTERCORE", "solid", 10.0,   {(26, 56): 1.0}),
    ("INNERCORE", "solid", 13.0,   {(26, 56): 1.0}),
    ("AIR",       "gas",    0.001225, {(7, 14): 0.7562, (8, 16): 0.2438}),
]


def _build_materials_xml() -> str:
    """Generate GDML <materials> content from the isotope and material tables."""
    lines = []
    for Z, A, mass in _ISOTOPES:
        iso = f"env_Z{Z}_A{A}"
        el = f"env_el_Z{Z}_A{A}"
        lines.append(f'    <isotope N="{A}" Z="{Z}" name="{iso}">')
        lines.append(f'      <atom unit="g/mole" value="{mass}"/>')
        lines.append(f'    </isotope>')
        lines.append(f'    <element name="{el}">')
        lines.append(f'      <fraction n="1.0" ref="{iso}"/>')
        lines.append(f'    </element>')
    for name, state, density, comp in _MATERIALS:
        lines.append(f'    <material name="{name}" state="{state}">')
        lines.append(f'      <D value="{density}" unit="g/cm3"/>')
        for (Z, A), frac in comp.items():
            lines.append(f'      <fraction n="{frac}" ref="env_el_Z{Z}_A{A}"/>')
        lines.append(f'    </material>')
    return "\n".join(lines)


# Ground-level elevation in BNB coordinates (meters above the BNB target).
# Fermilab grade is approximately 7.62 m above the BNB beam axis.
_GRADE_Y_BNB = 7.62

# Glacial till thickness below grade (meters). The Quaternary glacial
# deposits at Fermilab are 18-30 m thick; we use 20 m as a representative
# value. Below this is Silurian dolomite bedrock.
_TILL_THICKNESS = 20.0


def build_composite(
    abs_dir: str,
    sources: list[dict[str, Any]],
    cache_name: str = "composite.gdml",
) -> str:
    """Build a stub GDML that references source files via <file> elements.

    Each source dict has keys:
      file     : relative path to GDML file
      prefix   : name prefix (used for physvol naming only)
      position : (x, y, z) in meters (BNB frame)
      rotation : (rx, ry, rz) GDML Euler angles in radians, or None
      unwrap   : if True, set as_assembly="true" on the <file> element

    Returns the path to the written stub GDML.
    """
    cache_path = os.path.join(abs_dir, cache_name)

    _ensure_gdml_files(abs_dir, sources)

    physvols = []
    for spec in sources:
        if not spec.get("file"):
            continue
        prefix = spec["prefix"]
        pos = spec["position"]
        rot = spec.get("rotation")
        as_asm = ' as_assembly="true"' if spec.get("unwrap", False) else ""

        pv = []
        pv.append(f'      <physvol name="pv_{prefix}">')
        pv.append(f'        <file name="{spec["file"]}"{as_asm}/>')
        pv.append(f'        <position unit="m" x="{pos[0]:.6f}" y="{pos[1]:.6f}" z="{pos[2]:.6f}"/>')
        if rot is not None:
            pv.append(f'        <rotation unit="rad" x="{rot[0]:.10f}" y="{rot[1]:.10f}" z="{rot[2]:.10f}"/>')
        pv.append('      </physvol>')
        physvols.append("\n".join(pv))

    source_physvols = "\n".join(physvols)

    # GDML <box> x/y/z are half-widths
    # So this gives a full 1600 m x 400 m x 1800 m box
    box_half_x = 800.0
    box_half_y = 200.0
    box_half_z = 900.0
    margin = 10.0

    bedrock_y = _GRADE_Y_BNB - _TILL_THICKNESS

    atmo_height = box_half_y - _GRADE_Y_BNB
    atmo_center_y = _GRADE_Y_BNB + atmo_height / 2.0
    atmo_half_height = atmo_height / 2.0

    bedrock_height = box_half_y + bedrock_y
    bedrock_center_y = bedrock_y - bedrock_height / 2.0
    bedrock_half_height = bedrock_height / 2.0

    child_half_x = box_half_x - margin
    child_half_z = box_half_z - margin

    materials_xml = _build_materials_xml()

    stub = f"""\
<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
{materials_xml}
  </materials>
  <solids>
    <box name="sol_site_volume" lunit="m" x="{box_half_x}" y="{box_half_y}" z="{box_half_z}"/>
    <box name="sol_atmosphere" lunit="m" x="{child_half_x}" y="{atmo_half_height}" z="{child_half_z}"/>
    <box name="sol_dolomite_bedrock" lunit="m" x="{child_half_x}" y="{bedrock_half_height}" z="{child_half_z}"/>
  </solids>
  <structure>
    <volume name="vol_atmosphere">
      <materialref ref="env_Air"/>
      <solidref ref="sol_atmosphere"/>
    </volume>
    <volume name="vol_dolomite_bedrock">
      <materialref ref="env_Dolomite"/>
      <solidref ref="sol_dolomite_bedrock"/>
    </volume>
    <volume name="vol_site_geology">
      <materialref ref="env_GlacialTill"/>
      <solidref ref="sol_site_volume"/>
      <physvol name="pv_atmosphere">
        <volumeref ref="vol_atmosphere"/>
        <position unit="m" x="0" y="{atmo_center_y:.4f}" z="0"/>
      </physvol>
      <physvol name="pv_dolomite_bedrock">
        <volumeref ref="vol_dolomite_bedrock"/>
        <position unit="m" x="0" y="{bedrock_center_y:.4f}" z="0"/>
      </physvol>
{source_physvols}
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="vol_site_geology"/>
  </setup>
</gdml>
"""

    tmp_path = cache_path + ".tmp"
    try:
        with open(tmp_path, "w", encoding="utf-8") as f:
            f.write(stub)
        os.replace(tmp_path, cache_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise

    return cache_path
