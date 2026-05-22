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


_ENV_MATERIALS = """\
    <isotope N="1" Z="1" name="env_Z1_A1">
      <atom unit="g/mole" value="1.008"/>
    </isotope>
    <element name="env_el_Z1_A1">
      <fraction n="1.0" ref="env_Z1_A1"/>
    </element>
    <isotope N="12" Z="6" name="env_Z6_A12">
      <atom unit="g/mole" value="12"/>
    </isotope>
    <element name="env_el_Z6_A12">
      <fraction n="1.0" ref="env_Z6_A12"/>
    </element>
    <isotope N="14" Z="7" name="env_Z7_A14">
      <atom unit="g/mole" value="14"/>
    </isotope>
    <element name="env_el_Z7_A14">
      <fraction n="1.0" ref="env_Z7_A14"/>
    </element>
    <isotope N="16" Z="8" name="env_Z8_A16">
      <atom unit="g/mole" value="16"/>
    </isotope>
    <element name="env_el_Z8_A16">
      <fraction n="1.0" ref="env_Z8_A16"/>
    </element>
    <isotope N="24" Z="12" name="env_Z12_A24">
      <atom unit="g/mole" value="24"/>
    </isotope>
    <element name="env_el_Z12_A24">
      <fraction n="1.0" ref="env_Z12_A24"/>
    </element>
    <isotope N="27" Z="13" name="env_Z13_A27">
      <atom unit="g/mole" value="26.982"/>
    </isotope>
    <element name="env_el_Z13_A27">
      <fraction n="1.0" ref="env_Z13_A27"/>
    </element>
    <isotope N="28" Z="14" name="env_Z14_A28">
      <atom unit="g/mole" value="28.086"/>
    </isotope>
    <element name="env_el_Z14_A28">
      <fraction n="1.0" ref="env_Z14_A28"/>
    </element>
    <isotope N="40" Z="20" name="env_Z20_A40">
      <atom unit="g/mole" value="40"/>
    </isotope>
    <element name="env_el_Z20_A40">
      <fraction n="1.0" ref="env_Z20_A40"/>
    </element>
    <isotope N="56" Z="26" name="env_Z26_A56">
      <atom unit="g/mole" value="55.845"/>
    </isotope>
    <element name="env_el_Z26_A56">
      <fraction n="1.0" ref="env_Z26_A56"/>
    </element>
    <material name="env_GlacialTill" state="solid">
      <D value="1.9" unit="g/cm3"/>
      <fraction n="0.303830" ref="env_el_Z14_A28"/>
      <fraction n="0.534889" ref="env_el_Z8_A16"/>
      <fraction n="0.079387" ref="env_el_Z13_A27"/>
      <fraction n="0.034971" ref="env_el_Z26_A56"/>
      <fraction n="0.035734" ref="env_el_Z20_A40"/>
      <fraction n="0.011190" ref="env_el_Z1_A1"/>
    </material>
    <material name="env_Dolomite" state="solid">
      <D value="2.6" unit="g/cm3"/>
      <fraction n="0.2171" ref="env_el_Z20_A40"/>
      <fraction n="0.1318" ref="env_el_Z12_A24"/>
      <fraction n="0.1304" ref="env_el_Z6_A12"/>
      <fraction n="0.5207" ref="env_el_Z8_A16"/>
    </material>
    <material name="env_Air" state="gas">
      <D value="0.001225" unit="g/cm3"/>
      <fraction n="0.7562" ref="env_el_Z7_A14"/>
      <fraction n="0.2438" ref="env_el_Z8_A16"/>
    </material>"""


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

    site_half_y = 100.0
    site_half_x = 400.0
    site_half_z = 900.0
    margin = 10.0

    bedrock_y = _GRADE_Y_BNB - _TILL_THICKNESS

    atmo_height = site_half_y - _GRADE_Y_BNB
    atmo_center_y = _GRADE_Y_BNB + atmo_height / 2.0

    bedrock_height = site_half_y + bedrock_y
    bedrock_center_y = bedrock_y - bedrock_height / 2.0

    child_x = 2.0 * (site_half_x - margin)
    child_z = 2.0 * (site_half_z - margin)

    stub = f"""\
<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
{_ENV_MATERIALS}
  </materials>
  <solids>
    <box name="sol_site_volume" lunit="m" x="{2*site_half_x}" y="{2*site_half_y}" z="{2*site_half_z}"/>
    <box name="sol_atmosphere" lunit="m" x="{child_x}" y="{atmo_height}" z="{child_z}"/>
    <box name="sol_dolomite_bedrock" lunit="m" x="{child_x}" y="{bedrock_height}" z="{child_z}"/>
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
        <position unit="m" x="0" y="{atmo_center_y:.4f}" z="300"/>
      </physvol>
      <physvol name="pv_dolomite_bedrock">
        <volumeref ref="vol_dolomite_bedrock"/>
        <position unit="m" x="0" y="{bedrock_center_y:.4f}" z="300"/>
      </physvol>
{source_physvols}
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="vol_site_geology"/>
  </setup>
</gdml>
"""

    with open(cache_path, "w") as f:
        f.write(stub)

    return cache_path
