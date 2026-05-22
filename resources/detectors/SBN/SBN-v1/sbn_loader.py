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
    <isotope N="40" Z="20" name="env_Z20_A40">
      <atom unit="g/mole" value="40"/>
    </isotope>
    <element name="env_el_Z20_A40">
      <fraction n="1.0" ref="env_Z20_A40"/>
    </element>
    <isotope N="24" Z="12" name="env_Z12_A24">
      <atom unit="g/mole" value="24"/>
    </isotope>
    <element name="env_el_Z12_A24">
      <fraction n="1.0" ref="env_Z12_A24"/>
    </element>
    <isotope N="12" Z="6" name="env_Z6_A12">
      <atom unit="g/mole" value="12"/>
    </isotope>
    <element name="env_el_Z6_A12">
      <fraction n="1.0" ref="env_Z6_A12"/>
    </element>
    <material name="env_Dolomite" state="solid">
      <D value="2.24" unit="g/cm3"/>
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

    GRADE_Y_BNB = 7.62

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

    stub = f"""\
<?xml version="1.0"?>
<gdml>
  <define/>
  <materials>
{_ENV_MATERIALS}
  </materials>
  <solids>
    <box name="sol_composite_world" lunit="m" x="400" y="100" z="800"/>
    <box name="sol_atmosphere" lunit="m" x="800" y="100" z="1600"/>
  </solids>
  <structure>
    <volume name="vol_atmosphere">
      <materialref ref="env_Air"/>
      <solidref ref="sol_atmosphere"/>
    </volume>
    <volume name="volCompositeWorld">
      <materialref ref="env_Dolomite"/>
      <solidref ref="sol_composite_world"/>
      <physvol name="pv_atmosphere">
        <volumeref ref="vol_atmosphere"/>
        <position unit="m" x="0" y="{GRADE_Y_BNB + 50.0:.2f}" z="300"/>
      </physvol>
{source_physvols}
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="volCompositeWorld"/>
  </setup>
</gdml>
"""

    with open(cache_path, "w") as f:
        f.write(stub)

    return cache_path
