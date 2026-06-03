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


# ----------------------------------------------------------------------
# MiniBooNE detector + vault enclosure geometry
#
# Spherical carbon-steel tank (inner radius 6.096 m, 40 ft diameter) of
# Marcol 7 mineral oil, with an opaque optical barrier at r = 5.746 m
# (574.6 cm) splitting the oil into an inner signal region and a 35 cm
# veto shell. The tank hangs on supports inside a vertical cylindrical
# CONCRETE VAULT (air-filled), with an electronics room above and a dirt
# overburden berm on top -- per the MiniBooNE enclosure cross-section.
#
# Detector values verified against the primary MiniBooNE papers:
#   - inner radius 6.096 m, optical barrier 574.6 cm, veto 35 cm:
#     NIM A 599 (2009) 28 [arXiv:0806.4201], Sec. 1.3.
#   - Marcol 7 oil, rho = 0.845 g/cm3 (CH2, H:C ~2), n_D = 1.4684:
#     NIM A 599, Sec. 3 / Table 1.
#   - 1280 inner + 240 veto 8-inch PMTs (optical detail, not modelled).
#   - 541 m baseline / 1.9 m beam offset: flux PRD 79 (2009) 072002
#     [arXiv:0806.1449]; refined to (0, 1.896, 541.34) m by the G4BNB
#     survey (bsim::Location, NuBeamOutput.cc:136).
#
# The vault inner diameter is the PUBLISHED 45 ft / 13.716 m and the dirt
# overburden is the PUBLISHED >= 3 m (NIM A 599 Sec. 1.4), both trusted
# over the diagram (whose estimates were ~15.26 m and ~2.5 m). All other
# enclosure dimensions (heights, wall/slab thicknesses, clearance, room,
# cap) have no published value, so they are scaled off the enclosure
# cross-section diagram (anchor: detector sphere = 40 ft = 12.192 m =
# 913 px; 1 px = 1.3354 cm) and are ESTIMATES.
#
# Modelling ASSUMPTIONS: steel-shell thickness ~1 cm (back-computed from
# the ~37 t shell mass; not a published value); concrete rho = 2.30,
# dirt rho = 1.90; detector support legs and the side overflow tank are
# omitted (negligible mass); the berm is approximated by the dirt world
# block rather than a sloped frustum.
# ----------------------------------------------------------------------
_MB_BARRIER_RADIUS_M = 5.746    # optical barrier (inner signal / veto)
_MB_TANK_INNER_M = 6.096        # tank inner radius / oil outer (20 ft)
_MB_STEEL_THICKNESS_M = 0.010   # carbon-steel shell (back-computed)
_MB_STEEL_OUTER_M = _MB_TANK_INNER_M + _MB_STEEL_THICKNESS_M
_MB_OIL_DENSITY = 0.845         # Marcol 7, NIM A 599 Table 1

# Enclosure dimensions from the diagram (pixels -> metres).
_MB_PX_M = 12.192 / 913.0       # 40 ft sphere = 913 px
def _mbpx(px):
    return px * _MB_PX_M

_MB_VAULT_AIR_R = 13.716 / 2.0        # vault inner radius: 45 ft, NIM A 599 Sec. 1.4
_MB_VAULT_AIR_H = _mbpx(979)          # vault air cavity height (diagram estimate)
_MB_VAULT_WALL_T = _mbpx(35)          # vault wall thickness (diagram estimate)
_MB_SLAB_R = _MB_VAULT_AIR_R + _MB_VAULT_WALL_T  # slabs flush with outer wall
_MB_LSLAB_H = _mbpx(72)               # lower slab (floor) thickness
_MB_USLAB_H = _mbpx(39)               # upper slab thickness
_MB_TANK_CLEAR = _mbpx(76)            # sphere bottom above lower-slab top
_MB_CAP_R = _mbpx(183) / 2.0          # oil chimney radius
_MB_CAP_H = _mbpx(80)                 # oil chimney height
_MB_ROOM_AIR_H = _mbpx(270)           # electronics room air height
_MB_ROOM_WALL_T = _mbpx(26)           # electronics room wall thickness
_MB_ROOF_H = _mbpx(49)                # room roof slab thickness
_MB_BERM_H = 3.0                      # dirt overburden above roof: >= 3 m, NIM A 599 Sec. 1.4
_MB_VAULT_WALL_OUTER = _MB_VAULT_AIR_R + _MB_VAULT_WALL_T
_MB_ROOM_WALL_OUTER = _MB_VAULT_AIR_R + _MB_ROOM_WALL_T

# Vertical levels, sphere centre at y = 0, +y up (= BNB up axis).
_MB_Y_LSLAB_TOP = -_MB_TANK_INNER_M - _MB_TANK_CLEAR
_MB_Y_LSLAB_BOT = _MB_Y_LSLAB_TOP - _MB_LSLAB_H
_MB_Y_CAV_TOP = _MB_Y_LSLAB_TOP + _MB_VAULT_AIR_H
_MB_Y_GRADE = _MB_Y_CAV_TOP + _MB_USLAB_H          # upper-slab top = grade
_MB_Y_ROOM_TOP = _MB_Y_GRADE + _MB_ROOM_AIR_H
_MB_Y_ROOF_TOP = _MB_Y_ROOM_TOP + _MB_ROOF_H
_MB_Y_BERM_TOP = _MB_Y_ROOF_TOP + _MB_BERM_H
_MB_WORLD_HALF = 12.5                               # world dirt box half-x/z


def _build_miniboone_gdml():
    """Assemble the MiniBooNE enclosure GDML from the scaled dimensions.

    All cylindrical volumes are GDML <tube>s whose axis (local z) is rotated
    90 deg about x so it points along +y (vertical). They are placed flat in
    a dirt world block; the detector oil/steel/cap are emitted last so they
    take precedence in their overlaps with the vault air cavity.
    """
    def tube(name, rmin, rmax, h):
        return (f'    <tube name="{name}" lunit="m" aunit="deg" rmin="{rmin:.4f}" '
                f'rmax="{rmax:.4f}" z="{h / 2.0:.4f}" startphi="0" deltaphi="360"/>')

    def vpv(pv, vol, yc):
        return (f'      <physvol name="{pv}">\n'
                f'        <volumeref ref="{vol}"/>\n'
                f'        <position unit="m" x="0" y="{yc:.4f}" z="0"/>\n'
                f'        <rotation unit="deg" x="90" y="0" z="0"/>\n'
                f'      </physvol>')

    def opv(pv, vol):
        return f'      <physvol name="{pv}"><volumeref ref="{vol}"/></physvol>'

    cav_yc = 0.5 * (_MB_Y_LSLAB_TOP + _MB_Y_CAV_TOP)
    room_yc = _MB_Y_GRADE + 0.5 * _MB_ROOM_AIR_H
    world_half_y = max(abs(_MB_Y_LSLAB_BOT), abs(_MB_Y_BERM_TOP)) + 0.1

    solids = "\n".join([
        f'    <box name="mb_world" lunit="m" x="{_MB_WORLD_HALF}" y="{world_half_y:.4f}" z="{_MB_WORLD_HALF}"/>',
        f'    <orb name="mb_inner_oil" lunit="m" r="{_MB_BARRIER_RADIUS_M}"/>',
        f'    <sphere name="mb_veto_oil" lunit="m" aunit="deg" rmin="{_MB_BARRIER_RADIUS_M}" rmax="{_MB_TANK_INNER_M}" startphi="0" deltaphi="360" starttheta="0" deltatheta="180"/>',
        f'    <sphere name="mb_steel" lunit="m" aunit="deg" rmin="{_MB_TANK_INNER_M}" rmax="{_MB_STEEL_OUTER_M}" startphi="0" deltaphi="360" starttheta="0" deltatheta="180"/>',
        tube("mb_lslab", 0.0, _MB_SLAB_R, _MB_LSLAB_H),
        tube("mb_uslab", 0.0, _MB_SLAB_R, _MB_USLAB_H),
        tube("mb_vwall", _MB_VAULT_AIR_R, _MB_VAULT_WALL_OUTER, _MB_VAULT_AIR_H),
        tube("mb_vair", 0.0, _MB_VAULT_AIR_R, _MB_VAULT_AIR_H),
        tube("mb_rair", 0.0, _MB_VAULT_AIR_R, _MB_ROOM_AIR_H),
        tube("mb_rwall", _MB_VAULT_AIR_R, _MB_ROOM_WALL_OUTER, _MB_ROOM_AIR_H),
        tube("mb_roof", 0.0, _MB_SLAB_R, _MB_ROOF_H),
        tube("mb_cap", 0.0, _MB_CAP_R, _MB_CAP_H),
    ])

    vols = "\n".join([
        '    <volume name="vol_mb_inner_oil"><materialref ref="MINERAL_OIL"/><solidref ref="mb_inner_oil"/></volume>',
        '    <volume name="vol_mb_veto_oil"><materialref ref="MINERAL_OIL"/><solidref ref="mb_veto_oil"/></volume>',
        '    <volume name="vol_mb_steel"><materialref ref="MB_CARBON_STEEL"/><solidref ref="mb_steel"/></volume>',
        '    <volume name="vol_mb_cap"><materialref ref="MINERAL_OIL"/><solidref ref="mb_cap"/></volume>',
        '    <volume name="vol_mb_lslab"><materialref ref="MB_CONCRETE"/><solidref ref="mb_lslab"/></volume>',
        '    <volume name="vol_mb_uslab"><materialref ref="MB_CONCRETE"/><solidref ref="mb_uslab"/></volume>',
        '    <volume name="vol_mb_vwall"><materialref ref="MB_CONCRETE"/><solidref ref="mb_vwall"/></volume>',
        '    <volume name="vol_mb_vair"><materialref ref="MB_AIR"/><solidref ref="mb_vair"/></volume>',
        '    <volume name="vol_mb_rair"><materialref ref="MB_AIR"/><solidref ref="mb_rair"/></volume>',
        '    <volume name="vol_mb_rwall"><materialref ref="MB_CONCRETE"/><solidref ref="mb_rwall"/></volume>',
        '    <volume name="vol_mb_roof"><materialref ref="MB_CONCRETE"/><solidref ref="mb_roof"/></volume>',
    ])

    # Order matters: structural concrete/air first, detector oil/steel last.
    pvs = "\n".join([
        vpv("pv_mb_lslab", "vol_mb_lslab", _MB_Y_LSLAB_BOT + 0.5 * _MB_LSLAB_H),
        vpv("pv_mb_uslab", "vol_mb_uslab", _MB_Y_CAV_TOP + 0.5 * _MB_USLAB_H),
        vpv("pv_mb_vwall", "vol_mb_vwall", cav_yc),
        vpv("pv_mb_vair", "vol_mb_vair", cav_yc),
        vpv("pv_mb_rair", "vol_mb_rair", room_yc),
        vpv("pv_mb_rwall", "vol_mb_rwall", room_yc),
        vpv("pv_mb_roof", "vol_mb_roof", _MB_Y_ROOM_TOP + 0.5 * _MB_ROOF_H),
        opv("pv_mb_inner_oil", "vol_mb_inner_oil"),
        opv("pv_mb_veto_oil", "vol_mb_veto_oil"),
        opv("pv_mb_steel", "vol_mb_steel"),
        vpv("pv_mb_cap", "vol_mb_cap", _MB_Y_CAV_TOP + 0.5 * _MB_CAP_H),
    ])

    return f"""<?xml version="1.0" encoding="UTF-8"?>
<gdml>
  <define/>
  <materials>
    <isotope N="1" Z="1" name="mb_H1"><atom unit="g/mole" value="1.008"/></isotope>
    <element name="mb_H"><fraction n="1.0" ref="mb_H1"/></element>
    <isotope N="12" Z="6" name="mb_C12"><atom unit="g/mole" value="12.0"/></isotope>
    <element name="mb_C"><fraction n="1.0" ref="mb_C12"/></element>
    <isotope N="14" Z="7" name="mb_N14"><atom unit="g/mole" value="14.0"/></isotope>
    <element name="mb_N"><fraction n="1.0" ref="mb_N14"/></element>
    <isotope N="16" Z="8" name="mb_O16"><atom unit="g/mole" value="16.0"/></isotope>
    <element name="mb_O"><fraction n="1.0" ref="mb_O16"/></element>
    <isotope N="23" Z="11" name="mb_Na23"><atom unit="g/mole" value="22.99"/></isotope>
    <element name="mb_Na"><fraction n="1.0" ref="mb_Na23"/></element>
    <isotope N="27" Z="13" name="mb_Al27"><atom unit="g/mole" value="26.98"/></isotope>
    <element name="mb_Al"><fraction n="1.0" ref="mb_Al27"/></element>
    <isotope N="28" Z="14" name="mb_Si28"><atom unit="g/mole" value="28.09"/></isotope>
    <element name="mb_Si"><fraction n="1.0" ref="mb_Si28"/></element>
    <isotope N="40" Z="20" name="mb_Ca40"><atom unit="g/mole" value="40.08"/></isotope>
    <element name="mb_Ca"><fraction n="1.0" ref="mb_Ca40"/></element>
    <isotope N="56" Z="26" name="mb_Fe56"><atom unit="g/mole" value="55.845"/></isotope>
    <element name="mb_Fe"><fraction n="1.0" ref="mb_Fe56"/></element>
    <material name="MINERAL_OIL" state="liquid">
      <D unit="g/cm3" value="{_MB_OIL_DENSITY}"/>
      <fraction n="0.1437" ref="mb_H"/>
      <fraction n="0.8563" ref="mb_C"/>
    </material>
    <material name="MB_CARBON_STEEL" state="solid">
      <D unit="g/cm3" value="7.86"/>
      <fraction n="0.99" ref="mb_Fe"/>
      <fraction n="0.01" ref="mb_C"/>
    </material>
    <material name="MB_AIR" state="gas">
      <D unit="g/cm3" value="0.001225"/>
      <fraction n="0.7562" ref="mb_N"/>
      <fraction n="0.2438" ref="mb_O"/>
    </material>
    <material name="MB_CONCRETE" state="solid">
      <D unit="g/cm3" value="2.30"/>
      <fraction n="0.52" ref="mb_O"/>
      <fraction n="0.325" ref="mb_Si"/>
      <fraction n="0.06" ref="mb_Ca"/>
      <fraction n="0.015" ref="mb_Na"/>
      <fraction n="0.04" ref="mb_Fe"/>
      <fraction n="0.04" ref="mb_Al"/>
    </material>
    <material name="MB_DIRT" state="solid">
      <D unit="g/cm3" value="1.90"/>
      <fraction n="0.5326" ref="mb_O"/>
      <fraction n="0.4674" ref="mb_Si"/>
    </material>
  </materials>
  <solids>
{solids}
  </solids>
  <structure>
{vols}
    <volume name="vol_mb_world">
      <materialref ref="MB_DIRT"/>
      <solidref ref="mb_world"/>
{pvs}
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="vol_mb_world"/>
  </setup>
</gdml>
"""


_MINIBOONE_GDML = _build_miniboone_gdml()


def ensure_miniboone_gdml(abs_dir: str,
                          filename: str = "gdml/miniboone_tank.gdml") -> str:
    """Write the MiniBooNE placeholder tank GDML if not already present.

    The file is generated locally (there is no remote URL), so this must
    run before ``_ensure_gdml_files`` sees the MiniBooNE source spec.
    Returns the relative *filename* (matching the _DETECTOR_SPECS entry).
    """
    path = os.path.join(abs_dir, filename)
    if not os.path.isfile(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + ".tmp"
        try:
            with open(tmp, "w", encoding="utf-8") as f:
                f.write(_MINIBOONE_GDML)
            os.replace(tmp, path)
        except Exception:
            if os.path.exists(tmp):
                os.remove(tmp)
            raise
    return filename


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

    # GDML <box> x/y/z are FULL widths (GDMLSolidParser ParseBox). The
    # placement math below works in half-extents, so the box attributes are
    # written as 2x these values. Site volume: 1600 m x 400 m x 1800 m.
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
    <box name="sol_site_volume" lunit="m" x="{2.0 * box_half_x}" y="{2.0 * box_half_y}" z="{2.0 * box_half_z}"/>
    <box name="sol_atmosphere" lunit="m" x="{2.0 * child_half_x}" y="{2.0 * atmo_half_height}" z="{2.0 * child_half_z}"/>
    <box name="sol_dolomite_bedrock" lunit="m" x="{2.0 * child_half_x}" y="{2.0 * bedrock_half_height}" z="{2.0 * child_half_z}"/>
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
