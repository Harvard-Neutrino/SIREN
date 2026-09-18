"""Rebuild the bundled CCM GDML from facility.gdml and geometry.json.

Uses only the Python standard library. Run this file after editing the inputs;
normal detector loading reads the prebuilt GDML without writing any cache.
"""

import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET


def build(directory=None):
    root_dir = Path(directory) if directory is not None else Path(__file__).parent
    config = json.loads((root_dir / "geometry.json").read_text())
    facility = (root_dir / "facility.gdml").read_bytes()
    if hashlib.sha256(facility).hexdigest() != config["facility_sha256"]:
        raise ValueError("facility.gdml does not match the recorded input hash")
    ET.register_namespace("xsi", "http://www.w3.org/2001/XMLSchema-instance")
    root = ET.fromstring(facility)
    materials, solids, structure = (root.find(key) for key in
                                   ("materials", "solids", "structure"))
    world_ref = root.find("setup/world").attrib["ref"]
    world = next(v for v in structure if v.get("name") == world_ref)
    if any(e.get("name", "").startswith("ccm_v3_") for e in root.iter()):
        raise ValueError("facility input already contains a CCM-v3 insert")

    def add(parent, tag, **attributes):
        return ET.SubElement(parent, tag, {k: str(v) for k, v in attributes.items()})

    for code, mass in config["isotope_molar_masses_g_mol"].items():
        pdg = int(code)
        z, a = (pdg // 10000) % 1000, (pdg // 10) % 1000
        name = "ccm_v3_isotope_" + code
        isotope = add(materials, "isotope", name=name, Z=z, N=a)
        add(isotope, "atom", value=mass, unit="g/mole")
        element = add(materials, "element", name="ccm_v3_element_" + code)
        add(element, "fraction", n=1, ref=name)

    densities = {layer["material"]: layer["density_g_cm3"]
                 for layer in config["layers"]}
    densities["AIR"] = config["vacuum_density_g_cm3"]
    for name, spec in config["materials"].items():
        material = add(materials, "material", name="ccm_v3_" + name,
                       state="gas" if name == "AIR" else
                       "liquid" if name == "ARGON" else "solid")
        add(material, "D", value=densities[name], unit="g/cm3")
        for code, fraction in spec["fractions"].items():
            add(material, "fraction", n=fraction, ref="ccm_v3_element_" + code)

    # Nest the six cylinders: each daughter displaces its mother's material.
    # Reverse declaration order keeps every daughter defined before its parent.
    daughter = None
    for layer in reversed(config["layers"]):
        name = layer["name"]
        solid = "ccm_v3_" + name + "_solid"
        add(solids, "tube", name=solid, rmin=0, rmax=layer["radius_m"],
            z=layer["height_m"], startphi=0, deltaphi=360, aunit="deg", lunit="m")
        volume = add(structure, "volume", name=name)
        add(volume, "materialref", ref="ccm_v3_" + layer["material"])
        add(volume, "solidref", ref=solid)
        if daughter is not None:
            pv = add(volume, "physvol", name=daughter)
            add(pv, "volumeref", ref=daughter)
        daughter = name
    pv = add(world, "physvol", name=daughter)
    add(pv, "volumeref", ref=daughter)
    x, y, z = config["placement"]["center_facility_m"]
    add(pv, "position", name="ccm_v3_center", unit="m", x=x, y=y, z=z)
    # GDML physical-volume rotations are passive; DetectorRotation is active.
    yaw = config["placement"]["detector_to_facility_yaw_rad"]
    add(pv, "rotation", name="ccm_v3_rotation", unit="rad", x=0, y=0, z=-yaw)
    # The world must be declared after the new daughter volumes for Geant4.
    structure.remove(world)
    structure.append(world)
    ET.indent(root, space="  ")
    output = ET.tostring(root, encoding="utf-8", xml_declaration=True) + b"\n"
    (root_dir / "ccm-facility.gdml").write_bytes(output)

    # Compatibility metadata for utilities.get_fiducial_volume. This is not a
    # transport model; intentionally no materials.dat accompanies this file.
    fiducial = config["fiducial"]
    (root_dir / "densities.dat").write_text(
        "# Fiducial metadata only. Load transport with load_detector('CCM-v3').\n"
        "# The complete transport model is ccm-facility.gdml.\n"
        f"detector {x:.17g} {y:.17g} {z:.17g} 0 0 {yaw:.17g}\n"
        "fiducial detector_coords cylinder 0 0 0 0 0 0 "
        f"{fiducial['radius_m']:.17g} 0 {fiducial['height_m']:.17g}\n")
    return root_dir / "ccm-facility.gdml"


if __name__ == "__main__":
    print(build())
