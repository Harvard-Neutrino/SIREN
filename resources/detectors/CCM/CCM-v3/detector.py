"""Corrected target_sim facility with CCM's six-layer cylindrical detector.

Coordinates supplied to DetectorModel queries are detector-local metres.
The standalone GDML instead uses the target_sim facility/world frame.
See README.md for the provisional placement and retained approximations.
"""

import json
import math
from pathlib import Path

_DIRECTORY = Path(__file__).resolve().parent


def geometry_description():
    """Return a fresh copy of the geometry inputs and placement assumptions."""
    return json.loads((_DIRECTORY / "geometry.json").read_text())


def gdml_path():
    """Return the complete, directly viewable facility-plus-detector GDML."""
    return str(_DIRECTORY / "ccm-facility.gdml")


def fiducial_volume():
    """Return the legacy active-envelope cylinder in detector coordinates."""
    from siren.geometry import Cylinder
    spec = geometry_description()["fiducial"]
    return Cylinder(spec["radius_m"], 0.0, spec["height_m"])


def load_detector():
    """Load CCM-v3 with a strict GDML import and its detector-local frame."""
    from siren.detector import DetectorModel, GeometryPosition
    from siren.math import Matrix3D, Quaternion, Vector3D

    placement = geometry_description()["placement"]
    model = DetectorModel()
    model.LoadGDML(gdml_path(), True)
    model.DetectorOrigin = GeometryPosition(Vector3D(*placement["center_facility_m"]))
    yaw = placement["detector_to_facility_yaw_rad"]
    c, s = math.cos(yaw), math.sin(yaw)
    rotation = Quaternion()
    rotation.SetMatrix(Matrix3D(c, -s, 0, s, c, 0, 0, 0, 1))
    model.DetectorRotation = rotation
    return model
