"""
SBN geometry: coordinate frames and rigid transforms.

Minimal port of the transform graph from g4bnb/sbn_geometry.py into SIREN.
Provides composable frame-to-frame transforms for the Fermilab SBN program
detectors (ICARUS, SBND, MicroBooNE) and beamlines (BNB, NuMI).

World coordinate system (BNB frame):
  Origin : BNB beryllium target at MI-12
  x      : horizontal west
  y      : up
  z      : along BNB beam axis (approximately north)
  Units  : meters

References:
  ICARUS-NuMI rotation matrix and translation vector from the production
  GNuMIFlux.xml in icaruscode (public):
    https://github.com/SBNSoftware/icaruscode/blob/develop/fcl/gen/numi/GNuMIFlux.xml
  Confirmed by SBN DocDB 22998-v2 (A. Aduszkiewicz, 2022-10-26; restricted).

  Detector positions in BNB frame from G4BNB bsim::Location tuples (public):
    https://github.com/SBNSoftware/G4BNB
  Confirmed by SBN DocDB 22998-v2 for ICARUS z = 600 m.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class Frame:
    name: str
    description: str
    origin: str
    x_axis: str
    y_axis: str
    z_axis: str


@dataclass(frozen=True)
class Transform:
    """Rigid transform: r_target = R @ r_source + t (meters)."""
    source: str
    target: str
    R: np.ndarray
    t: np.ndarray
    reference: str = ""

    def __post_init__(self):
        object.__setattr__(self, "R", np.array(self.R, dtype=float))
        object.__setattr__(self, "t", np.array(self.t, dtype=float))
        self.R.flags.writeable = False
        self.t.flags.writeable = False

    def apply(self, r: np.ndarray) -> np.ndarray:
        return self.R @ np.asarray(r, dtype=float) + self.t

    def inverse(self) -> Transform:
        Rt = self.R.T
        return Transform(self.target, self.source, Rt, -Rt @ self.t,
                         self.reference + " (inverted)")

    def compose(self, other: Transform) -> Transform:
        return Transform(self.source, other.target,
                         other.R @ self.R, other.R @ self.t + other.t,
                         f"({self.reference}) then ({other.reference})")

    @classmethod
    def identity(cls, name: str) -> Transform:
        return cls(name, name, np.eye(3), np.zeros(3), "identity")

    @classmethod
    def translation(cls, src: str, tgt: str, t: Sequence[float],
                    reference: str = "") -> Transform:
        return cls(src, tgt, np.eye(3), np.asarray(t, dtype=float), reference)


class FrameGraph:
    def __init__(self) -> None:
        self.frames: dict[str, Frame] = {}
        self.edges: list[Transform] = []

    def add_frame(self, frame: Frame) -> None:
        self.frames[frame.name] = frame

    def add_transform(self, xform: Transform) -> None:
        self.edges.append(xform)

    def compose(self, src: str, dst: str) -> Transform:
        if src == dst:
            return Transform.identity(src)
        adjacency: dict[str, list[Transform]] = {}
        for e in self.edges:
            adjacency.setdefault(e.source, []).append(e)
            adjacency.setdefault(e.target, []).append(e.inverse())
        queue: deque[tuple[str, list[Transform]]] = deque([(src, [])])
        visited = {src}
        while queue:
            current, partial = queue.popleft()
            for edge in adjacency.get(current, []):
                if edge.target in visited:
                    continue
                new_partial = partial + [edge]
                if edge.target == dst:
                    out = new_partial[0]
                    for nxt in new_partial[1:]:
                        out = out.compose(nxt)
                    return out
                visited.add(edge.target)
                queue.append((edge.target, new_partial))
        raise ValueError(f"No transform path: {src} -> {dst}")

    def convert(self, r: Iterable[float], src: str, dst: str) -> np.ndarray:
        return self.compose(src, dst).apply(np.asarray(r, dtype=float))


# ======================================================================
# ICARUS LArSoft -> NuMI rotation and translation
#
# Source: icaruscode GNuMIFlux.xml <beamdir>/<beampos> elements
#   https://github.com/SBNSoftware/icaruscode/blob/develop/fcl/gen/numi/GNuMIFlux.xml
# Confirmed by SBN DocDB 22998-v2 Eq. (3) (restricted).
# ======================================================================
R_ICARUS_TO_NUMI = np.array([
    [ 0.921035925,  0.0,         -0.389477631],
    [ 0.022715103,  0.998297825,  0.053716629],
    [ 0.388814672, -0.058321970,  0.919468161],
])
T_ICARUS_IN_NUMI_M = np.array([4.503730, 80.153901, 795.112945])

# Z-position of the BNB target center in the BooNE GDML (SAND) frame,
# traced through the volume hierarchy:
#   SAND -> CAVE (+426.8 mm) -> CAVI (0) -> PILE (~0) -> PICV (-426.8 mm) -> TARG (+390.525 mm)
# Total: 390.525 mm = 0.390525 m
BNB_TARGET_Z_SAND_M = 0.390525


# ======================================================================
# Build the frame graph
# ======================================================================
def _build_graph() -> FrameGraph:
    g = FrameGraph()

    g.add_frame(Frame("BNB",
        "Booster Neutrino Beam global frame (g4bnb bsim::Location).",
        "BNB beryllium target at MI-12",
        "horizontal west", "up", "along BNB beam axis (north)"))

    # The NuMI frame origin is MCZERO (Horn 1 upstream face). The GDML
    # world volume for the NuMI beamline must be centered at this same
    # point so that the physvol placement T_numi.apply([0,0,0]) is correct.
    g.add_frame(Frame("NuMI",
        "NuMI beam frame (GNuMIFlux.xml convention in icaruscode).",
        "NuMI Horn 1 upstream face (MCZERO)",
        "horizontal southwest", "up",
        "along NuMI beam axis (northwest, slightly downward)"))

    g.add_frame(Frame("ICARUS_LArSoft",
        "ICARUS detector LArSoft World frame (axes = BNB axes).",
        "center of ICARUS TPC system",
        "horizontal west", "up", "along BNB beam axis (north)"))

    g.add_frame(Frame("SBND_LArSoft",
        "SBND detector LArSoft World frame (axes assumed = BNB axes).",
        "SBND LArSoft World origin",
        "horizontal west", "up", "along BNB beam axis"))

    g.add_frame(Frame("MicroBooNE_LArSoft",
        "MicroBooNE detector LArSoft World frame (axes assumed = BNB axes).",
        "MicroBooNE LArSoft World origin",
        "horizontal west", "up", "along BNB beam axis"))

    # ICARUS LArSoft -> NuMI (from icaruscode GNuMIFlux.xml, confirmed by DocDB 22998)
    g.add_transform(Transform(
        "ICARUS_LArSoft", "NuMI", R_ICARUS_TO_NUMI, T_ICARUS_IN_NUMI_M,
        "icaruscode GNuMIFlux.xml beamdir/beampos; confirmed by SBN DocDB 22998-v2"))

    # ICARUS LArSoft -> BNB: pure translation (0, 0, 600) m
    # G4BNB bsim::Location for ICARUS; confirmed by DocDB 22998.
    g.add_transform(Transform.translation(
        "ICARUS_LArSoft", "BNB", [0.0, 0.0, 600.0],
        "G4BNB ICARUS Location (0, 0, 60000) cm; confirmed by SBN DocDB 22998-v2"))

    # SBND LArSoft -> BNB (G4BNB bsim::Location)
    g.add_transform(Transform.translation(
        "SBND_LArSoft", "BNB", [0.7378, 0.0, 110.0],
        "G4BNB SBND Location (73.78, 0, 11000) cm"))

    # MicroBooNE LArSoft -> BNB (G4BNB bsim::Location)
    g.add_transform(Transform.translation(
        "MicroBooNE_LArSoft", "BNB", [0.0, 0.0, 470.0],
        "G4BNB MicroBooNE Location (0, 0, 47000) cm"))

    return g


GRAPH: FrameGraph = _build_graph()


# ======================================================================
# Detector active-volume data
# ======================================================================
@dataclass(frozen=True)
class Detector:
    name: str
    native_frame: str
    center_native: np.ndarray
    active_min: np.ndarray
    active_max: np.ndarray

    def __post_init__(self):
        for field in ("center_native", "active_min", "active_max"):
            arr = np.array(getattr(self, field), dtype=float)
            arr.flags.writeable = False
            object.__setattr__(self, field, arr)

    def active_size(self) -> np.ndarray:
        return self.active_max - self.active_min


def _build_detectors() -> dict[str, Detector]:
    detectors: dict[str, Detector] = {}

    _ic_y = -0.202
    _ic_hh = 1.58
    _ic_hl = 8.975
    detectors["ICARUS"] = Detector(
        "ICARUS", "ICARUS_LArSoft",
        np.array([0.0, _ic_y, 0.0]),
        np.array([-3.60, _ic_y - _ic_hh, -_ic_hl]),
        np.array([+3.60, _ic_y + _ic_hh, +_ic_hl]))

    for tag, xc in (("ICARUS_C0", -2.10215), ("ICARUS_C1", +2.10215)):
        detectors[tag] = Detector(
            tag, "ICARUS_LArSoft",
            np.array([xc, _ic_y, 0.0]),
            np.array([xc - 1.50, _ic_y - _ic_hh, -_ic_hl]),
            np.array([xc + 1.50, _ic_y + _ic_hh, +_ic_hl]))

    # SBND LArSoft World origin is at the cathode plane (x=0).
    # The LAr volume is symmetric in x (two drift volumes), but offset
    # in y and z relative to the LArSoft origin.
    #
    # LAr extent measured from sbnd_v02_06.gdml (1cm resolution):
    #   x: -5.19 to +5.19 m  (cathode at x=0, anodes at ~x=+/-2m,
    #                          inactive LAr fills the rest of the cryostat)
    #   y: -3.50 to +2.32 m  (center at y=-0.59)
    #   z: -1.39 to +7.22 m  (center at z=+2.92)
    #
    # The center_native below is the geometric center of the LAr volume,
    # which is offset from the LArSoft origin. This is used as the
    # DetectorCoordinates origin so that injectors aiming at (0,0,0)
    # hit the center of the detector.
    detectors["SBND"] = Detector(
        "SBND", "SBND_LArSoft",
        np.array([0.0, -0.59, 2.92]),
        np.array([-5.19, -3.50, -1.39]),
        np.array([+5.19, +2.32, +7.22]))

    detectors["MicroBooNE"] = Detector(
        "MicroBooNE", "MicroBooNE_LArSoft",
        np.array([0.0, 0.0, 0.0]),
        np.array([-1.281750, -1.165, -5.184]),
        np.array([+1.281750, +1.165, +5.184]))

    detectors["MiniBooNE"] = Detector(
        "MiniBooNE", "BNB",
        np.array([0.0, 1.89614, 541.34]),
        np.array([0.0, 1.89614, 541.34]) - 6.1,
        np.array([0.0, 1.89614, 541.34]) + 6.1)

    return detectors


DETECTORS: dict[str, Detector] = _build_detectors()


# ======================================================================
# Public API
# ======================================================================
def convert(r: Iterable[float], src: str, dst: str) -> np.ndarray:
    """Convert a position (meters) from frame src to frame dst."""
    return GRAPH.convert(r, src, dst)


def transform(src: str, dst: str) -> Transform:
    """Return the composed Transform from src to dst."""
    return GRAPH.compose(src, dst)


def detector_center(name: str, frame: str) -> np.ndarray:
    """Active-volume geometric center in the given frame (meters)."""
    d = DETECTORS[name]
    return GRAPH.convert(d.center_native, d.native_frame, frame)


def detector_origin(name: str, frame: str) -> np.ndarray:
    """Native-frame origin of a detector, expressed in the given frame (meters)."""
    d = DETECTORS[name]
    return GRAPH.convert(np.zeros(3), d.native_frame, frame)


def detector_transform(name: str, frame: str) -> Transform:
    """Full rigid transform from a detector's native frame to the given frame."""
    d = DETECTORS[name]
    return GRAPH.compose(d.native_frame, frame)



def gdml_rotation_angles(R: np.ndarray):
    """Decompose rotation matrix to GDML Euler angles (rx, ry, rz) in radians.

    GDML convention is extrinsic XYZ (static frame): R = Rz(rz) @ Ry(ry) @ Rx(rx).
    This matches the SIREN GDML parser which calls QFromXYZs(rx, ry, rz).
    """
    from siren.math import Quaternion, Matrix3D
    q = Quaternion()
    q.SetMatrix(Matrix3D(*R.flatten()))
    return q.GetEulerAnglesXYZs()
