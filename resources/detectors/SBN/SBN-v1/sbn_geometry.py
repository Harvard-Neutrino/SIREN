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

BNB_TARGET_Z_SAND_M = 0.81733


# ======================================================================
# Build the frame graph
# ======================================================================
def _build_graph() -> FrameGraph:
    g = FrameGraph()

    g.add_frame(Frame("BNB",
        "Booster Neutrino Beam global frame (g4bnb bsim::Location).",
        "BNB beryllium target at MI-12",
        "horizontal west", "up", "along BNB beam axis (north)"))

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

    detectors["SBND"] = Detector(
        "SBND", "SBND_LArSoft",
        np.array([0.0, 0.0, 0.0]),
        np.array([-2.0, -2.0, -2.5]),
        np.array([+2.0, +2.0, +2.5]))

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


def quaternion_from_matrix(R: np.ndarray) -> tuple[float, float, float, float]:
    """Convert a 3x3 rotation matrix to quaternion (x, y, z, w).

    Uses the Shepperd method for numerical stability across all rotations."""
    trace = float(R[0, 0] + R[1, 1] + R[2, 2])
    if trace > 0:
        s = 0.5 / math.sqrt(trace + 1.0)
        w = 0.25 / s
        x = (R[2, 1] - R[1, 2]) * s
        y = (R[0, 2] - R[2, 0]) * s
        z = (R[1, 0] - R[0, 1]) * s
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        s = 2.0 * math.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2])
        w = (R[2, 1] - R[1, 2]) / s
        x = 0.25 * s
        y = (R[0, 1] + R[1, 0]) / s
        z = (R[0, 2] + R[2, 0]) / s
    elif R[1, 1] > R[2, 2]:
        s = 2.0 * math.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2])
        w = (R[0, 2] - R[2, 0]) / s
        x = (R[0, 1] + R[1, 0]) / s
        y = 0.25 * s
        z = (R[1, 2] + R[2, 1]) / s
    else:
        s = 2.0 * math.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1])
        w = (R[1, 0] - R[0, 1]) / s
        x = (R[0, 2] + R[2, 0]) / s
        y = (R[1, 2] + R[2, 1]) / s
        z = 0.25 * s
    return float(x), float(y), float(z), float(w)


def gdml_rotation_angles(R: np.ndarray):
    """Decompose rotation matrix to GDML Euler angles (rx, ry, rz) in radians.
    Convention: R = Rz(rz) @ Ry(ry) @ Rx(rx)."""
    sy = float(np.clip(R[0, 2], -1.0, 1.0))
    ry = math.asin(sy)
    cy = math.cos(ry)
    if abs(cy) > 1e-10:
        rx = math.atan2(-R[1, 2], R[2, 2])
        rz = math.atan2(-R[0, 1], R[0, 0])
    else:
        rx = math.atan2(R[1, 0], R[1, 1])
        rz = 0.0
    return rx, ry, rz
