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
    https://github.com/SBNSoftware/G4BNB/blob/HEAD/src/NuBeamOutput.cc
  The ICARUS z = 600 m on-axis baseline is published (peer-reviewed): ICARUS
  Collaboration, Eur. Phys. J. C 83 (2023) 467 [arXiv:2301.08634] ("exposed
  at a 600 m baseline distance to the ... BNB"), and the SBN Proposal
  [arXiv:1503.01520] Sec. I.2 + Table 1 ("600 m from the BNB target"; all
  three detectors "on-axis in the BNB"). Also confirmed by SBN DocDB
  22998-v2 (restricted). 600 m is a round as-designed value: G4BNB encodes
  exactly 60000 cm with no survey-refinement comment (contrast MiniBooNE's
  54134 cm "refined by Zarko"); a finer surveyed value, if one exists, is in
  the Fermilab Nov-2024 BNB site survey (noted in arXiv:2501.06323), not public.
  SBND x = +73.78 cm off-axis offset confirmed by SBN DocDB 20891 (flux
  config H); see the SBND_LArSoft -> BNB edge below for the full provenance.
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

# The BNB frame origin coincides with the SAND world volume center
# in the BooNE GDML. G4BNB bsim::Location detector positions and
# dk2nu flux ray origins are all expressed relative to this point.
# The BNB target sits ~39 cm downstream (z = +0.3905 m in SAND).
BNB_TARGET_Z_SAND_M = 0.0

T_MiniBooNE_local = np.array([0.0, 1.89614, 541.34])

# ======================================================================
# Build the frame graph
# ======================================================================
def _build_graph() -> FrameGraph:
    g = FrameGraph()

    g.add_frame(Frame("BNB",
        "Booster Neutrino Beam global frame (g4bnb bsim::Location).",
        "G4BNB SAND world volume center (target is ~39 cm downstream)",
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

    g.add_frame(Frame("MiniBooNE_local",
        "MiniBooNE detector local frame (axes = BNB axes; origin at tank center).",
        "MiniBooNE spherical-tank center",
        "horizontal west", "up", "along BNB beam axis"))

    # ICARUS LArSoft -> NuMI (from icaruscode GNuMIFlux.xml, confirmed by DocDB 22998)
    g.add_transform(Transform(
        "ICARUS_LArSoft", "NuMI", R_ICARUS_TO_NUMI, T_ICARUS_IN_NUMI_M,
        "icaruscode GNuMIFlux.xml beamdir/beampos; confirmed by SBN DocDB 22998-v2"))

    # ICARUS LArSoft -> BNB: pure translation (0, 0, 600) m.
    # From G4BNB NuBeamOutput.cc:143:
    #   bsim::Location(0., 0., 60000., "T600")  [cm]
    # i.e. on-axis (x = 0, y = 0), 600 m downstream. This is the round
    # as-designed baseline, NOT a survey: the T600 tuple carries no
    # refinement comment, unlike MiniBooNE's 54134 cm ("refined by Zarko")
    # on NuBeamOutput.cc:136. The 600 m on-axis placement is published
    # (peer-reviewed): ICARUS Collaboration, Eur. Phys. J. C 83 (2023) 467
    # [arXiv:2301.08634] ("exposed at a 600 m baseline distance to the BNB"),
    # and the SBN Proposal [arXiv:1503.01520] Sec. I.2 + Table 1 ("600 m from
    # the BNB target"; all three SBN detectors "on-axis in the BNB"); the SBN
    # program overview [arXiv:2203.05814] repeats "Icarus at 600 m".
    # The (0, 0, 600) places the ICARUS LArSoft origin (TPC-system center) on
    # the beam axis; the active-LAr center offset is carried separately in
    # DETECTORS["ICARUS"] (_ic_y = -0.202 m; cryostats at x = +/-2.1 m).
    g.add_transform(Transform.translation(
        "ICARUS_LArSoft", "BNB", [0.0, 0.0, 600.0],
        "G4BNB T600 Location (0, 0, 60000) cm [NuBeamOutput.cc:143]; round "
        "as-designed 600 m on-axis baseline published in ICARUS Coll. EPJC 83 "
        "(2023) 467 [arXiv:2301.08634] and SBN Proposal [arXiv:1503.01520] "
        "Sec. I.2; confirmed by SBN DocDB 22998-v2"))

    # SBND LArSoft -> BNB (G4BNB bsim::Location).
    # From G4BNB NuBeamOutput.cc:145:
    #   bsim::Location(73.78, 0., 11000., "SBND")  [cm]
    # i.e. off-axis (x = +0.7378 m), on-axis vertically (y = 0), 110 m
    # downstream. This is the SBND coordinate origin expressed in the BNB
    # beam frame; we store it with that sign (detector-in-beam-frame), which
    # matches G4BNB and sbncode bnb_kaon_sbnd.fcl ("... in detector in beam
    # frame" = [73.78, 0, 11000]). The SBND-side view of the same offset is
    # the beam x-center at -73.78 cm in detector coordinates (sbncode
    # run_fluxreader_sbnd.fcl XShift = -73.78; SBND flux-files wiki) -- the
    # inverse of this edge, not a discrepancy.
    #
    # Provenance: the 73.78 cm horizontal offset is the as-designed SBND
    # position; the authoritative reference is SBN-DocDB 20891. It superseded
    # two earlier wrong values -- a uBooNE-frame ~1.3 m, then 45.7 cm in flux
    # configs F/G. The corrected -73.78 cm landed in flux config H (produced
    # by Z. Pavlovic, sbndcode PR #95, ~2021); see the SBND flux-files wiki
    # and Release Notes 09.21.00. G4BNB's OWN default stayed on-axis
    # (0, 0, 11000) until NuBeamOutput.cc:145 was fixed in 2025 (G4BNB commit
    # bdd5f29e). The 11000 cm = 110 m target-to-TPC baseline is current (an
    # early flux config A used 100 m).
    # Caveats: the offset is to the SBND LArSoft origin (cathode plane), not
    # the active-LAr center (carried separately in DETECTORS["SBND"]); and
    # 73.78 cm is a design value (DocDB 20891), not a published survey like
    # MiniBooNE's 1.9 m.
    g.add_transform(Transform.translation(
        "SBND_LArSoft", "BNB", [0.7378, 0.0, 110.0],
        "G4BNB SBND Location (73.78, 0, 11000) cm [NuBeamOutput.cc:145]; "
        "as-designed off-axis offset per SBN-DocDB 20891 (corrected from "
        "45.7 cm to -73.78 cm in flux config H, Z. Pavlovic, sbndcode PR #95)"))

    # MicroBooNE LArSoft -> BNB (G4BNB bsim::Location)
    g.add_transform(Transform.translation(
        "MicroBooNE_LArSoft", "BNB", [0.0, 0.0, 470.0],
        "G4BNB MicroBooNE Location (0, 0, 47000) cm"))

    # MiniBooNE tank center -> BNB (G4BNB bsim::Location).
    # From G4BNB NuBeamOutput.cc:136:
    #   bsim::Location(0., 189.614, 54134.0, "MiniBooNE")  [cm]
    # i.e. on-axis (x = 0), 1.89614 m vertical, 541.34 m downstream.
    #
    # Source: both the vertical offset and the baseline are PUBLISHED in the
    # MiniBooNE flux paper, PRD 79 (2009) 072002 [arXiv:0806.1449], Sec. II:
    # "The axis of the beam, defined by the center of the decay pipe, is
    # displaced vertically from the center of the MiniBooNE detector by 1.9
    # meters" (plus the 541 m target-to-detector-center baseline). The exact
    # 189.614 / 54134 cm digits are the survey-precision form of that published
    # 1.9 m / 541 m -- no finer value is in print. They were carried from the
    # MiniBooNE beam MC (NuBeamOutput.cc:136 note "Refined by Zarko" =
    # Z. Pavlovic, MiniBooNE flux physicist). The underlying positions come
    # from the Fermilab Alignment & Metrology survey (methodology in Oshinowo,
    # FERMILAB-CONF-02-425); detailed BNB geometry is in the 8 GeV Beam TDR
    # (the PRD's ref [6]).
    # Sign: stored as +y (BNB +y = up, per the BNB frame above); the PRD gives
    # only the magnitude, so the up/down sense is not independently pinned here.
    g.add_transform(Transform.translation(
        "MiniBooNE_local", "BNB", T_MiniBooNE_local,
        "1.9 m vertical offset + 541 m baseline published in PRD 79 (2009) "
        "072002 Sec. II [arXiv:0806.1449]; 189.614/54134 cm is the survey-"
        "precision form from G4BNB NuBeamOutput.cc:136 (refined by Z. Pavlovic)"))

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

    # MiniBooNE: spherical mineral-oil tank, inner radius 6.1 m. The native
    # frame is MiniBooNE_local (origin at the tank center); the tank position
    # in the BNB frame is carried by the MiniBooNE_local -> BNB edge above, so
    # center_native is the origin and the active volume is the tank bounding
    # box. This matches how the LArSoft detectors are handled (GDML in the
    # native frame, placed via the frame transform).
    _mb_tank_r = 6.1
    detectors["MiniBooNE"] = Detector(
        "MiniBooNE", "MiniBooNE_local",
        np.array([0.0, 0.0, 0.0]),
        np.array([-_mb_tank_r, -_mb_tank_r, -_mb_tank_r]),
        np.array([+_mb_tank_r, +_mb_tank_r, +_mb_tank_r]))

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
