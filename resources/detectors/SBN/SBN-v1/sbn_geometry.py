"""
SBN geometry: coordinate frames and rigid transforms.

Provides composable frame-to-frame transforms for the Fermilab SBN program
detectors (ICARUS, SBND, MicroBooNE) and beamlines (BNB, NuMI).

World coordinate system (BNB frame):
  Origin : G4BNB SAND world volume center (the Be target sits ~39 cm
           downstream, at z = +0.3905 m; bsim::Location positions and
           dk2nu flux rays use this same SAND-center reference)
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

# z position (BNB frame) at which the BooNE GDML SAND world center is placed.
# The BNB frame origin is DEFINED as the SAND world center -- the reference
# that G4BNB bsim::Location detector positions and dk2nu flux ray origins
# use -- so no shift is applied. The physical Be target sits ~39 cm
# downstream of this origin (z = +0.3905 m); positions in this frame are
# NOT target-relative.
SAND_CENTER_Z_BNB_M = 0.0

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

    # MiniBooNE tank center -> BNB (G4BNB bsim::Location, refined survey).
    # From G4BNB NuBeamOutput.cc:136:
    #   bsim::Location(0., 189.614, 54134.0, "MiniBooNE")  [cm]
    # i.e. on-axis (x = 0), +1.89614 m vertical, 541.34 m downstream.
    # PRD 79 (2009) 072002 quotes the nominal 541 m baseline; this tuple
    # is the beam-group's refined value (Z. Pavlovic).
    g.add_transform(Transform.translation(
        "MiniBooNE_local", "BNB", T_MiniBooNE_local,
        "G4BNB MiniBooNE bsim::Location (0, 189.614, 54134) cm; NuBeamOutput.cc:136"))

    # ==================================================================
    # DUNE / LBNF beamline + near detectors
    #
    # Authoritative beam<->detector transforms from the DUNE GNuMIFlux
    # configs -- the same kind of dk2nu/GENIE flux-driver artifact as the
    # icaruscode GNuMIFlux.xml used for NuMI above. GNuMIFlux <beampos>
    # "( user ) = ( beam )" gives one physical point in both frames;
    # <beamdir> is the user->beam rotation (here a single rotation about x
    # = the beam declination). So r_beam = Rx(theta) @ r_user + t.
    #
    # NOTE: these anchor the DUNE detectors to their *beams*. The 2x2
    # (ProtoDUNE-ND) is in the NuMI beam (MINOS/MINERvA hall) -> reaches
    # BNB via the existing NuMI edge. The full DUNE_ND is in the LBNF beam;
    # since no single flux artifact crosses experiments, the LBNF->BNB site
    # tie (LBNF vs BNB/NuMI target offset + azimuth) is instead derived below
    # from the Fermilab MI-10/MI-60 survey bridged through NuMI.
    # ==================================================================
    def _Rx(theta):
        c, s = math.cos(theta), math.sin(theta)
        return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])

    # LBNF beam frame: MCZERO at the target / Horn 1, +z along the
    # (horizontal-built) beam axis, +y up. Matches the g4lbnf Tunnel world.
    g.add_frame(Frame("LBNF",
        "LBNF beam frame (g4lbnf Tunnel world; MCZERO at target).",
        "LBNF target / Horn 1 (MCZERO)",
        "horizontal", "up", "along LBNF beam axis (downstream, built horizontal)"))

    # DUNE ND geometry frame (edep-sim / dunendggd world); level hall, the
    # beam enters dipping 101 mrad.
    g.add_frame(Frame("DUNE_ND",
        "DUNE near detector geometry frame (edep-sim / dunendggd world).",
        "DUNE ND hall reference", "beam-left", "up", "downstream (level)"))

    # DUNE_ND -> LBNF. GNuMIFlux DUNEND: beamdir Rx(-0.101); beampos
    # user(0, 0.05387, 6.66) m = beam(0, 0, 562.1179) m. Places the ND on
    # the beam axis ~555 m downstream (~0.73 m below); the target sits
    # ~552.6 m upstream and ~56.7 m above the ND (562 * 0.101 = 56.8 m check).
    _R_nd_lbnf = _Rx(-0.101)
    _t_nd_lbnf = (np.array([0.0, 0.0, 562.1179])
                  - _R_nd_lbnf @ np.array([0.0, 0.05387, 6.66]))
    g.add_transform(Transform(
        "DUNE_ND", "LBNF", _R_nd_lbnf, _t_nd_lbnf,
        "DUNE/ND_Production toolbox/generator/GNuMIFlux.xml param_set DUNEND "
        "(beamdir rotation x -0.101 rad; beampos (0,0.05387,6.66)=(0,0,562.1179) m)"))

    # 2x2 ProtoDUNE-ND: in the NuMI beam in the MINOS/MINERvA near hall
    # (NOT LBNF). Ties to NuMI -> BNB (existing), so the 2x2 prototype is
    # placeable in the BNB world. GNuMIFlux ProtoDUNE-ND: beamdir
    # Rx(-0.0582977560) (= NuMI 58.3 mrad declination); beampos
    # (0,0,0)=(0,0,1036.48837) m (2x2 origin on the NuMI axis, 1036.49 m
    # downstream of the NuMI target; matches g4lbne zdetNear[1]=1036.49).
    g.add_frame(Frame("ProtoDUNE_ND_2x2",
        "DUNE 2x2 ArgonCube prototype frame (MINOS/MINERvA near hall, NuMI beam).",
        "2x2 detector center", "beam-left", "up", "downstream (level)"))
    g.add_transform(Transform(
        "ProtoDUNE_ND_2x2", "NuMI", _Rx(-0.0582977560),
        np.array([0.0, 0.0, 1036.48837]),
        "DUNE/2x2_sim run-genie/flux/GNuMIFlux.xml param_set ProtoDUNE-ND "
        "(beamdir rotation x -0.0582977560 rad; beampos (0,0,0)=(0,0,1036.48837) m)"))

    # LBNF -> BNB: closes the DUNE system into the SBN/BNB world.
    #
    # Rotation: the LBNF beam frame (g4lbnf Tunnel world -- +z = beam axis,
    # built horizontal; +y up; +x = up x beam) pointed along the real LBNF
    # beam: true azimuth 287.79 deg (toward SURF) and 101 mrad downward.
    # Columns of R are the LBNF axes expressed in BNB (x=west, y=up, z=north).
    #
    # Translation: the LBNF target (MCZero = g4lbnf origin) in BNB. Solved from
    # the Fermilab-survey MI-10 / MI-60 straight-section centers (FSCS, ft):
    #   MI-10 = (E 99673.5943, N 97303.9630),  MI-60 = (E 101513.8611, N 97129.4878).
    # Projected to ENU (FSCS Y-axis geodetic azimuth alpha = 38 16 48.0143") and
    # bridged into BNB through NuMI: MI-60 + 350 m lands on the graph's NuMI
    # point; MI-10 + 328 m -> the LBNF target. Validated -- MI-10 tangent + 7.2
    # deg extraction = 288.1 deg ~ 287.79; the survey->BNB chain reproduces the
    # NuMI matrix to ~1 deg. Horizontal good to ~10-15 m. VERTICAL: the LBNF
    # MCZero is at grade. The SBN model takes the Fermilab site grade (minus
    # berm/overburden) to be the MiniBooNE computer-room floor, already encoded
    # in sbn_loader: _FNAL_SITE_GRADE_Y (6.4733 m, tank-center frame) + the
    # tank-center->BNB shift _MB_Y_BNB (1.89614 m) = 8.369 m in the BNB frame
    # (so the BNB beam axis is ~8.4 m below grade). Hence y_up = 8.369 m.
    # (Uniform-grade approximation; the actual LBNF-site grade is ~757.5 ft
    # DUSAF = MI-10 + ~42 ft per CDR Vol3 p107, but we use the model grade.)
    # Hardcoded -- keep in sync with sbn_loader._FNAL_SITE_GRADE_Y.
    _lb_b, _lb_d = math.radians(287.79), 0.101
    _cb, _sb = math.cos(_lb_b), math.sin(_lb_b)
    _cd, _sd = math.cos(_lb_d), math.sin(_lb_d)
    _z_lbnf = np.array([-_cd * _sb, -_sd, _cd * _cb])   # LBNF +z (dipping beam) in BNB
    _x_lbnf = np.array([_cb, 0.0, _sb])                 # LBNF +x (horiz perp = up x beam)
    _y_lbnf = np.cross(_z_lbnf, _x_lbnf)                # LBNF +y (tilted up)
    g.add_transform(Transform(
        "LBNF", "BNB",
        np.column_stack([_x_lbnf, _y_lbnf, _z_lbnf]),
        np.array([261.5, 8.369, 35.8]),                 # (x_west, y_up=site grade, z_north) m
        "LBNF beam (az 287.79 deg, 101 mrad down) placed in BNB via the Fermilab "
        "MI-10/MI-60 FSCS survey centers bridged through NuMI; horizontal ~10-15 m, "
        "y_up=8.369 m: LBNF MCZero at site grade (MiniBooNE room floor; "
        "sbn_loader _FNAL_SITE_GRADE_Y + _MB_Y_BNB)"))

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

    # DUNE near-detector hall. Native frame == the dunendggd/edep-sim hall
    # "world" frame, which is exactly the DUNE_ND frame already in the graph, so
    # placement in BNB (DUNE_ND -> LBNF -> BNB) is exact. center_native is the
    # hall-world origin and the active box is a nominal ND-LAr extent, not the
    # exact per-detector active-LAr center/extent.
    detectors["DUNE_ND"] = Detector(
        "DUNE_ND", "DUNE_ND",
        np.array([0.0, 0.0, 0.0]),
        np.array([-3.5, -1.7, -2.5]),
        np.array([+3.5, +1.7, +2.5]))

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
