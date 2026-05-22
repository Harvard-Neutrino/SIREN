"""Tests for SBN geometry coordinate transforms and rotation utilities."""
from __future__ import annotations

import importlib.util
import math
import os
import sys

import numpy as np
import pytest

_SBN_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "resources", "detectors",
    "SBN", "SBN-v1")


def test_detector_loader_imports_without_preloaded_siblings():
    """Importing detector.py must load dataclass sibling modules itself."""
    module_names = ("sbn_detector_import_smoke", "sbn_geometry", "sbn_loader")
    previous_modules = {name: sys.modules.pop(name, None) for name in module_names}
    try:
        spec = importlib.util.spec_from_file_location(
            "sbn_detector_import_smoke",
            os.path.join(_SBN_DIR, "detector.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        assert callable(mod.load_detector)
        assert "ICARUS" in mod._DETECTOR_SPECS
        assert mod.geo.DETECTORS["ICARUS"].name == "ICARUS"
    finally:
        for name in module_names:
            sys.modules.pop(name, None)
        for name, module in previous_modules.items():
            if module is not None:
                sys.modules[name] = module


def test_sbn_resource_loader_imports_from_source_tree(monkeypatch):
    """The public resource loader path must import the SBN detector loader."""
    from siren import _util

    resources_root = os.path.abspath(os.path.join(_SBN_DIR, "..", "..", ".."))
    module_names = ("siren-detector-SBN", "sbn_geometry", "sbn_loader")
    previous_modules = {name: sys.modules.pop(name, None) for name in module_names}
    monkeypatch.setattr(_util, "resource_package_dir", lambda: resources_root)

    try:
        loader = _util.get_resource_loader("detector", "SBN")
        assert callable(loader)
        assert "ICARUS" in loader._DETECTOR_SPECS
        assert "SBND" in loader._DETECTOR_SPECS
    finally:
        for name in module_names:
            sys.modules.pop(name, None)
        for name, module in previous_modules.items():
            if module is not None:
                sys.modules[name] = module


@pytest.fixture(scope="module")
def geo():
    spec = importlib.util.spec_from_file_location(
        "sbn_geometry", os.path.join(_SBN_DIR, "sbn_geometry.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules["sbn_geometry"] = mod
    spec.loader.exec_module(mod)
    return mod


# ======================================================================
# gdml_rotation_angles: round-trip via matrix reconstruction
# ======================================================================

def _reconstruct_matrix(rx, ry, rz):
    """Build R = Rz(rz) @ Ry(ry) @ Rx(rx) (extrinsic XYZ / GDML convention)."""
    cx, sx = math.cos(rx), math.sin(rx)
    cy, sy = math.cos(ry), math.sin(ry)
    cz, sz = math.cos(rz), math.sin(rz)
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


class TestGDMLRotationAngles:

    def test_identity(self, geo):
        rx, ry, rz = geo.gdml_rotation_angles(np.eye(3))
        assert abs(rx) < 1e-15
        assert abs(ry) < 1e-15
        assert abs(rz) < 1e-15

    def test_single_axis_x(self, geo):
        angle = 0.7
        R = _reconstruct_matrix(angle, 0, 0)
        rx, ry, rz = geo.gdml_rotation_angles(R)
        assert abs(rx - angle) < 1e-14
        assert abs(ry) < 1e-14
        assert abs(rz) < 1e-14

    def test_single_axis_y(self, geo):
        angle = -0.4
        R = _reconstruct_matrix(0, angle, 0)
        rx, ry, rz = geo.gdml_rotation_angles(R)
        assert abs(rx) < 1e-14
        assert abs(ry - angle) < 1e-14
        assert abs(rz) < 1e-14

    def test_single_axis_z(self, geo):
        angle = 1.2
        R = _reconstruct_matrix(0, 0, angle)
        rx, ry, rz = geo.gdml_rotation_angles(R)
        assert abs(rx) < 1e-14
        assert abs(ry) < 1e-14
        assert abs(rz - angle) < 1e-14

    def test_combined_rotation(self, geo):
        R = _reconstruct_matrix(0.3, -0.5, 0.8)
        rx, ry, rz = geo.gdml_rotation_angles(R)
        R_recon = _reconstruct_matrix(rx, ry, rz)
        np.testing.assert_allclose(R, R_recon, atol=1e-14)

    def test_90_degree_z(self, geo):
        R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
        rx, ry, rz = geo.gdml_rotation_angles(R)
        R_recon = _reconstruct_matrix(rx, ry, rz)
        np.testing.assert_allclose(R, R_recon, atol=1e-14)

    def test_icarus_rotation(self, geo):
        """The real ICARUS->NuMI rotation matrix must round-trip."""
        R = geo.R_ICARUS_TO_NUMI
        rx, ry, rz = geo.gdml_rotation_angles(R)
        R_recon = _reconstruct_matrix(rx, ry, rz)
        np.testing.assert_allclose(R, R_recon, atol=1e-9)

    def test_numi_to_bnb_transform(self, geo):
        T = geo.transform("NuMI", "BNB")
        rx, ry, rz = geo.gdml_rotation_angles(T.R)
        R_recon = _reconstruct_matrix(rx, ry, rz)
        np.testing.assert_allclose(T.R, R_recon, atol=1e-9)

    @pytest.mark.parametrize("seed", range(10))
    def test_random_rotations(self, geo, seed):
        rng = np.random.RandomState(seed)
        M = rng.randn(3, 3)
        Q, _ = np.linalg.qr(M)
        if np.linalg.det(Q) < 0:
            Q[:, 0] *= -1
        rx, ry, rz = geo.gdml_rotation_angles(Q)
        R_recon = _reconstruct_matrix(rx, ry, rz)
        np.testing.assert_allclose(Q, R_recon, atol=1e-14)

    def test_near_gimbal_lock(self, geo):
        R = _reconstruct_matrix(0.3, math.pi / 2 - 1e-8, 0.5)
        rx, ry, rz = geo.gdml_rotation_angles(R)
        R_recon = _reconstruct_matrix(rx, ry, rz)
        np.testing.assert_allclose(R, R_recon, atol=1e-7)


# ======================================================================
# gdml_rotation_angles + quaternion: match C++ GDML parser path
# ======================================================================

class TestGDMLRotationMatchesCpp:
    """Verify that angles from gdml_rotation_angles, when fed to the C++
    QFromXYZs (via Quaternion.SetEulerAnglesXYZs), produce a quaternion
    matching quaternion_from_matrix applied to the same rotation matrix."""

    @staticmethod
    def _quaternion_dot(q_py, q_cpp):
        """Absolute dot product of two quaternions (1.0 = same rotation)."""
        from siren.math import Quaternion as Q
        v = np.array([q_py[0], q_py[1], q_py[2], q_py[3]])
        w = np.array([q_cpp.X, q_cpp.Y, q_cpp.Z, q_cpp.W])
        return abs(np.dot(v, w))

    def _check(self, geo, R, atol=1e-10):
        from siren.math import Quaternion
        rx, ry, rz = geo.gdml_rotation_angles(R)
        q_py = geo.quaternion_from_matrix(R)
        q_cpp = Quaternion()
        q_cpp.SetEulerAnglesXYZs(rx, ry, rz)
        dot = self._quaternion_dot(q_py, q_cpp)
        assert abs(dot - 1.0) < atol, (
            f"Quaternion mismatch: dot={dot}, "
            f"py=({q_py}), cpp=({q_cpp.X},{q_cpp.Y},{q_cpp.Z},{q_cpp.W})")

    def test_identity(self, geo):
        self._check(geo, np.eye(3))

    def test_icarus_rotation(self, geo):
        self._check(geo, geo.R_ICARUS_TO_NUMI)

    def test_numi_to_bnb(self, geo):
        T = geo.transform("NuMI", "BNB")
        self._check(geo, T.R)

    def test_combined_45deg(self, geo):
        self._check(geo, _reconstruct_matrix(
            math.pi / 4, math.pi / 4, math.pi / 4))

    @pytest.mark.parametrize("seed", range(10))
    def test_random_rotations(self, geo, seed):
        rng = np.random.RandomState(seed)
        M = rng.randn(3, 3)
        Q, _ = np.linalg.qr(M)
        if np.linalg.det(Q) < 0:
            Q[:, 0] *= -1
        self._check(geo, Q)


# ======================================================================
# quaternion_from_matrix
# ======================================================================

class TestQuaternionFromMatrix:

    def test_identity(self, geo):
        qx, qy, qz, qw = geo.quaternion_from_matrix(np.eye(3))
        assert abs(qx) < 1e-15
        assert abs(qy) < 1e-15
        assert abs(qz) < 1e-15
        assert abs(qw - 1.0) < 1e-15

    def test_180_about_z(self, geo):
        R = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]], dtype=float)
        qx, qy, qz, qw = geo.quaternion_from_matrix(R)
        assert abs(qz) - 1.0 < 1e-14 or abs(qz) + 1.0 < 1e-14
        assert abs(qw) < 1e-14
        assert abs(qx) < 1e-14
        assert abs(qy) < 1e-14

    def test_unit_quaternion(self, geo):
        """Output must always be a unit quaternion."""
        rng = np.random.RandomState(99)
        for _ in range(20):
            M = rng.randn(3, 3)
            Q, _ = np.linalg.qr(M)
            if np.linalg.det(Q) < 0:
                Q[:, 0] *= -1
            qx, qy, qz, qw = geo.quaternion_from_matrix(Q)
            norm = math.sqrt(qx**2 + qy**2 + qz**2 + qw**2)
            assert abs(norm - 1.0) < 1e-14

    def test_roundtrip_with_rotation_vector(self, geo):
        """Quaternion -> matrix -> quaternion should be consistent."""
        from siren.math import Quaternion
        q_orig = Quaternion()
        q_orig.SetEulerAnglesXYZs(0.3, -0.7, 1.1)
        # Get the rotation matrix from the C++ quaternion
        # by reconstructing from the known angles
        R = _reconstruct_matrix(0.3, -0.7, 1.1)
        qx, qy, qz, qw = geo.quaternion_from_matrix(R)
        q_py = np.array([qx, qy, qz, qw])
        q_cx = np.array([q_orig.X, q_orig.Y, q_orig.Z, q_orig.W])
        dot = abs(np.dot(q_py, q_cx))
        assert abs(dot - 1.0) < 1e-12


# ======================================================================
# Frame graph transforms
# ======================================================================

class TestFrameGraph:

    def test_identity_transform(self, geo):
        T = geo.transform("BNB", "BNB")
        np.testing.assert_allclose(T.R, np.eye(3), atol=1e-15)
        np.testing.assert_allclose(T.t, np.zeros(3), atol=1e-15)

    def test_inverse_roundtrip(self, geo):
        T = geo.transform("ICARUS_LArSoft", "BNB")
        T_inv = T.inverse()
        T_rt = T.compose(T_inv)
        np.testing.assert_allclose(T_rt.R, np.eye(3), atol=1e-12)
        np.testing.assert_allclose(T_rt.t, np.zeros(3), atol=1e-12)

    def test_icarus_at_600m(self, geo):
        """ICARUS is at z=600m in BNB frame (pure translation)."""
        T = geo.transform("ICARUS_LArSoft", "BNB")
        np.testing.assert_allclose(T.t, [0.0, 0.0, 600.0], atol=1e-10)
        np.testing.assert_allclose(T.R, np.eye(3), atol=1e-15)

    def test_sbnd_position(self, geo):
        T = geo.transform("SBND_LArSoft", "BNB")
        np.testing.assert_allclose(T.t, [0.7378, 0.0, 110.0], atol=1e-10)

    def test_convert_origin(self, geo):
        origin_bnb = geo.convert([0, 0, 0], "ICARUS_LArSoft", "BNB")
        np.testing.assert_allclose(origin_bnb, [0.0, 0.0, 600.0], atol=1e-10)

    def test_detector_center(self, geo):
        c = geo.detector_center("ICARUS", "BNB")
        assert abs(c[2] - 600.0) < 1.0

    def test_numi_transform_has_rotation(self, geo):
        T = geo.transform("NuMI", "BNB")
        assert not np.allclose(T.R, np.eye(3), atol=1e-3)


# ======================================================================
# GDML placement convention: passive rotation (M^T applied)
# ======================================================================

class TestGDMLPlacementConvention:
    """Verify that gdml_rotation_angles produces angles for the GDML passive
    convention: SIREN applies r_parent = M^T @ r_child + position.

    For a child volume in frame S placed in parent frame P, the GDML rotation
    matrix M must satisfy M^T = R_child_to_parent, i.e. M = R_parent_to_child.
    """

    def test_numi_placement_roundtrip(self, geo):
        """NuMI GDML placed in BNB frame: M^T should equal R(NuMI->BNB)."""
        from siren.math import Quaternion
        T = geo.transform("NuMI", "BNB")
        # GDML angles should encode M = R(NuMI->BNB)^T = R(BNB->NuMI)
        rx, ry, rz = geo.gdml_rotation_angles(T.R.T)
        M = _reconstruct_matrix(rx, ry, rz)
        # SIREN applies M^T, which should equal R(NuMI->BNB)
        np.testing.assert_allclose(M.T, T.R, atol=1e-9)

    def test_identity_placement(self, geo):
        """Pure translation (ICARUS/SBND) should give zero angles."""
        R_identity = np.eye(3)
        rx, ry, rz = geo.gdml_rotation_angles(R_identity.T)
        assert abs(rx) < 1e-15
        assert abs(ry) < 1e-15
        assert abs(rz) < 1e-15

    @pytest.mark.parametrize("seed", range(5))
    def test_random_placement(self, geo, seed):
        """For any rotation R(child->parent), GDML angles of R^T should
        reconstruct M such that M^T = R."""
        rng = np.random.RandomState(seed + 100)
        M = rng.randn(3, 3)
        R_child_to_parent, _ = np.linalg.qr(M)
        if np.linalg.det(R_child_to_parent) < 0:
            R_child_to_parent[:, 0] *= -1
        rx, ry, rz = geo.gdml_rotation_angles(R_child_to_parent.T)
        M_gdml = _reconstruct_matrix(rx, ry, rz)
        np.testing.assert_allclose(M_gdml.T, R_child_to_parent, atol=1e-14)
