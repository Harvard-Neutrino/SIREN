"""siren.dk2nu coordinate frames: dk2nu files written in beam coordinate
systems other than the detector model's geometry frame.

dk2nu_to_primary_distribution takes the file's frame through the ``frame``
argument: None/"geometry" (the model's geometry frame, the default),
"detector" (detector-local, no model needed), or a rigid transform given
as a FrameTransform, a (rotation, translation) pair, or any object with
rotation/translation (or R/t) attributes. These tests exercise the
mechanism with synthetic single-row data, model-free through the
detector-frame target, and check that a real frame-graph transform (the
SBN NuMI-to-BNB survey) flows through the generic interface.
"""

import os
from types import SimpleNamespace

import numpy as np
import pytest

siren = pytest.importorskip("siren")

from siren import dk2nu
from siren import _util


M_PI = 0.13957039

# Rotation by +90 degrees about y: file z -> x, file x -> -z.
ROT_Y90 = np.array([[0.0, 0.0, 1.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, 0.0, 0.0]])


def _one_row_data():
    """A single pi+ row in file coordinates (positions in cm)."""
    return {
        "ptype": np.array([211]),
        "E": np.array([2.0]),
        "px": np.array([0.1]),
        "py": np.array([0.2]),
        "pz": np.array([1.5]),
        "vx": np.array([100.0]),
        "vy": np.array([-50.0]),
        "vz": np.array([3000.0]),
        "nimpwt": np.array([0.5]),
        "ntype": np.array([14]),
        "ndecay": np.array([13]),
        "t0": np.array([12.5]),
        "pol_x": np.array([0.6]),
        "pol_y": np.array([0.0]),
        "pol_z": np.array([0.8]),
        "pot": 2.0,
    }


def _sample_row(dist):
    """Sample the (single-row) distribution and finalize the record."""
    rand = siren.utilities.SIREN_random(1)
    record = siren.dataclasses.PrimaryDistributionRecord(
        siren.dataclasses.Particle.ParticleType.PiPlus)
    dist.Sample(rand, None, None, record)
    finalized = siren.dataclasses.InteractionRecord()
    record.finalize(finalized)
    return record, finalized


def test_detector_frame_is_a_passthrough():
    """frame="detector": cm -> m only, no model consulted (None is fine)."""
    dist = dk2nu.dk2nu_to_primary_distribution(
        _one_row_data(), detector_model=None, frame="detector")
    record, finalized = _sample_row(dist)

    assert record.interaction_vertex == pytest.approx([1.0, -0.5, 30.0])
    assert record.three_momentum == pytest.approx([0.1, 0.2, 1.5])
    assert record.initial_time == pytest.approx(12.5)
    p_sq = 0.1 ** 2 + 0.2 ** 2 + 1.5 ** 2
    assert record.energy == pytest.approx(np.sqrt(p_sq + M_PI ** 2))
    params = dict(finalized.interaction_parameters)
    assert params["weight"] == pytest.approx(0.25)  # nimpwt / pot
    assert (params["pol_x"], params["pol_y"], params["pol_z"]) \
        == pytest.approx((0.6, 0.0, 0.8))


def test_rigid_transform_rotates_everything_consistently():
    """A rotated, translated file frame: positions pick up rotation and
    translation; momenta and spin axes rotate; scalars are untouched."""
    translation = np.array([1.0, 2.0, 3.0])
    dist = dk2nu.dk2nu_to_primary_distribution(
        _one_row_data(), detector_model=None,
        frame=dk2nu.FrameTransform(ROT_Y90, translation, target="detector"))
    record, finalized = _sample_row(dist)

    assert record.interaction_vertex == pytest.approx(
        ROT_Y90 @ [1.0, -0.5, 30.0] + translation)
    assert record.three_momentum == pytest.approx(ROT_Y90 @ [0.1, 0.2, 1.5])
    assert record.initial_time == pytest.approx(12.5)
    p_sq = 0.1 ** 2 + 0.2 ** 2 + 1.5 ** 2
    assert record.energy == pytest.approx(np.sqrt(p_sq + M_PI ** 2))
    params = dict(finalized.interaction_parameters)
    assert (params["pol_x"], params["pol_y"], params["pol_z"]) \
        == pytest.approx(tuple(ROT_Y90 @ [0.6, 0.0, 0.8]))


def test_frame_argument_forms_normalize():
    """Tuple and duck-typed forms produce the same transform, targeting
    the geometry frame."""
    t = np.array([4.0, 5.0, 6.0])
    forms = (
        (ROT_Y90, t),
        SimpleNamespace(rotation=ROT_Y90, translation=t),
        SimpleNamespace(R=ROT_Y90, t=t),
    )
    for form in forms:
        xform = dk2nu._as_frame_transform(form)
        assert xform.target == "geometry"
        assert np.allclose(xform.rotation, ROT_Y90)
        assert np.allclose(xform.translation, t)

    identity = dk2nu._as_frame_transform(None)
    assert identity.target == "geometry"
    assert np.allclose(identity.rotation, np.eye(3))
    assert np.allclose(identity.translation, 0.0)
    assert dk2nu._as_frame_transform("detector").target == "detector"


def test_invalid_frames_fail_loud():
    with pytest.raises(siren.utilities.ConfigurationError):
        dk2nu.FrameTransform(np.diag([1.0, 1.0, -1.0]))  # reflection
    with pytest.raises(siren.utilities.ConfigurationError):
        dk2nu.FrameTransform(np.full((3, 3), 0.5))  # not orthonormal
    with pytest.raises(siren.utilities.ConfigurationError):
        dk2nu.FrameTransform(np.eye(4))  # wrong shape
    with pytest.raises(siren.utilities.ConfigurationError):
        dk2nu.FrameTransform(translation=[1.0, 2.0])  # wrong length
    with pytest.raises(siren.utilities.ConfigurationError):
        dk2nu.FrameTransform(target="lab")  # unknown target
    with pytest.raises(siren.utilities.ConfigurationError):
        dk2nu._as_frame_transform("bogus")
    with pytest.raises(siren.utilities.ConfigurationError):
        dk2nu._as_frame_transform(object())


def test_geometry_frame_requires_a_detector_model():
    """The default frame needs the model for the geometry-to-detector
    conversion; failing loud beats an AttributeError deep in the loop."""
    with pytest.raises(siren.utilities.ConfigurationError):
        dk2nu.dk2nu_to_primary_distribution(
            _one_row_data(), detector_model=None)


def test_frame_graph_transform_flows_through():
    """A detector resource's frame-graph transform (here the SBN NuMI to
    BNB survey) is an R/t-carrying object, so it feeds ``frame`` as-is.
    The FrameTransform validation doubles as a check that the survey
    rotation is a proper rotation."""
    try:
        base = _util.get_detector_model_path("SBN")
    except Exception:
        pytest.skip("SBN detector resource not installed")
    geo_py = os.path.join(base, "sbn_geometry.py")
    if not os.path.isfile(geo_py):
        pytest.skip("sbn_geometry.py not present")
    geo = _util.load_module("test_dk2nu_sbn_geometry", geo_py)

    survey = geo.transform("NuMI", "BNB")
    xform = dk2nu._as_frame_transform(survey)
    assert xform.target == "geometry"
    assert np.allclose(xform.rotation @ xform.rotation.T, np.eye(3))

    # The same transform retargeted to detector coordinates runs the full
    # builder without a detector model: a decay at the NuMI origin lands
    # at the frame's surveyed position, not at the file coordinates.
    data = _one_row_data()
    for k in ("vx", "vy", "vz"):
        data[k] = np.array([0.0])
    dist = dk2nu.dk2nu_to_primary_distribution(
        data, detector_model=None,
        frame=dk2nu.FrameTransform(survey.R, survey.t, target="detector"))
    record, _ = _sample_row(dist)
    assert record.interaction_vertex == pytest.approx(survey.t)
    assert record.three_momentum == pytest.approx(
        np.asarray(survey.R) @ [0.1, 0.2, 1.5])
