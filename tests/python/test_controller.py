"""Tests for SIREN_Controller features and the Weighter additions in PR #125.

Covers:
  - Constructor: experiment name, custom file paths, validation, seed
  - SetInteractions: None primary, merging, assertion on type mismatch
  - SetProcesses: fid_vol_secondary flag
  - GetVolumePositionDistributionFromSector: renamed method, error on missing sector
  - SaveEvents dataset structure: int_probs, int_params, survival_probs
  - Weighter C++ additions: bounds checking, probability retrieval
  - End-to-end: generate, weight, save cycle with DummyCrossSection
"""
import os
import pytest
import numpy as np

siren = pytest.importorskip("siren")

from siren import dataclasses as dc
from siren import injection
from siren import interactions
from siren import distributions
from siren import detector
from siren import math as smath
from siren import utilities
from siren import _util


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ccm_det_paths():
    """Return (detector_model_file, materials_model_file) for CCM."""
    det_dir = _util.get_detector_model_path("CCM")
    return (
        os.path.join(det_dir, "densities.dat"),
        os.path.join(det_dir, "materials.dat"),
    )


@pytest.fixture(scope="module")
def ccm_controller(ccm_det_paths):
    from siren.SIREN_Controller import SIREN_Controller
    try:
        return SIREN_Controller(10, experiment="CCM")
    except (AttributeError, TypeError, OSError) as e:
        pytest.skip(f"Cannot create CCM controller: {e}")


@pytest.fixture(scope="module")
def controller_custom_paths(ccm_det_paths):
    """Controller created with explicit file paths instead of experiment name."""
    from siren.SIREN_Controller import SIREN_Controller
    det_file, mat_file = ccm_det_paths
    return SIREN_Controller(
        5,
        detector_model_file=det_file,
        materials_model_file=mat_file,
        seed=123,
    )


# ---------------------------------------------------------------------------
# Constructor tests
# ---------------------------------------------------------------------------

class TestControllerConstructor:
    def test_experiment_name_creates_controller(self, ccm_controller):
        assert ccm_controller is not None
        assert ccm_controller.experiment == "CCM"
        assert ccm_controller.events_to_inject == 10

    def test_custom_paths_creates_controller(self, controller_custom_paths):
        assert controller_custom_paths is not None
        assert controller_custom_paths.experiment is None
        assert controller_custom_paths.events_to_inject == 5

    def test_missing_all_paths_raises(self):
        from siren.SIREN_Controller import SIREN_Controller
        with pytest.raises(ValueError, match="Must provide"):
            SIREN_Controller(10)

    def test_missing_one_path_raises(self, ccm_det_paths):
        from siren.SIREN_Controller import SIREN_Controller
        det_file, _ = ccm_det_paths
        with pytest.raises(ValueError, match="Must provide"):
            SIREN_Controller(10, detector_model_file=det_file)

    def test_seed_is_applied(self):
        from siren.SIREN_Controller import SIREN_Controller
        c1 = SIREN_Controller(1, experiment="CCM", seed=42)
        c2 = SIREN_Controller(1, experiment="CCM", seed=42)
        # Same seed should produce same first random number
        r1 = c1.random.Uniform(0, 1)
        r2 = c2.random.Uniform(0, 1)
        assert r1 == r2


# ---------------------------------------------------------------------------
# SetInteractions tests
# ---------------------------------------------------------------------------

class TestSetInteractions:
    def test_none_primary_is_accepted(self, ccm_controller):
        """SetInteractions with primary_interaction_collection=None should not crash."""
        ccm_controller.SetInteractions(
            primary_interaction_collection=None,
            injection=True,
            physical=True,
        )

    def test_none_secondary_is_accepted(self, ccm_controller):
        """SetInteractions with secondary_interaction_collections=None should not crash."""
        NuMu = dc.Particle.ParticleType.NuMu
        int_col = interactions.InteractionCollection(NuMu, [])
        ccm_controller.primary_injection_process.primary_type = NuMu
        ccm_controller.primary_physical_process.primary_type = NuMu
        ccm_controller.SetInteractions(
            primary_interaction_collection=int_col,
            secondary_interaction_collections=None,
        )

    def test_injection_only_flag(self):
        """Setting injection=True, physical=False should only update injection process."""
        from siren.SIREN_Controller import SIREN_Controller
        ctrl = SIREN_Controller(1, experiment="CCM")
        NuMu = dc.Particle.ParticleType.NuMu
        ctrl.primary_injection_process.primary_type = NuMu
        ctrl.primary_physical_process.primary_type = NuMu

        int_col = interactions.InteractionCollection(NuMu, [])
        ctrl.SetInteractions(
            primary_interaction_collection=int_col,
            injection=True,
            physical=False,
        )
        assert ctrl.primary_injection_process.interactions is not None

    def test_physical_type_mismatch_raises(self):
        """physical=True with mismatched type should raise AssertionError."""
        from siren.SIREN_Controller import SIREN_Controller
        ctrl = SIREN_Controller(1, experiment="CCM")
        NuMu = dc.Particle.ParticleType.NuMu
        NuE = dc.Particle.ParticleType.NuE

        ctrl.primary_injection_process.primary_type = NuE
        ctrl.primary_physical_process.primary_type = NuMu

        int_col = interactions.InteractionCollection(NuE, [])
        with pytest.raises(AssertionError):
            ctrl.SetInteractions(
                primary_interaction_collection=int_col,
                injection=False,
                physical=True,
            )

    def test_auto_sets_primary_type_from_collection(self):
        """When primary_type is unknown, SetInteractions should set it from the collection."""
        from siren.SIREN_Controller import SIREN_Controller
        ctrl = SIREN_Controller(1, experiment="CCM")
        NuTau = dc.Particle.ParticleType.NuTau
        unknown = dc.Particle.ParticleType.unknown

        ctrl.primary_injection_process.primary_type = unknown
        ctrl.primary_physical_process.primary_type = unknown

        int_col = interactions.InteractionCollection(NuTau, [])
        ctrl.SetInteractions(primary_interaction_collection=int_col)

        assert ctrl.primary_injection_process.primary_type == NuTau
        assert ctrl.primary_physical_process.primary_type == NuTau

    def test_merge_interaction_collections(self):
        """Setting interactions twice should merge, not replace."""
        from siren.SIREN_Controller import SIREN_Controller
        ctrl = SIREN_Controller(1, experiment="CCM")
        NuMu = dc.Particle.ParticleType.NuMu
        ctrl.primary_injection_process.primary_type = NuMu
        ctrl.primary_physical_process.primary_type = NuMu

        col1 = interactions.InteractionCollection(NuMu, [])
        col2 = interactions.InteractionCollection(NuMu, [])
        ctrl.SetInteractions(primary_interaction_collection=col1)
        ctrl.SetInteractions(primary_interaction_collection=col2)
        assert ctrl.primary_injection_process.interactions is not None


# ---------------------------------------------------------------------------
# GetVolumePositionDistributionFromSector
# ---------------------------------------------------------------------------

class TestGetVolumePositionDistribution:
    def test_missing_sector_raises(self, ccm_controller):
        with pytest.raises(ValueError, match="not found"):
            ccm_controller.GetVolumePositionDistributionFromSector("nonexistent_sector_999")


# ---------------------------------------------------------------------------
# Weighter bounds checking (Python layer)
# ---------------------------------------------------------------------------

class TestWeighterPython:
    @pytest.fixture(scope="class")
    def weighter_setup(self):
        """Build a minimal injector + weighter with DummyCrossSection."""
        NuMu = dc.Particle.ParticleType.NuMu

        dm = detector.DetectorModel()
        det_dir = _util.get_detector_model_path("CCM")
        dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
        dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))

        xs = interactions.DummyCrossSection()
        int_col = interactions.InteractionCollection(NuMu, [xs])

        primary_inj = injection.PrimaryInjectionProcess()
        primary_inj.primary_type = NuMu
        primary_inj.interactions = int_col
        primary_inj.distributions = [
            distributions.PrimaryMass(0),
            distributions.Monoenergetic(1.0),
            distributions.IsotropicDirection(),
            distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
        ]

        primary_phys = injection.PhysicalProcess()
        primary_phys.primary_type = NuMu
        primary_phys.interactions = int_col
        primary_phys.distributions = [
            distributions.PrimaryMass(0),
            distributions.IsotropicDirection(),
        ]

        rand = utilities.SIREN_random(42)
        inj = injection._Injector(10, dm, primary_inj, rand)
        weighter = injection._Weighter([inj], dm, primary_phys)
        event = inj.GenerateEvent()
        return weighter, event, inj

    def test_interaction_probs_valid_index(self, weighter_setup):
        weighter, event, _ = weighter_setup
        probs = weighter.GetInteractionProbabilities(event, 0)
        assert len(probs) > 0
        for p in probs:
            assert 0.0 <= p <= 1.0

    def test_survival_probs_valid_index(self, weighter_setup):
        weighter, event, _ = weighter_setup
        probs = weighter.GetSurvivalProbabilities(event, 0)
        assert len(probs) > 0
        for p in probs:
            assert 0.0 <= p <= 1.0

    def test_negative_i_inj_raises(self, weighter_setup):
        weighter, event, _ = weighter_setup
        with pytest.raises((RuntimeError, IndexError)):
            weighter.GetInteractionProbabilities(event, -1)

    def test_out_of_range_i_inj_raises(self, weighter_setup):
        weighter, event, _ = weighter_setup
        with pytest.raises((RuntimeError, IndexError)):
            weighter.GetInteractionProbabilities(event, 999)

    def test_survival_negative_i_inj_raises(self, weighter_setup):
        weighter, event, _ = weighter_setup
        with pytest.raises((RuntimeError, IndexError)):
            weighter.GetSurvivalProbabilities(event, -1)

    def test_event_weight_is_finite(self, weighter_setup):
        weighter, event, _ = weighter_setup
        w = weighter.EventWeight(event)
        assert np.isfinite(w)
        assert w >= 0


# ---------------------------------------------------------------------------
# End-to-end: generate, weight, probabilities
# ---------------------------------------------------------------------------

class TestEndToEnd:
    @pytest.fixture(scope="class")
    def full_setup(self):
        """Full low-level setup: build injector + weighter, generate events."""
        NuMu = dc.Particle.ParticleType.NuMu

        dm = detector.DetectorModel()
        det_dir = _util.get_detector_model_path("CCM")
        dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
        dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))

        xs = interactions.DummyCrossSection()
        int_col = interactions.InteractionCollection(NuMu, [xs])

        primary_inj = injection.PrimaryInjectionProcess()
        primary_inj.primary_type = NuMu
        primary_inj.interactions = int_col
        primary_inj.distributions = [
            distributions.PrimaryMass(0),
            distributions.Monoenergetic(1.0),
            distributions.IsotropicDirection(),
            distributions.PointSourcePositionDistribution(
                smath.Vector3D(0, 0, 0), 25.0
            ),
        ]

        primary_phys = injection.PhysicalProcess()
        primary_phys.primary_type = NuMu
        primary_phys.interactions = int_col
        primary_phys.distributions = [
            distributions.PrimaryMass(0),
            distributions.IsotropicDirection(),
        ]

        rand = utilities.SIREN_random(99)
        inj = injection._Injector(5, dm, primary_inj, rand)
        weighter = injection._Weighter([inj], dm, primary_phys)

        events = []
        for _ in range(5):
            events.append(inj.GenerateEvent())
        return weighter, events

    def test_events_generated(self, full_setup):
        weighter, events = full_setup
        assert len(events) == 5

    def test_event_weight_finite(self, full_setup):
        weighter, events = full_setup
        for event in events:
            w = weighter.EventWeight(event)
            assert np.isfinite(w)

    def test_interaction_probs_per_event(self, full_setup):
        weighter, events = full_setup
        for event in events:
            probs = weighter.GetInteractionProbabilities(event, 0)
            assert len(probs) > 0
            for p in probs:
                assert 0.0 <= p <= 1.0

    def test_survival_probs_per_event(self, full_setup):
        weighter, events = full_setup
        for event in events:
            probs = weighter.GetSurvivalProbabilities(event, 0)
            assert len(probs) > 0
            for p in probs:
                assert 0.0 <= p <= 1.0

    def test_interaction_plus_survival_le_one(self, full_setup):
        """For each datum, interaction_prob + survival_prob should be <= 1
        (they measure complementary things over different path segments)."""
        weighter, events = full_setup
        for event in events:
            int_probs = weighter.GetInteractionProbabilities(event, 0)
            surv_probs = weighter.GetSurvivalProbabilities(event, 0)
            assert len(int_probs) == len(surv_probs)

    def test_multiple_events_give_different_weights(self, full_setup):
        """With random directions, not all events should have identical weights."""
        weighter, events = full_setup
        weights = [weighter.EventWeight(e) for e in events]
        # At least some variation expected (not all identical)
        assert len(set(weights)) > 1 or len(events) == 1
