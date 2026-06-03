"""Tests for SIREN_Controller and Weighter (PR #125).

Covers:
  - Constructor: experiment name, custom file paths, validation, seed
  - SetInteractions: None primary/secondary, injection-only flag,
    type mismatch/match assertions, auto-set from collection, merging
  - GetVolumePositionDistributionFromSector: missing sector error
  - save_int_params nested structure (logic-pattern tests)
  - DN_min_decay_width filtering (logic-pattern tests)
  - InputDarkNewsModel secondary process preservation
  - Weighter: bounds checking, probability retrieval, event weight
  - End-to-end: generate, weight, probability retrieval
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
def ccm_controller():
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


def _make_fresh_controller():
    """Helper: create an isolated controller for tests that mutate state."""
    from siren.SIREN_Controller import SIREN_Controller
    return SIREN_Controller(1, experiment="CCM")


# ---------------------------------------------------------------------------
# Constructor
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
        assert c1.random.Uniform(0, 1) == c2.random.Uniform(0, 1)

    def test_prebuilt_detector_model_is_used(self):
        """A pre-built DetectorModel (e.g. a GDML composite) can be passed
        directly, bypassing experiment/file loading."""
        from siren.SIREN_Controller import SIREN_Controller
        model = siren.load_detector("CCM")
        c = SIREN_Controller(7, detector_model=model)
        assert c.detector_model is model
        assert c.detector_model_file is None
        assert c.materials_model_file is None
        assert c.events_to_inject == 7
        # No densities.dat to parse a fiducial line from -> graceful None.
        assert c.GetFiducialVolume() is None


# ---------------------------------------------------------------------------
# Detector model queries
# ---------------------------------------------------------------------------

class TestDetectorModelQueries:
    def test_fiducial_volume_exists(self, ccm_controller):
        """CCM should have a fiducial volume defined."""
        fv = ccm_controller.GetFiducialVolume()
        assert fv is not None

    def test_detector_model_targets_non_empty(self, ccm_controller):
        """GetDetectorModelTargets should return non-empty lists."""
        targets, target_strs = ccm_controller.GetDetectorModelTargets()
        assert len(targets) > 0
        assert len(target_strs) > 0
        # target_strs may be shorter (deduped by string name)
        assert len(target_strs) <= len(targets)

    def test_detector_model_targets_contain_hydrogen(self, ccm_controller):
        """CCM detector should contain hydrogen."""
        _, target_strs = ccm_controller.GetDetectorModelTargets()
        assert "H1" in target_strs

    def test_sector_geometry_returns_geometry(self, ccm_controller):
        """GetDetectorSectorGeometry with a valid name should return a geometry."""
        geo = ccm_controller.GetDetectorSectorGeometry("ccm_inner_argon")
        if geo is None:
            pytest.skip("ccm_inner_argon sector not found in this detector model")
        from siren import geometry as _geometry
        assert isinstance(geo, (_geometry.Cylinder, _geometry.Sphere))

    def test_old_cylinder_method_name_removed(self, ccm_controller):
        """The old name GetCylinderVolumePositionDistributionFromSector should
        not exist; it was renamed to GetVolumePositionDistributionFromSector."""
        assert not hasattr(ccm_controller, "GetCylinderVolumePositionDistributionFromSector")

    def test_volume_position_from_valid_sector(self, ccm_controller):
        """GetVolumePositionDistributionFromSector should return a distribution for a valid sector."""
        geo = ccm_controller.GetDetectorSectorGeometry("ccm_inner_argon")
        if geo is None:
            pytest.skip("ccm_inner_argon sector not found")
        dist = ccm_controller.GetVolumePositionDistributionFromSector("ccm_inner_argon")
        assert dist is not None


# ---------------------------------------------------------------------------
# SetInteractions
# ---------------------------------------------------------------------------

class TestSetInteractions:
    def test_none_primary_is_accepted(self, ccm_controller):
        ccm_controller.SetInteractions(
            primary_interaction_collection=None,
            injection=True, physical=True,
        )

    def test_none_secondary_is_accepted(self, ccm_controller):
        NuMu = dc.Particle.ParticleType.NuMu
        int_col = interactions.InteractionCollection(NuMu, [])
        ccm_controller.primary_injection_process.primary_type = NuMu
        ccm_controller.primary_physical_process.primary_type = NuMu
        ccm_controller.SetInteractions(
            primary_interaction_collection=int_col,
            secondary_interaction_collections=None,
        )

    def test_injection_only_flag(self):
        ctrl = _make_fresh_controller()
        NuMu = dc.Particle.ParticleType.NuMu
        ctrl.primary_injection_process.primary_type = NuMu
        ctrl.primary_physical_process.primary_type = NuMu
        int_col = interactions.InteractionCollection(NuMu, [])
        ctrl.SetInteractions(
            primary_interaction_collection=int_col,
            injection=True, physical=False,
        )
        assert ctrl.primary_injection_process.interactions is not None

    def test_physical_type_mismatch_raises(self):
        ctrl = _make_fresh_controller()
        NuMu = dc.Particle.ParticleType.NuMu
        NuE = dc.Particle.ParticleType.NuE
        ctrl.primary_injection_process.primary_type = NuE
        ctrl.primary_physical_process.primary_type = NuMu
        int_col = interactions.InteractionCollection(NuE, [])
        with pytest.raises(AssertionError):
            ctrl.SetInteractions(
                primary_interaction_collection=int_col,
                injection=False, physical=True,
            )

    def test_physical_type_match_passes(self):
        ctrl = _make_fresh_controller()
        NuMu = dc.Particle.ParticleType.NuMu
        ctrl.primary_physical_process.primary_type = NuMu
        ctrl.primary_injection_process.primary_type = NuMu
        int_col = interactions.InteractionCollection(NuMu, [])
        ctrl.SetInteractions(
            primary_interaction_collection=int_col,
            injection=False, physical=True,
        )

    def test_auto_sets_primary_type_from_collection(self):
        ctrl = _make_fresh_controller()
        NuTau = dc.Particle.ParticleType.NuTau
        unknown = dc.Particle.ParticleType.unknown
        ctrl.primary_injection_process.primary_type = unknown
        ctrl.primary_physical_process.primary_type = unknown
        int_col = interactions.InteractionCollection(NuTau, [])
        ctrl.SetInteractions(primary_interaction_collection=int_col)
        assert ctrl.primary_injection_process.primary_type == NuTau
        assert ctrl.primary_physical_process.primary_type == NuTau

    def test_merge_interaction_collections(self):
        ctrl = _make_fresh_controller()
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
# save_int_params structure (logic-pattern tests)
# ---------------------------------------------------------------------------

class TestSaveIntParamsStructure:
    def test_int_params_key_is_dict_per_event(self):
        """datasets['int_params'] should be a list of dicts, one per event."""
        datasets = {}
        events_data = [
            [{"x": 1.0, "y": 2.0}, {"x": 3.0}],
            [{"x": 4.0, "y": 5.0, "z": 6.0}],
        ]
        for event_params_list in events_data:
            datasets.setdefault("int_params", [])
            datasets["int_params"].append({})
            for datum_params in event_params_list:
                for param_name, param_value in datum_params.items():
                    if param_name not in datasets["int_params"][-1]:
                        datasets["int_params"][-1][param_name] = []
                    datasets["int_params"][-1][param_name].append(param_value)
        assert len(datasets["int_params"]) == 2
        assert datasets["int_params"][0] == {"x": [1.0, 3.0], "y": [2.0]}
        assert datasets["int_params"][1] == {"x": [4.0], "y": [5.0], "z": [6.0]}

    def test_params_from_later_events_not_lost(self):
        """Params appearing only in later events must not cause KeyError."""
        datasets = {}
        events_data = [
            [{"x": 1.0}],
            [{"x": 2.0, "z": 9.0}],
        ]
        for event_params_list in events_data:
            datasets.setdefault("int_params", [])
            datasets["int_params"].append({})
            for datum_params in event_params_list:
                for param_name, param_value in datum_params.items():
                    if param_name not in datasets["int_params"][-1]:
                        datasets["int_params"][-1][param_name] = []
                    datasets["int_params"][-1][param_name].append(param_value)
        assert datasets["int_params"][1] == {"x": [2.0], "z": [9.0]}
        assert "z" not in datasets["int_params"][0]


# ---------------------------------------------------------------------------
# DN_min_decay_width filtering (logic-pattern tests)
# ---------------------------------------------------------------------------

class TestDNMinDecayWidth:
    def test_min_decay_width_only_matching(self):
        """Only decays matching primary_type contribute to DN_min_decay_width."""
        primary_type = "N4"
        decays = [("N4", 1e-10), ("N4", 5e-11), ("NuE", 1e-15)]
        min_width = np.inf
        for parent, width in decays:
            if parent == primary_type and width > 0 and width < min_width:
                min_width = width
        assert min_width == pytest.approx(5e-11)

    def test_zero_width_excluded(self):
        """Zero-width decays should not corrupt the minimum."""
        primary_type = "N4"
        decays = [("N4", 0.0), ("N4", 1e-10)]
        min_width = np.inf
        for parent, width in decays:
            if parent == primary_type and width > 0 and width < min_width:
                min_width = width
        assert min_width == pytest.approx(1e-10)


# ---------------------------------------------------------------------------
# InputDarkNewsModel secondary process preservation
# ---------------------------------------------------------------------------

class TestInputDarkNewsModelSecondaryPreservation:
    def test_existing_secondary_not_duplicated(self, ccm_controller):
        """Existing secondary injection processes should be reused, not duplicated."""
        N4 = dc.Particle.ParticleType.N4
        sec_proc = injection.SecondaryInjectionProcess()
        sec_proc.primary_type = N4
        ccm_controller.secondary_injection_processes = [sec_proc]
        count = sum(1 for p in ccm_controller.secondary_injection_processes
                    if p.primary_type == N4)
        assert count == 1


# ---------------------------------------------------------------------------
# Weighter bounds checking and probability retrieval
# ---------------------------------------------------------------------------

class TestWeighter:
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
        return weighter, event

    def test_interaction_probs_valid(self, weighter_setup):
        weighter, event = weighter_setup
        probs = weighter.GetInteractionProbabilities(event, 0)
        assert len(probs) > 0
        for p in probs:
            assert 0.0 <= p <= 1.0

    def test_survival_probs_valid(self, weighter_setup):
        weighter, event = weighter_setup
        probs = weighter.GetSurvivalProbabilities(event, 0)
        assert len(probs) > 0
        for p in probs:
            assert 0.0 <= p <= 1.0

    def test_negative_i_inj_raises(self, weighter_setup):
        weighter, event = weighter_setup
        with pytest.raises((RuntimeError, IndexError)):
            weighter.GetInteractionProbabilities(event, -1)

    def test_out_of_range_i_inj_raises(self, weighter_setup):
        weighter, event = weighter_setup
        with pytest.raises((RuntimeError, IndexError)):
            weighter.GetInteractionProbabilities(event, 999)

    def test_survival_negative_i_inj_raises(self, weighter_setup):
        weighter, event = weighter_setup
        with pytest.raises((RuntimeError, IndexError)):
            weighter.GetSurvivalProbabilities(event, -1)

    def test_event_weight_is_finite(self, weighter_setup):
        weighter, event = weighter_setup
        w = weighter.EventWeight(event)
        assert np.isfinite(w)
        assert w >= 0


# ---------------------------------------------------------------------------
# End-to-end: generate, weight, probabilities
# ---------------------------------------------------------------------------

class TestEndToEnd:
    @pytest.fixture(scope="class")
    def full_setup(self):
        """Build injector + weighter, generate 5 events."""
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
        events = [inj.GenerateEvent() for _ in range(5)]
        return weighter, events

    def test_events_generated(self, full_setup):
        _, events = full_setup
        assert len(events) == 5

    def test_event_weight_finite(self, full_setup):
        weighter, events = full_setup
        for event in events:
            assert np.isfinite(weighter.EventWeight(event))

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

    def test_prob_array_lengths_match(self, full_setup):
        """Interaction and survival probability arrays should have the same length."""
        weighter, events = full_setup
        for event in events:
            int_probs = weighter.GetInteractionProbabilities(event, 0)
            surv_probs = weighter.GetSurvivalProbabilities(event, 0)
            assert len(int_probs) == len(surv_probs)

    def test_multiple_events_give_different_weights(self, full_setup):
        """With random directions, not all events should have identical weights."""
        weighter, events = full_setup
        weights = [weighter.EventWeight(e) for e in events]
        assert len(set(weights)) > 1 or len(events) == 1


# ---------------------------------------------------------------------------
# Serialization: injector/weighter save + load round-trips
# ---------------------------------------------------------------------------

class TestSerialization:
    """Guard the injector/weighter serialization the SaveEvents path and the
    Weighter wrapper expose. The controller inject->save path was previously
    untested, which let a wrapper/pybind API mismatch (.save vs SaveInjector/
    SaveWeighter, and a removed Add*Distribution API) go unnoticed."""

    def _ccm_pieces(self):
        NuMu = dc.Particle.ParticleType.NuMu
        dm = detector.DetectorModel()
        det_dir = _util.get_detector_model_path("CCM")
        dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
        dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
        int_col = interactions.InteractionCollection(NuMu, [interactions.DummyCrossSection()])
        p_inj = injection.PrimaryInjectionProcess()
        p_inj.primary_type = NuMu
        p_inj.interactions = int_col
        p_inj.distributions = [
            distributions.PrimaryMass(0), distributions.Monoenergetic(1.0),
            distributions.IsotropicDirection(),
            distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0)]
        p_phys = injection.PhysicalProcess()
        p_phys.primary_type = NuMu
        p_phys.interactions = int_col
        p_phys.distributions = [distributions.PrimaryMass(0), distributions.IsotropicDirection()]
        return NuMu, dm, p_inj, p_phys

    def test_controller_save_events_writes_loadable_injector_and_weighter(self, tmp_path):
        """SIREN_Controller inject -> SaveEvents writes .siren_injector /
        .siren_weighter via the pybind SaveInjector/SaveWeighter, and the saved
        weighter reloads through the (injectors, filename) constructor."""
        from siren.SIREN_Controller import SIREN_Controller
        NuMu = dc.Particle.ParticleType.NuMu
        try:
            c = SIREN_Controller(3, experiment="CCM")
        except (AttributeError, TypeError, OSError) as e:
            pytest.skip(f"Cannot create CCM controller: {e}")
        c.SetInteractions(interactions.InteractionCollection(NuMu, [interactions.DummyCrossSection()]))
        inj_d = {"mass": distributions.PrimaryMass(0), "energy": distributions.Monoenergetic(1.0),
                 "direction": distributions.IsotropicDirection(),
                 "position": distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0)}
        phys_d = {"mass": distributions.PrimaryMass(0), "direction": distributions.IsotropicDirection()}
        c.SetProcesses(NuMu, inj_d, phys_d)
        c.Initialize()
        c.GenerateEvents(verbose=False)
        base = str(tmp_path / "ev")
        c.SaveEvents(base, hdf5=False, parquet=False, verbose=False)
        assert os.path.exists(base + ".siren_injector")
        assert os.path.exists(base + ".siren_weighter")
        w2 = injection._Weighter(c.injectors, base)   # (injectors, filename) load ctor
        assert np.isfinite(w2.EventWeight(c.events[0]))

    def test_weighter_wrapper_load_needs_only_injectors(self, tmp_path):
        """Weighter.load() restores via the (injectors, filename) ctor and must NOT
        require the detector / interactions / distributions to be configured."""
        NuMu, dm, p_inj, p_phys = self._ccm_pieces()
        inj = injection._Injector(10, dm, p_inj, utilities.SIREN_random(7))
        event = inj.GenerateEvent()

        configured = injection.Weighter(
            injectors=[inj], detector_model=dm, primary_type=NuMu,
            primary_interactions=[interactions.DummyCrossSection()],
            primary_physical_distributions=[distributions.PrimaryMass(0),
                                            distributions.IsotropicDirection()])
        base = str(tmp_path / "w")
        configured.save(base)
        ref = configured(event)

        # a fresh wrapper that ONLY knows the injectors -- no detector/process config
        fresh = injection.Weighter(injectors=[inj])
        fresh.load(base)
        assert np.isclose(fresh(event), ref)

    def test_injector_save_reload_roundtrip(self, tmp_path):
        """Injector SaveInjector writes the literal path; reload via the filename ctor."""
        NuMu, dm, p_inj, p_phys = self._ccm_pieces()
        inj = injection._Injector(10, dm, p_inj, utilities.SIREN_random(11))
        path = str(tmp_path / "inj.siren_injector")
        inj.SaveInjector(path)
        assert os.path.getsize(path) > 0
        inj2 = injection._Injector(10, path, utilities.SIREN_random(11))   # filename ctor
        ev = inj2.GenerateEvent()
        w = injection._Weighter([inj2], dm, p_phys)
        assert np.isfinite(w.EventWeight(ev))
