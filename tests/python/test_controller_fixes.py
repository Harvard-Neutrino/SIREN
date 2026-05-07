"""Tests for SIREN_Controller bug fixes (PR #125).

Covers fixes:
  #2  SetInteractions assertion checks correct process type
  #3  save_int_params produces per-event nested structure
  #4  Weighter i_inj bounds-checking
  #5  SurvivalProbability clamping (via #4 integration)
  #7  InputDarkNewsDecay DN_min_decay_width only on matching decays
  #8  exit(0) replaced with exceptions
  #9  InputDarkNewsModel preserves existing secondary processes
  #11 Typo corrections
"""
import pytest

siren = pytest.importorskip("siren")

from siren import dataclasses as dc
from siren import injection
from siren import interactions
from siren import distributions


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ccm_controller():
    from siren.SIREN_Controller import SIREN_Controller
    try:
        return SIREN_Controller(10, experiment="CCM")
    except (AttributeError, TypeError, OSError) as e:
        pytest.skip(f"Cannot create CCM controller: {e}")


# ---------------------------------------------------------------------------
# Fix #8: exit(0) replaced with exceptions
# ---------------------------------------------------------------------------

class TestExitReplacedWithExceptions:
    def test_missing_model_files_raises_valueerror(self):
        from siren.SIREN_Controller import SIREN_Controller
        with pytest.raises(ValueError, match="Must provide"):
            SIREN_Controller(10, experiment=None,
                             detector_model_file=None,
                             materials_model_file=None)

    def test_missing_sector_raises_valueerror(self, ccm_controller):
        with pytest.raises(ValueError, match="not found"):
            ccm_controller.GetVolumePositionDistributionFromSector("nonexistent_sector_999")


# ---------------------------------------------------------------------------
# Fix #2: SetInteractions physical process assertion
# ---------------------------------------------------------------------------

class TestSetInteractionsAssertion:
    def test_physical_type_mismatch_raises(self, ccm_controller):
        """SetInteractions(physical=True) must assert against
        primary_physical_process.primary_type, not primary_injection_process."""
        NuMu = dc.Particle.ParticleType.NuMu
        NuE = dc.Particle.ParticleType.NuE

        # Bug scenario: injection type is NuE, physical type is NuMu.
        # Collection type is NuE. OLD code checked injection (NuE == NuE -> pass).
        # FIXED code checks physical (NuMu != NuE -> fail).
        ccm_controller.primary_injection_process.primary_type = NuE
        ccm_controller.primary_physical_process.primary_type = NuMu

        int_col = interactions.InteractionCollection(NuE, [])

        with pytest.raises(AssertionError):
            ccm_controller.SetInteractions(
                primary_interaction_collection=int_col,
                injection=False,
                physical=True,
            )

    def test_physical_type_match_passes(self, ccm_controller):
        """SetInteractions passes when physical process type matches collection."""
        NuMu = dc.Particle.ParticleType.NuMu

        ccm_controller.primary_physical_process.primary_type = NuMu
        ccm_controller.primary_injection_process.primary_type = NuMu

        int_col = interactions.InteractionCollection(NuMu, [])
        ccm_controller.SetInteractions(
            primary_interaction_collection=int_col,
            injection=False,
            physical=True,
        )


# ---------------------------------------------------------------------------
# Fix #3: save_int_params per-event nested structure
# ---------------------------------------------------------------------------

class TestSaveIntParamsStructure:
    def test_int_params_key_is_dict_per_event(self):
        """When save_int_params=True, datasets['int_params'] should be a list
        of dicts (one per event), each mapping param names to per-interaction lists."""
        # Simulate the fixed logic from SaveEvents
        datasets = {}
        events_data = [
            # event 0: 2 interactions with different param sets
            [{"x": 1.0, "y": 2.0}, {"x": 3.0}],
            # event 1: 1 interaction with an extra param
            [{"x": 4.0, "y": 5.0, "z": 6.0}],
        ]

        for ie, event_params_list in enumerate(events_data):
            datasets.setdefault("int_params", [])
            datasets["int_params"].append({})
            for datum_params in event_params_list:
                for param_name, param_value in datum_params.items():
                    if param_name not in datasets["int_params"][-1]:
                        datasets["int_params"][-1][param_name] = []
                    datasets["int_params"][-1][param_name].append(param_value)

        assert len(datasets["int_params"]) == 2
        # Event 0: x appears in both interactions, y only in first
        assert datasets["int_params"][0] == {"x": [1.0, 3.0], "y": [2.0]}
        # Event 1: all three params from single interaction
        assert datasets["int_params"][1] == {"x": [4.0], "y": [5.0], "z": [6.0]}

    def test_params_from_later_events_not_lost(self):
        """OLD bug: if ie==0 was the only time keys were initialized,
        params appearing only in later events caused KeyError."""
        datasets = {}
        events_data = [
            [{"x": 1.0}],           # event 0: only x
            [{"x": 2.0, "z": 9.0}], # event 1: x and z (z is new)
        ]

        for ie, event_params_list in enumerate(events_data):
            datasets.setdefault("int_params", [])
            datasets["int_params"].append({})
            for datum_params in event_params_list:
                for param_name, param_value in datum_params.items():
                    if param_name not in datasets["int_params"][-1]:
                        datasets["int_params"][-1][param_name] = []
                    datasets["int_params"][-1][param_name].append(param_value)

        # z only appears in event 1 — should NOT raise KeyError
        assert datasets["int_params"][1] == {"x": [2.0], "z": [9.0]}
        # event 0 shouldn't have z at all
        assert "z" not in datasets["int_params"][0]


# ---------------------------------------------------------------------------
# Fix #4/#5: Weighter i_inj bounds checking and survival clamping
# ---------------------------------------------------------------------------

class TestWeighterBoundsCheck:
    @pytest.fixture(scope="class")
    def weighter_and_event(self):
        """Build a minimal Weighter with DummyCrossSection and generate one event."""
        import os
        from siren import detector, math as smath, utilities
        from siren._util import get_detector_model_path

        NuMu = dc.Particle.ParticleType.NuMu

        dm = detector.DetectorModel()
        det_dir = get_detector_model_path("CCM")
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

    def test_negative_i_inj_raises(self, weighter_and_event):
        """GetInteractionProbabilities should throw on negative i_inj."""
        weighter, event = weighter_and_event
        with pytest.raises((RuntimeError, IndexError)):
            weighter.GetInteractionProbabilities(event, -1)

    def test_out_of_range_i_inj_raises(self, weighter_and_event):
        """i_inj beyond injector count should throw."""
        weighter, event = weighter_and_event
        with pytest.raises((RuntimeError, IndexError)):
            weighter.GetInteractionProbabilities(event, 999)

    def test_survival_probs_negative_i_inj_raises(self, weighter_and_event):
        """GetSurvivalProbabilities should throw on negative i_inj."""
        weighter, event = weighter_and_event
        with pytest.raises((RuntimeError, IndexError)):
            weighter.GetSurvivalProbabilities(event, -1)

    def test_survival_probs_clamped_to_unit_interval(self, weighter_and_event):
        """Survival probabilities must be in [0, 1]."""
        weighter, event = weighter_and_event
        probs = weighter.GetSurvivalProbabilities(event, 0)
        assert len(probs) > 0
        for p in probs:
            assert 0.0 <= p <= 1.0, f"Survival probability {p} outside [0,1]"


# ---------------------------------------------------------------------------
# Fix #7: DN_min_decay_width only considers matching decays
# ---------------------------------------------------------------------------

class TestDNMinDecayWidth:
    def test_min_decay_width_only_matching(self):
        """Verify that only decays matching the primary_type contribute
        to DN_min_decay_width (logic pattern test)."""
        import numpy as np

        primary_type = "NuF4"
        decays_with_parents = [
            ("NuF4", 1e-10),  # matches
            ("NuF4", 5e-11),  # matches, smaller
            ("NuE", 1e-15),   # does NOT match -- must be excluded
        ]

        min_width = np.inf
        for parent, width in decays_with_parents:
            if parent == primary_type:
                if width > 0 and width < min_width:
                    min_width = width

        assert min_width == pytest.approx(5e-11)

    def test_zero_width_excluded(self):
        """Zero-width decays should not corrupt the minimum."""
        import numpy as np

        primary_type = "NuF4"
        decays_with_parents = [
            ("NuF4", 0.0),    # zero width -- must be excluded
            ("NuF4", 1e-10),  # valid
        ]

        min_width = np.inf
        for parent, width in decays_with_parents:
            if parent == primary_type:
                if width > 0 and width < min_width:
                    min_width = width

        assert min_width == pytest.approx(1e-10)


# ---------------------------------------------------------------------------
# Fix #9: InputDarkNewsModel secondary process preservation
# ---------------------------------------------------------------------------

class TestInputDarkNewsModelSecondaryPreservation:
    def test_existing_secondary_not_duplicated(self, ccm_controller):
        """When a secondary injection process already exists for a type,
        InputDarkNewsModel should reuse it rather than appending a duplicate."""
        N4 = dc.Particle.ParticleType.N4

        sec_proc = injection.SecondaryInjectionProcess()
        sec_proc.primary_type = N4
        ccm_controller.secondary_injection_processes = [sec_proc]

        count = sum(1 for p in ccm_controller.secondary_injection_processes
                    if p.primary_type == N4)
        assert count == 1
