"""Layer-2 Injector: generate() success counting, shortfall policy, reset,
save guard, positional-Box warning, and legacy __iter__ semantics.

Data-free where possible; anything needing detector data files skips cleanly.
"""

import os
import warnings

import pytest

import siren
from siren import dataclasses as dc
from siren import injection
from siren import interactions
from siren import distributions
from siren import detector
from siren import math as smath
from siren import utilities
from siren import _util

_NuMu = dc.Particle.ParticleType.NuMu


def _skip_unless_ccm_data():
    try:
        det_dir = _util.get_detector_model_path("CCM")
    except ValueError as e:
        pytest.skip("CCM detector model path unavailable: {}".format(e))
    missing = [p for p in (os.path.join(det_dir, "materials.dat"),
                           os.path.join(det_dir, "densities.dat"))
               if not os.path.exists(p)]
    if missing:
        pytest.skip("CCM detector data missing: " + ", ".join(missing))
    return det_dir


def _load_ccm_detector():
    dm = detector.DetectorModel()
    det_dir = _util.get_detector_model_path("CCM")
    dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
    dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    return dm


def _distributions(max_distance):
    return [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(
            smath.Vector3D(0, 0, 0), max_distance),
    ]


def _working_injector(events=5, seed=11, max_distance=25.0, max_attempts=None):
    _skip_unless_ccm_data()
    dm = _load_ccm_detector()
    inj = injection.Injector(
        detector=dm,
        primary=None,
        number_of_events=events,
        seed=seed,
        primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_distributions(max_distance),
        max_attempts=max_attempts,
    )
    return inj


def _forced_failure_injector(events=5, seed=11, max_attempts=None):
    # A bare detector plus max_distance=0 makes every path miss its target,
    # so every attempt fails with NoTargetsOnPath -- no data files needed.
    dm = detector.DetectorModel()
    inj = injection.Injector(
        detector=dm,
        number_of_events=events,
        seed=seed,
        primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_distributions(0.0),
        max_attempts=max_attempts,
    )
    return inj


# ------------------------------------------------------------------ #
#  generate() counts successes                                        #
# ------------------------------------------------------------------ #

def test_generate_counts_successes_not_attempts():
    """generate(N) returns exactly N non-empty trees."""
    inj = _working_injector(events=10)
    trees = inj.generate(10, on_shortfall="raise")
    assert len(trees) == 10
    assert all(len(t.tree) > 0 for t in trees)


def test_generate_default_events_uses_configured_count():
    """generate() with no argument uses the configured event count."""
    inj = _working_injector(events=4)
    trees = inj.generate(on_shortfall="raise")
    assert len(trees) == 4


# ------------------------------------------------------------------ #
#  on_shortfall triple                                                #
# ------------------------------------------------------------------ #

def test_on_shortfall_raise():
    """A shortfall with on_shortfall='raise' raises InjectionShortfall."""
    inj = _forced_failure_injector(events=3, max_attempts=20)
    with pytest.raises(siren.errors.InjectionShortfall) as exc:
        inj.generate(3, on_shortfall="raise")
    assert exc.value.report is not None
    assert exc.value.report.successes == 0


def test_on_shortfall_warn():
    """A shortfall with on_shortfall='warn' warns and returns what succeeded."""
    inj = _forced_failure_injector(events=3, max_attempts=20)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        trees = inj.generate(3, on_shortfall="warn")
    assert len(trees) == 0
    assert any(isinstance(x.message, siren.errors.InjectionShortfall) for x in w)


def test_on_shortfall_ignore():
    """A shortfall with on_shortfall='ignore' is silent."""
    inj = _forced_failure_injector(events=3, max_attempts=20)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        trees = inj.generate(3, on_shortfall="ignore")
    assert len(trees) == 0


# ------------------------------------------------------------------ #
#  min_efficiency early abort                                         #
# ------------------------------------------------------------------ #

def test_min_efficiency_early_abort_raises_with_report():
    """Efficiency below min_efficiency aborts with an InjectionShortfall report."""
    inj = _forced_failure_injector(events=100, max_attempts=100000)
    with pytest.raises(siren.errors.InjectionShortfall) as exc:
        inj.generate(100, on_shortfall="raise", min_efficiency=0.5)
    assert exc.value.report is not None
    assert exc.value.report.dominant() is not None


# ------------------------------------------------------------------ #
#  reset / reuse                                                      #
# ------------------------------------------------------------------ #

def test_reset_allows_reuse():
    """reset() clears counters so a second generate() starts fresh."""
    inj = _working_injector(events=5)
    first = inj.generate(5, on_shortfall="raise")
    assert len(first) == 5
    inj.reset()
    assert inj.injected_events == 0
    second = inj.generate(5, on_shortfall="raise")
    assert len(second) == 5


# ------------------------------------------------------------------ #
#  typo'd secondary dict key raises                                   #
# ------------------------------------------------------------------ #

def test_typo_secondary_key_raises():
    """A secondary key present in one dict but not the interactions dict raises."""
    _skip_unless_ccm_data()
    dm = _load_ccm_detector()
    xs = interactions.DummyCrossSection()
    inj = injection.Injector(
        detector=dm,
        number_of_events=2,
        seed=1,
        primary_type=_NuMu,
        primary_interactions=[xs],
        primary_injection_distributions=_distributions(25.0),
        secondary_interactions={_NuMu: [interactions.DummyCrossSection()]},
        secondary_injection_distributions={_NuMu: [distributions.SecondaryPhysicalVertexDistribution()]},
        secondary_weighting_modes={dc.Particle.ParticleType.NuE: injection.VertexWeightingMode.Fixed()},
    )
    with pytest.raises((siren.errors.ConfigurationError, ValueError)):
        inj.generate(2, on_shortfall="ignore")


# ------------------------------------------------------------------ #
#  positional Box warns but still constructs                          #
# ------------------------------------------------------------------ #

def test_positional_box_warns_and_constructs():
    """Box(x, y, z) emits a DeprecationWarning and still builds a box."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        box = siren.geometry.Box(1.0, 2.0, 3.0)
    assert any(issubclass(x.category, DeprecationWarning) for x in w)
    assert box.X == 1.0 and box.Y == 2.0 and box.Z == 3.0


def test_widths_center_box_places_at_center():
    """Box(widths=..., center=...) places the box at the given center."""
    box = siren.geometry.Box(widths=(1.0, 1.0, 1.0), center=(0.0, 0.0, 5.0))
    assert box.IsInside(siren.math.Vector3D(0, 0, 5))
    assert not box.IsInside(siren.math.Vector3D(0, 0, 0))


# ------------------------------------------------------------------ #
#  __iter__ legacy semantics                                          #
# ------------------------------------------------------------------ #

def test_iter_warns_and_yields_attempt_count():
    """__iter__ warns and yields exactly number_of_events raw attempts."""
    inj = _working_injector(events=4)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        events = list(inj)
    assert any(issubclass(x.category, DeprecationWarning) for x in w)
    assert len(events) == 4


# ------------------------------------------------------------------ #
#  save() round-trip guard                                            #
# ------------------------------------------------------------------ #

def test_save_on_fixed_mode_raises_not_serializable(tmp_path):
    """save() refuses a Fixed-mode process the engine would not round-trip."""
    _skip_unless_ccm_data()
    dm = _load_ccm_detector()
    inj = injection.Injector(
        detector=dm,
        number_of_events=2,
        seed=1,
        primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_distributions(25.0),
        primary_weighting_mode=injection.VertexWeightingMode.Fixed(),
    )
    with pytest.raises(siren.errors.NotSerializableError):
        inj.save(str(tmp_path / "inj"))


def test_report_is_rendered_after_failures():
    """report() renders an attrition table with a dominant reason."""
    inj = _forced_failure_injector(events=3, max_attempts=10)
    inj.generate(3, on_shortfall="ignore")
    report = inj.report()
    text = str(report)
    assert "InjectionReport" in text
    assert report.dominant() is not None
    assert report.dominant().reason_name == "NoTargetsOnPath"
