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
#  positional Box raises                                              #
# ------------------------------------------------------------------ #

def test_positional_box_raises():
    """Box(x, y, z) raises with the widths=/center= fix-it."""
    with pytest.raises(siren.errors.ConfigurationError, match="widths="):
        siren.geometry.Box(1.0, 2.0, 3.0)


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


# ------------------------------------------------------------------ #
#  Secondary trampoline serialization guard                           #
# ------------------------------------------------------------------ #

class _PySecondaryCrossSection(interactions.CrossSection):
    """Data-free Python-subclassed CrossSection (a non-serializable trampoline)."""

    def __init__(self):
        interactions.CrossSection.__init__(self)

    def TotalCrossSection(self, *args, **kwargs):
        return 1.0

    def DifferentialCrossSection(self, *args):
        return 1.0

    def InteractionThreshold(self, *args):
        return 0.0

    def SampleFinalState(self, *args):
        pass

    def GetPossibleTargets(self):
        return [siren.particles.Nucleon]

    def GetPossibleTargetsFromPrimary(self, primary_type):
        return [siren.particles.Nucleon]

    def GetPossiblePrimaries(self):
        return [_NuMu]

    def GetPossibleSignatures(self):
        sig = dc.InteractionSignature()
        sig.primary_type = _NuMu
        sig.target_type = siren.particles.Nucleon
        sig.secondary_types = [_NuMu, siren.particles.Nucleon]
        return [sig]

    def GetPossibleSignaturesFromParents(self, primary_type, target_type):
        return self.GetPossibleSignatures()

    def FinalStateProbability(self, *args):
        return 1.0

    def DensityVariables(self):
        return []

    def equal(self, other):
        return isinstance(other, _PySecondaryCrossSection)


def test_save_secondary_trampoline_raises_not_serializable(tmp_path):
    """save() refuses a Python-subclassed SECONDARY interaction (not just primary).

    A bare detector suffices: the guard reads the Python-side secondary
    interaction dict before the engine's own archiving is reached.
    """
    dm = detector.DetectorModel()
    inj = injection.Injector(
        detector=dm,
        number_of_events=2,
        seed=1,
        primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_distributions(0.0),
        secondary_interactions={_NuMu: [_PySecondaryCrossSection()]},
        secondary_injection_distributions={
            _NuMu: [distributions.SecondaryPhysicalVertexDistribution()]},
    )
    inj._build()
    with pytest.raises(siren.errors.NotSerializableError) as exc:
        inj.save(str(tmp_path / "inj"))
    assert any("secondary" in o for o in exc.value.offenders)


def test_save_stopping_condition_raises_not_serializable(tmp_path):
    """save() refuses an injector with a stopping condition (not archived)."""
    dm = detector.DetectorModel()
    inj = injection.Injector(
        detector=dm,
        number_of_events=2,
        seed=1,
        primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_distributions(0.0),
        stopping_condition=lambda datum, i: True,
    )
    inj._build()
    with pytest.raises(siren.errors.NotSerializableError) as exc:
        inj.save(str(tmp_path / "inj"))
    assert any("stopping condition" in o for o in exc.value.offenders)


# ------------------------------------------------------------------ #
#  Duplicate secondary Vertex particle type                           #
# ------------------------------------------------------------------ #

def test_duplicate_secondary_vertex_type_raises():
    """Two secondary vertices resolving to the same particle type is loud.

    Both would otherwise silently overwrite each other's interactions and
    distributions in the by-type maps.
    """
    from siren.vertex import Vertex

    pxs = interactions.DummyCrossSection()
    psig = pxs.GetPossibleSignatures()[0]
    primary = Vertex(psig.primary_type, pxs,
                     distributions=_distributions(0.0))

    def _sec(label):
        xs = interactions.DummyCrossSection()
        sig = xs.GetPossibleSignatures()[0]
        return Vertex(
            sig.primary_type, xs,
            distributions=[distributions.SecondaryPhysicalVertexDistribution()],
            label=label)

    va, vb = _sec("a"), _sec("b")
    # Both resolve to DummyCrossSection's signature primary type.
    assert va._resolved_particle == vb._resolved_particle

    inj = injection.Injector(
        detector=detector.DetectorModel(),
        primary=primary,
        secondaries=(va, vb),
        events=2,
        seed=1,
    )
    with pytest.raises(siren.errors.ConfigurationError, match="same particle type"):
        inj._build()


# ------------------------------------------------------------------ #
#  Phase-space probe propagates typed SIREN errors                    #
# ------------------------------------------------------------------ #

def _process_with_phase_space():
    xs = interactions.DummyCrossSection()
    sig = xs.GetPossibleSignatures()[0]
    ic = interactions.InteractionCollection(sig.primary_type, [xs])
    proc = injection.PrimaryInjectionProcess(sig.primary_type, ic)
    mc = injection.MultiChannelPhaseSpace()
    mc.channels = [injection.Isotropic2BodyChannel(0)]
    mc.weights = [1.0]
    proc.SetPhaseSpace(sig, mc)
    return proc


class _RaisingValidation:
    """Stub whose probe_channel_densities raises a chosen exception."""

    def __init__(self, exc):
        self._exc = exc

    def probe_channel_densities(self, mcps, sig, detector_model, rng):
        raise self._exc


def test_probe_reraises_typed_configuration_error():
    """A ConfigurationError from the probe propagates (it is a real fault)."""
    inj = injection.Injector(
        detector=detector.DetectorModel(),
        number_of_events=1,
        primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_distributions(0.0),
    )
    proc = _process_with_phase_space()
    stub = _RaisingValidation(
        siren.errors.ConfigurationError("bad model configuration"))
    rng = utilities.SIREN_random(1)
    with pytest.raises(siren.errors.ConfigurationError, match="bad model"):
        inj._probe_phase_spaces(stub, proc, [], rng)


def test_probe_swallows_plain_probe_inapplicability():
    """A bare ValueError/RuntimeError means the probe could not run: swallowed."""
    inj = injection.Injector(
        detector=detector.DetectorModel(),
        number_of_events=1,
        primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_distributions(0.0),
    )
    proc = _process_with_phase_space()
    rng = utilities.SIREN_random(1)
    # A plain ValueError (synthetic template unsolvable) is skipped, not raised.
    inj._probe_phase_spaces(
        _RaisingValidation(ValueError("unsolvable template")), proc, [], rng)
    # A plain RuntimeError that is not a typed SIREN error is likewise skipped.
    inj._probe_phase_spaces(
        _RaisingValidation(RuntimeError("engine hiccup")), proc, [], rng)
