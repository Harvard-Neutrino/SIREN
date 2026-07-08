"""Layer-2 Injector: generate() success counting, shortfall policy, reset,
save guard, positional-Box warning, and legacy __iter__ semantics.

Data-free where possible; anything needing detector data files skips cleanly.
"""

import os
import pickle
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


def _depth_ge_1(tree, datum, i):
    """Expand one chain level. Defined at module scope so pickle can carry it
    by reference (a lambda would not survive the pickle state tuple)."""
    return datum.depth(tree) >= 1


def _chain_injector(events=20, seed=7, max_attempts=8000):
    """A NuMu -> NuMu chain (DummyCrossSection at both vertices) whose stopping
    condition expands one secondary level, mirroring the raw _build_chain in
    test_density_breakdown. Uses CCM detector data; skips cleanly without it."""
    _skip_unless_ccm_data()
    dm = _load_ccm_detector()
    return injection.Injector(
        detector=dm,
        number_of_events=events,
        seed=seed,
        primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_distributions(25.0),
        secondary_interactions={_NuMu: [interactions.DummyCrossSection()]},
        secondary_injection_distributions={
            _NuMu: [distributions.SecondaryPhysicalVertexDistribution()]},
        stopping_condition=_depth_ge_1,
        max_attempts=max_attempts,
    )


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

def test_save_on_fixed_mode_succeeds(tmp_path):
    """save() archives a Fixed weighting mode with the process instead of
    refusing it; the reloaded injector carries the same mode."""
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
    path = str(tmp_path / "inj")
    inj.save(path)   # no NotSerializableError: the mode rides with the process
    loaded = injection.Injector.load(path)
    assert loaded.engine.GetPrimaryProcess().GetWeightingMode() == (
        injection.VertexWeightingMode.Fixed())


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
    interaction dict before the engine's own archiving is reached. A
    stopping_condition is supplied so _build() succeeds on its own terms;
    without one, secondaries with no expansion control are rejected before
    save() is ever reached.
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
        stopping_condition=lambda datum, i: True,
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
#  Secondaries declared with no expansion wiring                      #
# ------------------------------------------------------------------ #

def test_secondaries_without_expand_raises_at_build():
    """A spec-form chain that registers a secondary but declares no expand
    rule anywhere is rejected at build.

    Without an expand declaration (and no legacy stopping condition) the
    engine's default stopping condition prunes every secondary, silently
    collapsing the chain to primary-only trees. The expansion-wiring
    validation runs whenever a chain carries secondaries, so the
    misconfiguration raises a typed ConfigurationError.
    """
    from siren.vertex import Vertex

    pxs = interactions.DummyCrossSection()
    psig = pxs.GetPossibleSignatures()[0]
    primary = Vertex(psig.primary_type, pxs,
                     distributions=_distributions(0.0))
    # A self-scatter secondary of the same type as the primary/root, so it is
    # reachable-as-root and the only wiring fault is the absent expand rule
    # (case (c)); no expand is declared anywhere and no stopping condition is
    # supplied.
    secondary = Vertex(
        psig.primary_type, interactions.DummyCrossSection(),
        distributions=[distributions.SecondaryPhysicalVertexDistribution()])

    inj = injection.Injector(
        detector=detector.DetectorModel(),
        primary=primary,
        secondaries=(secondary,),
        events=2,
        seed=1,
    )
    with pytest.raises(siren.errors.ConfigurationError, match="expand"):
        inj._build()


def test_legacy_secondaries_without_stopping_condition_raises_at_build():
    """A legacy-dict chain with secondaries but no stopping_condition is
    rejected at build.

    The legacy secondary_interactions/secondary_injection_distributions dict
    form has no Vertex objects for validate_expansion_wiring to check, so this
    independent safety net is the only thing that catches it: without a
    stopping_condition, the engine's default stopping condition prunes every
    secondary, silently collapsing the chain to primary-only trees.
    """
    inj = injection.Injector(
        detector=detector.DetectorModel(),
        number_of_events=2,
        seed=1,
        primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_distributions(0.0),
        secondary_interactions={_NuMu: [interactions.DummyCrossSection()]},
        secondary_injection_distributions={
            _NuMu: [distributions.SecondaryPhysicalVertexDistribution()]},
    )
    with pytest.raises(siren.errors.ConfigurationError, match="stopping"):
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


# ------------------------------------------------------------------ #
#  Pickle round-trip: usable engine, RNG install, guard, stopping cond #
# ------------------------------------------------------------------ #

def test_pickle_roundtrip_installs_rng_data_free():
    """generate() after a pickle round-trip must not dereference a null RNG.

    Data-free: every attempt of a bare-detector, zero-distance injector still
    draws its primary variables from the RNG before missing its target, so a
    missing RNG install would crash here. It returns zero events, not a crash.
    """
    inj = _forced_failure_injector(events=3, seed=11, max_attempts=20)
    restored = pickle.loads(pickle.dumps(inj))
    trees = restored.generate(3, on_shortfall="ignore")
    assert trees == []


def test_pickle_roundtrip_plain_injector_generates():
    """A plain injector round-trips to a usable engine: the restored injector
    installs a fresh RNG (re-seeded from the stored seed) and generate()
    produces events instead of segfaulting on a null RNG."""
    inj = _working_injector(events=5, seed=11)
    restored = pickle.loads(pickle.dumps(inj))
    trees = restored.generate(2, on_shortfall="raise")
    assert len(trees) == 2
    assert all(len(t.tree) > 0 for t in trees)


def test_pickle_roundtrip_preserves_stopping_condition_on_engine():
    """A chain's stopping condition survives the round-trip and re-attaches to
    the engine, so secondaries are expanded (some tree has depth >= 1). Without
    the re-attach the engine default prunes every secondary and all trees
    collapse to a single vertex."""
    reference = _chain_injector(seed=7)
    ref_trees = reference.generate(20, on_shortfall="ignore")
    if not any(len(t.tree) >= 2 for t in ref_trees):
        pytest.skip("chain fixture produced no multi-vertex events to pin")

    inj = _chain_injector(seed=7)
    restored = pickle.loads(pickle.dumps(inj))
    trees = restored.generate(20, on_shortfall="ignore")
    assert any(len(t.tree) >= 2 for t in trees)


def test_pickle_roundtrips_phase_space_map():
    """Pickle preserves the facade and engine views of a phase-space map."""
    inj = _forced_failure_injector(events=2, seed=1)
    inj._build()
    xs = interactions.DummyCrossSection()
    sig = xs.GetPossibleSignatures()[0]
    mc = injection.MultiChannelPhaseSpace()
    mc.channels = [injection.Isotropic2BodyChannel(0)]
    mc.weights = [1.0]
    inj.primary_phase_spaces = {sig: mc}

    restored = pickle.loads(pickle.dumps(inj))
    restored_mc = restored.primary_phase_spaces[sig]
    engine_mc = restored.engine.GetPrimaryProcess().GetPhaseSpaceMap()[sig]
    assert list(restored_mc.weights) == [1.0]
    assert restored_mc is engine_mc


def test_pickle_allows_stopping_condition():
    """A stopping condition is not a pickle offender (the state tuple carries
    it), so pickling an injector whose only save()-offender is a stopping
    condition succeeds and the condition round-trips."""
    inj = injection.Injector(
        detector=detector.DetectorModel(),
        number_of_events=2,
        seed=1,
        primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_distributions(0.0),
        stopping_condition=_depth_ge_1,
    )
    inj._build()
    blob = pickle.dumps(inj)   # must not raise
    restored = pickle.loads(blob)
    assert restored.stopping_condition is _depth_ge_1


# ------------------------------------------------------------------ #
#  seed setter on a built engine                                      #
# ------------------------------------------------------------------ #

def test_seed_setter_after_build_reseeds_and_generates():
    """Setting .seed on a built injector must reseed generation, not raise.

    The raw engine exposes only SetRandom, so the setter installs a fresh
    engine seeded from the new value (streams restart, not resume). Data-free:
    the forced-failure fixture draws from the RNG on every attempt, so a
    broken setter or a missing engine would surface here."""
    inj = _forced_failure_injector(events=2, seed=1, max_attempts=20)
    inj._build()
    inj.seed = 42
    assert inj.seed == 42
    trees = inj.generate(2, on_shortfall="ignore")
    assert trees == []
    assert inj.injection_attempts > 0
