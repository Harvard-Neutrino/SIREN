"""Whole-event physical factors and final Python weight validation."""

import copy
from fractions import Fraction
import math
import pickle

import numpy as np
import pytest

import siren
from siren import distributions as d, injection
from siren.errors import NotSerializableError, WeightCalculationError
from test_weighter_inheritance import _injector
from test_vertex_spec import _KeepAliveCrossSection


def _kinematic_factor(tree):
    first, last = tree.tree[0].record, tree.tree[-1].record
    return 1.0 + first.primary_momentum[0] / (1.0 + last.primary_momentum[0])


class _Factor:
    def __init__(self, value):
        self.value = value
        self.calls = []

    def __call__(self, tree):
        self.calls.append(tree)
        return self.value


class _SlottedWeighter(injection.Weighter):
    __slots__ = ("marker",)


@pytest.fixture
def cascade():
    inj = _injector()
    trees = inj.generate(5, on_failure="raise", on_shortfall="raise")
    assert all(len(tree.tree) == 2 for tree in trees)
    return inj, trees


@pytest.mark.parametrize("form", ["spec", "legacy", "raw", "pooled"])
def test_factor_agrees_across_scalar_batch_and_explanation(cascade, form):
    inj, trees = cascade
    reference = injection.Weighter(inj)
    kwargs = {}
    args = (inj,)
    if form == "legacy":
        args = ()
        kwargs = dict(injectors=[inj], detector_model=reference.detector_model,
                      primary_type=reference.primary_type,
                      primary_interactions=reference.primary_interactions,
                      primary_physical_distributions=reference.primary_physical_distributions,
                      secondary_interactions=reference.secondary_interactions,
                      secondary_physical_distributions=reference.secondary_physical_distributions)
    elif form == "raw":
        args = (inj.engine,)
        kwargs = dict(primary_physical=reference.primary_physical_distributions,
                      secondary_physical=reference.secondary_physical_distributions)
    elif form == "pooled":
        second = _injector()
        second.generate(3, on_failure="raise", on_shortfall="raise")
        args = (inj, second)
    base = injection.Weighter(*args, **kwargs)
    calls = []

    def factor(tree):
        calls.append(tree)
        return _kinematic_factor(tree)

    corrected = injection.Weighter(*args, event_factor=factor, **kwargs)
    expected = base.weight_all(trees) * [_kinematic_factor(t) for t in trees]
    for evaluate in (lambda: [corrected(t) for t in trees],
                     lambda: [corrected.event_weight(t) for t in trees],
                     lambda: corrected.weight_all(iter(trees)),
                     lambda: [corrected.explain(t).total for t in trees]):
        calls.clear()
        np.testing.assert_allclose(evaluate(), expected, rtol=1e-12, atol=0)
        assert len(calls) == len(trees)
        assert all(a is b for a, b in zip(calls, trees))
    bd = corrected.explain(trees[0])
    original = base.explain(trees[0])
    assert bd.base_total == original.total
    assert bd.event_factor == _kinematic_factor(trees[0])
    assert bd.flags == []
    assert [(v.generation, v.physical, v.flags) for v in bd.vertices] == [
        (v.generation, v.physical, v.flags) for v in original.vertices]
    assert "event_factor=" in str(bd)
    assert bd.culprit() is None
    with pytest.warns(DeprecationWarning):
        assert corrected.breakdown(trees[0]).total == bd.total
    assert corrected.engine.EventWeight(trees[0]) == base(trees[0])


def test_default_identity_and_replacement_leave_base_weight_unchanged(cascade):
    inj, trees = cascade
    w = injection.Weighter(inj)
    assert w.event_factor is None
    base = w.weight_all(trees)
    engine = w.engine
    w.event_factor = _Factor(np.float64(1.0))
    np.testing.assert_array_equal(w.weight_all(trees), base)
    w.event_factor = _Factor(2.0)
    np.testing.assert_array_equal(w.weight_all(trees), 2 * base)
    w.event_factor = None
    np.testing.assert_array_equal(w.weight_all(trees), base)
    assert w.engine is engine
    assert w.explain(trees[0]).event_factor is None


@pytest.mark.parametrize("value", [2.0, {}, [lambda tree: 1.0]])
def test_factor_configuration_requires_one_callable(value):
    with pytest.raises(TypeError, match="event_factor"):
        injection.Weighter(event_factor=value)
    w = injection.Weighter(event_factor=_kinematic_factor)
    with pytest.raises(TypeError, match="event_factor"):
        w.event_factor = value
    assert w.event_factor is _kinematic_factor


@pytest.mark.parametrize("value", [-1.0, float("nan"), float("inf"), -float("inf"),
                                   1j, None, "1", [1.0], np.array([1.0]),
                                   pytest.param(Fraction(-1, 10**400), id="negative-underflow")])
def test_invalid_factor_raises_and_explanation_retains_base(cascade, value):
    inj, trees = cascade
    factor = _Factor(value)
    w = injection.Weighter(inj, event_factor=factor)
    base = injection.Weighter(inj)(trees[0])
    with pytest.raises(WeightCalculationError, match="event factor"):
        w(trees[0])
    assert len(factor.calls) == 1
    factor.calls.clear()
    with pytest.raises(WeightCalculationError, match="Event 0.*event factor") as exc:
        w.weight_all(trees)
    assert exc.value.event_index == 0
    assert exc.value.event is trees[0]
    assert len(factor.calls) == 1
    bd = w.explain(trees[0])
    assert bd.base_total == pytest.approx(base, rel=1e-12, abs=0)
    assert math.isnan(bd.total)
    assert any("event factor" in flag for flag in bd.flags)
    assert "unusable" in str(bd)
    assert bd.culprit() is None  # The failure belongs to the event factor.


@pytest.mark.parametrize("zero_base", [False, True])
@pytest.mark.parametrize("zero", [0.0, -0.0])
def test_legitimate_zero_is_usable_and_factor_is_evaluated(cascade, zero_base, zero):
    inj, trees = cascade
    if zero_base:
        inj.primary.physical = [d.NormalizationConstant(zero)]
    factor = _Factor(2.0 if zero_base else zero)
    w = injection.Weighter(inj, event_factor=factor)
    assert w(trees[0]) == 0.0
    bd = w.explain(trees[0])
    assert bd.total == 0.0
    assert bd.flags == []
    assert bd.culprit() is None
    assert "unusable" not in str(bd)
    assert len(factor.calls) == 2


@pytest.mark.parametrize("value", [-1.0, float("nan"),
                                   pytest.param(Fraction(-1, 10**400), id="negative-underflow")])
def test_zero_base_does_not_hide_an_invalid_factor(cascade, value):
    inj, trees = cascade
    inj.primary.physical = [d.NormalizationConstant(0.0)]
    w = injection.Weighter(inj, event_factor=_Factor(value))
    with pytest.raises(WeightCalculationError, match="event factor"):
        w(trees[0])
    bd = w.explain(trees[0])
    assert bd.base_total == 0
    assert math.isnan(bd.total) and bd.flags


def test_corrected_weight_overflow_is_rejected(cascade):
    inj, trees = cascade
    inj.primary.physical = [d.NormalizationConstant(1e100)]
    w = injection.Weighter(inj, event_factor=_Factor(1e308))
    assert 1 < injection.Weighter(inj)(trees[0]) < float("inf")
    with pytest.raises(WeightCalculationError, match="corrected event weight"):
        w(trees[0])
    bd = w.explain(trees[0])
    assert math.isfinite(bd.base_total)
    assert bd.event_factor == 1e308
    assert math.isnan(bd.total)
    assert any("corrected event weight" in f for f in bd.flags)


def test_invalid_native_weight_fails_before_factor_evaluation(cascade):
    inj, trees = cascade
    inj.primary.physical = [d.NormalizationConstant(-1)]
    factor = _Factor(0.0)
    w = injection.Weighter(inj, event_factor=factor)
    with pytest.raises(WeightCalculationError):
        w(trees[0])
    bd = w.explain(trees[0])
    assert math.isnan(bd.total)
    assert bd.culprit() is not None
    assert factor.calls == []


def test_callback_exceptions_propagate_without_retry(cascade):
    inj, trees = cascade
    calls = []
    failure = ValueError("amplitude calculation failed")

    def factor(tree):
        calls.append(tree)
        raise failure

    w = injection.Weighter(inj, event_factor=factor)
    for evaluate in (lambda: w(trees[0]), lambda: w.weight_all(trees),
                     lambda: w.explain(trees[0])):
        calls.clear()
        with pytest.raises(ValueError) as exc:
            evaluate()
        assert exc.value is failure and len(calls) == 1


@pytest.mark.parametrize("value", [-1.0, float("nan"), float("inf"), -float("inf"), [1.0],
                                   pytest.param(Fraction(-1, 10**400), id="negative-underflow")])
def test_batch_validates_subclass_results_and_identifies_event(cascade, value):
    _, trees = cascade

    class CustomWeighter(injection.Weighter):
        def __init__(self):
            super().__init__()
            self.calls = 0

        def __call__(self, tree):
            self.calls += 1
            return value if self.calls == 3 else 0.0

    w = CustomWeighter()
    with pytest.raises(WeightCalculationError, match="Event 2") as exc:
        w.weight_all(iter(trees))
    assert exc.value.event_index == 2
    assert exc.value.event is trees[2]
    assert w.calls == 3
    np.testing.assert_array_equal(w.weight_all([]), np.array([], dtype=float))


def test_generate_uses_corrected_weights_and_propagates_validation(cascade):
    inj, _ = cascade
    factor = _Factor(2.0)
    w = injection.Weighter(inj, event_factor=factor)
    result = siren.generate(inj, w, events=3, on_failure="raise", on_shortfall="raise")
    assert len(factor.calls) == 3
    np.testing.assert_allclose(result.weights,
                               2 * injection.Weighter(inj).weight_all(result.events),
                               rtol=1e-12, atol=0)
    factor.value = 0.0
    result = siren.generate(inj, w, events=1, on_failure="raise", on_shortfall="raise")
    assert result.weights == [0.0] and len(result.events) == 1
    factor.value = float("nan")
    with pytest.raises(WeightCalculationError, match="Event 0.*event factor"):
        siren.generate(inj, w, events=1, on_failure="raise", on_shortfall="raise")


@pytest.mark.parametrize("values", [[float("nan")], [-1.0], [float("inf")], [], [1.0, 2.0],
                                    pytest.param([Fraction(-1, 10**400)], id="negative-underflow")])
def test_generate_validates_custom_batch_results(cascade, values):
    inj, _ = cascade

    class CustomBatch:
        def weight_all(self, trees):
            return values

    with pytest.raises(WeightCalculationError):
        siren.generate(inj, CustomBatch(), events=1,
                       on_failure="raise", on_shortfall="raise")


@pytest.mark.parametrize("factor", [_kinematic_factor, lambda tree: 1.0, _Factor(1.0)])
@pytest.mark.parametrize("method", ["save", "pickle"])
def test_factor_serialization_is_explicitly_unsupported(tmp_path, factor, method):
    w = injection.Weighter(event_factor=factor)
    with pytest.raises(NotSerializableError, match="event_factor"):
        if method == "save":
            w.save(str(tmp_path / "weight"))
        else:
            pickle.dumps(w)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("offender", ["expansion", "primary_interaction", "primary_distribution",
                                      "secondary_interaction", "secondary_distribution"])
def test_save_checks_lazy_injector_state_after_compilation(tmp_path, offender):
    class PythonMass(d.PrimaryMass):
        pass

    class PythonPosition(d.SecondaryPhysicalVertexDistribution):
        pass

    chain = offender == "expansion" or offender.startswith("secondary")
    inj = _injector(chain=chain)
    if offender == "expansion":
        inj.stopping_condition = None
        for vertex in [inj.primary, *inj.secondaries]:
            vertex.expand = [siren.expand.depth_below(1)]
        expected = "stopping condition"
    else:
        vertex = inj.secondaries[0] if chain else inj.primary
        if offender.endswith("interaction"):
            vertex.interactions = [_KeepAliveCrossSection()]
            expected = "interaction.*Python subclass"
        else:
            vertex.distributions[0] = PythonPosition() if chain else PythonMass(0)
            expected = "distribution.*Python subclass"

    # Spec fields have not yet been compiled into the legacy injector fields.
    assert inj.primary_interactions == []
    assert inj.primary_injection_distributions == []
    w = injection.Weighter(
        injectors=[inj], detector_model=inj.detector_model,
        primary_type=siren.particles.NuMu,
        primary_interactions=[siren.interactions.DummyCrossSection()],
        secondary_interactions={siren.particles.NuMu: [siren.interactions.DummyCrossSection()]}
        if chain else {})
    with pytest.raises(NotSerializableError, match=expected):
        w.save(str(tmp_path / "weight"))
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("factor", [_kinematic_factor, lambda tree: 1.0, _Factor(2.0)])
def test_shallow_copy_preserves_factor_weights_and_subclass_state(cascade, tmp_path, factor):
    inj, trees = cascade
    w = _SlottedWeighter(inj, event_factor=factor)
    w.marker = {"owner": w, "values": [1]}
    w.alias = w.marker["values"]
    expected = w.weight_all(trees)
    clone = copy.copy(w)

    assert type(clone) is type(w) and clone is not w
    assert clone.alias is clone.marker["values"]
    assert clone.engine is w.engine
    assert clone.marker is w.marker
    assert clone.event_factor is w.event_factor
    np.testing.assert_allclose(clone.weight_all(trees), expected, rtol=1e-12, atol=0)
    clone.event_factor = _Factor(0.0)
    np.testing.assert_array_equal(clone.weight_all(trees), np.zeros(len(trees)))
    np.testing.assert_allclose(w.weight_all(trees), expected, rtol=1e-12, atol=0)
    with pytest.raises(NotSerializableError, match="event_factor"):
        pickle.dumps(clone)
    with pytest.raises(NotSerializableError, match="event_factor"):
        clone.save(str(tmp_path / "weight"))
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("kind", ["function", "lambda", "instance"])
def test_deepcopy_preserves_callback_and_python_subclass_state(cascade, kind):
    _, trees = cascade
    factor = {"function": _kinematic_factor, "lambda": lambda tree: 1.0,
              "instance": _Factor(2.0)}[kind]
    # Native objects retain their own copy limitations; pure Python state can
    # be copied before configuring a detector or injector.
    w = _SlottedWeighter(event_factor=factor)
    w.marker = {"owner": w, "values": [1]}
    w.alias = w.marker["values"]
    clone = copy.deepcopy(w)

    assert type(clone) is type(w) and clone is not w
    assert clone.marker is not w.marker
    assert clone.marker["owner"] is clone
    assert clone.alias is clone.marker["values"]
    clone.alias.append(2)
    assert w.alias == [1]
    assert clone.event_factor(trees[0]) == factor(trees[0])
    if kind == "instance":
        assert clone.event_factor is not factor
        clone.event_factor.value = 0.0
        assert factor.value == 2.0
    else:
        assert clone.event_factor is factor
    with pytest.raises(NotSerializableError, match="event_factor"):
        pickle.dumps(clone)


def test_loading_a_base_archive_requires_an_unconfigured_factor(tmp_path):
    inj = _injector(chain=False)
    trees = inj.generate(3, on_failure="raise", on_shortfall="raise")
    base = injection.Weighter(inj)
    filename = str(tmp_path / "weight")
    base.save(filename)
    w = injection.Weighter(event_factor=_kinematic_factor)
    with pytest.raises(NotSerializableError, match="event_factor"):
        w.load(filename)
    assert w.event_factor is _kinematic_factor
    w.event_factor = None
    w.load(filename)
    np.testing.assert_allclose(w.weight_all(trees), base.weight_all(trees), rtol=1e-12, atol=0)
    w.event_factor = _kinematic_factor
    np.testing.assert_allclose(w.weight_all(trees),
                               base.weight_all(trees) * [_kinematic_factor(t) for t in trees],
                               rtol=1e-12, atol=0)


def test_pickle_without_factor_and_pre_factor_state_remain_supported():
    restored = pickle.loads(pickle.dumps(injection.Weighter()))
    assert restored.event_factor is None
    del restored.__dict__["_Weighter__event_factor"]
    restored = pickle.loads(pickle.dumps(restored))
    assert restored.event_factor is None


def test_pickle_preserves_existing_subclass_slots():
    w = _SlottedWeighter()
    w.marker = "subclass state"
    restored = pickle.loads(pickle.dumps(w))
    assert restored.marker == w.marker and restored.event_factor is None
