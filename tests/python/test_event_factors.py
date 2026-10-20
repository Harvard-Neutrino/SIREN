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


@pytest.mark.parametrize("factor", [0.0, 1.0, 1.75])
def test_scalar_and_batch_physical_factor(cascade, factor):
    inj, trees = cascade
    base = injection.Weighter(inj).weight_all(trees)
    calls = _Factor(factor)
    weighted = injection.Weighter(inj, event_factor=calls)
    np.testing.assert_array_equal([weighted(t) for t in trees], base * factor)
    np.testing.assert_array_equal(weighted.weight_all(trees), base * factor)
    assert len(calls.calls) == 2 * len(trees)


@pytest.mark.parametrize("factor", [-1.0, float("nan"), float("inf")])
def test_invalid_factor_is_not_hidden_by_zero_physical_weight(cascade, factor):
    inj, trees = cascade
    weighted = injection.Weighter(
        inj, primary_physical=[d.NormalizationConstant(0)],
        event_factor=lambda tree: factor,
    )
    with pytest.raises(WeightCalculationError, match="event factor"):
        weighted(trees[0])






@pytest.mark.parametrize("value", [2.0, {}, [lambda tree: 1.0]])
def test_factor_configuration_requires_one_callable(value):
    with pytest.raises(TypeError, match="event_factor"):
        injection.Weighter(event_factor=value)
    w = injection.Weighter(event_factor=_kinematic_factor)
    with pytest.raises(TypeError, match="event_factor"):
        w.event_factor = value
    assert w.event_factor is _kinematic_factor














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
            # Pure declarative expansion is native and archived. An arbitrary
            # Python predicate still cannot be preserved by a native archive.
            vertex.continue_if = lambda tree, parent, i: True
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
