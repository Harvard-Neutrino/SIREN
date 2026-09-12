"""Vertex physical declarations reach weighting without changing generation."""

import gc
import pickle
import weakref

import numpy as np
import pytest

import siren
from siren import distributions as d, injection, interactions
from siren.errors import NotSerializableError
from test_injector_api import _distributions, _load_ccm_detector, _depth_ge_1
from test_vertex_spec import _KeepAliveCrossSection

PT = siren.dataclasses.ParticleType


def _injector(*, chain=True):
    primary = siren.Vertex(
        PT.NuMu, interactions.DummyCrossSection(),
        distributions=_distributions(25.0),
        physical=[d.NormalizationConstant(0.2)],
        physical_interactions=interactions.DummyCrossSection())
    secondary = siren.Vertex(
        PT.NuMu, interactions.DummyCrossSection(),
        position=d.SecondaryPhysicalVertexDistribution(),
        physical=[d.NormalizationConstant(0.7)],
        physical_interactions=interactions.DummyCrossSection())
    return injection.Injector(
        detector=_load_ccm_detector(), primary=primary,
        secondaries=[secondary] if chain else [],
        stopping_condition=_depth_ge_1 if chain else None,
        events=5, max_attempts=5000, seed=71)


def test_inherited_primary_and_secondary_weights_match_explicit_configuration():
    inj = _injector()
    p, s = inj.primary, inj.secondaries[0]
    inherited = injection.Weighter(inj)
    explicit = injection.Weighter(
        inj, primary_physical=p.physical,
        secondary_physical={PT.NuMu: s.physical},
        overrides={"primary_interactions": p.physical_interactions,
                   "secondary_interactions": {PT.NuMu: s.physical_interactions}})
    conditional = injection.Weighter(inj, primary_physical=[], secondary_physical={})
    trees = inj.generate(5, on_shortfall="raise", on_failure="raise")
    assert all(len(tree.tree) == 2 for tree in trees)
    expected = explicit.weight_all(trees)
    assert np.all(expected > 0)
    np.testing.assert_array_equal(inherited.weight_all(trees), expected)
    np.testing.assert_allclose(expected / conditional.weight_all(trees), 0.2 * 0.7,
                               rtol=1e-12, atol=0)
    for tree, weight in zip(trees, expected):
        assert inherited.explain(tree).total == pytest.approx(weight, rel=1e-12, abs=0)
    assert inherited.engine.GetPrimaryPhysicalProcess().interactions.GetCrossSections()[0] is p.physical_interactions[0]
    assert inherited.engine.GetSecondaryPhysicalProcesses()[0].interactions.GetCrossSections()[0] is s.physical_interactions[0]
    assert inj.engine.GetPrimaryProcess().interactions.GetCrossSections()[0] is p.interactions[0]
    assert inj.engine.GetSecondaryProcessMap()[PT.NuMu].interactions.GetCrossSections()[0] is s.interactions[0]


def test_explicit_overrides_replace_inheritance_including_empty_collections():
    inj = _injector()
    model = interactions.DummyCrossSection()
    w = injection.Weighter(inj, primary_physical=[], secondary_physical={},
                           overrides={"primary_interactions": [model],
                                      "secondary_interactions": {}})
    assert w.primary_physical_distributions == []
    assert w.secondary_physical_distributions == {}
    assert w.primary_interactions[0] is model
    assert w.secondary_interactions == {}


@pytest.mark.parametrize("raw", [False, True])
def test_legacy_and_raw_injectors_keep_sampling_model_defaults(raw):
    spec = _injector(chain=False)
    inj = injection.Injector(
        detector=spec.detector_model, events=5,
        primary_type=PT.NuMu, primary_interactions=spec.primary.interactions,
        primary_injection_distributions=spec.primary.distributions)
    w = injection.Weighter(inj.engine if raw else inj)
    assert w.primary_interactions[0] is spec.primary.interactions[0]
    assert w.primary_physical_distributions == []
    assert w.secondary_physical_distributions == {}


def test_none_inherits_sampling_models_and_first_injector_defines_pooled_target():
    first, second = _injector(), _injector()
    first.primary.physical_interactions = None
    first.secondaries[0].physical_interactions = None
    w = injection.Weighter(first, second)
    assert w.primary_interactions[0] is first.primary.interactions[0]
    assert w.secondary_interactions[PT.NuMu][0] is first.secondaries[0].interactions[0]
    assert w.primary_physical_distributions[0] is first.primary.physical[0]


def test_weighter_owns_copies_of_inherited_lists_and_python_models():
    inj = _injector()
    model = _KeepAliveCrossSection()
    ref = weakref.ref(model)
    inj.primary.physical_interactions = [model]
    inj.secondaries[0].physical_interactions = [model]
    w = injection.Weighter(inj)
    inj.primary.physical_interactions.clear()
    inj.primary.physical.clear()
    inj.secondaries[0].physical_interactions.clear()
    inj.secondaries[0].physical.clear()
    del model
    gc.collect()
    assert w.primary_interactions[0] is ref()
    assert w.secondary_interactions[PT.NuMu][0] is ref()
    assert ref().TotalCrossSection(None) == 1.0
    assert len(w.primary_physical_distributions) == 1
    assert len(w.secondary_physical_distributions[PT.NuMu]) == 1
    with pytest.raises(NotSerializableError, match="Python subclass"):
        w._guard_serializable(siren.errors)


@pytest.mark.parametrize("field", ["physical", "physical_interactions"])
@pytest.mark.parametrize("secondary", [False, True])
@pytest.mark.parametrize("archive", ["save", "pickle"])
def test_injector_refuses_to_drop_physical_declarations(tmp_path, field, secondary, archive):
    inj = _injector(chain=secondary)
    for vertex in [inj.primary, *inj.secondaries]:
        vertex.physical = []
        vertex.physical_interactions = None
    vertex = inj.secondaries[0] if secondary else inj.primary
    setattr(vertex, field, [d.NormalizationConstant(0.3)] if field == "physical"
            else [interactions.DummyCrossSection()])
    with pytest.raises(NotSerializableError, match="declares physical"):
        if archive == "save":
            inj.save(str(tmp_path / "injector"))
        else:
            pickle.dumps(inj)


def test_weighter_archive_preserves_inherited_physical_processes(tmp_path):
    inj = _injector(chain=False)
    trees = inj.generate(3, on_shortfall="raise", on_failure="raise")
    w = injection.Weighter(inj)
    expected = w.weight_all(trees)
    filename = str(tmp_path / "weighter")
    w.save(filename)
    restored = injection.Weighter()
    restored.load(filename)
    np.testing.assert_allclose(restored.weight_all(trees), expected, rtol=1e-12, atol=0)
