"""Review witnesses for declared chain topology."""
import pickle
import base64

import pytest
import siren

from siren._validation import validate_expansion_wiring
from siren.expand import VertexSpec, child, depth_below
from test_injector_api import _distributions, _load_ccm_detector


def _declared_chain(*, terminal_rules, terminal_predicate=None):
    particle = siren.particles.NuMu
    primary = siren.Vertex(
        particle, siren.interactions.DummyCrossSection(),
        distributions=_distributions(25.0), expand=(depth_below(2),),
    )
    secondary = siren.Vertex(
        particle, siren.interactions.DummyCrossSection(),
        position=siren.distributions.SecondaryPhysicalVertexDistribution(),
        expand=terminal_rules, continue_if=terminal_predicate,
    )
    return siren.injection.Injector(
        events=3, detector=_load_ccm_detector(), seed=7,
        primary=primary, secondaries=[secondary], max_attempts=8000,
    )


@pytest.mark.parametrize("round_trip", ["live", "pickle", "archive"])
def test_repeated_particle_uses_secondary_expansion_policy(round_trip, tmp_path):
    injector = _declared_chain(terminal_rules=(depth_below(0),))
    if round_trip == "pickle":
        injector = pickle.loads(pickle.dumps(injector))
    elif round_trip == "archive":
        path = str(tmp_path / "separate-policies")
        injector.save(path)
        injector = siren.injection.Injector.load(path)
    trees = injector.generate(3, on_shortfall="raise")
    assert [len(tree.tree) for tree in trees] == [2, 2, 2]


def test_repeated_particle_keeps_secondary_continue_if():
    injector = _declared_chain(
        terminal_rules=(depth_below(2),),
        terminal_predicate=lambda *args: False,
    )
    trees = injector.generate(3, on_shortfall="raise")
    assert [len(tree.tree) for tree in trees] == [2, 2, 2]


def test_empty_expansion_marks_reachable_terminal_vertex():
    injector = _declared_chain(terminal_rules=())
    trees = injector.generate(3, on_shortfall="raise")
    assert [len(tree.tree) for tree in trees] == [2, 2, 2]


def test_named_expansion_requires_an_advertised_daughter():
    pt = siren.particles
    vertices = [
        VertexSpec(pt.NuMu, (child(pt.N4),), None, [pt.NuMu, pt.Nucleon]),
        VertexSpec(pt.N4, (depth_below(0),), None, [pt.NuLight, pt.Gamma]),
    ]
    with pytest.raises(siren.errors.ConfigurationError, match="daughter|produce"):
        validate_expansion_wiring(vertices)


def test_secondary_cycle_must_be_reachable_from_primary():
    pt = siren.particles
    vertices = [
        VertexSpec(pt.N4, (depth_below(0),), None, [pt.NuLight, pt.Gamma]),
        VertexSpec(pt.NuLight, (child(pt.Gamma),), None, [pt.Gamma]),
        VertexSpec(pt.Gamma, (child(pt.NuLight),), None, [pt.NuLight]),
    ]
    with pytest.raises(siren.errors.ConfigurationError, match="unreachable"):
        validate_expansion_wiring(vertices)


def _veto(*args):
    return False


def test_picklable_continue_if_survives_facade_pickle():
    injector = _declared_chain(terminal_rules=(depth_below(2),), terminal_predicate=_veto)
    restored = pickle.loads(pickle.dumps(injector))
    assert [len(tree.tree) for tree in restored.generate(3, on_shortfall="raise")] == [2, 2, 2]


@pytest.mark.parametrize("bad", ["NuMu", object()])
def test_invalid_rule_is_rejected(bad):
    vertex = VertexSpec(siren.particles.NuMu, (bad,), None)
    with pytest.raises(siren.errors.ConfigurationError, match="rules"):
        validate_expansion_wiring([vertex])
    with pytest.raises(siren.errors.ConfigurationError, match="rules"):
        siren.expand.compile_expansion([vertex])


@pytest.mark.parametrize("rules,predicate", [
    ((child(siren.particles.NuMu),), None),
    ((depth_below(1),), None),
    ((), _veto),
])
def test_orphaned_primary_expansion_is_rejected(rules, predicate):
    injector = _declared_chain(terminal_rules=())
    with pytest.raises(siren.errors.ConfigurationError, match="registered"):
        siren.injection.Injector(
            events=1, detector=injector.detector_model,
            primary=siren.Vertex(
                siren.particles.NuMu, siren.interactions.DummyCrossSection(),
                distributions=_distributions(25.0), expand=rules, continue_if=predicate),
        ).engine


def test_real_version_zero_expansion_pickle_still_loads():
    encoded = (
        "gASVowAAAAAAAACMD3NpcmVuLmluamVjdGlvbpSMElNlY29uZGFyeUV4cGFuc2lvbpSTlCmB"
        "lIxwMDEwMDAwODAwMDAwMDAwMDAxMDAwMDAwMDAwMDAwMDAwZTAwMDAwMDAwMDAwMDAw"
        "MGUwMDAwMDAwMDAwMDAwMGZmZmZmZmZmZmZmZmZmZmYwMTAwMDAwMDAwMDAwMDAwMDAw"
        "MDAwMDAwMDAwMDAwMJSFlGIu"
    )
    rules = pickle.loads(base64.b64decode(encoded))
    pt = int(siren.particles.NuMu)
    assert rules == siren.injection.SecondaryExpansion([[pt, pt, -1, 1, 0]])


@pytest.mark.parametrize("index", [1, 100, 0.5])
def test_child_index_must_reach_the_declared_secondary(index):
    pt = siren.particles
    vertices = [
        VertexSpec(pt.N4, (child(pt.NuLight, index=index),), None,
                   [pt.NuLight, pt.Gamma], ((pt.NuLight, pt.Gamma),)),
        VertexSpec(pt.NuLight, (), None, [pt.Gamma]),
    ]
    with pytest.raises(siren.errors.ConfigurationError, match="unreachable"):
        validate_expansion_wiring(vertices)
