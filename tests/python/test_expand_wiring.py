"""Tests for siren.expand vocabulary and expansion-wiring validation."""

import pytest

siren = pytest.importorskip("siren")

from siren import dataclasses as dc
from siren.expand import child, depth_below, compile_expansion, VertexSpec
from siren._validation import (
    validate_expansion_wiring,
    check_expand_vs_legacy_stopping,
)
from siren.errors import ConfigurationError

PT = dc.ParticleType


def _make_record(primary_type, secondary_types):
    rec = dc.InteractionRecord()
    rec.signature.primary_type = primary_type
    rec.signature.target_type = PT.Decay
    rec.signature.secondary_types = list(secondary_types)
    rec.secondary_masses = [0.0 for _ in secondary_types]
    rec.secondary_momenta = [[0, 0, 0, 0] for _ in secondary_types]
    rec.secondary_helicities = [0 for _ in secondary_types]
    return rec


def _tree_with_root(primary_type, secondary_types):
    """Build a real InteractionTree with one root datum via add_entry."""
    tree = dc.InteractionTree()
    rec = _make_record(primary_type, secondary_types)
    parent = tree.add_entry(rec, None)
    return tree, parent


def test_listed_daughter_expands_unlisted_terminates():
    """A daughter named in expand continues; an unlisted daughter prunes."""
    tree, parent = _tree_with_root(PT.N4, [PT.NuLight, PT.Gamma])
    vertices = [
        VertexSpec(particle=PT.N4, expand=(child("NuLight"),), continue_if=None),
    ]
    stopping_condition = compile_expansion(vertices)

    assert stopping_condition(tree, parent, 0) is False
    assert stopping_condition(tree, parent, 1) is True


def test_child_index_restricts_match():
    """child(index=0) expands only the secondary at position 0."""
    tree, parent = _tree_with_root(PT.NuLight, [PT.NuLight, PT.NuLight])
    vertices = [
        VertexSpec(
            particle=PT.NuLight,
            expand=(child("NuLight", index=0),),
            continue_if=None,
        ),
    ]
    stopping_condition = compile_expansion(vertices)

    assert stopping_condition(tree, parent, 0) is False
    assert stopping_condition(tree, parent, 1) is True


def test_depth_below_expands_while_shallow():
    """depth_below(n) expands while parent depth < n, prunes at depth >= n."""
    tree, shallow_parent = _tree_with_root(PT.N4, [PT.NuLight])
    assert shallow_parent.depth(tree) == 0

    grandchild_rec = _make_record(PT.NuLight, [PT.NuLight])
    deep_parent = tree.add_entry(grandchild_rec, shallow_parent)
    assert deep_parent.depth(tree) == 1

    vertices = [
        VertexSpec(particle=PT.N4, expand=(depth_below(1),), continue_if=None),
        VertexSpec(particle=PT.NuLight, expand=(depth_below(1),), continue_if=None),
    ]
    stopping_condition = compile_expansion(vertices)

    assert stopping_condition(tree, shallow_parent, 0) is False
    assert stopping_condition(tree, deep_parent, 0) is True


def test_continue_if_can_veto_a_matched_child():
    """A matched child rule still prunes when continue_if returns False."""
    tree, parent = _tree_with_root(PT.N4, [PT.NuLight])
    vertices = [
        VertexSpec(
            particle=PT.N4,
            expand=(child("NuLight"),),
            continue_if=lambda tree, parent, i: False,
        ),
    ]
    stopping_condition = compile_expansion(vertices)

    assert stopping_condition(tree, parent, 0) is True


def test_unregistered_parent_type_prunes():
    """A parent whose primary_type has no registered vertex always prunes."""
    tree, parent = _tree_with_root(PT.MuMinus, [PT.NuLight])
    vertices = [
        VertexSpec(particle=PT.N4, expand=(child("NuLight"),), continue_if=None),
    ]
    stopping_condition = compile_expansion(vertices)

    assert stopping_condition(tree, parent, 0) is True


def test_wiring_error_unregistered_child_target():
    """expand names a particle with no registered Vertex -> ConfigurationError."""
    vertices = [
        VertexSpec(particle=PT.N4, expand=(child("NuLight"),), continue_if=None),
    ]
    with pytest.raises(ConfigurationError):
        validate_expansion_wiring(vertices)


def test_wiring_error_unreachable_vertex():
    """A registered Vertex named by no expand list -> ConfigurationError."""
    vertices = [
        VertexSpec(particle=PT.N4, expand=(), continue_if=None),
        VertexSpec(particle=PT.NuLight, expand=(), continue_if=None),
    ]
    with pytest.raises(ConfigurationError):
        validate_expansion_wiring(vertices)


def test_wiring_error_secondaries_without_expansion():
    """Secondaries present but no expand declaration at all -> ConfigurationError."""
    vertices = [
        VertexSpec(
            particle=PT.N4,
            expand=(),
            continue_if=None,
            secondary_types=[PT.NuLight, PT.Gamma],
        ),
    ]
    with pytest.raises(ConfigurationError):
        validate_expansion_wiring(vertices)


def test_legacy_stopping_condition_with_expand_is_hard_error():
    """Legacy stopping_condition plus expand/continue_if -> ConfigurationError."""
    with pytest.raises(ConfigurationError):
        check_expand_vs_legacy_stopping(True, True)

    # Either alone is fine.
    check_expand_vs_legacy_stopping(True, False)
    check_expand_vs_legacy_stopping(False, True)
    check_expand_vs_legacy_stopping(False, False)
