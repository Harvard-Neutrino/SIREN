"""siren.metrics.vertex_in selection predicate."""
import pytest

import siren
from siren import metrics


def _tree_at(x, y, z, primary=None, secondaries=None):
    tree = siren.dataclasses.InteractionTree()
    rec = siren.dataclasses.InteractionRecord()
    if primary is not None:
        rec.signature.primary_type = primary
    if secondaries is not None:
        rec.signature.secondary_types = secondaries
    rec.interaction_vertex = [x, y, z]
    tree.add_entry(rec, None)
    return tree


@pytest.fixture
def box():
    # A 2m cube centered at the origin.
    return siren.geometry.Box(widths=[2.0, 2.0, 2.0], center=[0.0, 0.0, 0.0])


def test_vertex_inside_returns_true(box):
    pred = metrics.vertex_in(box)
    assert pred(_tree_at(0.0, 0.0, 0.0)) is True


def test_vertex_outside_returns_false(box):
    pred = metrics.vertex_in(box)
    assert pred(_tree_at(100.0, 0.0, 0.0)) is False


def test_empty_tree_returns_false(box):
    pred = metrics.vertex_in(box)
    assert pred(siren.dataclasses.InteractionTree()) is False


def test_particle_filter_matches(box):
    pt = siren.particles.NuMu
    pred = metrics.vertex_in(box, particle="NuMu")
    inside_numu = _tree_at(0.0, 0.0, 0.0, primary=pt)
    assert pred(inside_numu) is True


def test_particle_filter_excludes_other_particle(box):
    pred = metrics.vertex_in(box, particle="NuMu")
    inside_nue = _tree_at(0.0, 0.0, 0.0, primary=siren.particles.NuE)
    # The vertex is inside the box but its particle is not NuMu.
    assert pred(inside_nue) is False


def test_particle_filter_matches_secondary(box):
    pred = metrics.vertex_in(box, particle="Gamma")
    tree = _tree_at(0.0, 0.0, 0.0, primary=siren.particles.NuMu,
                    secondaries=[siren.particles.NuE, siren.particles.Gamma])
    assert pred(tree) is True


def test_multi_vertex_only_child_inside(box):
    # Root vertex is far outside; a child vertex is inside the box.
    tree = siren.dataclasses.InteractionTree()
    root = siren.dataclasses.InteractionRecord()
    root.interaction_vertex = [100.0, 0.0, 0.0]
    root_datum = tree.add_entry(root, None)
    child = siren.dataclasses.InteractionRecord()
    child.interaction_vertex = [0.0, 0.0, 0.0]
    tree.add_entry(child, root_datum)

    pred = metrics.vertex_in(box)
    assert pred(tree) is True


def test_predicate_wraps_as_tune_metric(box):
    is_in = metrics.vertex_in(box)
    tree_in = _tree_at(0.0, 0.0, 0.0)
    tree_out = _tree_at(100.0, 0.0, 0.0)
    metric = lambda ev, w: w if is_in(ev) else 0.0
    assert metric(tree_in, 3.0) == 3.0
    assert metric(tree_out, 3.0) == 0.0
