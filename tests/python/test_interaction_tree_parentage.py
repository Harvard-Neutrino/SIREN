"""Regression tests for parent-interaction reconstruction in event flattening.

The hdf5/parquet flattening in ``siren._util.SaveEvents`` (and its twin in
``siren.SIREN_Controller.SIREN_Controller.SaveEvents``) records, for every
interaction, the index of the parent interaction that produced it. This used to
be reconstructed by comparing a datum's primary four-momentum against earlier
interactions' secondary momenta. That is O(n^2) per tree and, more importantly,
silently mis-links parentage whenever two secondaries share identical
four-momenta (a common, physical situation).

The fix reads parentage directly from the tree's parent/daughter edges via
``InteractionTreeDatum.parent`` -- the authoritative link the injector records
when it builds the tree. Both flattening functions now delegate to the shared
``siren._util.get_parent_indices`` helper. These tests build a small tree by
hand with two secondaries that share identical momenta and prove the new
reconstruction links each daughter to the correct parent where the old momentum
matching does not, end-to-end through both flattening call sites.
"""
import types

import numpy as np
import pytest

siren = pytest.importorskip("siren")


@pytest.fixture(scope="module")
def dc():
    from siren import dataclasses
    return dataclasses


def _make_record(dc, primary_type, primary_momentum, primary_id,
                 secondary_types, secondary_momenta, secondary_ids):
    record = dc.InteractionRecord()
    sig = record.signature
    sig.primary_type = primary_type
    sig.secondary_types = list(secondary_types)
    record.signature = sig
    record.primary_momentum = list(primary_momentum)      # [E, px, py, pz]
    record.primary_id = primary_id
    record.secondary_momenta = [list(p) for p in secondary_momenta]
    record.secondary_ids = list(secondary_ids)
    return record


def _build_degenerate_momentum_tree(dc):
    """Build a tree in which momentum matching is provably ambiguous.

    Layout (insertion order == index order, parents before children)::

        idx0  root      -> secondaries B (mom M), B2 (mom M)   # identical momenta
        idx1  primary B  (parent idx0) -> secondary C (mom M)
        idx2  primary B2 (parent idx0) -> secondary D (mom M)
        idx3  primary C  (parent idx1) -> secondary E (mom M)
        idx4  primary D  (parent idx2) -> secondary F (mom X)

    Almost every particle shares the same momentum ``M``, so matching a primary
    momentum against earlier secondary momenta cannot tell the branches apart.
    The true parentage (from the ``.parent`` edges) is ``[-1, 0, 0, 1, 2]``.

    Note idx4's true parent is index 2 while its depth is 2, so ``parent index``
    and ``depth - 1`` deliberately disagree there -- a reconstruction that merely
    tracked depth (rather than the actual parent pointer) would get idx4 wrong.
    """
    PT = dc.ParticleType

    M = [10.0, 0.0, 0.0, 10.0]   # degenerate momentum shared by B, B2, C, D, E
    X = [5.0, 0.0, 0.0, 5.0]

    # unique IDs for every secondary particle (2-arg ctor sets major/minor)
    sB = dc.ParticleID(1, 1)
    sB2 = dc.ParticleID(1, 2)
    sC = dc.ParticleID(1, 3)
    sD = dc.ParticleID(1, 4)
    sE = dc.ParticleID(1, 5)
    sF = dc.ParticleID(1, 6)

    root = _make_record(
        dc, PT.KPlus, [20.0, 0.0, 0.0, 20.0], dc.ParticleID(9, 9),
        [PT.NuMu, PT.NuMu], [M, M], [sB, sB2])
    # daughter B: its primary IS the root's first secondary (same id + momentum)
    child_B = _make_record(dc, PT.NuMu, M, sB, [PT.NuMu], [M], [sC])
    # daughter B2: its primary IS the root's second secondary (identical momentum)
    child_B2 = _make_record(dc, PT.NuMu, M, sB2, [PT.NuMu], [M], [sD])
    # grand-daughter C: primary IS child_B's secondary (still momentum M)
    child_C = _make_record(dc, PT.NuMu, M, sC, [PT.NuMu], [M], [sE])
    # grand-daughter D: primary IS child_B2's secondary (still momentum M)
    child_D = _make_record(dc, PT.NuMu, M, sD, [PT.NuMu], [X], [sF])

    tree = dc.InteractionTree()
    d0 = tree.add_entry(root, None)     # root: single-arg overload is not bound
    d1 = tree.add_entry(child_B, d0)
    d2 = tree.add_entry(child_B2, d0)
    _d3 = tree.add_entry(child_C, d1)
    _d4 = tree.add_entry(child_D, d2)

    expected_parents = [-1, 0, 0, 1, 2]
    return tree, expected_parents


def _momentum_parent_indices(entries):
    """Faithful reproduction of the OLD momentum-matching reconstruction.

    Mirrors the removed inline snippet from ``_util.SaveEvents`` /
    ``SIREN_Controller.SaveEvents``: for each datum, scan every earlier
    interaction's secondary momenta and append the index on a match, with the
    ``break`` only exiting the inner loop (so a single datum can match -- and be
    appended -- more than once). Kept here so the test documents exactly the
    behavior being regressed away from.
    """
    secondary_momenta = []
    parent_idx = []
    for datum in entries:
        primary_momentum = np.array(datum.record.primary_momentum, dtype=float)
        if datum.depth() == 0:
            parent_idx.append(-1)
        else:
            for _id in range(len(secondary_momenta)):
                for secondary_momentum in secondary_momenta[_id]:
                    if (primary_momentum == secondary_momentum).all():
                        parent_idx.append(_id)
                        break
        secondary_momenta.append(
            [np.array(p, dtype=float) for p in datum.record.secondary_momenta])
    return parent_idx


def test_ground_truth_edges(dc):
    """Sanity: the hand-built tree has the parent/daughter edges we expect."""
    tree, _expected = _build_degenerate_momentum_tree(dc)
    entries = list(tree.tree)
    assert [d.depth() for d in entries] == [0, 1, 1, 2, 2]
    assert entries[0].parent is None
    assert entries[1].parent is entries[0]
    assert entries[2].parent is entries[0]
    assert entries[3].parent is entries[1]
    assert entries[4].parent is entries[2]
    # the two root secondaries genuinely share identical four-momenta
    m0, m1 = entries[0].record.secondary_momenta
    assert list(m0) == list(m1)


def test_get_parent_indices_uses_tree_edges(dc):
    """The shipped helper reconstructs parentage from the tree edges."""
    from siren._util import get_parent_indices
    tree, expected = _build_degenerate_momentum_tree(dc)
    result = get_parent_indices(tree.tree)
    # one index per interaction, each pointing at the true parent
    assert result == expected
    assert len(result) == len(tree.tree)
    # each entry points at the *actual* parent object, not merely a same-depth
    # node (idx4's parent is index 2, not depth-1 == 1)
    for i, datum in enumerate(tree.tree):
        if datum.parent is None:
            assert result[i] == -1
        else:
            assert tree.tree[result[i]] is datum.parent
    # guard specifically against a depth-based stand-in
    depth_minus_one = [d.depth() - 1 for d in tree.tree]
    assert result != depth_minus_one


def test_momentum_matching_would_mislink(dc):
    """The old momentum-based reconstruction is wrong on this tree.

    It cannot produce the correct ``[-1, 0, 0, 1, 2]``: because the secondaries
    all carry momentum M, each deeper primary matches multiple earlier
    interactions' secondaries, so entries are mis-linked and appended more than
    once -- corrupting the per-event alignment.
    """
    from siren._util import get_parent_indices
    tree, expected = _build_degenerate_momentum_tree(dc)

    momentum_result = _momentum_parent_indices(tree.tree)
    correct_result = get_parent_indices(tree.tree)

    assert correct_result == expected
    # the momentum reconstruction disagrees with the truth ...
    assert momentum_result != expected
    # ... and specifically over-counts (more entries than interactions),
    # demonstrating the silent misalignment the fix removes.
    assert len(momentum_result) != len(tree.tree)


def _read_parquet_parent_idx(path):
    import awkward as ak
    return ak.from_parquet(path)["parent_idx"].tolist()


def test_util_saveevents_parquet_parent_idx(dc, tmp_path):
    """End-to-end: siren._util.SaveEvents writes the correct parent_idx column."""
    pytest.importorskip("pyarrow")
    from siren import _util

    tree, expected = _build_degenerate_momentum_tree(dc)
    out = str(tmp_path / "events_util")
    _util.SaveEvents(
        [tree],
        weighter=None,
        gen_times=[0.0],
        save_hdf5=False,
        save_parquet=True,
        save_siren_events=False,
        fid_vol=None,
        output_filename=out)

    assert _read_parquet_parent_idx(out + ".parquet") == [expected]


def test_controller_saveevents_parquet_parent_idx(dc, tmp_path):
    """End-to-end: SIREN_Controller.SaveEvents writes the correct parent_idx.

    SaveEvents only duck-types a handful of ``self`` attributes, so we drive the
    unbound method with a lightweight stand-in rather than standing up a full
    detector/injector -- this exercises the real (modified) controller call site,
    including its delegation to ``_util.get_parent_indices``.
    """
    pytest.importorskip("pyarrow")
    from siren.SIREN_Controller import SIREN_Controller

    tree, expected = _build_degenerate_momentum_tree(dc)
    out = str(tmp_path / "events_ctrl")
    stub_self = types.SimpleNamespace(
        events=[tree],
        gen_times=[0.0],
        global_times=[0.0],
        fid_vol=None,
        injector=types.SimpleNamespace(SaveInjector=lambda path: None),
    )
    SIREN_Controller.SaveEvents(
        stub_self, out, hdf5=False, parquet=True, siren_events=False,
        verbose=False)

    assert _read_parquet_parent_idx(out + ".parquet") == [expected]
