"""Results snapshot/normalization contract: immutability, merge, variance, bad()."""
import numpy as np
import pytest

import siren
from siren.Results import Results
from siren.errors import ConfigurationError

_NuMu = siren.dataclasses.ParticleType.NuMu


class _FakeInjector:
    """Minimal injector exposing the snapshot attributes Results reads."""

    def __init__(self, attempts, requested, injected, primary_type=_NuMu):
        self.injection_attempts = attempts
        self.number_of_events = requested
        self.injected_events = injected
        self.primary_type = primary_type
        self.engine = None


def _results(weights, injector, weighter=None):
    events = list(range(len(weights)))
    gen_times = [0.0] * len(weights)
    return Results(events, list(weights), gen_times, weighter, injector)


# --------------------------------------------------------------------------- #
# Snapshot immutability under injector reuse                                    #
# --------------------------------------------------------------------------- #

def test_weights_are_immutable():
    inj = _FakeInjector(attempts=10, requested=5, injected=5)
    r = _results([1.0, 2.0, 3.0], inj)
    with pytest.raises(ValueError):
        r.weights[0] = 99.0


def test_snapshot_does_not_drift_when_injector_reused():
    inj = _FakeInjector(attempts=10, requested=5, injected=5)
    r = _results([1.0, 2.0, 3.0], inj)
    assert r.attempts == 10
    assert r.requested == 5
    assert r.injected == 5
    before = np.array(r.weights, copy=True)

    # Reuse the injector for a "later run": mutate its live counters.
    inj.injection_attempts = 999
    inj.injected_events = 42

    # The earlier Results still reflects the run at its creation.
    assert r.attempts == 10
    assert r.injected == 5
    np.testing.assert_array_equal(np.array(r.weights), before)


# --------------------------------------------------------------------------- #
# bad() indexing                                                                #
# --------------------------------------------------------------------------- #

def test_bad_classifies_inf_nan_zero():
    inj = _FakeInjector(attempts=6, requested=6, injected=6)
    w = [1.0, float("inf"), float("nan"), 0.0, 2.0, -float("inf")]
    r = _results(w, inj)
    bad = r.bad()
    assert bad["inf"] == [1, 5]
    assert bad["nan"] == [2]
    assert bad["zero"] == [3]
    # The three index lists are disjoint.
    allbad = bad["inf"] + bad["nan"] + bad["zero"]
    assert len(allbad) == len(set(allbad))


def test_ess_over_finite_weights():
    inj = _FakeInjector(attempts=4, requested=4, injected=4)
    r = _results([1.0, 1.0, 1.0, 1.0], inj)
    # Equal weights -> ESS == N.
    assert r.ess() == pytest.approx(4.0)


def test_ess_excludes_inf_and_nan():
    inf = float("inf")
    nan = float("nan")
    inj = _FakeInjector(attempts=5, requested=5, injected=5)
    r = _results([1.0, 2.0, inf, nan, 3.0], inj)
    # ESS is computed over the finite weights [1, 2, 3] only.
    finite = [1.0, 2.0, 3.0]
    expected = (sum(finite) ** 2) / sum(x * x for x in finite)
    assert r.ess() == pytest.approx(expected)


def test_ess_all_zero_is_zero():
    inj = _FakeInjector(attempts=3, requested=3, injected=3)
    r = _results([0.0, 0.0, 0.0], inj)
    assert r.ess() == 0.0


def test_where_filters_and_inherits_counts():
    inj = _FakeInjector(attempts=20, requested=5, injected=5)
    r = _results([1.0, 0.0, 3.0, 0.0, 5.0], inj)
    kept = r.where(lambda e, w: w > 0.0)
    assert len(kept) == 3
    assert list(kept.weights) == pytest.approx([1.0, 3.0, 5.0])
    # A filtered view still describes the same generation run.
    assert kept.attempts == 20
    assert kept.injected == 5
    # The filtered weights remain immutable.
    with pytest.raises(ValueError):
        kept.weights[0] = 9.0


def test_slice_inherits_run_counts():
    inj = _FakeInjector(attempts=20, requested=5, injected=5)
    r = _results([1.0, 2.0, 3.0, 4.0, 5.0], inj)
    sub = r[1:3]
    assert len(sub) == 2
    assert sub.attempts == 20
    assert sub.injected == 5


# --------------------------------------------------------------------------- #
# merge() rescaling exactness                                                   #
# --------------------------------------------------------------------------- #

def test_merge_rescales_by_injected_fraction():
    """Pooled sum(w) equals the concatenated raw sum divided by n_total scaling.

    SIREN weights already carry a 1/n normalization (Weighter seeds the
    generation probability with InjectedEvents), so the run estimator is
    sum(w). Pooling K runs with w -> w * n_i/n_total makes sum over the pool
    equal the single-n_total-run estimator.
    """
    inj_a = _FakeInjector(attempts=200, requested=100, injected=100)
    inj_b = _FakeInjector(attempts=200, requested=100, injected=100)
    # Same config key (same requested + primary_type) so merge is allowed.
    wa = [2.0, 4.0, 6.0]
    wb = [1.0, 3.0]
    ra = _results(wa, inj_a)
    rb = _results(wb, inj_b)

    merged = Results.merge([ra, rb])

    n_total = 200
    expected = ([w * (100 / n_total) for w in wa]
                + [w * (100 / n_total) for w in wb])
    np.testing.assert_allclose(np.array(merged.weights), expected, rtol=1e-12)
    assert merged.injected == n_total
    assert merged.attempts == 400
    assert merged.requested == 200
    assert len(merged) == 5


def test_merge_estimator_equals_single_double_length_run():
    """Pooled sum(w) equals a single n_total-run's sum(w) estimator.

    SIREN weights carry a 1/n normalization, so a run of injected n has
    weights proportional to 1/n and its rate estimator is sum(w). Modeling the
    same physical events as either two n-runs or one 2n-run: the 2n-run's
    weights are half the n-run's (1/2n vs 1/n). merge()'s w * n_i/n_total
    rescaling must reproduce the 2n-run's per-event weights exactly, so the
    pooled estimator equals the single-double-length-run estimator.
    """
    n = 100
    # Physical per-event rate contributions for 6 events, as an n-run would
    # weight them (proportional to 1/n).
    raw_run_a = [0.020, 0.040, 0.060]
    raw_run_b = [0.010, 0.030, 0.050]
    inj_a = _FakeInjector(attempts=300, requested=n, injected=n)
    inj_b = _FakeInjector(attempts=300, requested=n, injected=n)
    ra = _results(raw_run_a, inj_a)
    rb = _results(raw_run_b, inj_b)
    merged = Results.merge([ra, rb])

    # A single 2n-run over the same 6 events weights each at half the n-run
    # value (1/2n vs 1/n).
    single_2n = [w * 0.5 for w in raw_run_a + raw_run_b]

    assert float(np.sum(merged.weights)) == pytest.approx(
        float(np.sum(single_2n)), rel=1e-12)
    np.testing.assert_allclose(np.array(merged.weights), single_2n, rtol=1e-12)


def test_merge_rejects_differing_config():
    inj_a = _FakeInjector(attempts=1, requested=100, injected=100)
    inj_b = _FakeInjector(attempts=1, requested=250, injected=100)  # different requested
    ra = _results([1.0], inj_a)
    rb = _results([1.0], inj_b)
    with pytest.raises(ConfigurationError):
        Results.merge([ra, rb])


def test_merge_rejects_empty():
    with pytest.raises(ConfigurationError):
        Results.merge([])


# --------------------------------------------------------------------------- #
# variance_report channel shares sum to 1 per vertex                            #
# --------------------------------------------------------------------------- #

class _StubEngine:
    """Engine stub exposing the two hooks variance_report reads.

    GetPhaseSpaces returns pre-seeded mixtures; AccumulateEventToMixtures is a
    no-op so the seeded kp_accumulator/kp_count survive the read.
    """

    def __init__(self, mixtures):
        self._mixtures = mixtures

    def GetPhaseSpaces(self):
        return self._mixtures

    def AccumulateEventToMixtures(self, event, weight, discount, recurse):
        return None


class _SeededMixture:
    """A mixture-like object with fixed KP accumulator state and channels."""

    def __init__(self, kp_accumulator, kp_count, n_channels):
        self.kp_accumulator = list(kp_accumulator)
        self.kp_count = kp_count
        self.channels = [object() for _ in range(n_channels)]

    def ResetAccumulators(self, recurse):
        # Preserve the seeded state so variance_report reads it back.
        return None


class _EngineInjector:
    """Injector-like object whose .engine is a stub engine."""

    def __init__(self, engine):
        self.engine = engine
        self.injection_attempts = 0
        self.number_of_events = 0
        self.injected_events = 0
        self.primary_type = _NuMu


def test_variance_report_shares_sum_to_one_per_vertex():
    # Two mixtures: one populated (W = [3,1] -> shares [0.75, 0.25]), one
    # degenerate (no accumulated events -> uniform shares).
    count = 100
    populated = _SeededMixture([3.0 * count, 1.0 * count], count, 2)
    degenerate = _SeededMixture([], 0, 3)
    engine = _StubEngine([populated, degenerate])
    inj = _EngineInjector(engine)

    r = _results([1.0, 2.0, 3.0], inj)
    report = r.variance_report()

    assert len(report.vertices) == 2
    for v in report.vertices:
        assert sum(v.shares) == pytest.approx(1.0, abs=1e-12)

    # The populated vertex reports W_i / sum(W) shares.
    assert report.vertices[0].degenerate is False
    assert report.vertices[0].shares == pytest.approx([0.75, 0.25], abs=1e-12)
    # The degenerate vertex reports uniform shares over its channel count.
    assert report.vertices[1].degenerate is True
    assert report.vertices[1].shares == pytest.approx([1 / 3, 1 / 3, 1 / 3], abs=1e-12)


def test_variance_report_empty_when_no_mixtures():
    inj = _FakeInjector(attempts=1, requested=1, injected=1)
    r = _results([1.0], inj)
    report = r.variance_report()
    assert report.vertices == []


# --------------------------------------------------------------------------- #
# save(pot=) forwarding                                                         #
# --------------------------------------------------------------------------- #

def _one_record_tree():
    """A minimal single-record InteractionTree for the real save path."""
    dc = siren.dataclasses
    PT = dc.ParticleType
    record = dc.InteractionRecord()
    sig = record.signature
    sig.primary_type = PT.NuMu
    sig.secondary_types = [PT.MuMinus, PT.PPlus]
    record.signature = sig
    record.primary_momentum = [10.0, 0.0, 0.0, 10.0]
    record.primary_id = dc.ParticleID(1, 1)
    record.secondary_momenta = [[5.0, 0.0, 0.0, 5.0], [5.0, 0.0, 0.0, 5.0]]
    record.secondary_ids = [dc.ParticleID(1, 2), dc.ParticleID(1, 3)]
    tree = dc.InteractionTree()
    tree.add_entry(record, None)
    return tree


def test_save_pot_written_to_hdf5(tmp_path):
    """Results.save(pot=...) records POT in the real HDF5 output.

    Exercises the real ``siren._util.SaveEvents`` (no monkeypatch) and reads the
    ``pot`` attribute back off the saved ``Events`` group.
    """
    h5py = pytest.importorskip("h5py")

    inj = _FakeInjector(attempts=1, requested=1, injected=1)
    tree = _one_record_tree()
    r = Results([tree], [1.0], [0.0], None, inj)

    out = str(tmp_path / "events_pot")
    r.save(out, pot=5e20, save_hdf5=True, save_parquet=False,
           save_siren_events=False)

    with h5py.File(out + ".hdf5", "r") as f:
        assert f["Events"].attrs["pot"] == pytest.approx(5e20)


def test_save_without_pot_writes_no_pot_attr(tmp_path):
    """With no pot given, the HDF5 output carries no pot attribute."""
    h5py = pytest.importorskip("h5py")

    inj = _FakeInjector(attempts=1, requested=1, injected=1)
    tree = _one_record_tree()
    r = Results([tree], [1.0], [0.0], None, inj)

    out = str(tmp_path / "events_nopot")
    r.save(out, save_hdf5=True, save_parquet=False, save_siren_events=False)

    with h5py.File(out + ".hdf5", "r") as f:
        assert "pot" not in f["Events"].attrs
