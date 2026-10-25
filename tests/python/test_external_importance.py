"""Physical row weights and independent proposal densities."""

import math

import numpy as np
import pytest

siren = pytest.importorskip("siren")

from siren import distributions

M_PI = 0.13957039


# --------------------------------------------------------------------------- #
# Distribution-level semantics of the "weight" column                          #
# --------------------------------------------------------------------------- #

def _table(weights, energies=None):
    keys = ["E", "px", "py", "pz", "x", "y", "z", "m", "weight"]
    rows = []
    for i, w in enumerate(weights):
        E = 1.0 + i if energies is None else energies[i]
        p = math.sqrt(max(E * E - M_PI * M_PI, 0.0))
        rows.append([E, 0.0, 0.0, p, 0.0, 0.0, float(i), M_PI, float(w)])
    return keys, rows


def _sample_rows(dist, n, seed=7):
    rand = siren.utilities.SIREN_random(seed)
    rows = []
    for _ in range(n):
        record = siren.dataclasses.PrimaryDistributionRecord(
            siren.dataclasses.Particle.ParticleType.PiPlus)
        dist.Sample(rand, None, None, record)
        ir = siren.dataclasses.InteractionRecord()
        record.finalize(ir)
        rows.append((int(round(
            ir.interaction_parameters["PrimaryExternalDistribution_row"])), ir))
    return rows


def test_weight_column_sets_normalization():
    weights = [1e-3, 1e-1, 1.0, 10.0]
    keys, rows = _table(weights)
    dist = distributions.PrimaryExternalDistribution(keys, rows)
    assert dist.IsNormalizationSet()
    assert dist.normalization == pytest.approx(sum(weights), rel=1e-12)


def test_rows_sample_uniformly_despite_weight_column():
    """The supplied parent list is the intended injection ensemble: the
    weight column must not reshape the row selection."""
    weights = [1e-3, 1e-1, 1.0, 10.0]
    keys, rows = _table(weights)
    dist = distributions.PrimaryExternalDistribution(keys, rows)
    n = 20000
    counts = np.zeros(len(weights))
    for row, ir in _sample_rows(dist, n):
        counts[row] += 1
        assert dist.GenerationProbability(None, None, ir) == 1.0
    frac = counts / n
    expected = np.full(len(weights), 1.0 / len(weights))
    sigma = np.sqrt(expected * (1 - expected) / n)
    assert np.all(np.abs(frac - expected) < 5 * sigma + 1e-4), frac


def test_physical_density_reports_row_weights():
    """The importance weights come back through PhysicalDensity: the
    physical row density N * w_i / sum(w) for the sampled row, while the
    sampling density stays flat, and the pair is declared non-cancelling."""
    weights = [1e-3, 1e-1, 1.0, 10.0]
    keys, rows = _table(weights)
    dist = distributions.PrimaryExternalDistribution(keys, rows)
    assert dist.PhysicalDensityDiffers()
    total = sum(weights)
    for row, ir in _sample_rows(dist, 50):
        assert dist.GenerationProbability(None, None, ir) == 1.0
        assert dist.PhysicalDensity(None, None, ir) == pytest.approx(
            len(weights) * weights[row] / total, rel=1e-12)


def test_explicit_bias_is_independent_of_importance_weights():
    """An explicit sampling bias reshapes the selection only: one shared
    instance reports the biased sampling density through
    GenerationProbability and the physical row density through
    PhysicalDensity, so the weight ratio de-biases exactly."""
    weights = [1e-3, 1e-1, 1.0, 10.0]
    bias = [2.0, 1.0, 0.5, 0.25]
    keys, rows = _table(weights)
    dist = distributions.PrimaryExternalDistribution(keys, rows, bias)
    total_b = sum(bias)
    total_w = sum(weights)
    counts = np.zeros(len(weights))
    n = 20000
    for row, ir in _sample_rows(dist, n):
        counts[row] += 1
        assert dist.GenerationProbability(None, None, ir) == pytest.approx(
            len(rows) * bias[row] / total_b, rel=1e-12)
        assert dist.PhysicalDensity(None, None, ir) == pytest.approx(
            len(rows) * weights[row] / total_w, rel=1e-12)
    frac = counts / n
    expected = np.asarray(bias) / total_b
    sigma = np.sqrt(expected * (1 - expected) / n)
    assert np.all(np.abs(frac - expected) < 5 * sigma + 1e-4), frac
    assert dist.normalization == pytest.approx(total_w, rel=1e-12)


def test_emin_filter_rederives_weights_and_normalization():
    weights = [1e-3, 1e-1, 1.0, 10.0]
    keys, rows = _table(weights)  # energies 1, 2, 3, 4
    dist = distributions.PrimaryExternalDistribution(keys, rows, [], 2.5)
    assert dist.GetPhysicalNumEvents() == 2
    assert dist.normalization == pytest.approx(1.0 + 10.0, rel=1e-12)
    counts = np.zeros(2)
    for row, ir in _sample_rows(dist, 2000):
        counts[row] += 1
        # Physical density over the two survivors, still de-biasing the
        # (uniform) selection.
        assert dist.PhysicalDensity(None, None, ir) == pytest.approx(
            2 * [1.0, 10.0][row] / 11.0, rel=1e-12)
    assert counts[0] / 2000 == pytest.approx(0.5, abs=0.05)


def test_negative_weight_column_fails_loud():
    keys, rows = _table([1.0, -0.5])
    with pytest.raises(RuntimeError, match="finite and non-negative"):
        distributions.PrimaryExternalDistribution(keys, rows)


def test_tables_without_weight_column_are_unchanged():
    keys = ["E", "px", "py", "pz", "x", "y", "z", "m"]
    rows = [[1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, M_PI],
            [2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, M_PI]]
    dist = distributions.PrimaryExternalDistribution(keys, rows)
    assert not dist.IsNormalizationSet()
    assert not dist.PhysicalDensityDiffers()
    counts = np.zeros(2)
    for row, ir in _sample_rows(dist, 2000):
        counts[row] += 1
        assert dist.GenerationProbability(None, None, ir) == 1.0
        assert dist.PhysicalDensity(None, None, ir) == 1.0
    assert abs(counts[0] / 2000 - 0.5) < 0.05
