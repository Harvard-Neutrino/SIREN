"""Segment mode of PrimaryExternalDistribution: opting in with
segment_column="length" turns each row into a track segment whose interaction
vertex is sampled along it.

The physics contract: with the process weighted by ExternalBounds(), the sum
of event weights equals sum_i w_i * [1 - exp(-int n sigma dl)] over each
segment, i.e. the physical primaries per exposure on the segment times the
probability that they interact along it. The uniform-along-segment proposal
is unbiased for any material profile along the segment, so a segment that
straddles a material boundary still reproduces the closed form.
"""

import math
import os

import numpy as np
import pytest

import siren
from siren import _util
from siren import channels
from siren import detector
from siren import distributions
from siren.Injector import Injector
from siren.Weighter import Weighter

Gamma = siren.dataclasses.Particle.ParticleType.Gamma
ALP = siren.dataclasses.Particle.ParticleType.ALP
NuMu = siren.dataclasses.Particle.ParticleType.NuMu
Nucleon = siren.dataclasses.Particle.ParticleType.Nucleon

SIGMA = 1.0e-30  # cm^2, constant


# --------------------------------------------------------------------------- #
# Distribution-level semantics                                                  #
# --------------------------------------------------------------------------- #

def _segment_table(rows, extra_keys=(), weights=None):
    """rows: list of (start xyz, direction xyz, E, length)."""
    keys = ["E", "px", "py", "pz", "x0", "y0", "z0", "m", "length"] + list(extra_keys)
    if weights is not None:
        keys.append("weight")
    out = []
    for i, (start, direction, E, L) in enumerate(rows):
        d = np.asarray(direction, dtype=float)
        d = d / np.linalg.norm(d)
        row = [E, *(E * d), *start, 0.0, L] + [0.0] * len(extra_keys)
        if weights is not None:
            row.append(weights[i])
        out.append(row)
    return keys, out


def _sample(dist, n, seed=3):
    rand = siren.utilities.SIREN_random(seed)
    records = []
    for _ in range(n):
        pr = siren.dataclasses.PrimaryDistributionRecord(Gamma)
        dist.Sample(rand, None, None, pr)
        ir = siren.dataclasses.InteractionRecord()
        pr.finalize(ir)
        records.append(ir)
    return records


def test_length_column_declares_vertex_and_longitudinal_density():
    keys, rows = _segment_table([((0, 0, 0), (0, 0, 1), 0.05, 0.2)])
    dist = distributions.PrimaryExternalDistribution(keys, rows, segment_column="length")
    DV = distributions.DistributionVariable
    assert DV.InteractionVertex in dist.SetVariables()
    assert DV.InitialPosition in dist.SetVariables()
    assert set(dist.DensityVariables()) == {"External", "PrimaryPositionLongitudinal"}
    assert dist.PhysicalDensityVariables() == ["External"]
    assert dist.ProvidesExternalBounds()
    point = distributions.PrimaryExternalDistribution(keys[:-1], [r[:-1] for r in rows])
    assert DV.InteractionVertex not in point.SetVariables()
    assert point.DensityVariables() == ["External"]
    assert not point.ProvidesExternalBounds()
    # Without the opt-in a column named "length" is ordinary metadata.
    metadata = distributions.PrimaryExternalDistribution(keys, rows)
    assert not metadata.ProvidesExternalBounds()
    assert DV.InteractionVertex not in metadata.SetVariables()
    assert metadata.DensityVariables() == ["External"]
    assert metadata.GetSegmentColumn() == ""
    metadata.SetSegmentColumn("length")
    assert metadata.ProvidesExternalBounds()
    assert metadata == dist
    metadata.SetSegmentColumn("")
    assert not metadata.ProvidesExternalBounds()


def test_vertex_is_uniform_along_segment_and_bounds_are_end_points():
    start, direction, L = (1.0, -2.0, 0.5), (1.0, 2.0, 2.0), 0.3
    keys, rows = _segment_table([(start, direction, 0.05, L)])
    dist = distributions.PrimaryExternalDistribution(keys, rows, segment_column="length")
    d = np.asarray(direction) / np.linalg.norm(direction)
    records = _sample(dist, 4000)
    fractions = []
    for ir in records:
        v = np.asarray(ir.interaction_vertex) - np.asarray(start)
        s = float(v @ d)
        assert np.allclose(v - s * d, 0.0, atol=1e-12)
        assert 0.0 <= s <= L
        assert np.allclose(ir.primary_initial_position, start)
        assert ir.interaction_parameters["length"] == pytest.approx(L)
        fractions.append(s / L)
        lo, hi = dist.InjectionBounds(None, None, ir)
        assert np.allclose([lo.GetX(), lo.GetY(), lo.GetZ()], start)
        assert np.allclose([hi.GetX(), hi.GetY(), hi.GetZ()], np.asarray(start) + L * d)
    counts, _ = np.histogram(fractions, bins=10, range=(0.0, 1.0))
    expected = len(records) / 10
    assert np.all(np.abs(counts - expected) < 5 * math.sqrt(expected))
    # Generation density: 1/L per metre for one row.
    assert dist.GenerationProbability(None, None, records[0]) == pytest.approx(1.0 / L)


def test_generation_probability_uses_sampling_weights_and_length():
    keys, rows = _segment_table([((0, 0, 0), (0, 0, 1), 0.05, 0.1),
                                 ((0, 0, 1), (0, 0, 1), 0.05, 0.4)])
    dist = distributions.PrimaryExternalDistribution(keys, rows, [3.0, 1.0], segment_column="length")
    records = _sample(dist, 2000)
    seen = {}
    for ir in records:
        i = int(round(ir.interaction_parameters["PrimaryExternalDistribution_row"]))
        seen[i] = dist.GenerationProbability(None, None, ir)
    assert seen[0] == pytest.approx(2 * 3.0 / 4.0 / 0.1)
    assert seen[1] == pytest.approx(2 * 1.0 / 4.0 / 0.4)
    # Without a weight column the physical row density follows the sampled
    # ensemble (legacy semantics) but never includes the 1/length factor.
    for ir in records[:50]:
        i = int(round(ir.interaction_parameters["PrimaryExternalDistribution_row"]))
        assert dist.PhysicalDensity(None, None, ir) == pytest.approx(2 * [3.0, 1.0][i] / 4.0)
    assert dist.PhysicalDensityDiffers()


@pytest.mark.parametrize("bad_keys,bad_rows,message", [
    (["E", "px", "py", "pz", "m", "length"], [[0.05, 0, 0, 0.05, 0, 0.1]], "requires x0, y0, z0"),
    (["E", "x0", "y0", "z0", "m", "length"], [[0.05, 0, 0, 0, 0, 0.1]], "requires px, py, pz"),
    (["E", "px", "py", "pz", "x0", "y0", "z0", "x", "y", "z", "m", "length"],
     [[0.05, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0.1]], "incompatible with x, y, z"),
    (["E", "px", "py", "pz", "x0", "y0", "z0", "m", "length"],
     [[0.05, 0, 0, 0.05, 0, 0, 0, 0, 0.0]], "finite and positive"),
    (["E", "px", "py", "pz", "x0", "y0", "z0", "m", "length"],
     [[0.05, 0, 0, 0.05, 0, 0, 0, 0, -0.1]], "finite and positive"),
    (["E", "px", "py", "pz", "x0", "y0", "z0", "m", "length", "length"],
     [[0.05, 0, 0, 0.05, 0, 0, 0, 0, 0.1, 0.1]], "Duplicate"),
    (["E", "px", "py", "pz", "x0", "y0", "z0", "m"],
     [[0.05, 0, 0, 0.05, 0, 0, 0, 0]], "not found"),
])
def test_segment_table_validation(bad_keys, bad_rows, message):
    with pytest.raises(RuntimeError, match=message):
        distributions.PrimaryExternalDistribution(bad_keys, bad_rows, segment_column="length")


def test_reserved_names_cannot_be_the_segment_column():
    keys, rows = _segment_table([((0, 0, 0), (0, 0, 1), 0.05, 0.2)])
    for name in ("E", "px", "x0"):
        with pytest.raises(RuntimeError, match="cannot be the segment length column"):
            distributions.PrimaryExternalDistribution(keys, rows, segment_column=name)


def test_generation_density_is_zero_off_own_support():
    # Two proposals over the same source point/direction/energy with lengths
    # 6 cm and 12 cm: the shorter table must report zero density for the
    # longer table's vertices beyond 6 cm, so pooled weighting stays exact.
    short_keys, short_rows = _segment_table([((0, 0, 0), (0, 0, 1), 0.05, 0.06)])
    long_keys, long_rows = _segment_table([((0, 0, 0), (0, 0, 1), 0.05, 0.12)])
    short = distributions.PrimaryExternalDistribution(short_keys, short_rows, segment_column="length")
    long = distributions.PrimaryExternalDistribution(long_keys, long_rows, segment_column="length")
    inside = outside = 0
    for ir in _sample(long, 400):
        z = ir.interaction_vertex[2]
        assert long.GenerationProbability(None, None, ir) == pytest.approx(1 / 0.12)
        if z <= 0.06:
            assert short.GenerationProbability(None, None, ir) == pytest.approx(1 / 0.06)
            inside += 1
        else:
            assert short.GenerationProbability(None, None, ir) == 0.0
            outside += 1
    assert inside > 100 and outside > 100
    # Transverse or upstream vertices, or a different start, are off support.
    ir = _sample(long, 1)[0]
    for vertex in ((0.001, 0.0, 0.03), (0.0, 0.0, -0.01), (0.0, 0.0, 0.13)):
        ir.interaction_vertex = vertex
        assert long.GenerationProbability(None, None, ir) == 0.0
    ir.interaction_vertex = (0.0, 0.0, 0.03)
    ir.primary_initial_position = (0.0, 0.0, 0.001)
    assert long.GenerationProbability(None, None, ir) == 0.0


def test_segment_table_round_trips_through_csv(tmp_path):
    keys, rows = _segment_table([((0.1, 0.2, 0.3), (0, 1, 0), 0.07, 0.25)], weights=[2.5])
    path = tmp_path / "segments.csv"
    path.write_text(",".join(keys) + "\n" + "\n".join(
        ",".join(repr(float(v)) for v in row) for row in rows) + "\n")
    plain = distributions.PrimaryExternalDistribution(str(path))
    assert not plain.ProvidesExternalBounds(), "a CSV length column is metadata unless opted in"
    dist = distributions.PrimaryExternalDistribution(str(path), segment_column="length")
    assert dist.ProvidesExternalBounds()
    assert distributions.DistributionVariable.InteractionVertex in dist.SetVariables()
    ir = _sample(dist, 1)[0]
    assert ir.interaction_parameters["length"] == pytest.approx(0.25)
    assert dist.GenerationProbability(None, None, ir) == pytest.approx(4.0)


# --------------------------------------------------------------------------- #
# Injector / Weighter witnesses on the CCM detector model                       #
# --------------------------------------------------------------------------- #

def _ccm():
    try:
        det_dir = _util.get_detector_model_path("CCM")
        dm = detector.DetectorModel()
        dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
        dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    except Exception as e:  # pragma: no cover - environment dependent
        pytest.skip(f"CCM detector model unavailable: {e}")
    return dm


@pytest.fixture(scope="module")
def ccm():
    return _ccm()


def _flat_xs(primary, target):
    return siren.interactions.TrivialCrossSection(SIGMA, [primary], [target])


def _density(dm, position, target):
    return dm.GetParticleDensity(detector.DetectorPosition(siren.math.Vector3D(*position)), target)


def _closed_form(dm, table_rows, weights, target):
    """sum_i w_i (1 - exp(-int n sigma dl)) with the detector model's densities."""
    total = 0.0
    for (start, direction, _E, L), w in zip(table_rows, weights):
        d = np.asarray(direction, float)
        d /= np.linalg.norm(d)
        p0 = detector.DetectorPosition(siren.math.Vector3D(*start))
        p1 = detector.DetectorPosition(siren.math.Vector3D(*(np.asarray(start) + L * d)))
        depth = dm.GetInteractionDepthInCGS(p0, p1, [target], [SIGMA], math.inf)
        total += w * (1.0 - math.exp(-depth))
    return total


def _run(dm, table_rows, weights, weighting, events, seed=11, sampling=None):
    keys, rows = _segment_table(table_rows, weights=weights)
    external = distributions.PrimaryExternalDistribution(
        keys, rows, [] if sampling is None else list(sampling), segment_column="length")
    primary = siren.Vertex(
        Gamma,
        _flat_xs(Gamma, Nucleon),
        distributions=[external],
        physical=[external],
        weighting=weighting,
    )
    injector = Injector(detector=dm, primary=primary, events=events, seed=seed)
    weighter = Weighter(injector, primary_physical=primary.physical)
    results = siren.generate(injector, weighter, events=events, on_shortfall="error")
    weights_out = np.array([w for _, w in results], dtype=float)
    vertices = np.array([ev.tree[0].record.interaction_vertex for ev, _ in results])
    return weights_out, vertices


# The CCM-v2 detector frame is centred on the detector: the Lujan target sits
# at x = -23 m, and geometry z maps to z + 0.65. The lower tungsten target is
# the cylinder r = 0.05 m, detector z in [0.26, 0.558] m; segments below stay
# well inside its homogeneous region.
TX = -23.0
TUNGSTEN_ROWS = [((TX, 0.0, 0.30), (0.0, 0.0, 1.0), 0.05, 0.12)]


def test_native_weight_reproduces_closed_form_homogeneous(ccm):
    weights = [7.0]
    n_w = _density(ccm, (TX, 0.0, 0.35), Nucleon)
    assert n_w > 0
    expected = _closed_form(ccm, TUNGSTEN_ROWS, weights, Nucleon)
    # Depth is tiny, so the closed form is w n sigma L to 1e-4 relative.
    n_sigma = n_w * SIGMA * 100.0  # per metre
    assert expected == pytest.approx(weights[0] * n_sigma * 0.12, rel=1e-4)
    events = 500
    out, vertices = _run(ccm, TUNGSTEN_ROWS, weights, siren.ExternalBounds(), events=events)
    z = vertices[:, 2]
    assert z.min() >= 0.30 - 1e-12 and z.max() <= 0.42 + 1e-12
    # Event weight = (w / N) * n sigma L * exp(-n sigma (z - z0)): the physical
    # rate on the segment times the survival to the sampled vertex, exactly.
    predicted = weights[0] / events * n_sigma * 0.12 * np.exp(-n_sigma * (z - 0.30))
    assert np.allclose(out, predicted, rtol=1e-9, atol=0.0)
    # Summed over events this estimates w (1 - exp(-D)) with relative spread
    # ~ D / sqrt(12 N) ~ 2e-6.
    assert out.sum() == pytest.approx(expected, rel=2e-5)
    counts, _ = np.histogram(z, bins=8, range=(0.30, 0.42))
    assert np.all(np.abs(counts - len(z) / 8) < 5 * math.sqrt(len(z) / 8))


def test_sampling_bias_leaves_totals_unchanged(ccm):
    rows = [((TX, 0.0, 0.29), (0.0, 0.0, 1.0), 0.05, 0.02),
            ((TX, 0.0, 0.35), (0.0, 0.0, 1.0), 0.07, 0.08)]
    weights = [3.0, 1.0]
    sampling = [1.0, 9.0]
    expected = _closed_form(ccm, rows, weights, Nucleon)
    n_sigma = _density(ccm, (TX, 0.0, 0.35), Nucleon) * SIGMA * 100.0
    events = 2000
    flat, _ = _run(ccm, rows, weights, siren.ExternalBounds(), events=events)
    biased, vertices = _run(ccm, rows, weights, siren.ExternalBounds(), events=events,
                            sampling=sampling, seed=5)
    row = (vertices[:, 2] >= 0.35).astype(int)
    # Biased rows are drawn ~9:1 and de-biased per event by (sum s / (N s_i)).
    assert 0.85 < row.mean() < 0.95
    z0 = np.where(row == 1, 0.35, 0.29)
    L = np.where(row == 1, 0.08, 0.02)
    w = np.asarray(weights)[row]
    s = np.asarray(sampling)[row]
    predicted = w * sum(sampling) / (events * s) * n_sigma * L * np.exp(-n_sigma * (vertices[:, 2] - z0))
    assert np.allclose(biased, predicted, rtol=1e-9, atol=0.0)
    for out in (flat, biased):
        # Binomial row-selection noise dominates: sigma(sum) ~ |a_0 - a_1| sqrt(N p q)
        err = out.std(ddof=1) * math.sqrt(len(out))
        assert abs(out.sum() - expected) < 4 * err
        assert abs(out.sum() - expected) / expected < 0.05


def test_segment_across_material_boundary_matches_closed_form(ccm):
    # Start inside the tungsten cylinder and run out through its side wall
    # into the beryllium reflector.
    rows = [((TX, 0.0, 0.41), (1.0, 0.0, 0.0), 0.05, 0.10)]
    weights = [1.0]
    n_in = _density(ccm, (TX, 0.0, 0.41), Nucleon)
    n_out = _density(ccm, (TX + 0.09, 0.0, 0.41), Nucleon)
    assert n_in != n_out, "witness needs a material change along the segment"
    expected = _closed_form(ccm, rows, weights, Nucleon)
    events = 4000
    out, vertices = _run(ccm, rows, weights, siren.ExternalBounds(), events=events, seed=17)
    x = vertices[:, 0] - TX
    assert x.max() <= 0.10 + 1e-12
    inside = x < 0.05
    assert inside.any() and (~inside).any()
    # Weights follow the local density n(x) sigma (survival factors are < 1e-4).
    assert np.allclose(out[inside] / out[~inside].mean(), n_in / n_out, rtol=2e-4)
    err = out.std(ddof=1) * math.sqrt(len(out))
    assert abs(out.sum() - expected) < 4 * err
    assert abs(out.sum() - expected) / expected < 0.05


def test_fixed_weighting_rejected_for_segment_tables(ccm):
    with pytest.raises(Exception, match="ExternalBounds"):
        _run(ccm, TUNGSTEN_ROWS, [1.0], siren.Fixed(), events=5)
    with pytest.raises(Exception, match="ExternalBounds"):
        _run(ccm, TUNGSTEN_ROWS, [1.0], siren.Propagated(), events=5)


def test_point_table_with_length_metadata_keeps_point_semantics(ccm):
    keys = ["E", "px", "py", "pz", "x", "y", "z", "m", "length", "weight"]
    rows = [[0.05, 0.0, 0.0, 0.05, TX, 0.0, 0.35, 0.0, 0.12, 1.0]]
    external = distributions.PrimaryExternalDistribution(keys, rows)
    assert not external.ProvidesExternalBounds()
    primary = siren.Vertex(Gamma, _flat_xs(Gamma, Nucleon), distributions=[external],
                           physical=[external], weighting=siren.Fixed())
    injector = Injector(detector=ccm, primary=primary, events=5, seed=1)
    weighter = Weighter(injector, primary_physical=primary.physical)
    results = siren.generate(injector, weighter, events=5, on_shortfall="error")
    assert len(results) == 5
    assert all(ev.tree[0].record.interaction_parameters["length"] == 0.12 for ev, _ in results)


def test_zero_target_segment_draw_is_a_permissible_miss(ccm):
    # Segment leaving tungsten into beryllium with a W183-only cross section:
    # draws in beryllium have no configured target and must count as ordinary
    # misses (retried under strict generation), not abort it.
    W183 = siren.dataclasses.Particle.ParticleType.W183Nucleus
    rows = [((TX, 0.0, 0.41), (1.0, 0.0, 0.0), 0.05, 0.10)]
    keys, table = _segment_table(rows, weights=[1.0])
    external = distributions.PrimaryExternalDistribution(keys, table, segment_column="length")
    primary = siren.Vertex(Gamma, _flat_xs(Gamma, W183), distributions=[external],
                           physical=[external], weighting=siren.ExternalBounds())
    injector = Injector(detector=ccm, primary=primary, events=200, seed=17)
    trees = injector.generate(events=200, on_failure="raise", on_shortfall="error")
    assert len(trees) == 200
    x = np.array([t.tree[0].record.interaction_vertex[0] - TX for t in trees])
    assert x.max() <= 0.05 + 1e-12
    report = injector.report()
    assert report.failures > 0
    assert all(bucket.reason_name == "NoTargetsOnPath" for bucket in report.by_vertex)


def test_pooled_proposals_with_different_lengths_are_unbiased(ccm):
    # Two injectors over the same source point, 6 cm and 12 cm long, pooled in
    # one Weighter with the 12 cm table as the physical ensemble. The shorter
    # proposal has zero density on the tail of the longer one; the mixture
    # weight must still reproduce the 12 cm closed form (review finding P2-1).
    from siren.Weighter import Weighter as _W
    short_rows = [((TX, 0.0, 0.30), (0.0, 0.0, 1.0), 0.05, 0.06)]
    long_rows = [((TX, 0.0, 0.30), (0.0, 0.0, 1.0), 0.05, 0.12)]
    dists = []
    injectors = []
    for rows, seed in ((short_rows, 11), (long_rows, 23)):
        keys, table = _segment_table(rows, weights=[1.0])
        external = distributions.PrimaryExternalDistribution(keys, table, segment_column="length")
        primary = siren.Vertex(Gamma, _flat_xs(Gamma, Nucleon), distributions=[external],
                               physical=[external], weighting=siren.ExternalBounds())
        injector = Injector(detector=ccm, primary=primary, events=1000, seed=seed)
        trees = injector.generate(events=1000, on_shortfall="error")
        dists.append(external)
        injectors.append((injector, trees))
    weighter = _W(injectors=[inj for inj, _ in injectors], detector_model=ccm, primary_type=Gamma,
                  primary_interactions=[_flat_xs(Gamma, Nucleon)],
                  primary_physical_distributions=[dists[1]])
    events = [t for _, trees in injectors for t in trees]
    weights = np.array([weighter(t) for t in events])
    expected = _closed_form(ccm, long_rows, [1.0], Nucleon)
    # The only noise is the binomial split of the 12 cm draws across the two
    # mixture regions: sigma(sum)/sum = sqrt(250) * 0.08 / (1000 * 0.12) ~ 1%.
    assert abs(weights.sum() - expected) / expected < 4 * 0.0105
    z = np.array([t.tree[0].record.interaction_vertex[2] for t in events]) - 0.30
    n_sigma = _density(ccm, (TX, 0.0, 0.35), Nucleon) * SIGMA * 100.0
    # Mixture density: (1000/0.06 + 1000/0.12) on [0, 6 cm], 1000/0.12 beyond.
    proposal = np.where(z <= 0.06, 1000 / 0.06 + 1000 / 0.12, 1000 / 0.12)
    predicted = n_sigma * np.exp(-n_sigma * z) / proposal
    assert np.allclose(weights, predicted, rtol=1e-9, atol=0.0)


def test_external_bounds_rejected_for_point_tables(ccm):
    keys = ["E", "px", "py", "pz", "x", "y", "z", "m", "weight"]
    rows = [[0.05, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 1.0]]
    external = distributions.PrimaryExternalDistribution(keys, rows)
    primary = siren.Vertex(Gamma, _flat_xs(Gamma, Nucleon), distributions=[external],
                           physical=[external], weighting=siren.ExternalBounds())
    injector = Injector(detector=ccm, primary=primary, events=5, seed=1)
    with pytest.raises(Exception, match="segment table"):
        injector.generate(events=5, on_shortfall="error")


def test_saved_events_reweight_without_the_table(ccm, tmp_path):
    weights = [2.0]
    expected = _closed_form(ccm, TUNGSTEN_ROWS, weights, Nucleon)
    keys, rows = _segment_table(TUNGSTEN_ROWS, weights=weights)
    external = distributions.PrimaryExternalDistribution(keys, rows, segment_column="length")
    primary = siren.Vertex(Gamma, _flat_xs(Gamma, Nucleon), distributions=[external],
                           physical=[external], weighting=siren.ExternalBounds())
    injector = Injector(detector=ccm, primary=primary, events=20, seed=23)
    weighter = Weighter(injector, primary_physical=primary.physical)
    results = siren.generate(injector, weighter, events=20, on_shortfall="error")
    path = str(tmp_path / "segments.siren_events")
    siren.dataclasses.SaveInteractionTrees([ev for ev, _ in results], path)
    loaded = siren.dataclasses.LoadInteractionTrees(path)
    reweighted = np.array([weighter(ev) for ev in loaded], dtype=float)
    assert np.allclose(reweighted, [w for _, w in results], rtol=1e-12, atol=0.0)
    assert reweighted.sum() == pytest.approx(expected, rel=1e-4)


# --------------------------------------------------------------------------- #
# Row layout of pooled tables                                                   #
# --------------------------------------------------------------------------- #

TWO_ROWS = [((TX, 0.0, 0.30), (0.0, 0.0, 1.0), 0.05, 0.12),
            ((TX, 0.01, 0.30), (0.0, 0.0, 1.0), 0.05, 0.12)]


def _segment_injector(ccm, rows, seed, sampling=None, events=20):
    keys, table = _segment_table(rows, weights=[1.0] * len(rows))
    external = distributions.PrimaryExternalDistribution(
        keys, table, [] if sampling is None else list(sampling), segment_column="length")
    primary = siren.Vertex(Gamma, _flat_xs(Gamma, Nucleon), distributions=[external],
                           physical=[external], weighting=siren.ExternalBounds())
    injector = Injector(detector=ccm, primary=primary, events=events, seed=seed)
    return external, injector


def _pool(ccm, injectors, physical):
    from siren.Weighter import Weighter as _W
    weighter = _W(injectors=injectors, detector_model=ccm, primary_type=Gamma,
                  primary_interactions=[_flat_xs(Gamma, Nucleon)],
                  primary_physical_distributions=[physical])
    # The wrapper builds the native weighter lazily; the layout check runs
    # when it is configured, before any event is weighted.
    weighter.engine
    return weighter


def test_row_layout_mismatch_reports_first_difference():
    keys, rows = _segment_table(TWO_ROWS, weights=[1.0, 2.0])
    reference = distributions.PrimaryExternalDistribution(keys, rows, segment_column="length")
    reversed_rows = distributions.PrimaryExternalDistribution(keys, rows[::-1], segment_column="length")
    assert reference.RowLayoutMismatch(reference) == ""
    assert "row 0" in reference.RowLayoutMismatch(reversed_rows)
    subset = distributions.PrimaryExternalDistribution(keys, rows[:1], segment_column="length")
    assert "row counts differ" in reference.RowLayoutMismatch(subset)
    # Lengths, physical weights, sampling weights and the column name may differ.
    renamed_keys = [k if k != "length" else "step_length" for k in keys]
    other_rows = [list(r) for r in rows]
    other_rows[0][keys.index("length")] = 0.06
    other_rows[1][keys.index("weight")] = 9.0
    other = distributions.PrimaryExternalDistribution(
        renamed_keys, other_rows, [1.0, 3.0], segment_column="step_length")
    assert reference.RowLayoutMismatch(other) == ""
    assert other.RowLayoutMismatch(reference) == ""


def test_reordered_tables_are_rejected_when_pooled(ccm):
    # Records cache a row index that each table reads as an index into
    # itself, so the same rows in a different order would be compared with
    # the wrong row and count only the generating proposal (a factor-of-two
    # bias in the rates). The pool is rejected at configuration instead.
    first, inj_a = _segment_injector(ccm, TWO_ROWS, 17)
    second, inj_b = _segment_injector(ccm, TWO_ROWS[::-1], 17)
    with pytest.raises(siren.errors.ConfigurationError, match="row layout"):
        _pool(ccm, [inj_a, inj_b], second)
    # The physical-side table is checked against the injection side too.
    with pytest.raises(siren.errors.ConfigurationError, match="row layout"):
        _pool(ccm, [inj_a], second)


@pytest.mark.parametrize("sampling", [None, "uniform"])
def test_subset_tables_are_rejected_when_pooled(ccm, sampling):
    # A one-row subset and the full two-row table report row densities
    # relative to different uniform-over-rows measures (an asymptotic 25%
    # excess if pooled); with explicit sampling weights the foreign row index
    # is out of range as well. Both are rejected at configuration.
    first, inj_a = _segment_injector(ccm, TWO_ROWS[:1], 17,
                                     sampling=None if sampling is None else [1.0])
    second, inj_b = _segment_injector(ccm, TWO_ROWS, 17,
                                      sampling=None if sampling is None else [1.0, 1.0])
    with pytest.raises(siren.errors.ConfigurationError, match="row counts differ"):
        _pool(ccm, [inj_a, inj_b], second)
    # Directly asked, the subset reports zero density (off its support) for
    # the full table's second row rather than a row-index error.
    trees = inj_b.generate(events=40, on_shortfall="error")
    foreign = [t.tree[0].record for t in trees
               if t.tree[0].record.interaction_parameters["PrimaryExternalDistribution_row"] == 1.0]
    assert foreign
    assert first.GenerationProbability(ccm, None, foreign[0]) == 0.0


def test_compatible_tables_pool_and_reproduce_the_closed_form(ccm):
    # Same primaries at the same indices with a biased sampler on one side
    # and different physical weights: the pool is accepted and unbiased.
    keys, table = _segment_table(TWO_ROWS, weights=[1.0, 1.0])
    reference = distributions.PrimaryExternalDistribution(keys, table, segment_column="length")
    injectors = []
    events = []
    for seed, sampling in ((17, None), (23, [1.0, 3.0])):
        _external, injector = _segment_injector(ccm, TWO_ROWS, seed, sampling=sampling, events=1000)
        events.extend(injector.generate(events=1000, on_shortfall="error"))
        injectors.append(injector)
    weighter = _pool(ccm, injectors, reference)
    weights = np.array([weighter(t) for t in events])
    expected = _closed_form(ccm, TWO_ROWS, [1.0, 1.0], Nucleon)
    assert abs(weights.sum() - expected) / expected < 0.05


# --------------------------------------------------------------------------- #
# Short segments                                                                #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("length", [2e-5, 9e-6, 1e-6])
def test_short_scattering_segments_weight_to_the_closed_form(ccm, length):
    # Paths at or below the detector model's 10 um direction threshold used
    # to integrate only the decay hazard, so a scattering-only segment had zero
    # depth and its normalized position density divided by zero. The depth is
    # now the local interaction density times the length there.
    rows = [((TX, 0.0, 0.30), (0.0, 0.0, 1.0), 0.05, length)]
    external, injector = _segment_injector(ccm, rows, 17, events=50)
    trees = injector.generate(events=50, on_shortfall="error")
    weighter = _pool(ccm, [injector], external)
    weights = np.array([weighter(t) for t in trees])
    n_sigma = _density(ccm, (TX, 0.0, 0.30), Nucleon) * SIGMA * 100.0
    record = trees[0].tree[0].record
    lo, hi = external.InjectionBounds(ccm, None, record)
    depth = ccm.GetInteractionDepthInCGS(detector.DetectorPosition(lo), detector.DetectorPosition(hi),
                                         [Nucleon], [SIGMA], math.inf)
    assert depth == pytest.approx(n_sigma * length, rel=1e-9)
    # Every event: (w / N) n sigma L exp(-n sigma s) with s < 10 um.
    assert np.allclose(weights, n_sigma * length / 50, rtol=1e-7, atol=0.0)
    assert weights.sum() == pytest.approx(_closed_form(ccm, rows, [1.0], Nucleon), rel=1e-6)


def test_support_tolerance_does_not_extend_short_segments():
    # A 0.1 nm proposal must not report support along a 1 nm proposal over
    # the same start and direction: the tolerance is representational (ulps
    # of the coordinates), not a fixed 1e-9 m floor.
    short_keys, short_rows = _segment_table([((0, 0, 0), (0, 0, 1), 0.05, 1e-10)])
    long_keys, long_rows = _segment_table([((0, 0, 0), (0, 0, 1), 0.05, 1e-9)])
    short = distributions.PrimaryExternalDistribution(short_keys, short_rows, segment_column="length")
    long = distributions.PrimaryExternalDistribution(long_keys, long_rows, segment_column="length")
    inside = outside = 0
    for ir in _sample(long, 1000):
        z = ir.interaction_vertex[2]
        assert long.GenerationProbability(None, None, ir) == pytest.approx(1e9)
        if z <= 1e-10:
            assert short.GenerationProbability(None, None, ir) == pytest.approx(1e10)
            inside += 1
        else:
            assert short.GenerationProbability(None, None, ir) == 0.0
            outside += 1
    assert inside > 50 and outside > 800
    # Far from the origin the tolerance scales with the coordinates and the
    # table's own records stay on support.
    keys, rows = _segment_table([((TX, 0.0, 0.30), (1.0, 2.0, 2.0), 0.05, 1e-7)])
    far = distributions.PrimaryExternalDistribution(keys, rows, segment_column="length")
    for ir in _sample(far, 200):
        assert far.GenerationProbability(None, None, ir) == pytest.approx(1e7)


def test_short_segment_across_material_boundary(ccm):
    # An 8 um path centred on the tungsten/beryllium boundary: both halves
    # are resolved by the sector integration (the path direction comes from
    # the trajectory, not from normalizing the tiny difference), so the depth
    # is the piecewise sum and does not depend on the direction of traversal.
    W183 = siren.dataclasses.Particle.ParticleType.W183Nucleus
    points = [(TX + 0.05 - 4e-6, 0.0, 0.41), (TX + 0.05 + 4e-6, 0.0, 0.41)]
    for target in (Nucleon, W183):
        local = [_density(ccm, p, target) * SIGMA * 100.0 for p in points]
        assert local[0] != local[1]
        expected = 0.5 * 8e-6 * sum(local)
        a, b = (detector.DetectorPosition(siren.math.Vector3D(*p)) for p in points)
        forward = ccm.GetInteractionDepthInCGS(a, b, [target], [SIGMA], math.inf)
        reverse = ccm.GetInteractionDepthInCGS(b, a, [target], [SIGMA], math.inf)
        assert forward == pytest.approx(expected, rel=1e-6)
        assert reverse == pytest.approx(expected, rel=1e-6)
    # A W183-only segment entering tungsten from beryllium weights to the
    # tungsten half's closed form.
    keys, table = _segment_table([(points[1], (-1.0, 0.0, 0.0), 0.05, 8e-6)], weights=[1.0])
    external = distributions.PrimaryExternalDistribution(keys, table, segment_column="length")
    vertex = siren.Vertex(Gamma, _flat_xs(Gamma, W183), distributions=[external],
                          physical=[external], weighting=siren.ExternalBounds())
    injector = Injector(detector=ccm, primary=vertex, events=40, seed=17)
    trees = injector.generate(events=40, on_failure="raise", on_shortfall="error")
    weighter = Weighter(injector, primary_physical=vertex.physical)
    weights = np.array([weighter(t) for t in trees])
    n_w = _density(ccm, points[0], W183) * SIGMA * 100.0
    assert np.all(np.array([t.tree[0].record.interaction_vertex[0] for t in trees]) <= TX + 0.05 + 1e-12)
    assert weights.sum() == pytest.approx(n_w * 4e-6, rel=1e-6)
