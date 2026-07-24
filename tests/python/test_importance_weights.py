"""Importance-weighted external tables through injection and weighting.

G4BNB-style dk2nu production applies importance reweighting: the recorded
rows carry nimpwt values spanning orders of magnitude, and the physical
ensemble is the nimpwt-weighted one. That biasing of the recorded parent
list is intentional, so injection samples the supplied rows UNIFORMLY
(or by an explicit sampling bias) and the importance weights only encode
the return to the physical distribution: the reader puts nimpwt/POT in
the table's "weight" column, PhysicalDensity reports the physical row
density on the physical side of the weight ratio (never cancelled
against the injection side), and the weight total is the distribution's
physical normalization.

Files produced with importance mode "None" (nimpwt identically 1) cannot
see either defect, which is how the original bugs survived the NuMI-file
validation: these tests pin the equivalence of an importance-weighted
table and its expanded equal-weight representation, so an
importance-blind regression fails loudly here.

The chain tests also pin the two companion contracts the beam-decays
workflow exposed:

  - a BoundedVertex secondary whose ray misses the fiducial volume fails
    the attempt (no silent fallback to the unrestricted path), and
  - event weights normalize by generation ATTEMPTS, so those failed
    attempts enter the estimator as legitimate zero-weight samples and
    the surviving weights do not inflate by 1/efficiency.

The single-row test checks the injection efficiency against a pure-numpy
computation of the boosted-decay geometric acceptance, and the rate
against the detector-directed channel, so the absolute normalization has
an anchor outside the channel machinery itself.

ASCII only. No network. Data-free cross sections; the chain tests need
only the in-repo CCM detector tables and skip when they are absent.
"""

import math
import os

import numpy as np
import pytest

siren = pytest.importorskip("siren")

from siren import _util
from siren import channels
from siren import detector
from siren import distributions
from siren import expand
from siren.Injector import Injector
from siren.Weighter import Weighter

M_PI = 0.13957039
M_MU = 0.1056583755


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


# --------------------------------------------------------------------------- #
# Chain-level pins (CCM detector, data-free cross sections)                    #
# --------------------------------------------------------------------------- #

def _skip_unless_ccm_data():
    try:
        det_dir = _util.get_detector_model_path("CCM")
    except ValueError as e:
        pytest.skip(f"CCM detector model path unavailable: {e}")
    missing = [p for p in (os.path.join(det_dir, "materials.dat"),
                           os.path.join(det_dir, "densities.dat"))
               if not os.path.exists(p)]
    if missing:
        pytest.skip("CCM detector data missing: " + ", ".join(missing))


def _load_ccm_detector():
    dm = detector.DetectorModel()
    det_dir = _util.get_detector_model_path("CCM")
    dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
    dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    return dm


BOX_CENTER = (0.0, 0.0, 0.0)
BOX_WIDTHS = (0.8, 0.8, 0.8)
POT = 1.0e4


def _flat_cross_section():
    return siren.interactions.TrivialCrossSection(
        {siren.particles.NuMu: ([0.001, 100.0], [1.0e-38, 1.0e-38])},
        [siren.particles.Nucleon])


def _pion_table(positions, momenta, nimpwts):
    """dk2nu-shaped in-memory table: weight = nimpwt/POT."""
    keys = ["E", "px", "py", "pz", "x", "y", "z", "m", "weight"]
    rows = []
    for pos, mom, w in zip(positions, momenta, nimpwts):
        p = np.asarray(mom, dtype=float)
        E = math.sqrt(float(p @ p) + M_PI * M_PI)
        rows.append([E, p[0], p[1], p[2],
                     pos[0], pos[1], pos[2], M_PI, w / POT])
    return distributions.PrimaryExternalDistribution(keys, rows)


def _run_chain(detector_model, external, kinematics, events, seed):
    decay = siren.resources.processes.BeamDecays.MesonTwoBodyLeptonicDecay(211)
    box = siren.geometry.Box(widths=BOX_WIDTHS, center=BOX_CENTER)
    parent = siren.Vertex(
        siren.dataclasses.ParticleType(211),
        decay,
        distributions=[external],
        physical=[external],
        weighting=siren.Fixed(),
        kinematics=kinematics(box) if callable(kinematics) else kinematics,
        expand=(expand.child("NuMu"),),
    )
    neutrino = siren.Vertex(
        "NuMu",
        _flat_cross_section(),
        position=siren.dist.BoundedVertex(box, np.inf),
        expand=(expand.depth_below(0),),
    )
    injector = Injector(detector=detector_model, primary=parent,
                        secondaries=(neutrino,), events=events, seed=seed)
    weighter = Weighter(injector, primary_physical=parent.physical)
    results = siren.generate(injector, weighter, events=events,
                             on_shortfall="error")
    weights = np.array([w for _, w in results], dtype=float)
    energies = np.array(
        [ev.tree[0].record.primary_momentum[0] for ev, _ in results])
    return weights, energies, injector.report()


def _sigma_of_sum(weights):
    return math.sqrt(float((weights ** 2).sum()))


def test_importance_weighted_and_expanded_tables_agree():
    """An importance-weighted table and its expanded equal-weight
    representation of the same physical ensemble yield the same rate and
    spectrum. Importance-mode-None files exercise only the expanded shape,
    which is why they masked the original defect; this pins the pair."""
    _skip_unless_ccm_data()
    dm = _load_ccm_detector()

    rng = np.random.default_rng(4)
    n_base = 40
    positions = [(0.1 * float(rng.uniform(-1, 1)),
                  0.1 * float(rng.uniform(-1, 1)),
                  -2.0 + 0.2 * float(rng.uniform(-1, 1)))
                 for _ in range(n_base)]
    momenta = [(0.02 * float(rng.uniform(-1, 1)),
                0.02 * float(rng.uniform(-1, 1)),
                float(rng.uniform(0.3, 1.5)))
               for _ in range(n_base)]
    # Importance weights spanning four orders of magnitude.
    nimpwts = 10.0 ** rng.uniform(-3, 1, size=n_base)

    weighted = _pion_table(positions, momenta, nimpwts)

    # Expanded representation: duplicate each row proportional to its
    # weight; every clone carries an equal share so the total is preserved.
    copies = np.maximum(1, np.round(
        nimpwts / nimpwts.sum() * 4000)).astype(int)
    total = nimpwts.sum()
    per_clone = total / copies.sum()
    exp_positions, exp_momenta, exp_weights = [], [], []
    for i in range(n_base):
        exp_positions.extend([positions[i]] * copies[i])
        exp_momenta.extend([momenta[i]] * copies[i])
        exp_weights.extend([per_clone] * copies[i])
    expanded = _pion_table(exp_positions, exp_momenta, exp_weights)

    events = 400
    toward = lambda box: channels.toward("NuMu", box)
    w_a, e_a, _ = _run_chain(dm, weighted, toward, events, seed=11)
    w_b, e_b, _ = _run_chain(dm, expanded, toward, events, seed=13)

    rate_a, rate_b = w_a.sum(), w_b.sum()
    sigma = math.hypot(_sigma_of_sum(w_a), _sigma_of_sum(w_b))
    assert abs(rate_a - rate_b) < 5 * sigma, (rate_a, rate_b, sigma)
    # The expanded-copy rounding distorts the ensemble by up to a few
    # percent on top of the MC error, so pin the spectra loosely.
    mean_a = float((w_a * e_a).sum() / rate_a)
    mean_b = float((w_b * e_b).sum() / rate_b)
    assert mean_a == pytest.approx(mean_b, rel=0.25), (mean_a, mean_b)


def _numpy_geometric_acceptance(position, momentum, n=200000, seed=5):
    """Boosted isotropic two-body decay: fraction of neutrino rays through
    the box, by direct numpy construction (no siren machinery)."""
    rng = np.random.default_rng(seed)
    p_parent = np.asarray(momentum, dtype=float)
    p_mag = float(np.linalg.norm(p_parent))
    E_parent = math.hypot(p_mag, M_PI)
    p_star = (M_PI ** 2 - M_MU ** 2) / (2.0 * M_PI)

    cos_t = rng.uniform(-1, 1, size=n)
    phi = rng.uniform(0, 2 * np.pi, size=n)
    sin_t = np.sqrt(1 - cos_t ** 2)
    p_rest = p_star * np.stack(
        [sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t], axis=1)

    beta = p_parent / E_parent
    gamma = E_parent / M_PI
    beta_mag = float(np.linalg.norm(beta))
    beta_hat = beta / beta_mag
    p_par = p_rest @ beta_hat
    E_rest = p_star  # massless neutrino
    E_lab = gamma * (E_rest + beta_mag * p_par)
    p_par_lab = gamma * (p_par + beta_mag * E_rest)
    p_lab = p_rest + np.outer(p_par_lab - p_par, beta_hat)

    lo = np.asarray(BOX_CENTER) - np.asarray(BOX_WIDTHS) / 2
    hi = np.asarray(BOX_CENTER) + np.asarray(BOX_WIDTHS) / 2
    p0 = np.asarray(position, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        t0 = (lo - p0) / p_lab
        t1 = (hi - p0) / p_lab
    tmin = np.nanmax(np.minimum(t0, t1), axis=1)
    tmax = np.nanmin(np.maximum(t0, t1), axis=1)
    hits = tmax > np.maximum(tmin, 0.0)
    return float(hits.mean()), E_lab


def test_single_row_physical_channel_closed_form():
    """Single-row anchor: with kinematics=channels.physical() the injection
    efficiency must equal the geometric acceptance of the physically
    distributed decay (fail-on-miss + honest 1/(4 pi) sampling), the event
    weights must normalize by attempts (not successes), and the resulting
    rate must agree with the detector-directed channel on the same row.
    Any silent re-aiming, silent fallback vertex placement, or
    success-based normalization moves one of these pins by a large
    factor."""
    _skip_unless_ccm_data()
    dm = _load_ccm_detector()

    position = (0.0, 0.0, -2.5)
    momentum = (0.0, 0.0, 0.5)
    nimpwt = 0.37
    external = _pion_table([position], [momentum], [nimpwt])
    assert external.normalization == pytest.approx(nimpwt / POT, rel=1e-12)

    acceptance, _ = _numpy_geometric_acceptance(position, momentum)

    events = 150
    w_p, _, report_p = _run_chain(
        dm, external, channels.physical(), events, seed=21)
    eff = report_p.successes / report_p.attempts
    sigma_eff = math.sqrt(acceptance * (1 - acceptance) / report_p.attempts)
    assert abs(eff - acceptance) < 5 * sigma_eff + 0.01, (eff, acceptance)
    assert report_p.attempts > 2 * report_p.successes  # misses really occur

    external_b = _pion_table([position], [momentum], [nimpwt])
    toward = lambda box: channels.toward("NuMu", box)
    w_t, _, _ = _run_chain(dm, external_b, toward, events, seed=23)

    rate_p, rate_t = w_p.sum(), w_t.sum()
    sigma = math.hypot(_sigma_of_sum(w_p), _sigma_of_sum(w_t))
    assert abs(rate_p - rate_t) < 5 * sigma + 0.15 * rate_t, (rate_p, rate_t)


def test_missed_fiducial_volume_fails_the_attempt():
    """A row aimed away from the box: every attempt must fail with
    NoPathThroughVolume rather than silently placing an out-of-volume
    vertex."""
    _skip_unless_ccm_data()
    dm = _load_ccm_detector()

    external = _pion_table([(0.0, 0.0, -2.5)], [(0.0, 0.0, -5.0)], [1.0])
    decay = siren.resources.processes.BeamDecays.MesonTwoBodyLeptonicDecay(211)
    box = siren.geometry.Box(widths=BOX_WIDTHS, center=BOX_CENTER)
    parent = siren.Vertex(
        siren.dataclasses.ParticleType(211), decay,
        distributions=[external], physical=[external],
        weighting=siren.Fixed(),
        kinematics=channels.toward("NuMu", box),
        expand=(expand.child("NuMu"),),
    )
    neutrino = siren.Vertex(
        "NuMu", _flat_cross_section(),
        position=siren.dist.BoundedVertex(box, np.inf),
        expand=(expand.depth_below(0),),
    )
    injector = Injector(detector=dm, primary=parent, secondaries=(neutrino,),
                        events=5, seed=3, max_attempts=200)
    trees = injector.generate(events=5, on_shortfall="ignore")
    report = injector.report()
    # A 5 GeV/c backward pion cannot send the neutrino into the box: the
    # boosted cone never reaches it, so nothing may succeed.
    assert len(trees) == 0
    assert report.attempts == 200
    reasons = {bucket.reason_name for bucket in report.by_vertex}
    assert any("NoPathThroughVolume" in r for r in reasons), reasons
