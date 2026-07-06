"""siren.check_closure -- the explicit, quantitative Sample==Density gauge.

A model's SampleFinalState draws kinematics; its FinalStateProbability names
the analytic density those kinematics are supposed to follow. The two must
agree: every weight in the framework is PhysicalProbability / GenerationProbability,
and if a model's own sampler disagrees with its own density, weights computed
against that model are wrong by a kinematics-dependent factor that no amount
of statistics averages away.

check_closure runs two complementary checks:

* An ABSOLUTE normalization check. Where the declared measure has a
  self-contained reference density (SolidAngleRest 2-body: uniform 1/(4*pi)
  per steradian over the sphere), it draws validation samples from that
  reference, evaluates FinalStateProbability on each, and estimates the
  integral of the density over the measure as E[f / g_ref], which must equal
  1. This is the ONLY check that catches a uniform scale error -- a density
  off by a constant factor everywhere -- because it fixes the absolute scale.

* A SHAPE (flatness) check. It bins f / g_empirical over a fallback
  coordinate; because g_empirical is a 1-D marginal carrying an unknown
  constant, each bin ratio is normalized by their weighted mean before testing
  for a flat profile at 1.0. This detects a cos(theta)-dependent mismatch
  between sampler and density but is BLIND to a uniform scale error by
  construction. Each declared DensityVariable's sampled mean is also compared
  to its density-weighted mean, and a frame heuristic flags a rest-frame /
  lab-frame swap for SolidAngle*-type measures.

This module never generates events for physics use; it exists only to certify
that a model's sampler and density are the same distribution.
"""

import math
import weakref

from . import dataclasses as _dataclasses
from . import injection as _injection
from . import utilities as _utilities
from . import _validation
from .errors import ClosureError

_MIN_SAMPLES = 200
_MIN_BIN_COUNT = 5
_N_BINS = 20

# Below this empirical cos(theta) spread the coordinate is single-valued: it
# carries no shape, so the flatness shape check has nothing to gauge and the
# empirical-range rebin would divide by ~0. Above it the rebin resolves the
# shape a cone-directed or forward-collapsed sampler actually populates.
_EMPIRICAL_MIN_WIDTH = 1e-9

# Below this spread a coordinate is treated as kinematically pinned for the
# purpose of CHOOSING which secondary to read: binning a micro-range
# amplifies Poisson noise into spurious moment shifts, so the gauge prefers
# a secondary whose angular spread is at least this wide.
_RESOLVABLE_MIN_WIDTH = 1e-3

# Result cache keyed on the model INSTANCE (via a WeakKeyDictionary) so that
# two objects of the same class with different internal state never collide,
# and a cached report dies with the model that produced it. Each instance maps
# to {call_params: ClosureReport}. A model that cannot be weakly referenced is
# simply not cached.
_CACHE = weakref.WeakKeyDictionary()


def _params_key(primary_energy, target, samples, seed, tol_sigma):
    """The per-instance sub-key: the call parameters other than the model."""
    return (primary_energy, target, samples, seed, tol_sigma)


def _cache_get(model, params):
    """Return a cached ClosureReport for (model instance, params), or None.

    Unhashable/unweakrefable models never hit the cache; identity is exact.
    """
    try:
        per_model = _CACHE.get(model)
    except TypeError:
        return None
    if per_model is None:
        return None
    return per_model.get(params)


def _cache_put(model, params, report):
    """Store a report under the model instance; skip silently if unweakrefable."""
    try:
        per_model = _CACHE.get(model)
        if per_model is None:
            per_model = {}
            _CACHE[model] = per_model
    except TypeError:
        return
    per_model[params] = report


def _is_mixture(obj):
    """True for a siren.channels.Mixture-like object (compile + validate)."""
    return hasattr(obj, "compile") and hasattr(obj, "validate")


def _is_multichannel(obj):
    return hasattr(obj, "ValidateChannelDensities")


class ClosureReport:
    """Result of a check_closure run.

    Two orthogonal checks feed ``ok``:

    * ``normalization`` -- an ABSOLUTE test that the model's density integrates
      to 1. It is populated only when the declared measure has a self-contained
      reference density the gauge can integrate against (currently a
      SolidAngleRest 2-body measure, whose reference is the uniform
      ``1/(4*pi)`` per steradian over the sphere). Validation samples are drawn
      from that reference, so the estimate ``E[f / g_ref]`` equals the integral
      of the model's FinalStateProbability over the reference measure and must
      be 1.0. This is the ONLY check that catches a uniform scale error in
      FinalStateProbability -- a density off by a constant factor everywhere.
      When no self-contained reference exists it is ``(nan, nan)`` and does not
      gate ``ok``.

    * ``flatness`` -- a SHAPE-only test. It bins ``f / g_empirical`` over a
      fallback coordinate and, because ``g_empirical`` is a 1-D marginal that
      carries an unknown constant, normalizes each bin ratio by their weighted
      mean before testing for a flat profile at 1.0. It therefore detects a
      cos(theta)-dependent (shape) mismatch between sampler and density but is
      blind to a uniform scale error BY CONSTRUCTION. Reported as
      ``(mean, stderr)`` of the self-normalized profile (mean is ~1.0 by
      construction; ``worst_region`` names the most deviant bin).

    Attributes
    ----------
    ok : bool
        True when every populated check passes within tol_sigma and no fatal
        flag was raised while probing the model or mixture. Specifically: the
        normalization estimate (when present) is within tol_sigma of 1.0, no
        flatness bin deviates from 1.0 by more than tol_sigma of its Poisson
        uncertainty, and no moment z-score exceeds tol_sigma.
    normalization : tuple
        (estimate, stderr) of the absolute integral E[f / g_ref] over the
        reference measure; expected 1.0. ``(nan, nan)`` when the declared
        measure has no self-contained reference density.
    flatness : tuple
        (mean, stderr) of the self-normalized bin ratios of f / g_empirical;
        a shape diagnostic only. ``(nan, nan)`` when no binned profile was
        computed.
    moment_z : dict
        {variable_name: z_score} of sampled vs density-weighted mean, one
        entry per DensityVariable (or per fallback coordinate).
    worst_region : str
        Human string naming the flatness bin with the largest normalized
        deviation, or "" if no binned diagnostic was computed.
    frame_check : str or None
        A note when the declared SolidAngle* frame looks wrong, else None.
    """

    def __init__(self, ok, normalization, flatness, moment_z, worst_region,
                 frame_check, notes=None):
        self.ok = bool(ok)
        self.normalization = tuple(normalization)
        self.flatness = tuple(flatness)
        self.moment_z = dict(moment_z)
        self.worst_region = worst_region
        self.frame_check = frame_check
        self.notes = list(notes) if notes is not None else []

    def raise_if_failed(self):
        """Raise ClosureError with the rendered report if ok is False."""
        if not self.ok:
            raise ClosureError(str(self))

    def __str__(self):
        n_est, n_err = self.normalization
        f_est, f_err = self.flatness
        lines = []
        lines.append("Closure report: %s" % ("PASS" if self.ok else "FAIL"))
        if math.isfinite(n_est):
            lines.append(
                "  normalization (absolute) E[f/g_ref] = %.4f +/- %.4f "
                "(expect 1.0)" % (n_est, n_err))
        else:
            lines.append(
                "  normalization (absolute): (not available for this measure; "
                "cannot detect a uniform scale error)")
        lines.append(
            "  flatness (shape only) = %.4f +/- %.4f "
            "(self-normalized bin ratios; blind to a uniform scale error)"
            % (f_est, f_err))
        if self.moment_z:
            lines.append("  moment z-scores:")
            for name in sorted(self.moment_z):
                lines.append("    %s: z=%.2f" % (name, self.moment_z[name]))
        else:
            lines.append("  moment z-scores: (none computed)")
        if self.worst_region:
            lines.append("  worst region: %s" % self.worst_region)
        if self.frame_check:
            lines.append("  frame check: %s" % self.frame_check)
        for note in self.notes:
            lines.append("  note: %s" % note)
        return "\n".join(lines)

    def __repr__(self):
        return ("<ClosureReport ok=%r normalization=%r flatness=%r>"
                % (self.ok, self.normalization, self.flatness))


def check_closure(model_or_mixture, *, primary_energy=None, target=None,
                   samples=2000, seed=0, tol_sigma=4.0):
    """Quantify whether a model's sampler matches its own analytic density.

    Parameters
    ----------
    model_or_mixture : object
        A Python model (DecayModel/CrossSectionModel subclass, or any object
        exposing GetPossibleSignatures/SampleFinalState/FinalStateProbability/
        DensityVariables/Measure/Topology), a siren.channels.Mixture, or a
        siren.injection.MultiChannelPhaseSpace.
    primary_energy : float, optional
        Parent/primary energy for the template record. Defaults to a small
        placeholder energy when not given.
    target : object, optional
        Target particle type, forwarded to detector-dependent mixture probes
        when available. Ignored for a plain model (decays have no target).
    samples : int
        Number of samples to draw for the model-path gauge.
    seed : int
        Seed for the validation RNG stream. Always its own stream, never a
        generation stream, so closure checks do not perturb injection.
    tol_sigma : float
        Tolerance in standard errors for every gated check: the absolute
        normalization estimate must fall within tol_sigma of 1.0, each flatness
        bin within tol_sigma of its Poisson uncertainty, and each moment
        z-score below tol_sigma, for ok to be True.

    Returns
    -------
    ClosureReport
    """
    params = _params_key(primary_energy, target, samples, seed, tol_sigma)
    cached = _cache_get(model_or_mixture, params)
    if cached is not None:
        return cached

    if _is_mixture(model_or_mixture) or _is_multichannel(model_or_mixture):
        report = _check_closure_mixture(
            model_or_mixture, primary_energy=primary_energy, target=target,
            samples=samples, seed=seed, tol_sigma=tol_sigma)
    else:
        report = _check_closure_model(
            model_or_mixture, primary_energy=primary_energy, target=target,
            samples=samples, seed=seed, tol_sigma=tol_sigma)

    _cache_put(model_or_mixture, params, report)
    return report


# ---------------------------------------------------------------------- #
#  Mixture / MultiChannelPhaseSpace path                                  #
# ---------------------------------------------------------------------- #

def _check_closure_mixture(obj, *, primary_energy, target, samples, seed,
                            tol_sigma):
    """Gauge a Mixture or MultiChannelPhaseSpace via ValidateChannelDensities.

    Kept defensive: any failure to probe (missing detector model, unresolved
    signature, engine error) is caught and reported as a note rather than
    propagated, since a mixture is frequently checked before a detector model
    or concrete signature exists.
    """
    notes = []
    ok = True

    # The detector-dependent ValidateChannelDensities probe needs a detector
    # model and validation RNG, which exist only inside Injector._build. The
    # standalone mixture gauge runs the detector-free Mixture.validate()
    # (ValidateChannelsDetailed + normalization + ConvertDensity probe) and
    # surfaces its result; the full density probe fires at build time.
    try:
        if hasattr(obj, "validate"):
            obj.validate()
            notes.append("mixture validate() passed; the detector-dependent "
                         "density probe runs at Injector build")
        elif _is_multichannel(obj):
            notes.append(
                "bare MultiChannelPhaseSpace has no detector model here; the "
                "density probe runs at Injector build")
        else:
            notes.append("object exposes no mixture validation surface")
    except Exception as exc:  # noqa: BLE001 -- defensive by contract
        ok = False
        notes.append("mixture validation raised %s: %s"
                      % (type(exc).__name__, exc))

    return ClosureReport(
        ok=ok,
        normalization=(float("nan"), float("nan")),
        flatness=(float("nan"), float("nan")),
        moment_z={},
        worst_region="",
        frame_check=None,
        notes=notes if notes else ["mixture path: no diagnostics raised"],
    )


# ---------------------------------------------------------------------- #
#  Model path                                                              #
# ---------------------------------------------------------------------- #

def _mass_or(default, ptype):
    """Particle mass, or a fallback for BSM types absent from the mass map."""
    try:
        return _dataclasses.GetParticleMass(ptype)
    except Exception:  # noqa: BLE001 -- BSM particles are not in the C++ map
        return default


def _make_template(model, primary_energy, target):
    signatures = model.GetPossibleSignatures()
    if not signatures:
        raise ClosureError(
            "%s.GetPossibleSignatures() returned no signatures; cannot "
            "build a template record" % type(model).__name__)
    signature = signatures[0]
    energy = primary_energy if primary_energy is not None else 0.05
    # A BSM primary (e.g. an HNL) is not in the mass map; fall back to a mass
    # comfortably below the energy so the parent boost stays well-defined.
    primary_mass = _mass_or(0.5 * energy, signature.primary_type)
    rec = _validation.build_template_record(
        signature, primary_mass=primary_mass, energy=energy)
    rec.secondary_masses = [_mass_or(0.001, t) for t in signature.secondary_types]
    if target is not None:
        rec.signature.target_type = target
        rec.target_mass = _mass_or(0.001, target)
    return rec, signature


def _blank_like(template):
    """A fresh InteractionRecord carrying the template's kinematics.

    InteractionRecord has no copy constructor binding, so the fields the gauge
    and FinalStateProbability read (signature, primary, secondary masses) are
    copied onto a default-constructed record; secondary momenta are filled by
    the CSDR finalize step.
    """
    rec = _dataclasses.InteractionRecord()
    rec.signature = template.signature
    rec.primary_mass = template.primary_mass
    rec.primary_momentum = template.primary_momentum
    rec.secondary_masses = list(template.secondary_masses)
    rec.secondary_momenta = [list(p) for p in template.secondary_momenta]
    rec.secondary_helicities = list(template.secondary_helicities)
    rec.interaction_vertex = template.interaction_vertex
    rec.primary_initial_position = template.primary_initial_position
    rec.target_mass = template.target_mass
    return rec


def _draw_samples(model, template, random, n):
    """Draw n (record, f_i) pairs by sampling and evaluating the density.

    Each draw builds a fresh CrossSectionDistributionRecord from the template
    so samples are independent and never contaminate each other's secondary
    state.
    """
    drawn = []
    for _ in range(n):
        csdr = _dataclasses.CrossSectionDistributionRecord(template)
        model.SampleFinalState(csdr, random)
        out_rec = _blank_like(template)
        csdr.finalize(out_rec)
        f_i = float(model.FinalStateProbability(out_rec))
        if not math.isfinite(f_i) or f_i < 0.0:
            f_i = 0.0
        drawn.append((out_rec, f_i))
    return drawn


def _reference_channel(model, signature):
    """A self-contained reference channel with a known analytic density, or None.

    Returns (channel, g_ref) where g_ref is the channel's density over its own
    measure -- the value the model's FinalStateProbability is integrated
    against. Only SolidAngleRest 2-body has such a reference here: the
    Isotropic2BodyChannel draws directions uniform in rest-frame solid angle,
    so g_ref = 1/(4*pi) per steradian, independent of the sampled point.
    """
    measure = model.Measure()
    n_finals = len(signature.secondary_types)
    if (n_finals == 2
            and measure.type == _injection.PhaseSpaceMeasureType.SolidAngleRest):
        return _injection.Isotropic2BodyChannel(0), 1.0 / (4.0 * math.pi)
    return None, None


def _estimate_absolute_normalization(model, template, signature, random, n):
    """Estimate the absolute integral of FinalStateProbability over the measure.

    Draws n samples from a reference channel whose density g_ref over the
    declared measure is known analytically (uniform 1/(4*pi) for SolidAngleRest
    2-body). For each drawn record it evaluates f = FinalStateProbability and
    forms f / g_ref; the Monte Carlo mean of that ratio is an unbiased estimate
    of the integral of f over the reference measure, which must equal 1.0 for a
    correctly normalized density. Because the samples come from the reference
    -- NOT the model's own sampler -- this estimate is sensitive to a uniform
    scale error in f (a density scaled by a constant everywhere returns that
    constant instead of 1.0), which the shape-only flatness check cannot see.

    Returns (estimate, stderr), or (nan, nan) when no self-contained reference
    density is available for the declared measure.
    """
    channel, g_ref = _reference_channel(model, signature)
    if channel is None:
        return float("nan"), float("nan")

    ratios = []
    for _ in range(n):
        out_rec = _blank_like(template)
        channel.Sample(random, None, out_rec)
        f_i = float(model.FinalStateProbability(out_rec))
        if not math.isfinite(f_i) or f_i < 0.0:
            f_i = 0.0
        ratios.append(f_i / g_ref)

    if not ratios:
        return float("nan"), float("nan")
    m = len(ratios)
    mean = sum(ratios) / m
    if m > 1:
        var = sum((r - mean) ** 2 for r in ratios) / (m - 1)
        stderr = math.sqrt(var / m)
    else:
        stderr = float("inf")
    return mean, stderr


def _boost_to_rest(momentum, mass):
    """Boost a lab-frame 4-momentum into the rest frame of a parent 4-momentum.

    momentum, mass describe the parent; returns a function mapping a
    secondary 4-momentum (list[4], E,px,py,pz) into the parent rest frame.
    Standard 4-vector boost along the parent's 3-velocity.
    """
    e_p, px_p, py_p, pz_p = momentum
    if mass <= 0.0 or e_p <= mass:
        return lambda p4: p4
    beta = [px_p / e_p, py_p / e_p, pz_p / e_p]
    beta2 = beta[0] ** 2 + beta[1] ** 2 + beta[2] ** 2
    if beta2 <= 0.0:
        return lambda p4: p4
    gamma = e_p / mass

    def boost(p4):
        e, px, py, pz = p4
        p_dot_beta = px * beta[0] + py * beta[1] + pz * beta[2]
        coeff = (gamma - 1.0) / beta2
        new_e = gamma * (e - p_dot_beta)
        new_px = px + coeff * p_dot_beta * beta[0] - gamma * beta[0] * e
        new_py = py + coeff * p_dot_beta * beta[1] - gamma * beta[1] * e
        new_pz = pz + coeff * p_dot_beta * beta[2] - gamma * beta[2] * e
        return [new_e, new_px, new_py, new_pz]

    return boost


def _costheta(p4):
    """cos(theta) of a 3-momentum relative to the z-axis; 0.0 if at rest."""
    _, px, py, pz = p4
    p = math.sqrt(px * px + py * py + pz * pz)
    if p <= 0.0:
        return 0.0
    return pz / p


def _chi_square_uniform(values, n_bins=_N_BINS):
    """Chi-square of a cos(theta)-like sample against a uniform [-1, 1] fit."""
    n = len(values)
    if n == 0:
        return float("inf")
    counts = [0] * n_bins
    for v in values:
        v = min(max(v, -1.0), 1.0)
        idx = min(int((v + 1.0) / 2.0 * n_bins), n_bins - 1)
        counts[idx] += 1
    expected = n / float(n_bins)
    return sum((c - expected) ** 2 / expected for c in counts)


def _frame_check(model, template, samples_xy):
    """Flag a declared SolidAngle*-type measure that fits the other frame better.

    samples_xy is a list of (out_rec, f_i) pairs already drawn; reused here so
    the frame heuristic costs no extra sampling.
    """
    measure = model.Measure()
    mtype = measure.type
    solid_angle_types = (
        _injection.PhaseSpaceMeasureType.SolidAngleRest,
        _injection.PhaseSpaceMeasureType.SolidAngleLab,
    )
    if mtype not in solid_angle_types:
        return None

    parent_p4 = template.primary_momentum
    parent_mass = template.primary_mass
    boost = _boost_to_rest(parent_p4, parent_mass)

    lab_cos = []
    rest_cos = []
    for out_rec, _f in samples_xy:
        if not out_rec.secondary_momenta:
            continue
        p4_lab = out_rec.secondary_momenta[0]
        lab_cos.append(_costheta(p4_lab))
        p4_rest = boost(p4_lab)
        rest_cos.append(_costheta(p4_rest))

    if len(lab_cos) < _MIN_SAMPLES // 4:
        return None

    chi2_lab = _chi_square_uniform(lab_cos)
    chi2_rest = _chi_square_uniform(rest_cos)

    declared_name = "SolidAngleRest" if mtype == _injection.PhaseSpaceMeasureType.SolidAngleRest \
        else "SolidAngleLab"

    if mtype == _injection.PhaseSpaceMeasureType.SolidAngleRest and chi2_lab < chi2_rest:
        return ("declared SolidAngleRest but sampled directions look isotropic "
                "in LAB (chi2 rest=%.1f vs lab=%.1f); check the boost"
                % (chi2_rest, chi2_lab))
    if mtype == _injection.PhaseSpaceMeasureType.SolidAngleLab and chi2_rest < chi2_lab:
        return ("declared SolidAngleLab but sampled directions look isotropic "
                "in REST frame (chi2 lab=%.1f vs rest=%.1f); check the boost"
                % (chi2_lab, chi2_rest))
    return None


def _coordinate_frame(model):
    """Which frame the empirical coordinate is read in.

    Scatter topologies (Scatter2to2/Scatter2to3) read the collision
    centre-of-mass frame: the secondaries are isotropic (or, for a forward
    process, symmetric) about the beam axis in the CM, so the CM cos(theta) is
    the coordinate against which the sampler/density ratio is flat. The lab
    frame collimates every secondary of a fixed-target scatter to cos(theta) ~ 1
    (a single bin, stderr inf), and the primary rest frame is the wrong frame
    entirely for a two-body collision. The CM boost is by the total 4-momentum
    primary + target-at-rest (see _cm_boost).

    Otherwise the frame follows the declared measure. SolidAngleRest densities
    are flat over rest-frame solid angle, so the coordinate is the parent
    rest-frame cos(theta); SolidAngleLab uses the lab frame. Recursive2Body
    factors the final state into a spectator plus a pair whose orientation is
    isotropic in the parent rest frame, so its coordinate is the rest-frame
    cos(theta) too. Other measures default to the lab frame with a descriptive
    label.
    """
    topology = model.Topology()
    if topology in (_injection.PhaseSpaceTopology.Scatter2to2,
                    _injection.PhaseSpaceTopology.Scatter2to3):
        return "cm"
    mtype = model.Measure().type
    if mtype in (_injection.PhaseSpaceMeasureType.SolidAngleRest,
                 _injection.PhaseSpaceMeasureType.Recursive2Body):
        return "rest"
    return "lab"


def _cm_boost(template):
    """Boost into the collision centre-of-mass frame for a scatter template.

    The CM frame is the rest frame of the total 4-momentum primary +
    target-at-rest; its invariant mass is sqrt(s). Falls back to the identity
    when no target mass is present (a decay template) so the caller can use one
    code path for every frame.
    """
    p_primary = list(template.primary_momentum)
    m_target = getattr(template, "target_mass", 0.0) or 0.0
    total = [p_primary[0] + m_target, p_primary[1], p_primary[2], p_primary[3]]
    s = total[0] ** 2 - (total[1] ** 2 + total[2] ** 2 + total[3] ** 2)
    if s <= 0.0:
        return lambda p4: p4
    return _boost_to_rest(total, math.sqrt(s))


def _coordinate_secondary_index(model):
    """Which secondary's direction is the fallback coordinate.

    A Recursive2Body final state factors into a spectator plus a pair whose
    orientation is isotropic in the parent (or CM) rest frame; the SPECTATOR
    carries that isotropic orientation, so its cos(theta) is the coordinate f/g
    is flat in. A pair member's individual direction convolves the pair
    orientation with the sub-decay angle and is NOT flat, so reading it would
    manufacture a shape mismatch that is not in the density. Every other measure
    reads the first secondary.
    """
    measure = model.Measure()
    if measure.type == _injection.PhaseSpaceMeasureType.Recursive2Body:
        spectator = getattr(measure, "spectator", 0)
        if spectator is not None:
            return int(spectator)
    return 0


def _extract_coordinate(out_rec, template, frame, boost, index=0):
    """cos(theta) of the selected secondary, in the frame the gauge selected.

    The single fallback coordinate for the flatness/moment checks. Read in the
    parent rest frame for a rest-frame measure and in the collision CM frame for
    a scatter topology (so f/g is flat for a correct sampler); in the lab frame
    otherwise. Any non-lab frame applies the supplied boost. index selects which
    secondary is read (the Recursive2Body spectator; see
    _coordinate_secondary_index).
    """
    if not out_rec.secondary_momenta or index >= len(out_rec.secondary_momenta):
        return None
    p4 = list(out_rec.secondary_momenta[index])
    if frame != "lab":
        p4 = boost(p4)
    return _costheta(p4)


def _bin_ratios(coords, valid_idx, f_values, lo, hi, n_bins=_N_BINS):
    """Histogram f/g over [lo, hi] and return the per-bin ratio profile.

    Returns (bin_ratio_info, ratios, ratio_weights): bin_ratio_info is a list of
    (bin_lo, bin_hi, ratio, count) for every bin with at least _MIN_BIN_COUNT
    samples, and ratios/ratio_weights are the parallel ratio and count lists.
    g_b is the sampled fraction in the bin divided by the bin width; the ratio
    is the mean FinalStateProbability in the bin over g_b. A bin below the count
    floor is dropped so a sparse tail does not manufacture a spurious deviation.
    """
    width = hi - lo
    if width <= 0.0:
        return [], [], []
    bin_width = width / n_bins
    bin_counts = [0] * n_bins
    bin_members = [[] for _ in range(n_bins)]
    for local_i, c in enumerate(coords):
        cc = min(max(c, lo), hi)
        idx = min(int((cc - lo) / width * n_bins), n_bins - 1)
        bin_counts[idx] += 1
        bin_members[idx].append(valid_idx[local_i])

    n_valid = len(coords)
    ratios = []
    ratio_weights = []
    bin_ratio_info = []
    for b in range(n_bins):
        count = bin_counts[b]
        if count < _MIN_BIN_COUNT or n_valid == 0:
            continue
        g_b = count / float(n_valid) / bin_width
        f_mean_b = sum(f_values[i] for i in bin_members[b]) / count
        if g_b <= 0.0:
            continue
        ratio = f_mean_b / g_b
        ratios.append(ratio)
        ratio_weights.append(count)
        b_lo = lo + b * bin_width
        b_hi = b_lo + bin_width
        bin_ratio_info.append((b_lo, b_hi, ratio, count))
    return bin_ratio_info, ratios, ratio_weights


def _check_closure_model(model, *, primary_energy, target, samples, seed,
                          tol_sigma):
    n = max(int(samples), _MIN_SAMPLES)
    random = _utilities.SIREN_random(int(seed) & 0x7FFFFFFF)

    template, signature = _make_template(model, primary_energy, target)

    drawn = _draw_samples(model, template, random, n)
    f_values = [f for _rec, f in drawn]

    notes = []
    n_positive = sum(1 for f in f_values if f > 0.0)
    if n_positive == 0:
        return ClosureReport(
            ok=False,
            normalization=(float("nan"), float("nan")),
            flatness=(float("nan"), float("nan")),
            moment_z={},
            worst_region="",
            frame_check=None,
            notes=["every drawn sample has FinalStateProbability <= 0; "
                   "sampler and density cannot be compared"],
        )

    # Absolute normalization: integrate FinalStateProbability against a
    # reference channel whose density over the declared measure is known. This
    # is the only diagnostic that catches a uniform scale error in the density;
    # populated only when a self-contained reference exists (SolidAngleRest
    # 2-body). A separate validation stream so it does not perturb the draw
    # already consumed above.
    norm_random = _utilities.SIREN_random((int(seed) & 0x7FFFFFFF) ^ 0x5A5A5A5A)
    norm_est, norm_err = _estimate_absolute_normalization(
        model, template, signature, norm_random, n)
    has_reference = math.isfinite(norm_est)

    # Empirical sampling density g_i via a 1-D histogram of the fallback
    # coordinate (cos theta of the first secondary, read in the measure's
    # declared frame): each sample's weight is f_i / g_i, where g_i is the
    # fraction of samples landing in its bin divided by the bin width. This is
    # well-defined regardless of what DensityVariables the model declares, since
    # every topology has at least one secondary whose direction can be read off
    # the record.
    frame = _coordinate_frame(model)
    if frame == "cm":
        boost = _cm_boost(template)
    else:
        boost = _boost_to_rest(template.primary_momentum, template.primary_mass)
    sec_index = _coordinate_secondary_index(model)

    def _collect(index):
        cs, vi = [], []
        for i, (out_rec, _f_i) in enumerate(drawn):
            c = _extract_coordinate(out_rec, template, frame, boost, index)
            if c is None:
                continue
            cs.append(c)
            vi.append(i)
        return cs, vi

    coords, valid_idx = _collect(sec_index)

    # The gauge bins the f/g ratio, so any secondary with resolvable angular
    # spread is a valid coordinate; a kinematically pinned one (e.g. a heavy
    # recoil nucleus collinear in the CM) carries no shape and would force the
    # micro-range rebin to amplify noise. When the preferred secondary is
    # degenerate, read the secondary with the largest spread instead.
    def _spread(cs):
        return (max(cs) - min(cs)) if cs else 0.0

    if _spread(coords) < _RESOLVABLE_MIN_WIDTH:
        n_sec = len(template.secondary_masses)
        best = (sec_index, coords, valid_idx, _spread(coords))
        for alt in range(n_sec):
            if alt == sec_index:
                continue
            cs, vi = _collect(alt)
            if _spread(cs) > best[3]:
                best = (alt, cs, vi, _spread(cs))
        if best[0] != sec_index and best[3] >= _RESOLVABLE_MIN_WIDTH:
            notes.append(
                "secondary %d cos(theta) is kinematically pinned "
                "(spread %.3g); reading secondary %d (spread %.3g) instead"
                % (sec_index, _spread(coords), best[0], best[3]))
            sec_index, coords, valid_idx = best[0], best[1], best[2]

    coord_name = "costheta_secondary%d_%s" % (sec_index, frame)

    if len(coords) < _MIN_SAMPLES:
        notes.append(
            "fewer than %d samples had an extractable coordinate; "
            "normalization/moment estimates are unreliable" % _MIN_SAMPLES)

    # Bin f/g over the coordinate. The primary binning spans the full [-1, 1]
    # cos(theta) domain. When a forward-collapsed or cone-directed sampler leaves
    # every sample in one [-1, 1] bin (fewer than two populated bins), the shape
    # is unresolved there, so rebin over the empirical coordinate range
    # [min(c), max(c)] to gauge the shape the sampler actually populates. The
    # fallback fires only when the primary binning is degenerate, so a healthy
    # model keeps its [-1, 1] bins unchanged. A zero-width empirical range is a
    # single-valued coordinate carrying no shape and is guarded below.
    bin_ratio_info, ratios, ratio_weights = _bin_ratios(
        coords, valid_idx, f_values, -1.0, 1.0)
    shape_unresolved = False
    if len(ratios) < 2 and coords:
        c_lo = min(coords)
        c_hi = max(coords)
        if c_hi - c_lo > _EMPIRICAL_MIN_WIDTH:
            bin_ratio_info, ratios, ratio_weights = _bin_ratios(
                coords, valid_idx, f_values, c_lo, c_hi)
            coord_name += "_empirical"
            notes.append(
                "cos(theta) collapsed to a single [-1, 1] bin; the flatness "
                "shape was gauged over the empirical range [%.4g, %.4g]"
                % (c_lo, c_hi))
        else:
            shape_unresolved = True
            notes.append(
                "cos(theta) is single-valued (empirical width %.2e <= %.0e); "
                "the coordinate carries no shape, so the flatness shape check "
                "is trivially satisfied" % (c_hi - c_lo, _EMPIRICAL_MIN_WIDTH))

    # Flatness (SHAPE-only) check. f_i / g_i must be FLAT across bins for the
    # sampler and density to describe the same distribution. Its absolute scale
    # carries a fixed measure/marginalization constant (f is a density over the
    # full measure; g is the 1-D cos(theta) marginal), so the gauge normalizes
    # each bin ratio by the weighted mean and tests that the normalized profile
    # sits at 1.0. A shape mismatch bends the profile away from flat, but a
    # UNIFORM scale error cancels in the self-normalization and is invisible
    # here -- the absolute normalization check above is what catches that.
    if ratios:
        total_w = sum(ratio_weights)
        scale = sum(r * w for r, w in zip(ratios, ratio_weights)) / total_w
        norm_ratios = [r / scale for r in ratios] if scale > 0.0 else ratios
        mean_ratio = sum(nr * w for nr, w in zip(norm_ratios, ratio_weights)) / total_w
        if len(norm_ratios) > 1:
            var = sum(w * (nr - mean_ratio) ** 2
                      for nr, w in zip(norm_ratios, ratio_weights)) / total_w
            stderr = math.sqrt(max(var, 0.0) / len(norm_ratios))
        else:
            stderr = float("inf")
    else:
        scale = float("nan")
        norm_ratios = []
        mean_ratio = float("nan")
        stderr = float("nan")

    # The worst region is the bin whose normalized ratio deviates most from 1.0
    # in units of its own Poisson uncertainty, matching the ok gate.
    worst_region = ""
    if bin_ratio_info and scale > 0.0:
        best = None
        for lo, hi, ratio, count in bin_ratio_info:
            norm = ratio / scale
            rel = 1.0 / math.sqrt(count)
            z = (norm - 1.0) / rel if rel > 0.0 else 0.0
            if best is None or abs(z) > abs(best[4]):
                best = (lo, hi, norm, count, z)
        if best is not None:
            lo, hi, norm, _count, z = best
            worst_region = ("%s in [%.2f, %.2f): ratio %.2f (z=%.1f)"
                             % (coord_name, lo, hi, norm, z))

    # Per-DensityVariable moment check: the coordinate's density mean vs its
    # sampled mean, built from the SAME flatness bins so it gauges shape at the
    # flatness resolution rather than a frame mismatch. The density mean
    # reweights each bin's occupancy by the bin's self-normalized flatness ratio
    # r_b = ratio_b / scale (the density-over-sampler shape there); the sampled
    # mean weights bins by occupancy alone. Their difference is driven purely by
    # r_b deviating from 1. The uncertainty on that difference propagates each
    # bin's Poisson relative error 1/sqrt(count_b) through the same weights, so
    # the z-score is a fixed-scale test statistic: when the sampler matches the
    # density r_b = 1 +/- 1/sqrt(count_b) is pure noise and the shift stays within
    # its own uncertainty (z ~ O(1) independent of sample count), while a real
    # shape mismatch drives r_b systematically off 1 and lifts z. Only one
    # fallback coordinate is available (see _extract_coordinate), so every
    # declared DensityVariable name is reported against it.
    density_vars = list(model.DensityVariables())
    moment_z = {}
    if bin_ratio_info and scale > 0.0 and not shape_unresolved:
        total_count = sum(count for _lo, _hi, _r, count in bin_ratio_info)
        if total_count > 0:
            centers = [0.5 * (lo + hi) for lo, hi, _r, _c in bin_ratio_info]
            occ = [count / total_count for _lo, _hi, _r, count in bin_ratio_info]
            r_b = [ratio / scale for _lo, _hi, ratio, _c in bin_ratio_info]
            dens_w = [o * r for o, r in zip(occ, r_b)]
            dens_norm = sum(dens_w)
            sample_mean = sum(c * o for c, o in zip(centers, occ))
            if dens_norm > 0.0:
                density_mean = sum(c * w for c, w in zip(centers, dens_w)) / dens_norm
                # Propagate each bin's Poisson error on r_b (relative 1/sqrt(count))
                # through the density-mean weights; the shift is measured against
                # this estimator uncertainty so the z does not grow with N.
                var = 0.0
                for (lo, hi, ratio, count), c, o in zip(bin_ratio_info, centers, occ):
                    r = ratio / scale
                    d_shift = o * r * (c - density_mean) / dens_norm
                    rel = 1.0 / math.sqrt(count)
                    var += (d_shift * rel) ** 2
                sigma = math.sqrt(var)
                z = (density_mean - sample_mean) / sigma if sigma > 0.0 else 0.0
                labels = density_vars if density_vars else [coord_name]
                for label in labels:
                    moment_z[label] = z

    frame_note = _frame_check(model, template, drawn)

    # The single-coordinate flatness/moment construction compares the mean
    # JOINT density in a bin against the bin's occupancy. Those agree under
    # perfect closure only when the joint density is constant along the
    # dimensions the coordinate integrates out: single-dof measures
    # (SolidAngleRest 2-body, MandelstamQ2 2->2) and decays whose overall
    # orientation is isotropic. A Scatter2to3 Recursive2Body density is
    # forward-peaked in the pair angle, so the per-bin ratio varies even when
    # sampler and density agree exactly -- the shape verdict has no
    # jurisdiction there and must not gate ok in either direction. The
    # mixture-level E[f/g] closure test is the authoritative statement for
    # such processes.
    shape_applicable = not (
        model.Topology() == _injection.PhaseSpaceTopology.Scatter2to3
        and model.Measure().type
            == _injection.PhaseSpaceMeasureType.Recursive2Body)
    if not shape_applicable:
        notes.append(
            "shape verdict not applicable: the joint density of a Scatter2to3 "
            "Recursive2Body process varies along integrated-out dimensions, so "
            "bin ratios are not flat even under perfect closure; certify this "
            "model with a mixture-level E[f/g] closure test")

    # The gauge passes only when every populated check passes:
    #  * the absolute normalization estimate (when a reference density exists)
    #    is within tol_sigma of 1.0 -- this catches a uniform scale error;
    #  * the shape profile is flat: no well-populated bin's self-normalized
    #    ratio deviates from 1.0 by more than tol_sigma of its own Poisson-scale
    #    uncertainty (two or more usable bins are required to see a shape);
    #  * no declared-variable moment z-score exceeds tol_sigma.
    # A single-valued coordinate carries no shape, so the flatness requirement is
    # trivially satisfied and does not gate ok (the normalization and moment
    # checks still apply).
    ok = ((not shape_applicable) or shape_unresolved
          or (math.isfinite(mean_ratio) and len(norm_ratios) >= 2))

    if has_reference:
        if not math.isfinite(norm_est):
            ok = False
        elif norm_err > 0.0:
            if abs(norm_est - 1.0) > tol_sigma * norm_err:
                ok = False
        elif abs(norm_est - 1.0) > 1e-9:
            # Zero stderr means f/g_ref was constant over the samples: the
            # estimate is exact, so it must equal 1 outright.
            ok = False
    else:
        notes.append(
            "no self-contained reference density for measure %r; the absolute "
            "normalization was not checked, so a uniform scale error in "
            "FinalStateProbability would not be caught (only shape is tested)"
            % (model.Measure(),))

    if ok and shape_applicable:
        for (lo, hi, ratio, count), w in zip(bin_ratio_info, ratio_weights):
            norm = ratio / scale if scale > 0.0 else ratio
            # Poisson relative error on the bin's occupancy sets the tolerance.
            rel = 1.0 / math.sqrt(count)
            if abs(norm - 1.0) > tol_sigma * rel:
                ok = False
                break
    if shape_applicable:
        for z in moment_z.values():
            if math.isfinite(z) and abs(z) > tol_sigma:
                ok = False

    return ClosureReport(
        ok=ok,
        normalization=(norm_est, norm_err),
        flatness=(mean_ratio, stderr),
        moment_z=moment_z,
        worst_region=worst_region,
        frame_check=frame_note,
        notes=notes,
    )
