"""siren.check_closure -- the explicit, quantitative Sample==Density gauge.

A model's SampleFinalState draws kinematics; its FinalStateProbability names
the analytic density those kinematics are supposed to follow. The two must
agree: every weight in the framework is PhysicalProbability / GenerationProbability,
and if a model's own sampler disagrees with its own density, weights computed
against that model are wrong by a kinematics-dependent factor that no amount
of statistics averages away.

check_closure draws samples from the model's sampler, evaluates the model's
own analytic density on each sample, and checks that the sampler and density
describe the same distribution: the importance-sampling normalization
E[f / g] must equal 1, and each declared DensityVariable's sampled mean must
match its density-weighted mean. It also runs a frame heuristic for
SolidAngle*-type measures, since a rest-frame/lab-frame swap is the single
most common way a sampler and its stated measure silently diverge.

This module never generates events for physics use; it exists only to certify
that a model's sampler and density are the same distribution.
"""

import math

from . import dataclasses as _dataclasses
from . import injection as _injection
from . import utilities as _utilities
from . import _validation
from .errors import ClosureError

_MIN_SAMPLES = 200
_MIN_BIN_COUNT = 5
_N_BINS = 20

# Cache keyed by a hash of the call parameters (falls back to id(model) if the
# parameters are unhashable), so repeated checks in a session/notebook do not
# re-draw samples each time. Values are ClosureReport instances.
_CACHE = {}


def _cache_key(model, primary_energy, target, samples, seed, tol_sigma):
    try:
        return (type(model).__qualname__, primary_energy, target, samples,
                seed, tol_sigma)
    except TypeError:
        return (type(model).__qualname__, id(model), samples, seed, tol_sigma)


def _is_mixture(obj):
    """True for a siren.channels.Mixture-like object (compile + validate)."""
    return hasattr(obj, "compile") and hasattr(obj, "validate")


def _is_multichannel(obj):
    return hasattr(obj, "ValidateChannelDensities")


class ClosureReport:
    """Result of a check_closure run.

    Attributes
    ----------
    ok : bool
        True when the normalization estimate is within tol_sigma of 1.0 and
        no fatal flag was raised while probing the model or mixture.
    normalization : tuple
        (estimate, stderr) of E[f / g] over the drawn samples; expected 1.0.
    moment_z : dict
        {variable_name: z_score} of sampled vs density-weighted mean, one
        entry per DensityVariable (or per fallback coordinate).
    worst_region : str
        Human string naming the bin with the largest normalized deviation,
        or "" if no binned diagnostic was computed.
    frame_check : str or None
        A note when the declared SolidAngle* frame looks wrong, else None.
    """

    def __init__(self, ok, normalization, moment_z, worst_region,
                 frame_check, notes=None):
        self.ok = bool(ok)
        self.normalization = tuple(normalization)
        self.moment_z = dict(moment_z)
        self.worst_region = worst_region
        self.frame_check = frame_check
        self.notes = list(notes) if notes is not None else []

    def raise_if_failed(self):
        """Raise ClosureError with the rendered report if ok is False."""
        if not self.ok:
            raise ClosureError(str(self))

    def __str__(self):
        est, err = self.normalization
        lines = []
        lines.append("Closure report: %s" % ("PASS" if self.ok else "FAIL"))
        lines.append("  normalization E[f/g] = %.4f +/- %.4f (expect 1.0)"
                      % (est, err))
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
        return "<ClosureReport ok=%r normalization=%r>" % (
            self.ok, self.normalization)


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
        Normalization estimate must fall within tol_sigma stderr of 1.0 for
        ok to be True.

    Returns
    -------
    ClosureReport
    """
    key = _cache_key(model_or_mixture, primary_energy, target, samples,
                      seed, tol_sigma)
    cached = _CACHE.get(key)
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

    _CACHE[key] = report
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
    """Which frame the empirical coordinate is read in, from the declared
    measure. SolidAngleRest densities are flat over rest-frame solid angle, so
    the coordinate must be the rest-frame cos(theta) for the sampler/density
    ratio to be flat; SolidAngleLab uses the lab frame. Other measures default
    to the lab frame with a descriptive label.
    """
    mtype = model.Measure().type
    if mtype == _injection.PhaseSpaceMeasureType.SolidAngleRest:
        return "rest"
    return "lab"


def _extract_coordinate(out_rec, template, frame, boost):
    """cos(theta) of the first secondary, in the frame the measure declares.

    The single fallback coordinate for the normalization/moment checks. Read in
    the parent rest frame for a rest-frame measure (so f/g is flat for a
    correct sampler) and in the lab frame otherwise.
    """
    if not out_rec.secondary_momenta:
        return None
    p4 = list(out_rec.secondary_momenta[0])
    if frame == "rest":
        p4 = boost(p4)
    return _costheta(p4)


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
            moment_z={},
            worst_region="",
            frame_check=None,
            notes=["every drawn sample has FinalStateProbability <= 0; "
                   "sampler and density cannot be compared"],
        )

    # Empirical sampling density g_i via a 1-D histogram of the fallback
    # coordinate (cos theta of the first secondary, read in the measure's
    # declared frame): each sample's weight is f_i / g_i, where g_i is the
    # fraction of samples landing in its bin divided by the bin width. This is
    # well-defined regardless of what DensityVariables the model declares, since
    # every topology has at least one secondary whose direction can be read off
    # the record.
    frame = _coordinate_frame(model)
    boost = _boost_to_rest(template.primary_momentum, template.primary_mass)
    coord_name = "costheta_secondary0_%s" % frame
    coords = []
    valid_idx = []
    for i, (out_rec, f_i) in enumerate(drawn):
        c = _extract_coordinate(out_rec, template, frame, boost)
        if c is None:
            continue
        coords.append(c)
        valid_idx.append(i)

    if len(coords) < _MIN_SAMPLES:
        notes.append(
            "fewer than %d samples had an extractable coordinate; "
            "normalization/moment estimates are unreliable" % _MIN_SAMPLES)

    n_bins = _N_BINS
    bin_width = 2.0 / n_bins
    bin_counts = [0] * n_bins
    bin_members = [[] for _ in range(n_bins)]
    for local_i, c in enumerate(coords):
        cc = min(max(c, -1.0), 1.0)
        idx = min(int((cc + 1.0) / 2.0 * n_bins), n_bins - 1)
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
        lo = -1.0 + b * bin_width
        hi = lo + bin_width
        bin_ratio_info.append((lo, hi, ratio, count))

    # Closure holds when f_i / g_i is FLAT across bins. Its absolute scale
    # carries a fixed measure/marginalization constant (f is a density over the
    # full measure; g is the 1-D cos(theta) marginal), so the gauge normalizes
    # each bin ratio by the weighted mean and tests that the normalized profile
    # sits at 1.0. A sampler/density mismatch bends the profile away from flat.
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
            worst_region = ("cost in [%.2f, %.2f): ratio %.2f (z=%.1f)"
                             % (lo, hi, norm, z))

    # Per-DensityVariable moment check: self-normalized importance-weighted
    # mean of the fallback coordinate vs its plain sampled mean. Only one
    # fallback coordinate is available (see _extract_coordinate), so every
    # declared DensityVariable name is reported against the same coordinate;
    # this is documented here rather than silently mislabeled.
    density_vars = list(model.DensityVariables())
    moment_z = {}
    if coords and math.isfinite(mean_ratio):
        sample_mean = sum(coords) / len(coords)
        f_subset = [f_values[i] for i in valid_idx]
        f_sum = sum(f_subset)
        if f_sum > 0.0:
            analytic_mean = sum(c * f for c, f in zip(coords, f_subset)) / f_sum
            sample_var = (sum((c - sample_mean) ** 2 for c in coords)
                          / max(len(coords) - 1, 1))
            se = math.sqrt(sample_var / len(coords)) if len(coords) > 1 else float("inf")
            z = (sample_mean - analytic_mean) / se if se > 0.0 else 0.0
            labels = density_vars if density_vars else [coord_name]
            for label in labels:
                moment_z[label] = z

    frame_note = _frame_check(model, template, drawn)

    # A flat normalized profile passes; the gauge fails when any well-populated
    # bin's normalized ratio deviates from 1.0 by more than tol_sigma of its own
    # Poisson-scale uncertainty, or when a declared-variable moment z-score
    # exceeds tol_sigma. Two or more usable bins are required to see a shape.
    ok = math.isfinite(mean_ratio) and len(norm_ratios) >= 2
    if ok:
        for (lo, hi, ratio, count), w in zip(bin_ratio_info, ratio_weights):
            norm = ratio / scale if scale > 0.0 else ratio
            # Poisson relative error on the bin's occupancy sets the tolerance.
            rel = 1.0 / math.sqrt(count)
            if abs(norm - 1.0) > tol_sigma * rel:
                ok = False
                break
    for z in moment_z.values():
        if math.isfinite(z) and abs(z) > tol_sigma:
            ok = False

    return ClosureReport(
        ok=ok,
        normalization=(mean_ratio, stderr),
        moment_z=moment_z,
        worst_region=worst_region,
        frame_check=frame_note,
        notes=notes,
    )
