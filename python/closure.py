"""Check model sampling against an independent physical-density reference.

Two-body SolidAngleRest decays have a built-in uniform reference. Its weighted
samples test absolute normalization and the joint (cos(theta), phi) shape.
Other measures and detector-dependent mixtures report incomplete coverage;
configuration checks alone do not certify their sampling densities.
"""

import math

import numpy as np
from scipy.special import gammaincc

from . import dataclasses as _dataclasses
from . import injection as _injection
from . import utilities as _utilities
from . import _validation
from .errors import ClosureError


class ClosureReport:
    """Closure diagnostics at the tested kinematics and statistical resolution.

    ``checks`` maps check names to ``passed``, ``failed``, or ``incomplete``.
    ``ok`` requires every check to pass; ``complete`` requires every check to
    have been performed. A failed check takes precedence over incomplete ones.
    ``normalization`` is the absolute density integral and its standard error.
    ``flatness`` summarizes sampled/reference bin-probability ratios, normalized
    to remove the overall density scale. ``moment_z`` names only coordinates
    actually measured. ``joint_shape`` is the joint histogram's chi-square
    statistic, degrees of freedom, and p-value, or None when not tested.
    Missing normalization/flatness diagnostics are ``(nan, nan)``.
    """

    def __init__(self, *, checks, normalization=(float('nan'), float('nan')),
                 flatness=(float('nan'), float('nan')), moment_z=None,
                 worst_region='', frame_check=None, notes=(), joint_shape=None):
        self.checks = dict(checks)
        self.normalization = tuple(normalization)
        self.flatness = tuple(flatness)
        self.moment_z = dict(moment_z or {})
        self.worst_region = worst_region
        self.frame_check = frame_check
        self.notes = list(notes)
        self.joint_shape = joint_shape

    @property
    def complete(self):
        return bool(self.checks) and 'incomplete' not in self.checks.values()

    @property
    def ok(self):
        return bool(self.checks) and all(v == 'passed' for v in self.checks.values())

    @property
    def status(self):
        if 'failed' in self.checks.values():
            return 'failed'
        return 'passed' if self.ok else 'incomplete'

    def raise_if_failed(self):
        """Raise ClosureError for failed or incomplete certification."""
        if not self.ok:
            raise ClosureError(str(self))

    def __str__(self):
        lines = ['Closure report: ' + self.status.upper()]
        lines.extend('  %s: %s' % item for item in self.checks.items())
        if math.isfinite(self.normalization[0]):
            lines.append('  normalization E[f/g_ref] = %.4f +/- %.4f (expect 1.0)'
                         % self.normalization)
        if math.isfinite(self.flatness[0]):
            lines.append('  flatness (shape only) = %.4f +/- %.4f' % self.flatness)
        lines.extend('  %s: z=%.2f' % item for item in self.moment_z.items())
        if self.joint_shape is not None:
            lines.append('  joint angular test: chi2=%.2f, df=%d, p=%.3g'
                         % self.joint_shape)
        if self.worst_region:
            lines.append('  worst region: ' + self.worst_region)
        if self.frame_check:
            lines.append('  frame check: ' + self.frame_check)
        lines.extend('  note: ' + note for note in self.notes)
        return '\n'.join(lines)

    def __repr__(self):
        return '<ClosureReport status=%r checks=%r>' % (self.status, self.checks)


def check_closure(model_or_mixture, *, primary_energy=None, target=None,
                  record=None, samples=2000, seed=0, tol_sigma=4.0):
    """Test density normalization and sampling shape with independent RNG streams.

    ``record`` supplies the initial kinematics and signature without being
    modified. Otherwise a synthetic template uses the first model signature,
    ``primary_energy`` (default 0.05 GeV), and ``target``. Unknown masses get
    synthetic defaults; supply a record to test specific model parameters.
    ``record`` cannot be combined with the energy/target shortcuts.

    ``samples`` draws are made from each of the model and reference samplers.
    ``tol_sigma`` bounds normalization, bin and moment deviations in standard
    errors, including reference Monte Carlo uncertainty. The joint histogram
    test uses the corresponding two-sided Gaussian tail probability. Tests cover the
    joint angular histogram at the reported resolution, not arbitrarily fine
    structure. Unsupported references or insufficient coverage are incomplete.
    Results are recomputed because model state can change between calls.
    """
    if int(samples) != samples or samples < 2:
        raise ValueError('samples must be an integer >= 2')
    if not math.isfinite(tol_sigma) or tol_sigma <= 0:
        raise ValueError('tol_sigma must be finite and positive')
    if record is not None and (primary_energy is not None or target is not None):
        raise ValueError('record cannot be combined with primary_energy or target')
    obj = model_or_mixture
    if hasattr(obj, 'ValidateChannelDensities') or (
            hasattr(obj, 'compile') and hasattr(obj, 'validate')):
        report = ClosureReport(checks={'configuration': 'incomplete',
                                      'normalization': 'incomplete', 'shape': 'incomplete'})
        if hasattr(obj, 'validate'):
            try:
                obj.validate()
                report.checks['configuration'] = 'passed'
            except Exception as exc:
                report.checks['configuration'] = 'failed'
                report.notes.append('%s: %s' % (type(exc).__name__, exc))
        report.notes.append('The detector-dependent density probe was not run; '
                            'configuration validation does not establish closure.')
        return report
    template = (_dataclasses.InteractionRecord(record) if record is not None
                else _make_template(obj, primary_energy, target)[0])
    return _check_closure_model(obj, template, int(samples), seed, tol_sigma)


def _density(model, record, label):
    value = float(model.FinalStateProbability(record))
    if not math.isfinite(value) or value < 0:
        raise ClosureError('%s: FinalStateProbability must be finite and '
                           'nonnegative, got %r' % (label, value))
    return value


def _reference_context_issue(model, template, drawn, densities):
    """Do sampler-written parameters change the density at fixed kinematics?

    The angular reference cannot reconstruct a model's cached coordinates.
    Compare the density on each sampled point with the same point carrying
    only template parameters before trusting density evaluations on references.
    Unused sampler bookkeeping is harmless; required cached values make this
    reference incomplete, even if the density supplies a plausible fallback.
    """
    parameters = template.interaction_parameters
    for i, (record, density) in enumerate(zip(drawn, densities)):
        if record.interaction_parameters == parameters:
            continue
        probe = _dataclasses.InteractionRecord(record)
        probe.interaction_parameters = parameters
        try:
            reference_value = _density(model, probe, 'density context probe')
        except Exception as exc:
            return ('Model sample %d requires sampler-written interaction parameters '
                    'unavailable to the reference (%s: %s).'
                    % (i, type(exc).__name__, exc))
        if not math.isclose(density, reference_value, rel_tol=1e-12, abs_tol=0.0):
            return ('Model sample %d changes density when sampler-written interaction '
                    'parameters are replaced by template parameters; the reference '
                    'cannot reconstruct these values.' % i)
    return None


def _prepare_reference_template(model, template):
    """Allocate reference output storage without modifying supplied kinematics."""
    record = _dataclasses.InteractionRecord(template)
    count = len(record.signature.secondary_types)
    if len(record.secondary_masses) != count:
        record.secondary_masses = _validation.model_secondary_masses(
            model, record.signature.secondary_types)
    if len(record.secondary_momenta) != count:
        record.secondary_momenta = [[0.0] * 4 for _ in range(count)]
    if len(record.secondary_helicities) != count:
        record.secondary_helicities = [0.0] * count
    if len(record.secondary_times) != count:
        record.secondary_times = [record.interaction_time] * count
    return record


def _check_closure_model(model, template, n, seed, tol_sigma):
    report = ClosureReport(checks={'sampling': 'incomplete',
                                  'normalization': 'incomplete', 'shape': 'incomplete'})
    random = _utilities.SIREN_random(int(seed) & 0x7fffffff)
    drawn = []
    densities = []
    try:
        for i in range(n):
            source = _dataclasses.InteractionRecord(template)
            csdr = _dataclasses.CrossSectionDistributionRecord(source)
            model.SampleFinalState(csdr, random)
            out = _dataclasses.InteractionRecord(template)
            csdr.finalize(out)
            density = _density(model, out, 'model sample %d' % i)
            if density == 0:
                raise ClosureError('model sample %d has zero physical density' % i)
            drawn.append(out)
            densities.append(density)
    except ClosureError as exc:
        report.checks['sampling'] = 'failed'
        report.worst_region = str(exc)
        return report
    report.checks['sampling'] = 'passed'

    frame = _coordinate_frame(model)
    index = _coordinate_secondary_index(model)
    try:
        boost = (_cm_boost(template) if frame == 'cm'
                 else _boost_to_rest(template.primary_momentum, template.primary_mass))
        actual = [_angles(r, frame, boost, index) for r in drawn]
    except ClosureError as exc:
        report.checks['sampling'] = 'failed'
        report.worst_region = str(exc)
        return report
    resolved = [point for point in actual if point is not None]
    if resolved:
        report.notes.append('Measured costheta_secondary%d_%s range [%.4g, %.4g].'
                            % (index, frame, min(c[0] for c in resolved),
                               max(c[0] for c in resolved)))

    channel, reference_density = _reference_channel(model, template.signature)
    if channel is None:
        report.notes.append('No independent reference for %r / %r; normalization '
                            'and joint shape were not checked.'
                            % (model.Topology(), model.Measure()))
        return report
    context_issue = _reference_context_issue(model, template, drawn, densities)
    if context_issue:
        report.notes.append(context_issue)
        report.notes.append('Normalization and shape need a density computed from '
                            'kinematics and configured template parameters.')
        return report
    reference_template = _prepare_reference_template(model, template)
    reference = []
    weights = []
    random = _utilities.SIREN_random((int(seed) & 0x7fffffff) ^ 0x5a5a5a5a)
    try:
        for i in range(n):
            out = _dataclasses.InteractionRecord(reference_template)
            channel.Sample(random, None, out)
            weights.append(_density(model, out, 'reference sample %d' % i)
                           / reference_density)
            reference.append(out)
    except ClosureError as exc:
        report.checks['normalization'] = 'failed'
        report.worst_region = str(exc)
        return report

    weights = np.asarray(weights)
    mean = float(weights.mean())
    stderr = float(weights.std(ddof=1) / math.sqrt(n))
    report.normalization = (mean, stderr)
    if not math.isfinite(mean) or not math.isfinite(stderr) or mean <= 0:
        report.checks['normalization'] = 'failed'
        report.notes.append('The reference density integral is invalid or zero.')
        return report
    report.checks['normalization'] = (
        'passed' if abs(mean - 1.0) <= max(tol_sigma * stderr, 1e-9) else 'failed')

    try:
        expected = [_angles(r, frame, boost, index) for r in reference]
    except ClosureError as exc:
        report.checks['shape'] = 'failed'
        report.worst_region = str(exc)
        return report
    if any(c is None for c in actual + expected):
        report.notes.append('An angular coordinate is unresolved; shape was not checked.')
        return report
    actual, expected = np.asarray(actual), np.asarray(expected)
    _compare_shape(report, actual, expected, weights, frame, index, tol_sigma)
    if report.checks['shape'] == 'failed' and np.ptp(weights) <= 1e-10 * mean:
        report.frame_check = _frame_check(model, template, drawn)
    return report


def _compare_shape(report, actual, reference, weights, frame, index, tol_sigma):
    """Compare joint angular bins and moments using independent reference weights."""
    n = len(actual)
    side = max(2, min(8, int(math.sqrt(n / 20))))
    bounds = [[-1, 1], [-math.pi, math.pi]]
    observed, edges = np.histogramdd(actual, bins=(side, side), range=bounds)
    target, _ = np.histogramdd(reference, bins=(side, side), range=bounds, weights=weights)
    squared, _ = np.histogramdd(reference, bins=(side, side), range=bounds, weights=weights**2)
    ref_count, _ = np.histogramdd(reference, bins=(side, side), range=bounds)
    probability = target / weights.sum()
    w2 = float(np.sum(weights**2))
    ref_var = (squared * (1-probability)**2 + (w2-squared) * probability**2)
    ref_var /= n * (n-1) * float(weights.mean())**2
    observed_probability = observed / n
    # Each independent sample estimates its own variance. Using the reference
    # probability twice understates the error when a reference bin fluctuates
    # low. With constant weights this is symmetric under exchanging samples.
    variance = observed_probability * (1-observed_probability) / (n-1) + ref_var
    deviation = observed_probability - probability
    z = np.divide(deviation, np.sqrt(variance), out=np.zeros_like(deviation), where=variance > 0)
    z[(variance == 0) & (deviation != 0)] = float('inf')
    adequate = (ref_count >= 5) & ((n * probability >= 5) | ((target == 0) & (observed == 0)))
    # A bin absent from the finite reference sample is unmeasured, not evidence
    # that its physical probability is zero.
    failed = bool(np.any((np.abs(z) > tol_sigma) & adequate))
    report.checks['shape'] = 'failed' if failed else ('passed' if adequate.all() else 'incomplete')
    if not adequate.all():
        report.notes.append('Some angular bins have insufficient reference coverage; '
                            'increase samples to resolve them.')
    report.notes.append('Joint angular shape tested in %d x %d (cos(theta), phi) bins.'
                        % (side, side))
    ratios = np.divide(observed/n, probability, out=np.zeros_like(probability), where=probability > 0)
    mean_ratio = float(np.sum(ratios * probability))
    error = math.sqrt(float(np.sum(probability * (ratios-mean_ratio)**2)) / (side*side))
    report.flatness = (mean_ratio, error)
    if adequate.any():
        i, j = np.unravel_index(np.argmax(np.where(adequate, np.abs(z), -1)), z.shape)
        report.worst_region = ('secondary%d_%s: cos(theta) in [%.3g, %.3g], '
                               'phi in [%.3g, %.3g], z=%.2f'
                               % (index, frame, edges[0][i], edges[0][i+1],
                                  edges[1][j], edges[1][j+1], z[i, j]))
    if adequate.all():
        # Broad correlations can move many bins without making any individual
        # bin a four-sigma outlier. Test the full probability vector with the
        # covariance of both independent samples, including correlations from
        # normalization of the weighted reference.
        try:
            report.joint_shape = _joint_shape_test(
                observed_probability.ravel(), probability.ravel(), squared.ravel(),
                n, float(weights.mean()))
        except np.linalg.LinAlgError:
            if report.checks['shape'] != 'failed':
                report.checks['shape'] = 'incomplete'
            report.notes.append('Joint angular covariance could not be resolved.')
        else:
            if report.joint_shape[2] < math.erfc(tol_sigma / math.sqrt(2)):
                report.checks['shape'] = 'failed'
    observed_moments = [actual[:, 0], np.cos(actual[:, 1]), np.sin(actual[:, 1])]
    reference_moments = [reference[:, 0], np.cos(reference[:, 1]), np.sin(reference[:, 1])]
    for name, obs, ref in zip(('costheta', 'cosphi', 'sinphi'), observed_moments, reference_moments):
        ref_mean = float(np.average(ref, weights=weights))
        ref_variance = float(np.sum((weights * (ref-ref_mean))**2)
                             / (n*(n-1)*weights.mean()**2))
        variance = float(obs.var(ddof=1)/n) + ref_variance
        difference = float(obs.mean()) - ref_mean
        z_value = difference/math.sqrt(variance) if variance > 0 else (0.0 if difference == 0 else float('inf'))
        report.moment_z['%s_secondary%d_%s' % (name, index, frame)] = z_value
        if n >= 200 and (not math.isfinite(z_value) or abs(z_value) > tol_sigma):
            report.checks['shape'] = 'failed'


def _joint_shape_test(observed, reference, squared_weights, n, mean_weight):
    """Wald test for two independent, normalized histogram estimates.

    The model contributes the multinomial sample-mean covariance. For the
    reference, the influence of draw j is w_j * (one_hot(bin_j) - reference)
    divided by the mean weight; summing its outer products gives the second
    term below. This retains bin correlations and importance-weight variance.
    Constant reference weights recover an exchange-symmetric two-sample test.
    """
    active = (observed + reference) > 0
    observed, reference = observed[active], reference[active]
    squared_weights = squared_weights[active]
    # Empty-probability bins add no degrees of freedom. The remaining vector
    # sums to one, so omit its last component to remove the redundant equation.
    degrees = len(observed) - 1
    if degrees <= 0:
        return (0.0, 0, 1.0)
    covariance = (np.diag(observed) - np.outer(observed, observed)) / (n-1)
    covariance += (
        np.diag(squared_weights)
        - np.outer(squared_weights, reference)
        - np.outer(reference, squared_weights)
        + squared_weights.sum() * np.outer(reference, reference)
    ) / (n * (n-1) * mean_weight**2)
    deviation = (observed - reference)[:-1]
    statistic = float(deviation @ np.linalg.solve(covariance[:-1, :-1], deviation))
    if not math.isfinite(statistic) or statistic < 0:
        raise np.linalg.LinAlgError('Invalid joint angular statistic')
    return (statistic, degrees, float(gammaincc(degrees / 2, statistic / 2)))


def _angles(record, frame, boost, index):
    if index >= len(record.secondary_momenta):
        return None
    momentum = list(record.secondary_momenta[index])
    if not all(math.isfinite(x) for x in momentum):
        raise ClosureError('Sampled momentum is not finite')
    if frame != 'lab':
        momentum = boost(momentum)
    _, x, y, z = momentum
    magnitude = math.hypot(x, y, z)
    if not math.isfinite(magnitude):
        raise ClosureError('Boosted momentum is not finite')
    if magnitude == 0:
        return None
    return (z/magnitude, math.atan2(y, x))


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
    if rec.signature.target_type != _dataclasses.ParticleType.Decay:
        rec.target_mass = _mass_or(0.001, rec.signature.target_type)
    return rec, signature


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
    if (model.Topology() == _injection.PhaseSpaceTopology.Decay2Body
            and n_finals == 2
            and measure.type == _injection.PhaseSpaceMeasureType.SolidAngleRest):
        return _injection.Isotropic2BodyChannel(0), 1.0 / (4.0 * math.pi)
    return None, None


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


def _chi_square_uniform(values, n_bins=20):
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


def _frame_check(model, template, samples):
    """Flag a declared SolidAngle*-type measure that fits the other frame better.

    Reuses the model samples without drawing additional events.
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
    for out_rec in samples:
        if not out_rec.secondary_momenta:
            continue
        p4_lab = out_rec.secondary_momenta[0]
        lab_cos.append(_costheta(p4_lab))
        p4_rest = boost(p4_lab)
        rest_cos.append(_costheta(p4_rest))

    if len(lab_cos) < 50:
        return None

    chi2_lab = _chi_square_uniform(lab_cos)
    chi2_rest = _chi_square_uniform(rest_cos)

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
    """Read declared lab angles in the lab, scattering in CM, and decays at rest."""
    measure = model.Measure().type
    if measure == _injection.PhaseSpaceMeasureType.SolidAngleLab:
        return "lab"
    if model.Topology() in (_injection.PhaseSpaceTopology.Scatter2to2,
                            _injection.PhaseSpaceTopology.Scatter2to3):
        return "cm"
    if measure in (_injection.PhaseSpaceMeasureType.SolidAngleRest,
                   _injection.PhaseSpaceMeasureType.Recursive2Body):
        return "rest"
    return "lab"


def _cm_boost(template):
    """Boost by the total primary plus stationary-target four-momentum."""
    momentum = list(template.primary_momentum)
    momentum[0] += template.target_mass
    mass_squared = momentum[0]**2 - sum(p*p for p in momentum[1:])
    if mass_squared <= 0 or not math.isfinite(mass_squared):
        raise ClosureError("The collision has no timelike centre-of-mass frame")
    return _boost_to_rest(momentum, math.sqrt(mass_squared))


def _coordinate_secondary_index(model):
    """Use the Recursive2Body spectator, otherwise the first secondary."""
    measure = model.Measure()
    if measure.type == _injection.PhaseSpaceMeasureType.Recursive2Body:
        return int(measure.spectator)
    return 0
