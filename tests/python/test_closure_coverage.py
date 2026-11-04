"""Independent-reference closure, incomplete coverage, and changing physics."""

import math

import numpy as np
import pytest
import siren
from siren import closure
from siren.errors import ClosureError


def template(model, boost=1.0):
    record = siren.dataclasses.InteractionRecord()
    record.signature = model.GetPossibleSignatures()[0]
    record.primary_mass = 0.02
    record.primary_momentum = [0.02 * boost, 0, 0, 0.02 * math.sqrt(boost**2 - 1)]
    n = len(record.signature.secondary_types)
    record.secondary_masses = [0.0] * n
    record.secondary_momenta = [[0.0] * 4 for _ in range(n)]
    record.secondary_helicities = [0.0] * n
    return record


def rest_angles(momentum, primary):
    """Independent longitudinal Lorentz boost for these test kinematics."""
    beta = primary[3] / primary[0]
    gamma = 1 / math.sqrt(1 - beta**2)
    energy, x, y, z = momentum
    z = gamma * (z - beta * energy)
    return z / math.hypot(x, y, z), math.atan2(y, x)


class AngularDecay(siren.DecayModel):
    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleRest()
    skew = 0.0
    correlation = 0.0
    scale = 1.0

    def total_width(self):
        return 1.0

    def differential_width(self, record):
        cost, phi = rest_angles(record.secondary_momenta[0], record.primary_momentum)
        return self.scale * (1 + self.skew * cost + self.correlation * cost * math.cos(phi)) / (4 * math.pi)

    def sample(self, record, random):
        bound = 1 + abs(self.skew) + abs(self.correlation)
        for _ in range(10000):
            self.sample_isotropic(record, random)
            momentum = record.get_secondary_particle_records()[0].four_momentum
            cost, phi = rest_angles(momentum, record.record.primary_momentum)
            if random.Uniform(0, bound) <= 1 + self.skew * cost + self.correlation * cost * math.cos(phi):
                return
        raise AssertionError("test sampler exhausted")


@pytest.mark.parametrize("skew,correlation", [(0.8, 0), (0, 0.8)])
def test_matching_anisotropic_density_and_sampler_pass(skew, correlation):
    model = AngularDecay()
    model.skew, model.correlation = skew, correlation
    report = siren.check_closure(model, record=template(model, boost=3), samples=8000, seed=19)
    assert report.ok and report.complete, str(report)
    assert report.normalization[0] == pytest.approx(1.0, abs=0.025)
    assert set(report.moment_z) == {name + "_secondary0_rest" for name in ("costheta", "cosphi", "sinphi")}


def test_joint_angle_mismatch_with_unchanged_marginal_means_fails():
    class WrongSampler(AngularDecay):
        correlation = 0.8

        def sample(self, record, random):
            self.sample_isotropic(record, random)

    model = WrongSampler()
    report = siren.check_closure(model, record=template(model), samples=10000, seed=11)
    assert report.checks['normalization'] == 'passed'
    assert report.checks['shape'] == 'failed'
    assert 'phi in' in report.worst_region


@pytest.mark.parametrize('seed', range(60))
def test_joint_correlation_mismatch_is_detected_at_default_sample_count(seed):
    class WrongSampler(AngularDecay):
        correlation = 0.8

        def sample(self, record, random):
            self.sample_isotropic(record, random)

    model = WrongSampler()
    # Keep the public default: the per-bin-only test missed 25/60 seeds here.
    report = siren.check_closure(model, record=template(model), seed=seed)
    assert report.checks['shape'] == 'failed', str(report)
    assert report.joint_shape is not None
    statistic, degrees, pvalue = report.joint_shape
    assert statistic > 0 and degrees == 63
    assert pvalue < math.erfc(4.0 / math.sqrt(2))
    assert 'joint angular test:' in str(report)


def test_missing_sampling_support_is_not_removed_by_empirical_rebinning():
    class Hemisphere(AngularDecay):
        def sample(self, record, random):
            while True:
                self.sample_isotropic(record, random)
                if record.get_secondary_particle_records()[0].four_momentum[3] > 0:
                    return

    model = Hemisphere()
    report = siren.check_closure(model, record=template(model), samples=2500, seed=4)
    assert report.checks['normalization'] == 'passed'
    assert report.checks['shape'] == 'failed'


@pytest.mark.parametrize('density', [-1.0, float('nan'), float('inf')])
def test_invalid_sample_density_is_a_failure(density):
    class Invalid(AngularDecay):
        def differential_width(self, record):
            return density

    model = Invalid()
    report = siren.check_closure(model, record=template(model), samples=200)
    assert report.status == 'failed'
    assert 'model sample 0' in report.worst_region
    with pytest.raises(ClosureError):
        report.raise_if_failed()


def test_invalid_reference_density_cannot_be_converted_to_zero():
    class InvalidReference(AngularDecay):
        def sample(self, record, random):
            while True:
                self.sample_isotropic(record, random)
                if record.get_secondary_particle_records()[0].four_momentum[3] > 0:
                    return

        def differential_width(self, record):
            return 1 / (2 * math.pi) if record.secondary_momenta[0][3] > 0 else float('nan')

    model = InvalidReference()
    report = siren.check_closure(model, record=template(model), samples=200)
    assert report.checks['normalization'] == 'failed'
    assert 'reference sample' in report.worst_region


def test_current_model_state_is_always_rechecked():
    model = AngularDecay()
    record = template(model)
    assert siren.check_closure(model, record=record, samples=2500, seed=2).ok
    model.scale = 2.0
    report = siren.check_closure(model, record=record, samples=2500, seed=2)
    assert report.checks['normalization'] == 'failed'
    assert report.normalization[0] == pytest.approx(2.0)


def test_sparse_reference_coverage_is_incomplete():
    model = AngularDecay()
    report = siren.check_closure(model, record=template(model), samples=3, seed=0)
    assert report.status == 'incomplete', str(report)
    assert not report.ok and not report.complete
    assert report.worst_region == ''
    assert report.joint_shape is None
    with pytest.raises(ClosureError):
        report.raise_if_failed()


def test_unresolved_angular_coordinate_is_incomplete():
    class AtRest(AngularDecay):
        def sample(self, record, random):
            for secondary in record.get_secondary_particle_records():
                secondary.four_momentum = [0.01, 0.0, 0.0, 0.0]

        def differential_width(self, record):
            return 1 / (4 * math.pi)

    model = AtRest()
    report = siren.check_closure(model, record=template(model), samples=200)
    assert report.status == 'incomplete'
    assert report.checks['normalization'] == 'passed'
    assert report.checks['shape'] == 'incomplete'


@pytest.mark.parametrize('scatter', [False, True])
def test_three_body_coverage_is_incomplete(scatter):
    class ThreeBody(AngularDecay):
        daughters = ("NuLight", "Gamma", "NuMu")
        measure = siren.Measure.Recursive2Body(2, 0, 1)

        def Topology(self):
            return siren.Topology.Scatter2to3 if scatter else siren.Topology.Decay3Body

        def sample(self, record, random):
            # Three massless particles, equally spaced in the rest-frame plane.
            for i, secondary in enumerate(record.get_secondary_particle_records()):
                phi = 2 * math.pi * i / 3
                secondary.four_momentum = [0.02/3, 0.02/3*math.cos(phi), 0.02/3*math.sin(phi), 0]

        def differential_width(self, record):
            return 1.0

    model = ThreeBody()
    report = siren.check_closure(model, record=template(model), samples=200)
    assert report.status == 'incomplete'
    assert report.checks['normalization'] == report.checks['shape'] == 'incomplete'
    assert any('secondary2_' + ('cm' if scatter else 'rest') in note for note in report.notes)


def test_mixture_configuration_success_is_not_density_certification():
    class Unprobed:
        def ValidateChannelDensities(self, *args):
            raise AssertionError('no detector context supplied')

    bare = siren.check_closure(Unprobed())
    assert bare.status == 'incomplete'
    mixture = siren.check_closure(siren.channels.isotropic(0) + siren.channels.isotropic(1))
    assert mixture.checks['configuration'] == 'passed'
    assert mixture.status == 'incomplete' and not mixture.ok


def test_cm_frame_and_explicit_lab_frame():
    record = siren.dataclasses.InteractionRecord()
    record.primary_momentum = [2.0, 0.0, 0.0, 2.0]
    record.primary_mass = 0.0
    record.target_mass = 1.0
    transform = closure._cm_boost(record)
    total = transform([3.0, 0.0, 0.0, 2.0])
    assert total == pytest.approx([math.sqrt(5), 0, 0, 0], abs=1e-14)

    class Scatter(AngularDecay):
        def Topology(self):
            return siren.Topology.Scatter2to2

    assert closure._coordinate_frame(Scatter()) == 'cm'
    model = Scatter()
    model.measure = siren.Measure.SolidAngleLab()
    assert closure._coordinate_frame(model) == 'lab'


def test_explicit_record_is_preserved():
    model = AngularDecay()
    record = template(model, boost=2)
    before = list(record.primary_momentum), [list(p) for p in record.secondary_momenta]
    siren.check_closure(model, record=record, samples=1000, seed=9)
    assert before == (list(record.primary_momentum), [list(p) for p in record.secondary_momenta])
    with pytest.raises(ValueError, match='cannot be combined'):
        siren.check_closure(model, record=record, primary_energy=1.0)


def test_lab_sampler_for_rest_density_fails_with_frame_diagnostic():
    class LabSampler(AngularDecay):
        def sample(self, record, random):
            cost = random.Uniform(-1, 1)
            phi = random.Uniform(-math.pi, math.pi)
            sine = math.sqrt(1-cost*cost)
            for i, secondary in enumerate(record.get_secondary_particle_records()):
                sign = 1 if i == 0 else -1
                secondary.four_momentum = [0.01, sign*0.01*sine*math.cos(phi), sign*0.01*sine*math.sin(phi), sign*0.01*cost]

    model = LabSampler()
    report = siren.check_closure(model, record=template(model, boost=3), samples=2500, seed=3)
    assert report.checks['shape'] == 'failed'
    assert report.frame_check and 'LAB' in report.frame_check


def test_record_copy_preserves_helicities_times_and_parameters():
    model = AngularDecay()
    source = template(model, boost=2)
    source.primary_helicity = -0.5
    source.target_helicity = 0.5
    source.primary_initial_time = 3.0
    source.interaction_time = 7.0
    source.secondary_times = [8.0, 9.0]
    source.secondary_helicities = [-0.5, 0.5]
    source.interaction_parameters = {"polarization": 0.75}
    clone = siren.dataclasses.InteractionRecord(source)
    for field in ("signature", "primary_helicity", "target_helicity",
                  "primary_initial_time", "interaction_time", "secondary_times",
                  "secondary_helicities", "interaction_parameters"):
        assert getattr(clone, field) == getattr(source, field)
    clone.signature.secondary_types = [siren.particles.Gamma, siren.particles.Gamma]
    clone.interaction_parameters = {"polarization": 0.0}
    assert source.interaction_parameters == {"polarization": 0.75}
    assert list(source.signature.secondary_types) != list(clone.signature.secondary_types)


def test_report_does_not_claim_unmeasured_density_variables():
    class ExtraVariables(AngularDecay):
        def density_variables(self):
            return ["s_pair", "cos_theta_sub", "phi_sub"]

    model = ExtraVariables()
    report = siren.check_closure(model, record=template(model), samples=2000, seed=7)
    assert report.ok, str(report)
    assert set(report.moment_z) == {"costheta_secondary0_rest",
                                    "cosphi_secondary0_rest", "sinphi_secondary0_rest"}


@pytest.mark.parametrize('correlation', [False, True])
def test_sampler_cached_density_cannot_certify_an_isotropic_mismatch(correlation):
    class CachedAngle(AngularDecay):
        def sample(self, record, random):
            self.sample_isotropic(record, random)
            momentum = record.get_secondary_particle_records()[0].four_momentum
            cost, _ = rest_angles(momentum, record.record.primary_momentum)
            record.interaction_parameters = {'costheta': cost}

        def differential_width(self, record):
            # The plain cos(theta) case is the exact regression witness:
            # the parent rejected it at z=-41.17; 4371232 falsely passed it.
            cost = record.interaction_parameters.get('costheta', 0.0)
            _, phi = rest_angles(record.secondary_momenta[0], record.primary_momentum)
            return (1 + .8 * cost * (math.cos(phi) if correlation else 1)) / (4*math.pi)

    report = siren.check_closure(CachedAngle(), samples=8000, seed=19)
    assert report.status == 'incomplete', str(report)
    assert report.checks['sampling'] == 'passed'
    assert report.checks['normalization'] == report.checks['shape'] == 'incomplete'
    assert any('sampler-written interaction parameters' in note for note in report.notes)
    with pytest.raises(ClosureError):
        report.raise_if_failed()


def test_required_sampler_parameter_is_incomplete_even_without_a_fallback():
    class CachedAngle(AngularDecay):
        def sample(self, record, random):
            self.sample_isotropic(record, random)
            momentum = record.get_secondary_particle_records()[0].four_momentum
            cost, _ = rest_angles(momentum, record.record.primary_momentum)
            record.interaction_parameters = {'costheta': cost}

        def differential_width(self, record):
            return (1 + .8 * record.interaction_parameters['costheta']) / (4*math.pi)

    report = siren.check_closure(CachedAngle(), samples=200, seed=19)
    assert report.status == 'incomplete', str(report)
    assert any('KeyError' in note for note in report.notes)


@pytest.mark.parametrize('matching_sampler', [False, True])
def test_cached_density_with_a_kinematic_fallback_can_be_checked(matching_sampler):
    class CachedWithFallback(AngularDecay):
        skew = 0.8

        def sample(self, record, random):
            if matching_sampler:
                super().sample(record, random)
            else:
                self.sample_isotropic(record, random)
            momentum = record.get_secondary_particle_records()[0].four_momentum
            cost, _ = rest_angles(momentum, record.record.primary_momentum)
            record.interaction_parameters = {'costheta': cost}

        def differential_width(self, record):
            parameters = record.interaction_parameters
            cost = (parameters['costheta'] if 'costheta' in parameters else
                    rest_angles(record.secondary_momenta[0], record.primary_momentum)[0])
            return (1 + self.skew * cost) / (4*math.pi)

    report = siren.check_closure(CachedWithFallback(), samples=8000, seed=19)
    assert report.complete, str(report)
    assert report.ok == matching_sampler, str(report)
    assert report.checks['shape'] == ('passed' if matching_sampler else 'failed')


def test_unused_sampler_parameters_and_configured_density_parameters_pass():
    class Bookkeeping(AngularDecay):
        def sample(self, record, random):
            super().sample(record, random)
            parameters = dict(record.interaction_parameters)
            parameters['unused'] = random.Uniform(0, 1)
            record.interaction_parameters = parameters

        def differential_width(self, record):
            return record.interaction_parameters['scale'] * super().differential_width(record)

    model = Bookkeeping()
    source = template(model)
    source.interaction_parameters = {'scale': 1.0}
    report = siren.check_closure(model, record=source, samples=2000, seed=7)
    assert report.ok, str(report)
    assert source.interaction_parameters == {'scale': 1.0}


@pytest.mark.parametrize('supplied_masses', [False, True])
def test_initial_record_needs_no_secondary_output_storage(supplied_masses):
    class Massive(AngularDecay):
        def SecondaryMasses(self, types):
            return [0.003, 0.004]

        def differential_width(self, record):
            expected = [0.005, 0.006] if supplied_masses else [0.003, 0.004]
            assert list(record.secondary_masses) == expected
            return 1 / (4*math.pi)

    model = Massive()
    source = siren.dataclasses.InteractionRecord()
    source.signature = model.GetPossibleSignatures()[0]
    source.primary_mass = 0.02
    source.primary_momentum = [0.02, 0, 0, 0]
    if supplied_masses:
        source.secondary_masses = [0.005, 0.006]
    report = siren.check_closure(model, record=source, samples=2000, seed=7)
    assert report.ok, str(report)
    assert list(source.secondary_momenta) == []
    assert list(source.secondary_masses) == ([0.005, 0.006] if supplied_masses else [])


def test_uniform_shape_comparison_is_symmetric_under_sample_exchange(monkeypatch):
    compare = closure._compare_shape
    captured = {}

    def capture(report, actual, reference, weights, frame, index, tol_sigma):
        captured.update(actual=actual, reference=reference, weights=weights)
        compare(report, actual, reference, weights, frame, index, tol_sigma)

    monkeypatch.setattr(closure, '_compare_shape', capture)
    from test_closure_gauge import GoodIsoDecay
    forward = siren.check_closure(GoodIsoDecay(), samples=2000, seed=6)
    backward = siren.ClosureReport(checks={'shape': 'incomplete'})
    compare(backward, captured['reference'], captured['actual'], captured['weights'],
            'rest', 0, 4.0)
    assert forward.ok, str(forward)
    assert backward.ok, str(backward)
    forward_region, forward_z = forward.worst_region.rsplit('z=', 1)
    backward_region, backward_z = backward.worst_region.rsplit('z=', 1)
    assert forward_region == backward_region
    assert float(forward_z) == -float(backward_z)
    assert forward.joint_shape == pytest.approx(backward.joint_shape)
    for name, value in forward.moment_z.items():
        assert value == pytest.approx(-backward.moment_z[name])


@pytest.mark.parametrize('correlation', [0.0, 0.8])
def test_shape_statistics_have_a_calibrated_false_rejection_rate(correlation):
    # Independent draws isolate statistical calibration from the native sampler.
    # At four sigma across 68 checks, 2% is a generous rejection bound; the old
    # reference-only denominator exceeded it for the uniform RNG sequence.
    random = np.random.default_rng(20260913)
    failures = 0
    for _ in range(500):
        actual = random.uniform([-1, -math.pi], [1, math.pi], size=(2000, 2))
        if correlation:
            # Rejection-sample the same angular law as the weighted reference.
            accepted = []
            count = 0
            while count < 2000:
                weight = 1 + correlation * actual[:, 0] * np.cos(actual[:, 1])
                batch = actual[random.uniform(0, 1 + correlation, len(actual)) < weight]
                accepted.append(batch)
                count += len(batch)
                actual = random.uniform([-1, -math.pi], [1, math.pi], size=(2000, 2))
            actual = np.concatenate(accepted)[:2000]
        reference = random.uniform([-1, -math.pi], [1, math.pi], size=(2000, 2))
        weights = 1 + correlation * reference[:, 0] * np.cos(reference[:, 1])
        report = siren.ClosureReport(checks={'shape': 'incomplete'})
        closure._compare_shape(report, actual, reference, weights, 'rest', 0, 4.0)
        if not report.complete:
            # A low-density bin can lack reference coverage even in a correct
            # model. Keep it incomplete and count it against this bound too.
            assert report.joint_shape is None
            assert any('insufficient reference coverage' in note for note in report.notes)
        failures += not report.ok
    assert failures <= 10
