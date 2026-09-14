"""Independent kinematic and normalization checks for the public utilities."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

import siren
from siren import three_body as tb


class FlatDecay:
    m_M, m_l, m_phi = 1.0, 0.2, 0.15
    E_nu_max = (m_M**2 - (m_l + m_phi)**2) / (2 * m_M)

    def _E_phi_limits(self, energy):
        return tb.dalitz_band(self.m_M, self.m_l, self.m_phi, energy)

    def _matel_sq(self, energy, other_energy):
        return 1.0


def _record(physics, momenta, parent):
    record = siren.dataclasses.InteractionRecord()
    signature = record.signature
    signature.primary_type = siren.particles.N4
    signature.secondary_types = [siren.particles.NuMu, siren.particles.Electron,
                                 siren.particles.Gamma]
    record.signature = signature
    record.primary_mass = physics.m_M
    record.primary_momentum = parent
    record.secondary_momenta = momenta
    record.secondary_masses = [0, physics.m_l, physics.m_phi]
    return record


def _reference_width(physics):
    """Unit-|M|^2 width from an independent invariant-mass phase-space integral."""
    M, a, b = physics.m_M, physics.m_l, physics.m_phi

    def integrand(s):
        lam = (M*M - (a + math.sqrt(s))**2) * (M*M - (a - math.sqrt(s))**2)
        return math.sqrt(max(lam, 0)) * (s - b*b) / s

    return quad(integrand, b*b, (M-a)**2, epsabs=1e-11)[0] / (256 * math.pi**3 * M**3)


def test_width_matches_independent_invariant_phase_space_integral():
    physics = FlatDecay()
    assert tb.dalitz_width(physics) == pytest.approx(_reference_width(physics), rel=2e-5, abs=0)


@pytest.mark.parametrize('scale', [1e6, 1.0, 1e-6, 1e-12, 1e-19])
def test_width_accuracy_does_not_depend_on_amplitude_normalization(scale):
    """Small couplings must not relax quadrature's relative accuracy.

    Scale the independent reference analytically. Disable pytest's default
    absolute tolerance too: it would accept even zero for tiny widths.
    """
    class Scaled(FlatDecay):
        def _matel_sq(self, energy, other_energy):
            return scale

    expected = _reference_width(FlatDecay()) * scale
    assert tb.dalitz_width(Scaled()) == pytest.approx(expected, rel=2e-5, abs=0)


@pytest.mark.parametrize('scale', [1.0, 1e-20])
def test_peaked_width_matches_resolved_one_dimensional_reference(scale):
    class Peaked(FlatDecay):
        def _matel_sq(self, energy, other_energy):
            return scale / (1 + ((energy - 0.137) / 0.003)**2)

    physics = Peaked()

    # The amplitude is constant in E_phi, so integrate that coordinate
    # analytically and resolve the peak explicitly in the reference integral.
    def marginal(energy):
        low, high = physics._E_phi_limits(energy)
        return (high - low) / (1 + ((energy - 0.137) / 0.003)**2)

    area = quad(marginal, 0, physics.E_nu_max, points=[0.137],
                epsabs=0, epsrel=1e-10)[0]
    expected = scale * area / (64 * math.pi**3 * physics.m_M)
    assert tb.dalitz_width(physics) == pytest.approx(expected, rel=2e-6, abs=0)


def test_energy_sampler_matches_the_nonuniform_band_marginal():
    physics = FlatDecay()
    rng = siren.utilities.SIREN_random(179)
    bound = tb.find_max_weight(physics)
    energies = np.array([tb.sample_energies(physics, bound, rng) for _ in range(4000)])
    area = quad(lambda e: np.diff(physics._E_phi_limits(e))[0], 0, physics.E_nu_max)[0]
    for edge in np.linspace(0.1, 0.9, 9) * physics.E_nu_max:
        expected = quad(lambda e: np.diff(physics._E_phi_limits(e))[0], 0, edge)[0] / area
        observed = np.mean(energies[:, 0] < edge)
        assert abs(observed - expected) < 5 * math.sqrt(expected*(1-expected)/len(energies))


@pytest.mark.parametrize('momentum', [[0, 0, 0], [0.3, -0.5, 2.0], [10, 2, -3]])
def test_mass_shells_conservation_boost_and_density(momentum):
    physics = FlatDecay()
    rng = siren.utilities.SIREN_random(278)
    parent = np.array([math.sqrt(physics.m_M**2 + np.dot(momentum, momentum)), *momentum])
    width, bound = tb.dalitz_width(physics), tb.find_max_weight(physics)
    for _ in range(30):
        energies = tb.sample_energies(physics, bound, rng)
        rest = np.array(tb.build_rest_momenta(physics, *energies, rng))
        lab = np.array([tb.boost_four_vector(parent, p, physics.m_M) for p in rest])
        np.testing.assert_allclose(rest.sum(axis=0), [physics.m_M, 0, 0, 0], atol=1e-14)
        np.testing.assert_allclose(lab.sum(axis=0), parent, atol=1e-13)
        np.testing.assert_allclose(lab[:, 0]**2 - np.sum(lab[:, 1:]**2, axis=1),
                                   [0, physics.m_l**2, physics.m_phi**2], atol=1e-13)
        np.testing.assert_allclose([tb.boost_to_rest(parent, p) for p in lab], rest, atol=1e-13)
        stationary = _record(physics, rest, [physics.m_M, 0, 0, 0])
        moving = _record(physics, lab, parent)
        expected = tb.final_state_probability(physics, width, 14, 22, stationary)
        assert expected > 0
        assert tb.final_state_probability(physics, width, 14, 22, moving) == pytest.approx(expected, rel=2e-11, abs=0)
        # Density identifies the two pair daughters by PDG, not record order.
        signature = moving.signature
        signature.secondary_types = list(reversed(signature.secondary_types))
        moving.signature = signature
        moving.secondary_momenta = list(reversed(moving.secondary_momenta))
        moving.secondary_masses = list(reversed(moving.secondary_masses))
        assert tb.final_state_probability(physics, width, 14, 22, moving) == pytest.approx(expected, rel=2e-11, abs=0)


def test_rest_orientation_is_isotropic():
    physics = FlatDecay()
    rng = siren.utilities.SIREN_random(279)
    energy = physics.E_nu_max * 0.4
    other = sum(physics._E_phi_limits(energy)) / 2
    directions = np.array([tb.build_rest_momenta(physics, energy, other, rng)[0][1:]
                           / energy for _ in range(4000)])
    assert np.max(np.abs(directions.mean(axis=0))) < 0.04
    np.testing.assert_allclose(directions.T @ directions / len(directions), np.eye(3)/3, atol=0.025)


def test_muon_width_matches_the_analytic_mass_correction():
    beam = siren.resources.processes.BeamDecays
    r = (beam.M_E / beam.M_MU)**2
    f = 1 - 8*r + 8*r**3 - r**4 - 12*r*r*math.log(r)
    expected = beam.G_F**2 * beam.M_MU**5 * f / (192 * math.pi**3)
    assert beam.MuonThreeBodyDecay(13).total_width() == pytest.approx(expected, rel=1e-6, abs=0)


@pytest.mark.parametrize('charge', [-13, 13])
@pytest.mark.parametrize('polarization', [None, [0.2, 0.3, -0.7]])
def test_native_muon_sampling_preserves_mass_shells_and_density(charge, polarization):
    from test_authoring_bases import _template_record
    beam = siren.resources.processes.BeamDecays
    model = beam.MuonThreeBodyDecay(charge)
    rng = siren.utilities.SIREN_random(84)
    for _ in range(30):
        record = _template_record(model.GetPossibleSignatures()[0], energy=0.3, primary_mass=beam.M_MU)
        if polarization:
            record.interaction_parameters = dict(zip(['pol_x', 'pol_y', 'pol_z'], polarization))
        sample = siren.dataclasses.CrossSectionDistributionRecord(record)
        siren.interactions.Decay.SampleFinalState(model, sample, rng)
        sample.finalize(record)
        momenta = np.array(record.secondary_momenta)
        np.testing.assert_allclose(momenta.sum(axis=0), record.primary_momentum, atol=1e-14)
        np.testing.assert_allclose(momenta[:, 0]**2 - np.sum(momenta[:, 1:]**2, axis=1),
                                   [beam.M_E**2, 0, 0], atol=1e-14)
        assert model.FinalStateProbability(record) > 0


@pytest.mark.parametrize('bound', [0, -1, math.inf, math.nan])
def test_invalid_envelope_is_rejected(bound):
    with pytest.raises(ValueError, match='finite positive'):
        tb.sample_energies(FlatDecay(), bound, siren.utilities.SIREN_random(1))


class MidpointRandom:
    def Uniform(self, low, high):
        return (low + high) / 2


@pytest.mark.parametrize('weight', [-1, math.nan, math.inf])
def test_invalid_matrix_element_is_rejected(weight):
    physics = FlatDecay()
    physics._matel_sq = lambda *args: weight
    with pytest.raises(RuntimeError, match='invalid acceptance weight'):
        tb.sample_energies(physics, 1, MidpointRandom())


def test_envelope_overrun_is_rejected():
    with pytest.raises(RuntimeError, match='exceeds max_weight'):
        tb.sample_energies(FlatDecay(), 1e-10, MidpointRandom())


def test_exhaustion_does_not_return_an_unsampled_midpoint():
    physics = FlatDecay()
    physics._matel_sq = lambda *args: 0
    with pytest.raises(RuntimeError, match='exhausted 10000'):
        tb.sample_energies(physics, 1, MidpointRandom())


def test_zero_density_is_not_accepted_at_zero_uniform_draw():
    class ZeroRandom:
        def Uniform(self, low, high):
            return low

    physics = FlatDecay()
    physics._matel_sq = lambda *args: 0
    with pytest.raises(RuntimeError, match='exhausted 10000'):
        tb.sample_energies(physics, 1, ZeroRandom())


class UnusedRandom:
    def Uniform(self, low, high):
        pytest.fail('invalid inputs must be rejected before drawing randomness')


@pytest.mark.parametrize('masses', [
    (0, 0.2, 0.15), (-1, 0.2, 0.15),
    (1, -0.2, 0.15), (1, 0.2, -0.15),
    (1, 0.5, 0.5), (1, 0.7, 0.5), (1, 0, 0),
    (math.nan, 0.2, 0.15), (1, math.nan, 0.15), (1, 0.2, math.nan),
    (math.inf, 0.2, 0.15), (1, math.inf, 0.15), (1, 0.2, math.inf),
    (-math.inf, 0.2, 0.15),
])
@pytest.mark.parametrize('operation', [
    lambda d: tb.dalitz_band(d.m_M, d.m_l, d.m_phi, 0.2),
    lambda d: tb.dalitz_width(d),
    lambda d: tb.find_max_weight(d),
    lambda d: tb.sample_energies(d, 1, UnusedRandom()),
    lambda d: tb.build_rest_momenta(d, 0.2, 0.3, UnusedRandom()),
    lambda d: tb.final_state_probability(d, 0, 14, 22, None),
], ids=['band', 'width', 'envelope', 'energies', 'momenta', 'density'])
def test_unsupported_masses_fail_consistently(masses, operation):
    physics = FlatDecay()
    physics.m_M, physics.m_l, physics.m_phi = masses
    with pytest.raises(ValueError, match='masses'):
        operation(physics)


@pytest.mark.parametrize('energy', [math.nan, math.inf, -math.inf])
def test_band_rejects_nonfinite_energy(energy):
    with pytest.raises(ValueError, match='finite'):
        tb.dalitz_band(1, 0.2, 0.15, energy)


@pytest.mark.parametrize('energy', [-0.01, 0.5])
def test_band_returns_no_support_for_finite_outside_energy(energy):
    assert tb.dalitz_band(1, 0.2, 0.15, energy) == (None, None)


@pytest.mark.parametrize('energies', [
    (0.2, 0.7), (0.2, 0.1), (-0.1, 0.3), (0.5, 0.3),
    (math.nan, 0.3), (math.inf, 0.3), (-math.inf, 0.3),
    (0.2, math.nan), (0.2, math.inf), (0.2, -math.inf),
])
def test_momenta_reject_invalid_energies_before_sampling(energies):
    with pytest.raises(ValueError, match='energ'):
        tb.build_rest_momenta(FlatDecay(), *energies, UnusedRandom())


def _decimal_band(masses, fraction):
    """Resolve boundary energies independently at higher precision."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = 80
        M, a, b = map(Decimal.from_float, masses)
        energy = Decimal.from_float(fraction) * (M*M - (a+b)**2) / (2*M)
        energy = float(energy)
    return energy, *_decimal_band_at_energy(masses, energy)


def _decimal_band_at_energy(masses, energy):
    """Evaluate at the exact supplied float, not its unrounded precursor."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = 80
        M, a, b = map(Decimal.from_float, masses)
        limit = (M*M - (a+b)**2) / (2*M)
        energy = min(Decimal.from_float(energy), limit)
        s = M*M - 2*M*energy
        center = (M-energy) * (s+b*b-a*a) / (2*s)
        half_band = energy * max(Decimal(0), (s-(a+b)**2)*(s-(a-b)**2)).sqrt() / (2*s)
        return float(center-half_band), float(center+half_band)


@pytest.mark.parametrize('masses', [
    (1., 0.2, 0.15), (1., 0.2, 0.), (1., 0., 0.15),
    (0.13957039, 0.10565837, 0.017), (0.49368, 0.000511, 0.020),
    (1., 0.4, 0.6-1e-10), (1e-6, 2e-7, 1.5e-7), (1e6, 2e5, 1.5e5),
])
@pytest.mark.parametrize('fraction', [0., 1e-12, 0.4, 1-1e-12, 1.])
def test_physical_boundaries_preserve_mass_shells(masses, fraction):
    physics = FlatDecay()
    physics.m_M, physics.m_l, physics.m_phi = masses
    e_nu, low, high = _decimal_band(masses, fraction)
    band = tb.dalitz_band(*masses, e_nu)
    assert band[0] is not None
    # Away from the precision fallback, small invariant-mass roundoff is
    # amplified in a nearly collapsed band. Also check its mass shells.
    if fraction in (0., 0.4, 1.):
        np.testing.assert_allclose(band, [low, high], rtol=0, atol=2e-14*masses[0])
    for e_phi in (low, (low+high)/2, high, *band):
        rest = np.array(tb.build_rest_momenta(
            physics, e_nu, e_phi, siren.utilities.SIREN_random(219))) / masses[0]
        np.testing.assert_allclose(rest.sum(axis=0), [1, 0, 0, 0], rtol=0, atol=3e-15)
        np.testing.assert_allclose(
            rest[:, 0]**2 - np.sum(rest[:, 1:]**2, axis=1),
            [0, (masses[1]/masses[0])**2, (masses[2]/masses[0])**2],
            rtol=0, atol=3e-14)
        assert np.all(rest[:, 0] >= 0)


@pytest.mark.parametrize('fraction', [0., 0.4, 1.])
@pytest.mark.parametrize('direction', [-math.inf, math.inf])
def test_boundary_roundoff_is_tolerated_but_material_violations_raise(fraction, direction):
    physics = FlatDecay()
    energy, low, high = _decimal_band((physics.m_M, physics.m_l, physics.m_phi), fraction)
    for e_phi in [low, high]:
        rest = np.array(tb.build_rest_momenta(
            physics, np.nextafter(energy, direction), np.nextafter(e_phi, direction),
            siren.utilities.SIREN_random(32)))
        np.testing.assert_allclose(rest.sum(axis=0), [1, 0, 0, 0], atol=1e-15)
        np.testing.assert_allclose(rest[:, 0]**2-np.sum(rest[:, 1:]**2, axis=1),
                                   [0, 0.04, 0.0225], rtol=0, atol=3e-14)
    # At a collapsed endpoint the mass-shell violation is quadratic in
    # this displacement; choose a violation well above floating roundoff.
    outside = low - 1e-5 if direction < 0 else high + 1e-5
    with pytest.raises(ValueError, match='energ'):
        tb.build_rest_momenta(physics, energy, outside, UnusedRandom())


def test_original_endpoint_support_witness():
    energy = 0.49999999999999484
    expected = _decimal_band_at_energy((1., 1e-7, 0.), energy)
    actual = tb.dalitz_band(1., 1e-7, 0., energy)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-17)
    assert actual[1] > 0.015


@pytest.mark.parametrize('masses', [
    (1., 1e-7, 0.), (1., 0., 1e-7), (1., 1e-7, 2e-7), (1., .2, .15),
])
def test_adjacent_interior_energies_keep_resolved_band_width(masses):
    energy, _, _ = _decimal_band(masses, 1.)
    for _ in range(33):
        low, high = _decimal_band_at_energy(masses, energy)
        actual = tb.dalitz_band(*masses, energy)
        np.testing.assert_allclose(actual, [low, high], rtol=0,
                                   atol=4*math.ulp(masses[0]))
        if high-low > 4*math.ulp(masses[0]):
            assert actual[1] > actual[0]
        energy = math.nextafter(energy, 0.)


@pytest.mark.parametrize('dtype', [np.float32, np.float64, np.asarray])
@pytest.mark.parametrize('masses', [(1., .2, .15), (1., .4, .59999)])
def test_numpy_masses_use_the_same_precision_as_python_floats(dtype, masses):
    physics, promoted = FlatDecay(), FlatDecay()
    values = tuple(dtype(m) for m in masses)
    physics.m_M, physics.m_l, physics.m_phi = values
    exact = tuple(float(m) for m in values)
    promoted.m_M, promoted.m_l, promoted.m_phi = exact
    for fraction in (0., .4, 1.):
        energy, low, high = _decimal_band(exact, fraction)
        assert tb.dalitz_band(*values, energy) == tb.dalitz_band(*exact, energy)
        for e_phi in (low, (low+high)/2, high):
            actual = np.array(tb.build_rest_momenta(
                physics, energy, e_phi, siren.utilities.SIREN_random(219)))
            expected = np.array(tb.build_rest_momenta(
                promoted, energy, e_phi, siren.utilities.SIREN_random(219)))
            np.testing.assert_array_equal(actual, expected)
            np.testing.assert_allclose(actual[:, 0]**2-np.sum(actual[:, 1:]**2, axis=1),
                                       [0, exact[1]**2, exact[2]**2], rtol=0, atol=3e-14)


@pytest.mark.parametrize('dtype', [np.float32, np.float64, np.asarray])
def test_numpy_energy_scalars_are_promoted_before_arithmetic(dtype):
    energies = dtype(.2), dtype(.4)
    actual = tb.build_rest_momenta(FlatDecay(), *energies, siren.utilities.SIREN_random(219))
    expected = tb.build_rest_momenta(FlatDecay(), *map(float, energies),
                                     siren.utilities.SIREN_random(219))
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize('energy', [.4015659375, .43875*.914, .43875*.915])
def test_returned_band_edges_preserve_shells_at_phi_turning_point(energy):
    physics = FlatDecay()
    expected = _decimal_band_at_energy((1., .2, .15), energy)
    band = tb.dalitz_band(1., .2, .15, energy)
    np.testing.assert_allclose(band, expected, rtol=0, atol=2e-16)
    for e_phi in band:
        p = np.array(tb.build_rest_momenta(
            physics, energy, e_phi, siren.utilities.SIREN_random(219)))
        np.testing.assert_allclose(p[:, 0]**2-np.sum(p[:, 1:]**2, axis=1),
                                   [0, .04, .0225], rtol=0, atol=3e-14)


def test_turning_point_repair_does_not_relax_the_invariant_bound():
    # This was the old band helper's output. Its shell error exceeds the
    # tolerance, so accepting it unchanged would conceal the band defect.
    with pytest.raises(ValueError, match='energ'):
        tb.build_rest_momenta(FlatDecay(), .4015659375, .15000013605343188,
                              UnusedRandom())
    energy = .4*FlatDecay.E_nu_max
    low, high = _decimal_band_at_energy((1., .2, .15), energy)
    for outside in (low-1e-12, high+1e-12):
        with pytest.raises(ValueError, match='energ'):
            tb.build_rest_momenta(FlatDecay(), energy, outside, UnusedRandom())
