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
