"""Directed proposal contracts independent of proposal implementation details."""
import math
import pickle
from unittest.mock import patch

import numpy as np
import pytest
import siren
from siren import directed_tiling, channels


def box(center=(0, 0, 0), width=10):
    return siren.geometry.Box(siren.geometry.Placement(siren.math.Vector3D(*center)),
                              width, width, width)


def test_containment_does_not_drop_a_thin_protrusion():
    a, b = box(width=2), box(center=(-.0001, 0, 0), width=2)
    point = siren.math.Vector3D(.99995, 0, 0)
    assert a.IsInside(point) and not b.IsInside(point)
    assert directed_tiling.dedup_contained([a, b]) == [a, b]


@pytest.mark.parametrize('n', [0, -1, True, 2.5, (1, 0, 2), (1, 2), (1, 2.5, 3)])
def test_grid_requires_positive_integer_subdivisions(n):
    with pytest.raises(ValueError, match='positive integers'):
        directed_tiling.grid_cells(box(), n)


def test_directed_scattering_retains_named_daughter():
    spec = siren.Directed('O16Nucleus', box(),
                          variable=siren.injection.ScatteringVariable.Q2)
    with patch.object(channels, 'scatter_toward', return_value='sentinel') as factory:
        assert spec._directed_channel() == 'sentinel'
    assert factory.call_args.kwargs['daughter'] == 'O16Nucleus'


def record3(gamma=1):
    r = siren.dataclasses.InteractionRecord()
    sig = siren.dataclasses.InteractionSignature()
    sig.primary_type = siren.particles.N4
    sig.target_type = siren.dataclasses.ParticleType.Decay
    sig.secondary_types = [siren.particles.NuLight, siren.particles.EMinus,
                           siren.particles.EPlus]
    r.signature = sig
    r.primary_mass = 1.
    r.primary_momentum = [gamma, 0., 0., math.sqrt(gamma * gamma - 1)]
    r.secondary_masses = [.02, .03, .04]
    r.secondary_momenta = [[0.] * 4] * 3
    r.secondary_helicities = [0.] * 3
    r.interaction_vertex = [0.] * 3
    return r


@pytest.mark.parametrize('recursive', [False, True])
@pytest.mark.parametrize('mass_mode', ['Uniform', 'BreitWigner', 'PowerLaw', 'Tabulated'])
def test_three_body_pickle_continues_sampling_and_density(recursive, mass_mode):
    inj = siren.injection
    cfg = dict(factorization=inj.ThreeBodyMode.Recursive if recursive else inj.ThreeBodyMode.Direct,
               target=box(), directed_index=1 if recursive else 0,
               mass_mode=getattr(inj.InvariantMassMode, mass_mode),
               resonance_mass=.4, resonance_width=.1, power_law_nu=.7,
               mass_cdf_nodes=[0., .25, 1.], mass_cdf_values=[0., .1, 1.])
    if recursive:
        cfg.update(spectator_index=0, pair_first_index=1, pair_second_index=2)
    channel = inj.DetectorDirected3BodyChannel(**cfg)
    loaded = pickle.loads(pickle.dumps(channel))
    a, b = record3(20), record3(20)
    ra, rb = siren.utilities.SIREN_random(42), siren.utilities.SIREN_random(42)
    for _ in range(10):
        channel.Sample(ra, None, a)
        loaded.Sample(rb, None, b)
        np.testing.assert_array_equal(a.secondary_momenta, b.secondary_momenta)
        assert channel.Density(None, a) == loaded.Density(None, b) > 0
        np.testing.assert_allclose(np.sum(a.secondary_momenta, axis=0),
                                   a.primary_momentum, rtol=2e-14, atol=2e-14)

@pytest.mark.parametrize('name,args', [
    ('UniformMapping', (1., 1.)), ('UniformMapping', (0., math.inf)),
    ('UniformMapping', (math.nan, 1.)),
    ('BreitWignerMapping', (1., 0., 0., 2.)),
    ('BreitWignerMapping', (1., .1, 2., 1.)),
    ('BreitWignerMapping', (1., math.inf, 0., 2.)),
    ('LogMapping', (1., math.inf)),
    ('ExponentialMapping', (1., -1000., -999.)),
    ('ExponentialMapping', (math.nan, 0., 1.)),
    ('GaussianMapping', (0., math.inf, 0., 1.)),
])
def test_fixed_maps_reject_invalid_or_unrepresentable_normalization(name, args):
    with pytest.raises((ValueError, RuntimeError)):
        getattr(siren.injection, name)(*args)

@pytest.mark.parametrize('recursive', [False, True])
@pytest.mark.parametrize('gamma', [1., 100., 1e4, 1e8, 1e12])
def test_pair_density_uses_known_masses_under_large_boost(recursive, gamma):
    inj = siren.injection
    cfg = dict(factorization=inj.ThreeBodyMode.Recursive if recursive else inj.ThreeBodyMode.Direct,
               target=box(), directed_index=1 if recursive else 0,
               mass_mode=inj.InvariantMassMode.BreitWigner,
               resonance_mass=.4, resonance_width=.1)
    if recursive:
        cfg.update(spectator_index=0, pair_first_index=1, pair_second_index=2)
    channel = inj.DetectorDirected3BodyChannel(**cfg)
    r, random = record3(gamma), siren.utilities.SIREN_random(1907)
    mg = .4 * .1
    norm = math.atan(((1-.02)**2-.4**2)/mg) - math.atan(((.03+.04)**2-.4**2)/mg)
    for _ in range(100):
        channel.Sample(random, None, r)
        sampled = r.interaction_parameters['phase_space_s_pair' if recursive else 'phase_space_s_X']
        oracle = mg / (norm * ((sampled-.4**2)**2+mg**2)) / (4*math.pi)**2
        # The independent analytic density uses the actual sampled pair invariant.
        assert channel.Density(None, r) == pytest.approx(oracle, rel=2e-10, abs=0)


def test_direct_pair_mass_sugar_does_not_require_a_spectator():
    directed = siren.Directed('NuLight', box(), pair_mass=channels.PairMass.uniform())
    with patch.object(channels, 'toward_3body', return_value='sentinel') as factory:
        assert directed._directed_channel() == 'sentinel'
    assert factory.call_args.kwargs['strategy'] == 'direct'


@pytest.mark.parametrize('fraction', [-1., 1.1, math.nan, math.inf])
def test_directed_fraction_cannot_silently_drop_fallback(fraction):
    with pytest.raises(siren.utilities.ConfigurationError):
        channels.toward('Gamma', box(), fraction=fraction)


def test_conflicting_directed_factories_fail_explicitly():
    with pytest.raises(siren.utilities.ConfigurationError, match='cannot be combined'):
        siren.Directed('Gamma', box(), tiling=2, variable=siren.injection.ScatteringVariable.Q2)

@pytest.mark.parametrize('case', ['two_mode', 'scatter_index', 'scatter_variable',
                                  'nan_mediator', 'topology', 'spectator', 'factorization'])
def test_invalid_directed_configuration_is_rejected_at_construction(case):
    j = siren.injection
    common = dict(factorization=j.ThreeBodyMode.Recursive, target=box(),
                  spectator_index=0, pair_first_index=1, pair_second_index=2)
    factories = {
        'two_mode': lambda: j.DetectorDirected2BodyChannel(box(), 0, j.DirectedMode(33)),
        'scatter_index': lambda: j.DetectorDirectedScatteringChannel(box(), directed_index=2),
        'scatter_variable': lambda: j.DetectorDirectedScatteringChannel(box(), variable=j.ScatteringVariable(33)),
        'nan_mediator': lambda: j.DetectorDirectedScatteringChannel(box(), q2_mode=j.ScatteringQ2Mode.Propagator, mediator_mass=math.nan),
        'topology': lambda: j.DetectorDirected3BodyChannel(**common, topology=j.PhaseSpaceTopology.Decay2Body),
        'spectator': lambda: j.DetectorDirected3BodyChannel(**common, directed_index=0),
        'factorization': lambda: j.DetectorDirected3BodyChannel(target=box(), factorization=j.ThreeBodyMode(33)),
    }
    with pytest.raises((ValueError, RuntimeError)):
        factories[case]()


@pytest.mark.parametrize('bounds', [(1., 1.), (2., 1.), (0., math.inf)])
def test_propagator_map_requires_a_normalizable_interval(bounds):
    with pytest.raises(ValueError):
        siren.injection.PropagatorMapping(.1, *bounds)
