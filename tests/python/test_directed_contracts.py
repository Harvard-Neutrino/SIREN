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

@pytest.mark.parametrize('gamma', [1001., 1e4, 1e8])
@pytest.mark.parametrize('oblique', [False, True])
def test_large_boost_geometry_sampling_uses_declared_full_sphere_fallback(gamma, oblique):
    r = record3(gamma)
    sig = r.signature
    sig.secondary_types = [siren.particles.NuLight, siren.particles.Gamma]
    r.signature = sig
    r.secondary_masses = [0., 0.]
    r.secondary_momenta = [[0.] * 4] * 2
    r.secondary_helicities = [0.] * 2
    if oblique:
        p = math.sqrt(gamma * gamma - 1.)
        r.primary_momentum = [gamma, *[p*v/math.sqrt(14) for v in (1,2,3)]]
    channel = siren.injection.DetectorDirected2BodyChannel(box((0,0,100), 10))
    random = siren.utilities.SIREN_random(9)
    for _ in range(30):
        channel.Sample(random, None, r)
        assert not channel.DirectingActive(r)
        assert channel.Density(None, r) == pytest.approx(1/(4*math.pi), rel=1e-14)


@pytest.mark.parametrize('index', [0, 1, 2])
def test_three_body_off_shell_daughters_have_zero_proposal_support(index):
    j = siren.injection
    channel = j.DetectorDirected3BodyChannel(factorization=j.ThreeBodyMode.Direct, target=box())
    r = record3(10)
    channel.Sample(siren.utilities.SIREN_random(91), None, r)
    momenta = list(r.secondary_momenta)
    momenta[index] = list(momenta[index])
    momenta[index][0] *= 1.01
    r.secondary_momenta = momenta
    assert channel.Density(None, r) == 0
    assert not channel.DirectingActive(r)


def test_python_q2_parameters_reach_native_constructor():
    j = siren.injection
    r = record3()
    sig = r.signature
    sig.primary_type = siren.particles.NuMu
    sig.target_type = siren.dataclasses.ParticleType.O16Nucleus
    sig.secondary_types = [siren.particles.NuMu, siren.dataclasses.ParticleType.O16Nucleus]
    for kwargs in [dict(q2_mode=j.ScatteringQ2Mode.Propagator, mediator_mass=.2),
                   dict(q2_mode=j.ScatteringQ2Mode.Tabulated,
                        q2_cdf_nodes=[0.,1.,10.], q2_cdf_values=[0.,.5,1.])]:
        spec = siren.Directed('NuMu', box(), fraction=1., volume=1000., **kwargs)
        native = spec._directed_channel()._build(sig)
        assert isinstance(native, j.DetectorDirectedScatteringChannel)
        assert pickle.loads(pickle.dumps(native)).Measure() == native.Measure()
    three = channels.toward_3body('NuLight', box(), strategy='direct', volume=1000.)
    assert isinstance(three._build(record3().signature), j.DetectorDirected3BodyChannel)

@pytest.mark.parametrize('gamma', [100., 1000.])
def test_active_cone_density_matches_massless_aberration(gamma):
    r = record3(gamma)
    sig = r.signature
    sig.secondary_types = [siren.particles.NuLight, siren.particles.Gamma]
    r.signature = sig
    r.secondary_masses = [0., 0.]
    r.secondary_momenta = [[0.] * 4] * 2
    r.secondary_helicities = [0.] * 2
    channel = siren.injection.DetectorDirected2BodyChannel(
        box((0,0,100), 10), 0, siren.injection.DirectedMode.Cone)
    random = siren.utilities.SIREN_random(102)
    beta = np.sqrt(np.longdouble(1)-np.longdouble(1)/gamma**2)
    omega = 2*math.pi * 75/(100*(100+math.sqrt(10000-75)))
    for _ in range(100):
        channel.Sample(random, None, r)
        p = np.asarray(r.secondary_momenta[0][1:], dtype=np.longdouble)
        c = p[2]/np.sqrt(np.dot(p,p))
        oracle = float(gamma**2 * (1-beta*c)**2 / omega)
        assert channel.DirectingActive(r)
        assert channel.Density(None, r) == pytest.approx(oracle, rel=3e-8, abs=0)


def scattering_channel_and_record(qmode):
    j = siren.injection
    kwargs = {'Geometry': {}, 'Propagator': {'mediator_mass': .2},
              'Tabulated': {'q2_cdf_nodes': [0., 1., 10.],
                            'q2_cdf_values': [0., .5, 1.]}}[qmode]
    channel = j.DetectorDirectedScatteringChannel(
        box((0, 0, 100)), q2_mode=getattr(j.ScatteringQ2Mode, qmode), **kwargs)
    r = siren.dataclasses.InteractionRecord()
    sig = siren.dataclasses.InteractionSignature()
    sig.primary_type = siren.particles.NuMu
    sig.target_type = siren.particles.O16Nucleus
    sig.secondary_types = [siren.particles.NuMu, siren.particles.O16Nucleus]
    r.signature = sig
    r.primary_mass = .1
    r.primary_momentum = [math.hypot(.1, .99), 0., 0., .99]
    r.target_mass = .938
    r.secondary_masses = [.1, .938]
    r.secondary_momenta = [[0.] * 4] * 2
    r.secondary_helicities = [0., 0.]
    r.interaction_vertex = [0., 0., 0.]
    return channel, r


@pytest.mark.parametrize('qmode', ['Geometry', 'Propagator', 'Tabulated'])
@pytest.mark.parametrize('bad_input', ['off_shell', 'zero_target', 'nan_target',
                                      'negative_mass', 'nan_momentum', 'zero_energy'])
def test_scattering_invalid_initial_state_has_zero_support(qmode, bad_input):
    channel, r = scattering_channel_and_record(qmode)
    if bad_input == 'off_shell':
        r.primary_momentum = [1., 0., 0., .99]
    elif bad_input == 'zero_target':
        r.target_mass = 0.
    elif bad_input == 'nan_target':
        r.target_mass = math.nan
    elif bad_input == 'negative_mass':
        r.primary_mass = -.1
    elif bad_input == 'nan_momentum':
        r.primary_momentum = [1., 0., 0., math.nan]
    else:
        r.primary_mass = 0.
        r.primary_momentum = [0., 0., 0., 0.]
    assert channel.Density(None, r) == 0.
    assert not channel.DirectingActive(r)
    with pytest.raises(siren.utilities.InjectionFailure):
        channel.Sample(siren.utilities.SIREN_random(19), None, r)


@pytest.mark.parametrize('qmode', ['Geometry', 'Propagator', 'Tabulated'])
def test_scattering_retains_supported_rounded_source(qmode):
    channel, r = scattering_channel_and_record(qmode)
    raw = list(r.primary_momentum)
    raw[0] *= 1 + 1e-6
    r.primary_momentum = raw
    channel.Sample(siren.utilities.SIREN_random(19), None, r)
    assert channel.Density(None, r) > 0.
    np.testing.assert_array_equal(r.primary_momentum, raw)
