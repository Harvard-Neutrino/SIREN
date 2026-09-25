"""Constrained measure and independent sequential-reference checks."""

import math
import pickle

import numpy as np
import pytest
import siren


def record(gamma=2.0, mass=0.01):
    r = siren.dataclasses.InteractionRecord()
    pt = siren.dataclasses.ParticleType
    r.signature.primary_type = pt.Pi0
    r.signature.target_type = pt.Decay
    r.signature.secondary_types = [pt.Gamma, pt.N4, pt.N5]
    r.primary_mass = 0.1349768
    p = r.primary_mass * math.sqrt(gamma * gamma - 1)
    r.primary_momentum = [r.primary_mass * gamma, 0.0, 0.0, p]
    r.secondary_masses = [0.0, mass, mass]
    r.secondary_momenta = [[0.0] * 4 for _ in range(3)]
    r.secondary_helicities = [0.0] * 3
    r.interaction_vertex = [0.0] * 3
    r.primary_initial_position = [0.0] * 3
    return r


def target():
    return siren.geometry.Sphere(
        siren.geometry.Placement(siren.math.Vector3D(15, 0, 20)), 8.0, 0.0
    )


def boost(p, frame, inverse=False):
    b = np.asarray(frame[1:]) / frame[0] * (-1 if inverse else 1)
    b2 = b @ b
    if b2 == 0:
        return np.asarray(p)
    g = 1 / math.sqrt(1 - b2)
    dot = b @ p[1:]
    return np.r_[g * (p[0] + dot), p[1:] + ((g - 1) * dot / b2 + g * p[0]) * b]


def direction(rng):
    u = rng.uniform(-1, 1)
    phi = rng.uniform(-math.pi, math.pi)
    return np.array(
        [math.sqrt(1 - u * u) * math.cos(phi), math.sqrt(1 - u * u) * math.sin(phi), u]
    )


def independent_sequential(r, rng, ma):
    """Uniform four angular coordinates, no native sampler/kinematics helpers."""
    M = r.primary_mass
    m = r.secondary_masses[1]
    pa = (M * M - ma * ma) / (2 * M)
    a = np.r_[(M * M + ma * ma) / (2 * M), pa * direction(rng)]
    d = boost(np.r_[ma / 2, math.sqrt(ma * ma / 4 - m * m) * direction(rng)], a)
    s = np.r_[M, 0, 0, 0] - a
    r.secondary_momenta = [boost(v, r.primary_momentum).tolist() for v in [s, d, a - d]]


@pytest.mark.parametrize(
    "ma,gamma", [(0.03, 1.0), (0.03, 4.0), (0.09, 2.0), (0.12, 2.0)]
)
def test_conservation_shells_density_mixture_and_archive(ma, gamma):
    r = record(gamma, ma / 3)
    channel = siren.injection.OnShellCascadeChannel(target(), ma, 0.4)
    copy = pickle.loads(pickle.dumps(channel))
    q1 = siren.injection.OnShellCascadeChannel(
        target(), ma, 0.4, first_daughter_probability=1
    )
    q2 = siren.injection.OnShellCascadeChannel(
        target(), ma, 0.4, first_daughter_probability=0
    )
    rng = siren.utilities.SIREN_random(715)
    for _ in range(500):
        channel.Sample(rng, None, r)
        p = np.array(r.secondary_momenta)
        np.testing.assert_allclose(
            p.sum(axis=0), r.primary_momentum, rtol=1e-12, atol=1e-13
        )
        np.testing.assert_allclose(
            p[:, 0] ** 2 - np.sum(p[:, 1:] ** 2, axis=1),
            np.array(r.secondary_masses) ** 2,
            atol=1e-13,
        )
        pair = p[1] + p[2]
        assert pair[0] ** 2 - pair[1:] @ pair[1:] == pytest.approx(ma * ma, rel=1e-10)
        params = r.interaction_parameters
        np.testing.assert_allclose(
            [params["cascade_pair_" + k] for k in ["energy", "px", "py", "pz"]],
            pair,
            rtol=1e-12,
            atol=1e-13,
        )
        q = channel.Density(None, r)
        assert q > 0
        assert q == pytest.approx(
            (q1.Density(None, r) + q2.Density(None, r)) / 2, rel=1e-13
        )
        assert q == pytest.approx(copy.Density(None, r), rel=1e-13)
    mixture = siren.injection.MultiChannelPhaseSpace([channel, copy], [0.4, 0.6])
    assert pickle.loads(pickle.dumps(mixture)).Density(None, r) == pytest.approx(q)


@pytest.mark.parametrize("kappa", [-1.0, 0.0, 5 / 13])
def test_weighted_moments_and_independent_density_integral(kappa):
    ma = 0.06
    r = record(2.0, ma / 3)
    channel = siren.injection.OnShellCascadeChannel(target(), ma, kappa)
    rng = siren.utilities.SIREN_random(29712)
    values = []
    M = r.primary_mass
    ea = (M * M + ma * ma) / (2 * M)
    pa = (M * M - ma * ma) / (2 * M)
    pd = math.sqrt(ma * ma / 4 - (ma / 3) ** 2)
    for _ in range(12000):
        channel.Sample(rng, None, r)
        p = np.array(r.secondary_momenta)
        a = boost(p[1] + p[2], r.primary_momentum, True)
        d = boost(p[1], r.primary_momentum, True)
        c = (d[0] - ea / 2) / (pa * pd / ma)
        physical = (1 + kappa * c * c) / (16 * math.pi**2 * (1 + kappa / 3))
        w = physical / channel.Density(None, r)
        values.append([w, w * c * c, w * (a[3] / pa) ** 2])
    values = np.array(values)
    truth = np.array([1, (1 / 3 + kappa / 5) / (1 + kappa / 3), 1 / 3])
    assert np.all(
        abs(values.mean(axis=0) - truth)
        < 7 * values.std(axis=0, ddof=1) / math.sqrt(len(values))
    )
    prng = np.random.default_rng(153)
    ratios = []
    for _ in range(18000):
        independent_sequential(r, prng, ma)
        ratios.append(channel.Density(None, r) * 16 * math.pi**2)
    ratios = np.array(ratios)
    assert abs(ratios.mean() - 1) < 7 * ratios.std(ddof=1) / math.sqrt(len(ratios))


def test_measure_constraints_reject_continuous_and_other_surfaces():
    inj = siren.injection
    m = inj.PhaseSpaceMeasure.OnShellCascade(0.03)
    assert m != inj.PhaseSpaceMeasure.OnShellCascade(0.06)
    for other in [
        inj.PhaseSpaceMeasure.Recursive2Body(),
        inj.PhaseSpaceMeasure.OnShellCascade(0.06),
        inj.PhaseSpaceMeasure.OnShellCascade(0.03, 1, 0, 2),
    ]:
        assert not inj.PhaseSpaceDensityConvertible(
            inj.PhaseSpaceTopology.Decay3Body, m, other
        )
    a = inj.OnShellCascadeChannel(target(), 0.03)
    b = inj.OnShellCascadeChannel(target(), 0.06)
    with pytest.raises((ValueError, RuntimeError)):
        inj.MultiChannelPhaseSpace([a, b], [0.5, 0.5])
    with pytest.raises((ValueError, RuntimeError)):
        inj.OnShellCascadeChannel(target(), 0.03, orientation_weights=[0, 0.5, 0.5])
    r = record()
    r.secondary_masses = [0, 0.01, 0.011]
    with pytest.raises((ValueError, RuntimeError)):
        a.Sample(siren.utilities.SIREN_random(1), None, r)


@pytest.mark.parametrize('parent_mass,energy,pair_mass', [(.1349768,10.,.001),(.547862,20.,.003),(.1349768,60.,.01)])
@pytest.mark.parametrize('weights', [[1,0,0],[.1,.9,0],[.1,0,.9]])
def test_high_boost_own_samples_have_positive_density(parent_mass, energy, pair_mass, weights):
    r = record(energy/parent_mass, pair_mass/3)
    r.primary_mass = parent_mass
    r.primary_momentum = [energy,0,0,math.sqrt(energy**2-parent_mass**2)]
    channel = siren.injection.OnShellCascadeChannel(target(),pair_mass,0.4,weights)
    rng = siren.utilities.SIREN_random(9204)
    for _ in range(400):
        channel.Sample(rng,None,r)
        assert channel.Density(None,r)>0
        p = np.array(r.secondary_momenta, dtype=np.longdouble)
        np.testing.assert_allclose(p.sum(axis=0), r.primary_momentum, rtol=2e-12, atol=2e-13)
        pair = p[1]+p[2]
        error = abs(pair[0]**2-pair[1:]@pair[1:]-pair_mass**2)
        assert error < 128*np.finfo(float).eps*float(pair@pair)+1e-10*pair_mass**2


@pytest.mark.parametrize('energy,pair_mass,angle', [(40.,.001,2.9),(135.,.001,2.9),(135.,.003,2.9),(60.,.001,6.)])
@pytest.mark.parametrize('weights', [[1,0,0],[.1,0,.9],[.1,.45,.45]])
def test_off_axis_high_boost_own_samples_have_positive_density(energy, pair_mass, angle, weights):
    # A target off the parent's axis selects pairs emitted backward in the parent
    # frame. Their lab vectors carry boost rounding of order eps*gamma*E*, far above
    # eps times their own small lab energy.
    M = .1349768
    r = record(energy/M, pair_mass*.3)
    r.primary_momentum = [energy,0,0,math.sqrt(energy**2-M**2)]
    shape = siren.geometry.Sphere(siren.geometry.Placement(
        siren.math.Vector3D(100*math.tan(math.radians(angle)),0,100)),2.,0.)
    channel = siren.injection.OnShellCascadeChannel(shape,pair_mass,0.,weights)
    rng = siren.utilities.SIREN_random(5150)
    for _ in range(1000):
        channel.Sample(rng,None,r)
        assert channel.Density(None,r)>0


@pytest.mark.parametrize('gamma,relative', [(2.,1e-6),(1000.,1e-2)])
def test_pair_mass_check_still_rejects_other_pair_masses(gamma, relative):
    r = record(gamma,.003)
    sampler = siren.injection.OnShellCascadeChannel(target(),.01*(1+relative),0.,[1,0,0])
    channel = siren.injection.OnShellCascadeChannel(target(),.01,0.,[1,0,0])
    rng = siren.utilities.SIREN_random(77)
    for _ in range(200):
        sampler.Sample(rng,None,r)
        assert sampler.Density(None,r)>0
        assert channel.InternalDensity(r)==0
        assert channel.Density(None,r)==0


@pytest.mark.parametrize('weights', [[1,0,0],[.1,.9,0],[.1,0,.9]])
def test_tabulated_parent_keeps_its_energy_and_gross_mismatch_fails(weights):
    r=record(8)
    r.primary_momentum=np.asarray(r.primary_momentum,dtype=np.float32).astype(float).tolist()
    before=list(r.primary_momentum)
    onshell=[math.hypot(r.primary_mass,math.hypot(*before[1:])),*before[1:]]
    assert onshell[0]!=before[0]
    channel=siren.injection.OnShellCascadeChannel(target(),.03,.4,weights)
    rng=siren.utilities.SIREN_random(882)
    channel.Sample(rng,None,r)
    # The parent belongs to the upstream vertex, so its tabulated energy is kept;
    # the daughters conserve the on-shell four-momentum (hypot(m,|p|), p).
    assert list(r.primary_momentum)==before
    assert 'decay_input_energy' not in r.interaction_parameters
    assert channel.Density(None,r)>0
    np.testing.assert_allclose(np.sum(r.secondary_momenta,axis=0),onshell,atol=2e-13)
    r.primary_momentum=[r.primary_momentum[0]*1.01,*r.primary_momentum[1:]]
    assert channel.Density(None,r)==0
    with pytest.raises(siren.utilities.InjectionFailure,match='mass shell'):
        channel.Sample(rng,None,r)


@pytest.mark.parametrize('proposal', [True, False])
def test_parent_at_source_emin_keeps_its_energy_in_a_run(proposal):
    # A row 1e-6 above its mass shell with the source's emin equal to its energy.
    # Writing the projected (lower) energy back into the record would move the
    # parent below emin and make the source density zero.
    pt = siren.dataclasses.ParticleType
    M, mA, m, pz = .1349768, .03, .01, .3
    e = math.hypot(M,pz)*(1+1e-6)
    rows = [[e,M,0.,0.,pz,0.,0.,0.,1.]]
    source = siren.distributions.PrimaryExternalDistribution(
        ['E','m','px','py','pz','x','y','z','weight'],rows,e)
    signature = siren.dataclasses.InteractionSignature()
    signature.primary_type = pt.Pi0
    signature.target_type = pt.Decay
    signature.secondary_types = [pt.Gamma,pt.N4,pt.N5]
    physical = siren.injection.OnShellCascadeChannel(target(),mA,0.,[1.,0.,0.])
    model = siren.injection.PhaseSpaceDecay(signature,[0.,m,m],1e-3*7.8e-9,7.8e-9,physical)
    kinematics = {'kinematics':siren.channels.on_shell_cascade(target(),mA,0.,(.1,.45,.45),.5)} if proposal else {}
    vertex = siren.Vertex('Pi0',model,distributions=[source],physical=[source],
                          weighting=siren.Fixed(),**kinematics)
    result = siren.Simulation(detector=siren.detector.DetectorModel(),primary=vertex,
                              events=20,seed=7).run(on_failure='raise',on_shortfall='raise')
    assert len(result.events)==20
    weights = np.asarray(result.weights)
    assert np.all(np.isfinite(weights)) and np.all(weights>0)
    for tree in result.events:
        rec = tree.tree[0].record
        assert rec.primary_momentum[0]==e
        assert 'decay_input_energy' not in rec.interaction_parameters
        assert np.sum(rec.secondary_momenta,axis=0)[0]==pytest.approx(math.hypot(M,pz),rel=1e-12)


def test_normalization_survives_repeated_archive_loads_bitwise():
    channel=siren.injection.OnShellCascadeChannel(target(),.03,.4,[.1,.1,.35])
    r=record()
    rng=siren.utilities.SIREN_random(191)
    restored=channel
    for _ in range(10): restored=pickle.loads(pickle.dumps(restored))
    assert pickle.dumps(restored)==pickle.dumps(channel)
    model=siren.injection.PhaseSpaceDecay(r.signature,r.secondary_masses,1.,1.,channel)
    assert model==pickle.loads(pickle.dumps(model))
    for _ in range(100):
        channel.Sample(rng,None,r)
        assert channel.Density(None,r)==restored.Density(None,r)


def test_public_conversion_rejects_distinct_constraints_even_at_zero():
    inj=siren.injection
    a=inj.PhaseSpaceMeasure.OnShellCascade(.03)
    channel=inj.OnShellCascadeChannel(target(),.03)
    mixture=inj.MultiChannelPhaseSpace([channel],[1.])
    r=record();channel.Sample(siren.utilities.SIREN_random(75),None,r)
    for b in [inj.PhaseSpaceMeasure.OnShellCascade(.06),inj.PhaseSpaceMeasure.OnShellCascade(.03,1,0,2)]:
        for q in [0.,1.]:
            with pytest.raises((ValueError,RuntimeError)):
                inj.ConvertDensity(q,a,b,inj.PhaseSpaceTopology.Decay3Body,r)
        with pytest.raises((ValueError,RuntimeError)):
            mixture.DensityIn(None,r,b)



@pytest.mark.parametrize('kappa',[-1.,0.,5/13])
def test_bounded_weight_closure_detects_six_percent_density_error(kappa):
    # 80% physical component gives p/q <= 1.25, bounding variance independently
    # of a pilot or a lucky seed. The older narrow proposal remains a stress
    # test, but no longer supplies the normalization acceptance threshold.
    r=record(2.,.02)
    channel=siren.injection.OnShellCascadeChannel(target(),.06,kappa,[.8,.1,.1])
    rng=siren.utilities.SIREN_random(4029)
    M=r.primary_mass; ea=(M*M+.06**2)/(2*M); pa=(M*M-.06**2)/(2*M)
    pd=math.sqrt(.06**2/4-.02**2)
    values=[]
    for _ in range(50000):
        channel.Sample(rng,None,r)
        p=np.array(r.secondary_momenta)
        c=(boost(p[1],r.primary_momentum,True)[0]-ea/2)/(pa*pd/.06)
        physical=(1+kappa*c*c)/(16*math.pi**2*(1+kappa/3))
        values.append(physical/channel.Density(None,r))
    values=np.array(values)
    assert values.max()<=1.25*(1+1e-10)
    assert values.std(ddof=1)/math.sqrt(len(values))<.002
    assert values.mean()==pytest.approx(1.,abs=.012)
    with pytest.raises(AssertionError):
        assert (values/.94).mean()==pytest.approx(1.,abs=.012)


def test_nonanalytic_volume_is_validated_once_and_preserved():
    # Two disjoint 6 m cubes: a composite with no analytic volume in SIREN.
    cube=lambda x: siren.geometry.Box(widths=[6.,6.,6.],center=[x,0.,20.])
    shape=siren.geometry.BooleanGeometry(siren.geometry.BooleanOperation.UNION,cube(12.),cube(19.))
    with pytest.raises(RuntimeError,match='volume'):
        siren.injection.OnShellCascadeChannel(shape,.03)
    volume=2*6.**3
    channel=siren.injection.OnShellCascadeChannel(shape,.03,volume=volume)
    restored=pickle.loads(pickle.dumps(channel))
    r=record();rng=siren.utilities.SIREN_random(23)
    for _ in range(100):
        channel.Sample(rng,None,r)
        assert channel.Density(None,r)>0
        assert channel.Density(None,r)==restored.Density(None,r)
    # No volume component needs no volume, even for nonanalytic geometry.
    siren.injection.OnShellCascadeChannel(shape,.03,orientation_weights=[.1,.9,0])
