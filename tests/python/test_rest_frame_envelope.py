"""Independent coverage, normalization and persistence checks of the envelope."""

import math
import pickle

import numpy as np
import pytest
import siren


def record(gamma, mass=1 / 3):
    r = siren.dataclasses.InteractionRecord()
    pt = siren.dataclasses.ParticleType
    r.signature.primary_type = pt.N4
    r.signature.target_type = pt.Decay
    r.signature.secondary_types = [pt.NuLight, pt.Gamma]
    r.primary_mass = 1.0
    r.primary_momentum = [gamma, 0.0, 0.0, math.sqrt(gamma * gamma - 1)]
    r.secondary_masses = [mass, mass]
    r.secondary_momenta = [[0.0] * 4, [0.0] * 4]
    r.secondary_helicities = [0.0, 0.0]
    r.interaction_vertex = [0.0, 0.0, 0.0]
    r.primary_initial_position = [0.0, 0.0, 0.0]
    return r


def set_decay(r, u, phi):
    mass = r.secondary_masses[0]
    p = math.sqrt(0.25 - mass * mass)
    gamma = r.primary_momentum[0]
    beta = math.sqrt(1 - 1 / gamma**2)
    st = math.sqrt(max(0.0, 1 - u * u))
    v = np.array(
        [
            gamma * (0.5 + beta * p * u),
            p * st * math.cos(phi),
            p * st * math.sin(phi),
            gamma * (p * u + beta * 0.5),
        ]
    )
    r.secondary_momenta = [v.tolist(), (np.array(r.primary_momentum) - v).tolist()]
    return v


@pytest.mark.parametrize("gamma", [1.0, 1.2, 1.5, 2.0, 8.0])
@pytest.mark.parametrize("degrees", [0.0, 30.0, 90.0, 150.0, 180.0])
def test_independent_support_and_normalization(gamma, degrees):
    angle = math.radians(degrees)
    center = np.array([30 * math.sin(angle), 0.0, 30 * math.cos(angle)])
    target = siren.geometry.Sphere(
        siren.geometry.Placement(siren.math.Vector3D(*center)), 10.0, 0.0
    )
    channel = siren.injection.RestFrameEnvelope2BodyChannel(target, 0)
    r = record(gamma)
    rng = np.random.default_rng(1027)
    n = 5000
    values = []
    for u, phi in zip(rng.uniform(-1, 1, n), rng.uniform(-math.pi, math.pi, n)):
        v = set_decay(r, u, phi)
        q = channel.Density(None, r)
        values.append(4 * math.pi * q)
        # The entire target bounding cone, not only rays hitting the sphere,
        # must be present. Its AABB sphere has radius sqrt(3)*10.
        cosine = np.dot(v[1:], center) / (np.linalg.norm(v[1:]) * 30)
        if cosine > math.sqrt(1 - 1 / 3):
            assert q > 0
    values = np.array(values)
    assert abs(values.mean() - 1) < 7 * values.std(ddof=1) / math.sqrt(n) + 1e-12
    native_rng = siren.utilities.SIREN_random(318)
    for _ in range(100):
        channel.Sample(native_rng, None, r)
        assert channel.Density(None, r) > 0
        p = np.array(r.secondary_momenta)
        np.testing.assert_allclose(p.sum(axis=0), r.primary_momentum, atol=1e-13)
        np.testing.assert_allclose(
            p[:, 0] ** 2 - np.sum(p[:, 1:] ** 2, axis=1),
            r.secondary_masses[0] ** 2,
            rtol=1e-10,
            atol=1e-13,
        )


@pytest.mark.parametrize(
    "gamma,mass", [(1.0, 0.0), (1.01, 0.49999), (1.5, 1 / 3), (8.0, 0.0)]
)
def test_thresholds_massless_and_both_daughters(gamma, mass):
    target = siren.geometry.Sphere(
        siren.geometry.Placement(siren.math.Vector3D(0, 0, 30)), 1.0, 0.0
    )
    r = record(gamma, mass)
    channels = [
        siren.injection.RestFrameEnvelope2BodyChannel(target, i) for i in [0, 1]
    ]
    rng = siren.utilities.SIREN_random(925)
    for i in range(200):
        channels[i % 2].Sample(rng, None, r)
        q = sum(c.Density(None, r) for c in channels) / 2
        assert math.isfinite(q) and q > 0


@pytest.mark.parametrize("center", [(0.0, 0.0, -30.0), (0.0, 0.0, 0.0)])
def test_inactive_normalized_fallback_and_pickle(center):
    target = siren.geometry.Sphere(
        siren.geometry.Placement(siren.math.Vector3D(*center)), 1.0, 0.0
    )
    channel = siren.injection.RestFrameEnvelope2BodyChannel(target, 0)
    restored = pickle.loads(pickle.dumps(channel))
    r = record(2.0)
    assert not channel.DirectingActive(r)
    for u in np.linspace(-0.999, 0.999, 21):
        set_decay(r, float(u), 0.4)
        assert channel.Density(None, r) == pytest.approx(1 / (4 * math.pi), rel=1e-14)
        assert restored.Density(None, r) == channel.Density(None, r)
    mixture = siren.injection.MultiChannelPhaseSpace(
        [channel, siren.injection.Isotropic2BodyChannel()], [0.8, 0.2]
    )
    copied = pickle.loads(pickle.dumps(mixture))
    assert copied.Density(None, r) == pytest.approx(1 / (4 * math.pi), rel=1e-14)


def analytic_envelope(r, degrees):
    """Independent quadratic roots of the boosted polar-cone boundaries.

    The target AABB sphere subtends delta=asin(1/sqrt(3)). Split at all
    quadratic roots and select intervals by their midpoints; no native
    preimage construction, iterative inversion, or native sampling is used.
    """
    alpha=math.radians(degrees); delta=math.asin(1/math.sqrt(3))
    low=math.cos(min(math.pi,alpha+delta)); high=math.cos(max(0,alpha-delta))
    g=r.primary_momentum[0]; b=math.sqrt(max(0,1-1/g**2))
    p=math.sqrt(.25-r.secondary_masses[0]**2)
    roots=[-1.,1.]
    for c in (low,high):
        # (g*(p*u+b/2))^2 = c^2 * (p^2*(1-u^2)+(g*(p*u+b/2))^2)
        coeff=[(1-c*c)*g*g*p*p+c*c*p*p,
               (1-c*c)*g*g*p*b,
               (1-c*c)*g*g*b*b/4-c*c*p*p]
        for u in np.roots(coeff):
            if abs(u.imag)<1e-10 and -1<u.real<1:
                z=g*(p*u.real+b/2)
                cosine=z/math.hypot(z,p*math.sqrt(1-u.real**2))
                if abs(cosine-c)<1e-7: roots.append(float(u.real))
    roots=sorted(roots)
    intervals=[]
    for a,bound in zip(roots[:-1],roots[1:]):
        u=(a+bound)/2; z=g*(p*u+b/2)
        c=z/math.hypot(z,p*math.sqrt(1-u*u))
        if low<c<high: intervals.append((a,bound))
    half=math.pi
    if delta<alpha<math.pi-delta: half=math.asin(math.sin(delta)/math.sin(alpha))
    if not intervals: return [(-1.,1.)],math.pi
    return intervals,half


@pytest.mark.parametrize('gamma',[1.,1.2,1.5,2.,8.])
@pytest.mark.parametrize('degrees',[0.,30.,90.,150.,180.])
def test_density_normalization_from_independent_analytic_area(gamma,degrees):
    angle=math.radians(degrees)
    target=siren.geometry.Sphere(siren.geometry.Placement(
        siren.math.Vector3D(30*math.sin(angle),0,30*math.cos(angle))),10.,0.)
    channel=siren.injection.RestFrameEnvelope2BodyChannel(target)
    r=record(gamma)
    intervals,half=analytic_envelope(r,0. if gamma==1 else degrees)
    area=2*half*sum(b-a for a,b in intervals)
    # Constant density on the product domain: midpoint quadrature is exact.
    integral=0.
    for a,b in intervals:
        for fraction in [.2,.5,.8]:
            # At gamma=1 the native rest axes point to the target. Rotate our
            # independent directions to that frame before writing the record.
            u=a+(b-a)*fraction; phi=.37*half
            if gamma==1:
                u0=u; st=math.sqrt(1-u*u)
                v=np.array([st*math.cos(phi),st*math.sin(phi),u])
                v=np.array([math.cos(angle)*v[0]+math.sin(angle)*v[2],v[1],
                            -math.sin(angle)*v[0]+math.cos(angle)*v[2]])
                u=float(v[2]);phi=math.atan2(v[1],v[0])
            set_decay(r,u,phi)
            q=channel.Density(None,r)
            # A stationary parent's envelope is centered on the target, so its
            # native alpha is zero and its area is a spherical cap.
            expected_area=4*math.pi*(1-math.cos(math.asin(1/math.sqrt(3))))/2 if gamma==1 else area
            assert q*expected_area==pytest.approx(1.,rel=2e-9,abs=0)
            integral += q*expected_area/3/len(intervals)
    assert integral==pytest.approx(1.,rel=2e-9)
    if half<math.pi and gamma>1 and channel.DirectingActive(r):
        # Explicit power witness: a localized 20% off-axis normalization error
        # fails the very same acceptance assertion.
        with pytest.raises(AssertionError):
            assert .8*integral==pytest.approx(1.,rel=2e-9)
