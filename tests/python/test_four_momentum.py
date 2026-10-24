"""Independent decimal Lorentz reference, value ownership and persistence."""
from decimal import Decimal, localcontext
import math
import pickle
import struct

import pytest

import siren
from siren.math import FourMomentum as P4, Vector3D


def decimal_boost(mass, momentum, frame_mass, frame_momentum, negative=False):
    with localcontext() as ctx:
        ctx.prec = 100
        m, fm = Decimal(mass), Decimal(frame_mass)
        p, fp = [Decimal(x) for x in momentum], [Decimal(x) for x in frame_momentum]
        energy = (m*m + sum(x*x for x in p)).sqrt()
        if negative:
            energy = -energy
        norm = sum(x*x for x in fp).sqrt()
        if norm == 0:
            return [float(energy), *momentum]
        gamma = (fm*fm + norm*norm).sqrt() / fm
        bg = norm / fm
        axis = [x / norm for x in fp]
        parallel = sum(x*y for x, y in zip(p, axis))
        return [float(gamma*energy + bg*parallel),
                *(float(x + n*((gamma-1)*parallel + bg*energy)) for x, n in zip(p, axis))]


@pytest.mark.parametrize('frame_scale', [0., 1., 1e3, 1e9])
@pytest.mark.parametrize('mass,momentum', [(0., [3., 4., 0.]), (0.001, [0., 0., 0.]), (0.4, [-0.3, 0.2, -0.7])])
@pytest.mark.parametrize('negative', [False, True])
def test_known_mass_boost_against_decimal(frame_scale, mass, momentum, negative):
    fp = [3*frame_scale, -4*frame_scale, 12*frame_scale]
    particle = P4.from_mass(mass, Vector3D(*momentum), negative)
    frame = P4.from_mass(2., Vector3D(*fp))
    result = particle.boost_to_lab(frame)
    reference = decimal_boost(mass, momentum, 2., fp, negative)
    gamma = math.hypot(2., *fp) / 2.
    # Absolute component conditioning, independent of a possibly tiny result.
    bound = 32*math.ulp(1.)*gamma*(abs(particle.energy()) + math.hypot(*momentum))
    assert result.components() == pytest.approx(reference, rel=0, abs=bound)
    assert result.mass() == mass
    assert result.has_known_mass()
    restored = pickle.loads(pickle.dumps(result))
    assert restored.components() == result.components()
    assert restored.mass() == mass
    assert restored.has_known_mass()
    if frame_scale <= 1e3:
        back = result.boost_to_rest(frame)
        assert back.components() == pytest.approx(particle.components(), rel=0, abs=8*gamma*bound)
        assert back.mass() == mass


@pytest.mark.parametrize('negative', [False, True])
def test_pickle_keeps_small_known_mass_and_raw_components(negative):
    physical = P4.from_mass(1e-3, Vector3D(0, 0, 1e9), negative)
    loaded = pickle.loads(pickle.dumps(physical))
    assert loaded.components() == physical.components()
    assert loaded.mass() == 1e-3
    assert loaded.has_known_mass()
    assert physical.mass_squared() == 0.
    raw = P4.from_components(physical.components())
    raw_loaded = pickle.loads(pickle.dumps(raw))
    assert raw_loaded.components() == raw.components()
    assert raw_loaded.mass() == 0.
    assert not raw_loaded.has_known_mass()
    for components in ([1., 2., 3., 4.], [-0., 0., -0., 0.], [math.nan, math.inf, -math.inf, 0.]):
        raw = P4.from_components(components)
        loaded = pickle.loads(pickle.dumps(raw))
        assert struct.pack('4d', *loaded.components()) == struct.pack('4d', *components)


def test_record_copy_setters_and_record_pickle_keep_array_contract():
    r = siren.dataclasses.InteractionRecord()
    r.primary_mass = 9.
    r.primary_momentum = [1., 2., 3., 4.]
    r.secondary_momenta = [[-2., 0., 0., 0.]]
    raw = r.get_primary_momentum()
    assert raw.mass_squared() == -28.
    doubled = raw + raw
    assert r.primary_momentum == [1., 2., 3., 4.]
    r.set_primary_momentum(doubled)
    assert r.primary_momentum == [2., 4., 6., 8.]
    assert r.primary_mass == 9.
    r.set_secondary_momentum(0, r.get_secondary_momentum(0) + r.get_secondary_momentum(0))
    assert r.secondary_momenta == [[-4., 0., 0., 0.]]
    with pytest.raises(IndexError):
        r.get_secondary_momentum(1)
    with pytest.raises(IndexError):
        r.set_secondary_momentum(1, raw)
    loaded = pickle.loads(pickle.dumps(r))
    assert loaded.primary_momentum == r.primary_momentum
    assert loaded.primary_mass == 9.
    assert loaded.get_primary_momentum().mass_squared() == -112.
    assert not loaded.get_primary_momentum().has_known_mass()


@pytest.mark.parametrize('components', [[1., 2., 0., 0.], [math.nan, 0., 0., 0.], [math.inf, 0., 0., 0.]])
def test_raw_nonparticle_operations_fail_deliberately(components):
    raw = P4.from_components(components)
    rest = P4.from_mass(2., Vector3D(0., 0., 0.))
    with pytest.raises(ValueError):
        raw.mass()
    with pytest.raises(ValueError):
        raw.boost_to_lab(rest)
    with pytest.raises(ValueError):
        raw.boost_to_rest(rest)
    with pytest.raises(ValueError):
        rest.boost_to_lab(raw)
    with pytest.raises(ValueError):
        P4.invariant_mass(raw, rest)


def test_legacy_rk_exports_remain_available():
    rk = siren.utilities.rk
    v = rk.geom3.Vector3(0., 0., 3.)
    p = rk.P4(v, 2.)
    assert p.m() == 2.
    assert p.e() == pytest.approx(math.sqrt(13.))
    assert rk.invMass(p, p) == pytest.approx(4.)


def test_raw_mass_residue_cannot_be_promoted_by_a_boost():
    raw = P4.from_components([1.0000001192092896, 0., 0., 1.])
    frame = P4.from_mass(0.135, Vector3D(0., 0., 100.))
    assert raw.mass() > 0.
    assert not raw.has_known_mass()
    for operation in (raw.boost_to_lab, raw.boost_to_rest):
        with pytest.raises(ValueError, match="explicitly known"):
            operation(frame)
    particle = P4.from_mass(0., raw.momentum())
    with pytest.raises(ValueError, match="explicitly known"):
        particle.boost_to_lab(P4.from_components(frame.components()))
    assert particle.boost_to_lab(frame).mass() == 0.
