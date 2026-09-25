"""Regressions for the 2026-09-24 follow-up review and its round-2 audit."""
import io
import math
import os
import pickle
import struct
import subprocess
import sys
import textwrap

import numpy as np
import pytest
import siren
from siren.source_importance import SourceImportanceTable

from test_on_shell_cascade import record, target

inj = siren.injection
pt = siren.dataclasses.ParticleType
M = .1349768


def float32_parent(gamma):
    """A pi0 row rounded to float32, as a tabulated source would store it."""
    return np.array([M*gamma, 0., 0., M*math.sqrt(gamma*gamma-1)],
                    np.float32).astype(float).tolist()


def two_body_record(parent, masses=(.03, 0.)):
    r = siren.dataclasses.InteractionRecord()
    r.signature.primary_type = pt.N4
    r.signature.target_type = pt.Decay
    r.signature.secondary_types = [pt.NuLight, pt.Gamma]
    r.primary_mass = M
    r.primary_momentum = list(parent)
    r.secondary_masses = list(masses)
    r.secondary_momenta = [[0.]*4 for _ in range(2)]
    r.secondary_helicities = [0., 0.]
    r.interaction_vertex = [0.]*3
    return r


@pytest.mark.parametrize('gamma', [300., 1000., 3000.])
def test_isotropic_sampler_uses_the_on_shell_parent_of_rounded_rows(gamma):
    # A float32 row at gamma 3000 has E == |p|: the raw four-vector is massless.
    parent = float32_parent(gamma)
    on_shell = [math.hypot(M, parent[3]), *parent[1:]]
    r = two_body_record(parent)
    channel = inj.Isotropic2BodyChannel(0)
    rng = siren.utilities.SIREN_random(74)
    for _ in range(200):
        r.primary_momentum = parent
        channel.Sample(rng, None, r)
        p = np.array(r.secondary_momenta)
        assert np.all(np.isfinite(p))
        assert r.primary_momentum == parent
        np.testing.assert_allclose(p.sum(axis=0), on_shell, rtol=1e-12, atol=1e-12)
        masses2 = p[:, 0]**2 - np.sum(p[:, 1:]**2, axis=1)
        np.testing.assert_allclose(masses2, [.03**2, 0.], atol=1e-6)


@pytest.mark.parametrize('gamma', [300., 1000., 3000.])
def test_directed_density_and_rest_lab_conversion_use_the_on_shell_parent(gamma):
    # Target inside the boosted cone, so no isotropic fallback hides a mismatch.
    shape = siren.geometry.Sphere(
        siren.geometry.Placement(siren.math.Vector3D(50/gamma, 0, 100)), 2/gamma, 0.)
    parent = float32_parent(gamma)
    on_shell = [math.hypot(M, parent[3]), *parent[1:]]
    r = two_body_record(parent)
    envelope = inj.RestFrameEnvelope2BodyChannel(shape, 0)
    directed = inj.DetectorDirected2BodyChannel(shape, 0)
    rest, lab = inj.PhaseSpaceMeasure.SolidAngleRest(), inj.PhaseSpaceMeasure.SolidAngleLab(0)
    rng = siren.utilities.SIREN_random(129)
    supported = 0
    for _ in range(300):
        r.primary_momentum = parent
        envelope.Sample(rng, None, r)
        values = []
        for four_momentum in (parent, on_shell):
            r.primary_momentum = four_momentum
            values.append((directed.Density(None, r), inj.ConvertDensity(
                1., rest, lab, inj.PhaseSpaceTopology.Decay2Body, r)))
        assert values[0] == values[1]
        supported += values[0][0] > 0
    assert supported > 0


def test_older_decay_channels_reject_parents_off_their_mass_shell():
    r = two_body_record([2*M*(1+1e-3), 0., 0., math.sqrt(3)*M])
    for channel in (inj.Isotropic2BodyChannel(0), inj.DetectorDirected2BodyChannel(target(), 0)):
        with pytest.raises(siren.utilities.InjectionFailure, match='mass shell'):
            channel.Sample(siren.utilities.SIREN_random(1), None, r)
        assert channel.Density(None, r) == 0


def test_point_source_decays_follow_the_exponential_flight_law():
    # The sampler's decay length used to come from a record without momentum,
    # which put every decay at the source.
    mass, energy, width, length = .1, 1., 2e-17, 25.
    momentum = math.sqrt(energy*energy - mass*mass)
    lam = momentum/mass * siren.utilities.Constants.hbarc / width
    mean = lam - length*math.exp(-length/lam)/(1 - math.exp(-length/lam))

    class Decay(siren.DecayModel):
        parent = 'N4'
        daughters = ('NuLight', 'Gamma')
        measure = siren.Measure.SolidAngleRest()
        def total_width(self): return width
        def differential_width(self, record): return width/(4*math.pi)
        def SecondaryMasses(self, types): return [0., 0.]
        def sample(self, record, random): self.sample_isotropic(record, random)

    d = siren.distributions
    distributions = [d.PrimaryMass(mass), d.Monoenergetic(energy),
                     d.FixedDirection(siren.math.Vector3D(0, 0, 1)),
                     d.PointSourcePositionDistribution(siren.math.Vector3D(0, 0, 0), length)]
    result = siren.Simulation(detector=siren.detector.DetectorModel(),
                              primary=siren.Vertex(pt.N4, Decay(), distributions=distributions),
                              events=4000, seed=17).run(on_failure='raise', on_shortfall='ignore')
    z = np.array([event.tree[0].record.interaction_vertex[2] for event in result.events])
    w = np.asarray(result.weights)
    assert np.mean(z < 1e-9) < .01
    # The truncated exponential has a standard deviation of about 7.2 m here.
    assert abs(z.mean() - mean) < .5
    assert abs(np.sum(w*z)/np.sum(w) - mean) < .5
    assert np.sum(w) == pytest.approx(1 - math.exp(-length/lam), rel=1e-6)


_AT_REST = textwrap.dedent('''
    import math, sys
    import siren
    inj, pt = siren.injection, siren.dataclasses.ParticleType
    detector = siren.detector.DetectorModel()
    materials = detector.Materials
    materials.AddMaterial('ARGON', {1000180400: 1.0})
    detector.Materials = materials
    sector = siren.detector.DetectorSector()
    sector.name, sector.level = 'argon', 0
    sector.material_id = materials.GetMaterialId('ARGON')
    sector.geo = siren.geometry.Sphere(siren.geometry.Placement(siren.math.Vector3D(0, 0, 0)), 100., 0.)
    sector.density = siren.detector.ConstantDensityDistribution(1.4)
    detector.AddSector(sector)
    xs = siren.interactions.TrivialCrossSection(1e-25, [pt.N4], [pt.Ar40Nucleus])
    source = siren.distributions.PrimaryExternalDistribution(
        ['E', 'm', 'px', 'py', 'pz', 'x', 'y', 'z'], [[.1, .1, 0., 0., 0., 0., 0., .3]])
    vertex = siren.Vertex('N4', [xs], distributions=[source], physical=[source], weighting=siren.Fixed())
    try:
        siren.Simulation(detector=detector, primary=vertex, events=5, seed=2).run(
            on_failure='raise', on_shortfall='raise')
        print('NO ERROR')
    except siren.errors.ConfigurationError as error:
        print('CONFIGURATION ERROR', error)
''')


def test_parent_at_rest_with_cross_sections_is_a_configuration_error():
    # Run in a subprocess: the old code crashed natively on this input.
    result = subprocess.run([sys.executable, '-c', _AT_REST], capture_output=True, text=True,
                            env=dict(os.environ), timeout=300)
    assert result.returncode == 0, result.stderr[-2000:]
    assert 'CONFIGURATION ERROR' in result.stdout and 'nonzero primary momentum' in result.stdout


def _two_boxes():
    """A composite target with no analytic volume in SIREN: two disjoint 6 m cubes."""
    cube = lambda x: siren.geometry.Box(widths=[6., 6., 6.], center=[x, 0., 20.])
    return siren.geometry.BooleanGeometry(siren.geometry.BooleanOperation.UNION, cube(12.), cube(19.)), 432.


def test_supplied_volume_must_agree_with_an_analytic_or_estimated_volume():
    sphere = siren.geometry.Sphere(siren.geometry.Placement(siren.math.Vector3D(15, 0, 20)), 8., 0.)
    exact = 4/3*math.pi*8**3
    with pytest.raises(RuntimeError, match='disagrees with the analytic volume'):
        inj.OnShellCascadeChannel(sphere, .03, volume=2*exact)
    with pytest.raises(RuntimeError, match='disagrees with the analytic volume'):
        inj.DetectorDirected2BodyChannel(sphere, 0, inj.DirectedMode.Volume, 2*exact)
    near = inj.OnShellCascadeChannel(sphere, .03, volume=exact*(1+1e-6))
    default = inj.OnShellCascadeChannel(sphere, .03)
    r = record()
    rng = siren.utilities.SIREN_random(3)
    for _ in range(50):
        default.Sample(rng, None, r)
        assert near.Density(None, r) == default.Density(None, r)
    # The ellipsoid now has an analytic volume, like the other primitives.
    ellipsoid = siren.geometry.Ellipsoid(
        siren.geometry.Placement(siren.math.Vector3D(15, 0, 20)), 8., 7., 6.)
    with pytest.raises(RuntimeError, match='disagrees with the analytic volume'):
        inj.OnShellCascadeChannel(ellipsoid, .03, volume=1.5*4*math.pi*8*7*6/3)
    inj.OnShellCascadeChannel(ellipsoid, .03)
    # A composite is checked against an independent estimate of the solid.
    composite, volume = _two_boxes()
    with pytest.raises(RuntimeError, match='volume'):
        inj.OnShellCascadeChannel(composite, .03)
    for wrong in (1.05*volume, .95*volume, 2*volume):
        # 2x also exceeds the bounding box, which is checked first.
        with pytest.raises(RuntimeError, match='chord-integration estimate|bounding-box'):
            inj.OnShellCascadeChannel(composite, .03, volume=wrong)
    inj.OnShellCascadeChannel(composite, .03, volume=volume)
    # SetVolume is checked the same way; a refused value leaves the channel as it was.
    parent = [M*1.5, 0., 0., M*math.sqrt(1.25)]
    r = two_body_record(parent)
    channel = inj.DetectorDirected2BodyChannel(sphere, 0, inj.DirectedMode.Volume)
    channel.Sample(siren.utilities.SIREN_random(5), None, r)
    density = channel.Density(None, r)
    assert density > 0
    with pytest.raises(RuntimeError, match='disagrees'):
        channel.SetVolume(2*exact)
    assert channel.Density(None, r) == density
    channel.SetVolume(exact*(1+1e-6))
    assert channel.Density(None, r) == density
    boxed = inj.DetectorDirected2BodyChannel(composite, 0, inj.DirectedMode.Volume, volume)
    with pytest.raises(RuntimeError, match='chord-integration estimate'):
        boxed.SetVolume(.9*volume)
    with pytest.raises(RuntimeError, match='bounding-box'):
        boxed.SetVolume(1.5*volume)


@pytest.mark.parametrize('name,make', [
    ('ellipsoid with z cuts', lambda: siren.geometry.Ellipsoid(8., 7., 6., -3., 5.)),
    ('elliptical tube', lambda: siren.geometry.EllipticalTube(3., 2., 4.)),
    ('cone segment', lambda: siren.geometry.Cone(1., 3., .5, 2., 6., .3, 4.)),
    ('torus segment', lambda: siren.geometry.Torus(5., 1.5, .5, .2, 5.)),
    ('Trd', lambda: siren.geometry.Trd(2., 1., 3., 1.5, 2.5)),
    ('parallelepiped', lambda: siren.geometry.Para(2., 1.5, 3., .3, .4, .7)),
    ('generic polycone cup', lambda: siren.geometry.GenericPolycone([0., 1., 1., .8, .8, 0.], [0., 0., 1., 1., .2, .2], .3, 2.1)),
    ('polycone tube', lambda: siren.geometry.Polycone([-.5, 0., .5], [.2, .1, .3], [.6, .9, .7], 1., 2.5)),
    ('closed mesh plate with a boss', lambda: _plate_with_boss(.3, .8, h=.2)),
])
def test_new_analytic_volumes_match_the_solids_inside_test(name, make):
    # Independent reference: the fraction of uniform bounding-box points that the
    # solid's own IsInside accepts.
    geometry = make()
    exact = inj.geometry_volume(geometry)
    box = geometry.GetWorldBoundingBox()
    low = np.array([box.min_corner.GetX(), box.min_corner.GetY(), box.min_corner.GetZ()])
    high = np.array([box.max_corner.GetX(), box.max_corner.GetY(), box.max_corner.GetZ()])
    points = np.random.default_rng(11).uniform(low, high, size=(100000, 3))
    inside = np.array([geometry.IsInside(siren.math.Vector3D(*point)) for point in points])
    fraction = inside.mean()
    volume = np.prod(high - low)
    se = volume*math.sqrt(fraction*(1 - fraction)/len(points))
    assert abs(exact - fraction*volume) < 5*se, (name, exact, fraction*volume, se)


def test_sparse_and_overfull_supplied_volumes():
    # Codex rounds 5-6: point counting accepted wrong volumes for thin shells,
    # refused a correct thin slab and a rotated shell, allowed more than the
    # box, and accepted a volume for an empty solid. Chord integration resolves
    # shells and slabs; a solid the rays cannot resolve is refused outright.
    B = siren.geometry.BooleanOperation
    box = lambda width: siren.geometry.Box(widths=[width]*3, center=[0., 0., 0.])
    make = lambda shape, volume: inj.DetectorDirected2BodyChannel(shape, 0, inj.DirectedMode.Volume, volume)
    h = 1e-4
    shell = siren.geometry.BooleanGeometry(B.SUBTRACTION, box(2.), box(2*(1 - h)))
    truth = 8*h*(3 - 3*h + h*h)
    make(shell, truth)
    with pytest.raises(RuntimeError, match='chord-integration estimate'):
        make(shell, truth/2)
    rotation = siren.math.Quaternion()
    rotation.SetEulerAnglesXYZs(0.7648539807834135, 0.7864798872285119, -0.12069603162830649)
    h = 7.55513980040491e-05
    rotated = siren.geometry.BooleanGeometry(
        siren.geometry.Placement(siren.math.Vector3D(0, 0, 0), rotation), B.SUBTRACTION, box(2.), box(2*(1 - h)))
    make(rotated, 8*h*(3 - 3*h + h*h))
    union = siren.geometry.BooleanGeometry(B.UNION, box(2.), box(1.))
    make(union, 8.)
    with pytest.raises(RuntimeError, match='bounding-box'):
        make(union, 8.08)
    # Codex round 6: a much thinner shell accepted half its volume, and a 1e-8 m
    # slab the point grid happened to hit once refused its correct volume.
    h = 1e-6
    thin = siren.geometry.BooleanGeometry(B.SUBTRACTION, box(2.), box(2*(1 - h)))
    make(thin, 8*h*(3 - 3*h + h*h))
    for factor in (.5, 2., 1.05, .95):
        with pytest.raises(RuntimeError, match='chord-integration estimate'):
            make(thin, factor*8*h*(3 - 3*h + h*h))
    x0, h = 0.00482307945804572, 1e-8
    left = siren.geometry.Box(widths=[(x0 - h/2) + 1.5, 3., 3.], center=[((x0 - h/2) - 1.5)/2, 0., 0.])
    right = siren.geometry.Box(widths=[1.5 - (x0 + h/2), 3., 3.], center=[((x0 + h/2) + 1.5)/2, 0., 0.])
    slab = siren.geometry.BooleanGeometry(B.SUBTRACTION, box(2.), siren.geometry.BooleanGeometry(B.UNION, left, right))
    make(slab, 4*h)
    with pytest.raises(RuntimeError, match='chord-integration estimate'):
        make(slab, 2*h)
    # A solid the rays cannot resolve is not accepted on trust.
    with pytest.raises(RuntimeError, match='cannot be checked'):
        make(siren.geometry.BooleanGeometry(B.SUBTRACTION, box(2.), box(2.)), 1e-5)


def _plate_with_boss(r, H, h=1e-8, as_triangles=False):
    # Codex's round-9 mesh: one closed, connected surface of a 2 m x 2 m plate h
    # thick with a central r x r boss rising to H, outward-oriented.
    xs, ys, triangles = [-1., -r, r, 1.], [-1., -r, r, 1.], []

    def quad(points, normal):
        a = np.array(points, dtype=float)
        if np.dot(np.cross(a[1] - a[0], a[2] - a[0]), normal) < 0:
            a = a[::-1]
        for ids in ((0, 1, 2), (0, 2, 3)):
            triangles.append([siren.math.Vector3D(*a[k]) for k in ids])
    for i in range(3):
        for j in range(3):
            top = H if i == 1 and j == 1 else h
            quad([(xs[i], ys[j], top), (xs[i + 1], ys[j], top), (xs[i + 1], ys[j + 1], top), (xs[i], ys[j + 1], top)], (0, 0, 1))
            quad([(xs[i], ys[j], 0), (xs[i + 1], ys[j], 0), (xs[i + 1], ys[j + 1], 0), (xs[i], ys[j + 1], 0)], (0, 0, -1))
    for j in range(3):
        for x, n in ((-1., -1), (1., 1)):
            quad([(x, ys[j], 0), (x, ys[j + 1], 0), (x, ys[j + 1], h), (x, ys[j], h)], (n, 0, 0))
    for i in range(3):
        for y, n in ((-1., -1), (1., 1)):
            quad([(xs[i], y, 0), (xs[i + 1], y, 0), (xs[i + 1], y, h), (xs[i], y, h)], (0, n, 0))
    for x, n in ((xs[1], -1), (xs[2], 1)):
        quad([(x, ys[1], h), (x, ys[2], h), (x, ys[2], H), (x, ys[1], H)], (n, 0, 0))
    for y, n in ((ys[1], -1), (ys[2], 1)):
        quad([(xs[1], y, h), (xs[2], y, h), (xs[2], y, H), (xs[1], y, H)], (0, n, 0))
    return triangles if as_triangles else siren.geometry.TriangularMesh(triangles)


def _cube_triangles(half, inward=False):
    triangles = []
    for axis in range(3):
        u, v = (axis + 1) % 3, (axis + 2) % 3
        for sign in (-1., 1.):
            quad = []
            for a, b in ((0, 0), (1, 0), (1, 1), (0, 1)):
                p = [0., 0., 0.]
                p[axis], p[u], p[v] = sign*half, half*(2*a - 1), half*(2*b - 1)
                quad.append(siren.math.Vector3D(*p))
            for ids in ((0, 1, 2), (0, 2, 3)):
                triangle = [quad[k] for k in ids]
                triangles.append(triangle[::-1] if (sign < 0) != inward else triangle)
    return triangles


def test_aligned_submicron_and_malformed_volume_triggers():
    # Codex round 7: three fixed ray directions could all run along the walls
    # of a shell rotated into their frame and agree on a third of its volume;
    # intersections merge crossings closer than 1e-9 m, so a thinner shell kept
    # only grazing chords and was confidently underestimated; and nested or
    # open meshes, whose inside test and crossings disagree, passed. The
    # estimate now averages 48 directions with their scatter as its error, and
    # refuses solids whose chords are at the intersection resolution or whose
    # crossings do not alternate entering/exiting.
    B = siren.geometry.BooleanOperation
    V = siren.math.Vector3D
    box = lambda width: siren.geometry.Box(widths=[width]*3, center=[0., 0., 0.])
    make = lambda shape, volume: inj.DetectorDirected2BodyChannel(shape, 0, inj.DirectedMode.Volume, volume)
    rotation = siren.math.Quaternion()
    rotation.SetEulerAnglesXYZs(.4, .7, .9)
    h = 1e-5
    aligned = siren.geometry.BooleanGeometry(
        siren.geometry.Placement(V(0, 0, 0), rotation), B.SUBTRACTION, box(2.), box(2*(1 - h)))
    truth = 8*h*(3 - 3*h + h*h)
    make(aligned, truth)
    with pytest.raises(RuntimeError, match='chord-integration estimate'):
        make(aligned, truth/3)
    h = 1e-10
    shell = siren.geometry.BooleanGeometry(B.SUBTRACTION, siren.geometry.Sphere(1., 0.), siren.geometry.Sphere(1. - h, 0.))
    with pytest.raises(RuntimeError, match='cannot be checked: .*shorter than'):
        make(shell, 4*math.pi*h)
    bar = siren.geometry.Box(widths=[.01, .01, 10.], center=[0., 0., 0.])
    make(siren.geometry.BooleanGeometry(B.UNION, bar, bar), 1e-3)
    mesh = siren.geometry.TriangularMesh
    make(mesh(_cube_triangles(1.)), 8.)
    make(mesh(_cube_triangles(1.) + _cube_triangles(.5, inward=True)), 7.)
    for triangles in (_cube_triangles(1.) + _cube_triangles(.5), _cube_triangles(1.)[2:]):
        with pytest.raises(RuntimeError, match='cannot be checked: .*(alternate|does not bound a solid)'):
            make(mesh(triangles), 8.)


def test_small_hidden_part_cannot_be_checked():
    # Codex round 8: a 1e-8 m sphere shell united with a 7 mm core sphere. The
    # rays met the shell everywhere and the core almost never, so the shell's
    # volume alone (8% of the truth) was accepted. A part met by too few rays
    # now bounds, by its bounding box, the volume the rays may have missed.
    B = siren.geometry.BooleanOperation
    h, r = 1e-8, .007
    shell = siren.geometry.BooleanGeometry(B.SUBTRACTION, siren.geometry.Sphere(1., 0.), siren.geometry.Sphere(1. - h, 0.))
    shell_volume = 4*math.pi*h*(1 - h + h*h/3)
    make = lambda shape, volume: inj.DetectorDirected2BodyChannel(shape, 0, inj.DirectedMode.Volume, volume)
    hidden = siren.geometry.BooleanGeometry(B.UNION, shell, siren.geometry.Sphere(r, 0.))
    for volume in (shell_volume, shell_volume + 4*math.pi*r**3/3):
        with pytest.raises(RuntimeError, match='cannot be checked: part of it may be missed'):
            make(hidden, volume)
    # A core the rays do meet is weighed and checked.
    visible = siren.geometry.BooleanGeometry(B.UNION, shell, siren.geometry.Sphere(.3, 0.))
    truth = shell_volume + 4*math.pi*.3**3/3
    make(visible, truth)
    with pytest.raises(RuntimeError, match='chord-integration estimate'):
        make(visible, truth/2)


def test_parts_the_rays_meet_but_do_not_resolve():
    # Codex round 9: a part met by many rays can still hide most of its volume:
    # a 1 cm cube left by intersecting broad boxes inside a 1e-8 m shell, a
    # generic polycone whose thin wall and floor surround a small boss, and one
    # closed mesh plate with a small boss. The polycone and the mesh now have
    # exact volumes; an intersection or subtraction node the rays barely meet
    # bounds what may be missed by its box.
    B = siren.geometry.BooleanOperation
    make = lambda shape, volume: inj.DetectorDirected2BodyChannel(shape, 0, inj.DirectedMode.Volume, volume)
    box = lambda lo, hi: siren.geometry.Box(widths=[b - a for a, b in zip(lo, hi)], center=[(a + b)/2 for a, b in zip(lo, hi)])
    h, r = 1e-8, .005
    shell_volume = 8*h*(3 - 3*h + h*h)
    shell = siren.geometry.BooleanGeometry(B.SUBTRACTION, box([-1]*3, [1]*3), box([-1 + h]*3, [1 - h]*3))
    core = siren.geometry.BooleanGeometry(B.INTERSECTION, box([-1]*3, [r]*3), box([-r]*3, [1]*3))
    for volume in (shell_volume, shell_volume + (2*r)**3):
        with pytest.raises(RuntimeError, match='cannot be checked: part of it may be missed: an intersection'):
            make(siren.geometry.BooleanGeometry(B.UNION, shell, core), volume)
    base = math.pi*h + math.pi*h*(2 - h)*(1 - h)
    cup = siren.geometry.GenericPolycone([0, 1, 1, 1 - h, 1 - h, r, r, 0], [0, 0, 1, 1, h, h, .1, .1])
    make(cup, base + math.pi*r*r*(.1 - h))
    with pytest.raises(RuntimeError, match='analytic volume'):
        make(cup, base)
    # A mesh's exact volume counts only where the rays resolve the mesh: here
    # they miss the boss, so neither volume can be checked (Codex round 10: a
    # closed mesh need not bound one solid, and rays cannot check where they
    # do not go). A chunky plate with a boss is resolved and uses its volume.
    plate = _plate_with_boss(r, .01)
    for volume in (4e-8 + (2*r)**2*(.01 - 1e-8), 4e-8):
        with pytest.raises(RuntimeError, match='cannot be checked'):
            make(plate, volume)
    chunky = _plate_with_boss(.3, .8, h=.2)
    make(chunky, 4*.2 + .36*.6)
    with pytest.raises(RuntimeError, match='analytic volume'):
        make(chunky, 4*.2)


def test_far_mesh_volume_and_retraced_profile():
    # Codex round 10: the divergence sum about the origin lost a 1 cm cube 3 km
    # away to cancellation (35% too large, the true volume refused), and a
    # generic polycone profile traced twice was given twice its volume.
    cube = siren.geometry.TriangularMesh(
        [[siren.math.Vector3D(*(3000. + .01*(c + .5) for c in (p.GetX(), p.GetY(), p.GetZ()))) for p in tri]
         for tri in _cube_triangles(.5)])
    assert inj.geometry_volume(cube) == pytest.approx(1e-6, rel=1e-9)
    inj.DetectorDirected2BodyChannel(cube, 0, inj.DirectedMode.Volume, 1e-6)
    twice = siren.geometry.GenericPolycone([0, 1, 1, 0]*2, [0, 0, 1, 1]*2)
    with pytest.raises(ValueError):
        inj.geometry_volume(twice)
    for volume in (math.pi, 2*math.pi):
        with pytest.raises(RuntimeError):
            inj.DetectorDirected2BodyChannel(twice, 0, inj.DirectedMode.Volume, volume)


def _box_triangles(lo, hi, inward=False):
    # A closed box surface, outward-oriented unless inward.
    triangles = []
    for axis in range(3):
        a, b = [k for k in range(3) if k != axis]
        for face in (0, 1):
            corners = []
            for u, v in ((0, 0), (1, 0), (1, 1), (0, 1)):
                p = [0., 0., 0.]
                p[axis], p[a], p[b] = (lo, hi)[face][axis], (lo, hi)[u][a], (lo, hi)[v][b]
                corners.append(np.array(p))
            normal = np.zeros(3)
            normal[axis] = 2*face - 1
            if (np.dot(np.cross(corners[1] - corners[0], corners[2] - corners[0]), normal) < 0) != inward:
                corners.reverse()
            triangles.extend([[siren.math.Vector3D(*corners[k]) for k in ids] for ids in ((0, 1, 2), (0, 2, 3))])
    return triangles


def test_mesh_shells_thin_layers_and_profile_precision():
    # Codex round 11. (1) A closed, vertex-connected mesh whose inward box cancels
    # its boss in the divergence sum; (2) a 10 nm plate, whose top 1 nm the
    # channel's inside test samples as outside; (3) a thin annulus far along z,
    # whose first moment lost 3.4%; (4) a touch hidden by rounding, and a
    # collinear retrace, treated as simple profiles; (5) a thin rotated mesh bar
    # whose divergence sum lost 1%; (6) a rotated slender bar refused as
    # unresolved because its world box is loose.
    from fractions import Fraction
    make = lambda shape, volume: inj.DetectorDirected2BodyChannel(shape, 0, inj.DirectedMode.Volume, volume)
    mesh = siren.geometry.TriangularMesh
    cancelling = _plate_with_boss(.002, .004, as_triangles=True) + _box_triangles([.002, .002, .002000005], [.006, .010, .004], inward=True)
    with pytest.raises(RuntimeError, match='does not bound a solid'):
        make(mesh(cancelling), 4e-8)
    with pytest.raises(RuntimeError, match='cannot be sampled consistently'):
        make(mesh(_box_triangles([-1, -1, 0], [1, 1, 1e-8])), 4e-8)
    make(mesh(_box_triangles([-1, -1, 0], [1, 1, 1e-6])), 4e-6)
    annulus = siren.geometry.GenericPolycone([1., 1.000001, 1.000001, 1.], [3000., 3000., 3000.000001, 3000.000001])
    rz = [(Fraction(r), Fraction(z)) for r, z in zip([1., 1.000001, 1.000001, 1.], [3000., 3000., 3000.000001, 3000.000001])]
    moment = sum((rz[i][0]*rz[(i + 1) % 4][1] - rz[(i + 1) % 4][0]*rz[i][1])*(rz[i][0] + rz[(i + 1) % 4][0]) for i in range(4))
    reference = float(abs(moment)/6)*2*math.pi
    assert inj.geometry_volume(annulus) == pytest.approx(reference, rel=1e-9)
    make(annulus, reference)
    for r, z in (([0.1, 0.2, 0.3, 0.4], [0.010000000000000002, 0.020000000000000004, 0.029999999999999995, 0.04000000000000001]),
                 ([1., 3., 5.], [0., 1., 2.])):
        with pytest.raises(ValueError):
            inj.geometry_volume(siren.geometry.GenericPolycone(r, z))
    c, s = math.cos, math.sin
    rotation = (np.array([[c(.9), -s(.9), 0], [s(.9), c(.9), 0], [0, 0, 1]])
                @ np.array([[c(.7), 0, s(.7)], [0, 1, 0], [-s(.7), 0, c(.7)]])
                @ np.array([[1, 0, 0], [0, c(.4), -s(.4)], [0, s(.4), c(.4)]]))
    bar = [[siren.math.Vector3D(*(rotation @ np.array([p.GetX(), p.GetY(), p.GetZ()]))) for p in tri]
           for tri in _box_triangles([-5., -5e-8, -5e-8], [5., 5e-8, 5e-8])]
    exact = Fraction(0)
    for tri in bar:
        a, b, e = [[Fraction(p.GetX()), Fraction(p.GetY()), Fraction(p.GetZ())] for p in tri]
        exact += a[0]*(b[1]*e[2] - b[2]*e[1]) - a[1]*(b[0]*e[2] - b[2]*e[0]) + a[2]*(b[0]*e[1] - b[1]*e[0])
    assert inj.geometry_volume(mesh(bar)) == pytest.approx(float(exact/6), rel=1e-8)
    quaternion = siren.math.Quaternion()
    quaternion.SetEulerAnglesXYZs(.4, .7, .9)
    slender = mesh(siren.geometry.Placement(siren.math.Vector3D(0, 0, 0), quaternion), _box_triangles([-5, -.025, -.025], [5, .025, .025]))
    assert inj.DetectorDirected2BodyChannel(slender, 0, inj.DirectedMode.Volume) is not None


def test_thin_ellipsoid_cap_volume_does_not_cancel():
    a = b = 100.
    c, h = 1e8, 1.
    cap = siren.geometry.Ellipsoid(a, b, c, c - h, c)
    assert inj.geometry_volume(cap) == pytest.approx(math.pi*a*b*h*h/c*(1 - h/(3*c)), rel=1e-12)
    bottom = siren.geometry.Ellipsoid(a, b, c, -c, -c + h)
    assert inj.geometry_volume(bottom) == pytest.approx(inj.geometry_volume(cap), rel=1e-12)


def test_cone_target_channels_round_trip():
    # The geometry Cone shared an include guard with the direction distribution
    # Cone, so its archive registration was skipped in some translation units.
    cone = siren.geometry.Cone(0., 3., 0., 1., 6.)
    channel = inj.DetectorDirected2BodyChannel(cone, 0, inj.DirectedMode.Volume)
    restored = pickle.loads(pickle.dumps(channel))
    r = two_body_record([M*1.5, 0., 0., M*math.sqrt(1.25)])
    channel.Sample(siren.utilities.SIREN_random(8), None, r)
    assert restored.Density(None, r) == channel.Density(None, r)


def _three_body_record(energy_offset):
    gamma = 1000.
    r = siren.dataclasses.InteractionRecord()
    r.signature.primary_type = pt.N4
    r.signature.target_type = pt.Decay
    r.signature.secondary_types = [pt.NuLight, pt.Gamma, pt.NuMu]
    r.primary_mass = M
    r.primary_momentum = [M*gamma*(1 + energy_offset), 0., 0., M*math.sqrt(gamma*gamma - 1)]
    r.secondary_masses = [.03, .01, .005]
    r.secondary_momenta = [[0.]*4 for _ in range(3)]
    r.secondary_helicities = [0.]*3
    r.interaction_vertex = [0.]*3
    return r


@pytest.mark.parametrize('mode', ['direct', 'recursive'])
def test_three_body_channel_uses_the_on_shell_parent(mode):
    # Codex round 5: a 1e-8 energy offset at gamma 1000 gave a 5 MeV daughter a
    # mass of 14 MeV, because the frames mixed the recorded energy and the mass.
    shape = siren.geometry.Sphere(100., 0.)
    channel = (inj.DetectorDirected3BodyChannel(inj.ThreeBodyMode.Direct, shape, 0, mode=inj.DirectedMode.Cone)
               if mode == 'direct' else
               inj.DetectorDirected3BodyChannel(inj.ThreeBodyMode.Recursive, shape, 1, 0, 1, 2,
                                                mode=inj.DirectedMode.Cone))
    rng = siren.utilities.SIREN_random(722)
    for _ in range(200):
        r = _three_body_record(1e-8)
        channel.Sample(rng, None, r)
        p = np.array(r.secondary_momenta)
        masses2 = p[:, 0]**2 - np.sum(p[:, 1:]**2, axis=1)
        np.testing.assert_allclose(masses2, np.array(r.secondary_masses)**2, atol=1e-6*M*M)
        assert channel.Density(None, r) > 0
    r = _three_body_record(1e-3)
    with pytest.raises(siren.utilities.InjectionFailure, match='mass shell'):
        channel.Sample(rng, None, r)
    assert channel.Density(None, r) == 0


def test_translated_box_volume_survives_bounding_box_rounding():
    # Codex: a 1 cm box 100 m from the origin has a bounding box whose volume
    # rounds 2.7e-12 below the exact one. The exact volume needs no box check.
    box = siren.geometry.Box(widths=[.01, .01, .01], center=[100., 100., 100.])
    for channel in (inj.DetectorDirected2BodyChannel(box, 0, inj.DirectedMode.Volume),
                    inj.DetectorDirected2BodyChannel(box, 0, inj.DirectedMode.Volume, 1e-6),
                    inj.OnShellCascadeChannel(box, .03, 0., [.1, 0., .9])):
        assert pickle.loads(pickle.dumps(channel)) is not None


def _with_archived_volume(channel, volume, factor):
    """The channel's pickle with its stored target volume multiplied by factor."""
    state = channel.__getstate__()[0]
    raw = bytearray(bytes.fromhex(state))
    offsets = [i for i in range(len(raw) - 7)
               if abs(struct.unpack('<d', raw[i:i+8])[0]/volume - 1) < 1e-12]
    assert len(offsets) == 1
    stored = struct.unpack('<d', raw[offsets[0]:offsets[0]+8])[0]
    raw[offsets[0]:offsets[0]+8] = struct.pack('<d', factor*stored)
    tampered = raw.hex().upper() if state.isupper() else raw.hex()
    payload = pickle.dumps(channel)
    assert payload.count(state.encode()) == 1
    return payload.replace(state.encode(), tampered.encode())


def test_archived_volume_is_checked_on_load():
    # An archive holding a volume the constructor would refuse (for example one
    # written before the check existed) must not load as a different density.
    sphere = siren.geometry.Sphere(siren.geometry.Placement(siren.math.Vector3D(15, 0, 20)), 8., 0.)
    exact = 4/3*math.pi*8**3
    for channel in (inj.DetectorDirected2BodyChannel(sphere, 0, inj.DirectedMode.Volume),
                    inj.OnShellCascadeChannel(sphere, .03, 0., [.1, 0., .9])):
        assert pickle.loads(_with_archived_volume(channel, exact, 1.)) is not None
        with pytest.raises(RuntimeError, match='Archived directed-channel target volume rejected'):
            pickle.loads(_with_archived_volume(channel, exact, 2.))
    # Cone mode does not use the volume, so its archive is not checked.
    cone = inj.DetectorDirected2BodyChannel(sphere, 0, inj.DirectedMode.Cone)
    assert pickle.loads(_with_archived_volume(cone, exact, 2.)) is not None


def test_measure_fields_are_read_only():
    measure = siren.Measure.OnShellCascade(.03)
    for field, value in [('type', inj.PhaseSpaceMeasureType.SolidAngleRest),
                         ('spectator', 1), ('pair_first', 2), ('pair_second', 1)]:
        with pytest.raises(AttributeError):
            setattr(measure, field, value)
    assert pickle.loads(pickle.dumps(measure)) == measure
    assert {measure: 'stored'}[siren.Measure.OnShellCascade(.03)] == 'stored'


def test_source_table_saves_and_loads_file_objects():
    metadata = {key: {} for key in ['source', 'physics', 'geometry', 'scoring']}
    table = SourceImportanceTable([1., 2.], [1., 3.], metadata)
    stream = io.BytesIO()
    table.save(stream)
    stream.seek(0)
    loaded = SourceImportanceTable.load(stream, metadata)
    np.testing.assert_array_equal(loaded.proposal, table.proposal)


def test_on_shell_rows_canonicalizes_rounded_tables_and_rejects_off_shell_rows():
    keys = ['E', 'm', 'px', 'py', 'pz', 'x', 'y', 'z', 'weight']
    rows = [[*float32_parent(g)[:1], M, 0., 0., float32_parent(g)[3], 0., 0., .3, 1.]
            for g in (2., 300., 3000.)]
    new_keys, new_rows = siren.dist.on_shell_rows(keys, rows)
    assert new_keys == keys + ['E_table']
    for old, new in zip(rows, new_rows):
        assert new[0] == math.hypot(M, math.sqrt(old[4]**2))
        assert new[1:9] == old[1:9] and new[9] == old[0]
    with pytest.raises(ValueError, match='mass shell'):
        siren.dist.on_shell_rows(keys, rows + [[2*M*(1+1e-3), M, 0., 0., math.sqrt(3)*M, 0., 0., .3, 1.]])
    # Repeated columns are ambiguous: the distribution reads the last one.
    with pytest.raises(ValueError, match='more than once'):
        siren.dist.on_shell_rows(keys + ['E'], [row + [99.] for row in rows])
    # |p|^2 overflows: the rebuilt energy is infinite, not a canonical row.
    with pytest.raises(ValueError, match='mass shell'):
        siren.dist.on_shell_rows(keys, [[1e200, 1., 1e200, 0., 0., 0., 0., .3, 1.]])
    with pytest.raises(ValueError, match='values for'):
        siren.dist.on_shell_rows(keys, [rows[0][:-1]])
    # The saved energy must not land in a column the distribution interprets.
    plain = ['E', 'm', 'px', 'py', 'pz']
    for name in ('weight', 't0', 'x0', 'PrimaryExternalDistribution_row', '', 3):
        with pytest.raises(ValueError, match='interprets|column name'):
            siren.dist.on_shell_rows(plain, [row[:5] for row in rows], keep_input_energy=name)



def _closure(sampler, density, reference, seed, n=250000, observable=None):
    """density(x)/reference(x) for x from sampler, with observable(x) if given."""
    r = record(2., .01)
    rng = siren.utilities.SIREN_random(seed)
    values = np.empty(n)
    seen = np.empty(n)
    for i in range(n):
        sampler.Sample(rng, None, r)
        values[i] = density(r)/reference(r)
        if observable:
            seen[i] = observable(r)
    return (values, seen) if observable else values


def _resolved(values):
    """Whether the mean of values is more than 3.5 standard errors from zero."""
    return abs(values.mean()) > 3.5*values.std(ddof=1)/math.sqrt(len(values))


def _directed_share(weights):
    return 2 if weights[2] else 1


@pytest.mark.parametrize('weights', [[.1, 0., .9], [.1, .9, 0.]])
def test_production_mixture_forward_closure_resolves_bias_and_mixing_errors(weights):
    # E_q[p/q] = 1 for the production mixture q (kappa = 0, so p is isotropic).
    # It resolves a 2% total bias and, as checked below, drawing the directed
    # component 6% less often than the density assumes. A density error confined
    # to the directed component moves it only by that component's share of the
    # physical probability (Codex: 1.3% for a 6% volume error); the reverse test
    # below covers that case.
    channel = inj.OnShellCascadeChannel(target(), .03, 0., weights)
    physical = 1/(16*math.pi**2)
    values = _closure(channel, lambda r: physical, lambda r: channel.Density(None, r), 2718)
    assert values.std(ddof=1)/math.sqrt(len(values)) < .02/3.5
    assert not _resolved(values - 1)
    # Drawing the directed component 6% less often than the density assumes.
    wrong = list(weights)
    wrong[_directed_share(weights)] *= .94
    wrong[0] = 1 - wrong[_directed_share(weights)]
    sampler = inj.OnShellCascadeChannel(target(), .03, 0., wrong)
    bad = _closure(sampler, lambda r: physical, lambda r: channel.Density(None, r), 2718, 50000)
    assert _resolved(bad - 1)


@pytest.mark.parametrize('weights', [[.1, 0., .9], [.1, .9, 0.]])
def test_production_mixture_density_integrates_to_one_component_by_component(weights):
    # Reverse closure: sample the physical law (the isotropic-only channel, exact
    # for kappa = 0) and average q/p, the integral of the mixture density. A
    # normalization error in one component's density moves it by that
    # component's weight times the error, whatever its share of the physical
    # probability. A shape error that keeps the integral is invisible to it. The
    # y-reflection moment below catches an error that moves the mean sign of
    # p_y (as Codex's mutation does); one moment cannot certify every shape, and
    # the angular distribution tests cover others.
    law = inj.OnShellCascadeChannel(target(), .03, 0., [1., 0., 0.])
    mixture = inj.OnShellCascadeChannel(target(), .03, 0., weights)
    values, side = _closure(law, lambda r: mixture.Density(None, r), lambda r: law.Density(None, r),
                            3141, observable=lambda r: math.copysign(1., r.secondary_momenta[1][2]))
    assert not _resolved(values - 1)
    # Target and parent are symmetric under y -> -y, so the integral of q times
    # the sign of daughter 1's p_y vanishes.
    assert not _resolved(values*side)
    # Deliberate errors on the same samples, with q/p = w_iso + w_d*q_d/p: the
    # directed density 6% low, or reshaped by 1 + 0.06*sign(p_y) (Codex's
    # integral-preserving mutation).
    directed = values - weights[0]
    assert _resolved(weights[0] + .94*directed - 1)
    assert _resolved((weights[0] + directed*(1 + .06*side))*side)
