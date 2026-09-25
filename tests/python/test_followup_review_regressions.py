"""Regressions for the 2026-09-24 follow-up review and its round-2 audit."""
import math
import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest
import siren

pt = siren.dataclasses.ParticleType


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


