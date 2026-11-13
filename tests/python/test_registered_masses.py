"""Registered masses reach authoring defaults through Python and native callers."""
import pytest
import siren
from siren import particles, models
from siren.errors import ConfigurationError

BASES = [
    siren.DecayModel,
    siren.CrossSectionModel,
    models.decay_model_base(siren.interactions.DarkNewsDecay),
    models.cross_section_model_base(siren.interactions.DarkNewsCrossSection),
]


@pytest.mark.parametrize('base', BASES)
def test_registered_and_standard_daughters(base):
    chi = particles.define('MassTestChi', 9010091, 0.125)
    daughters = [chi, particles.Electron, particles.NuMu]
    model = base()
    expected = [0.125, siren.dataclasses.GetParticleMass(particles.Electron), 0.0]
    assert model.SecondaryMasses(daughters) == expected
    root = (siren.interactions.Decay if isinstance(model, siren.interactions.Decay)
            else siren.interactions.CrossSection)
    assert root.SecondaryMasses(model, daughters) == expected
    assert particles.mass('MassTestChi') == 0.125
    assert particles.mass('Electron') == expected[1]


@pytest.mark.parametrize('base', BASES)
def test_explicit_mass_and_helicity_overrides(base):
    class Model(base):
        def SecondaryMasses(self, types):
            return [0.3 for _ in types]

        def SecondaryHelicities(self, record):
            return [0.5 for _ in record.signature.secondary_types]

    model = Model()
    record = siren.dataclasses.InteractionRecord()
    sig = record.signature
    sig.secondary_types = [siren.dataclasses.ParticleType(9010092)]
    record.signature = sig
    root = (siren.interactions.Decay if isinstance(model, siren.interactions.Decay)
            else siren.interactions.CrossSection)
    assert root.SecondaryMasses(model, sig.secondary_types) == [0.3]
    assert root.SecondaryHelicities(model, record) == [0.5]


@pytest.mark.parametrize('base', BASES)
def test_unknown_mass_has_actionable_error(base):
    with pytest.raises(ConfigurationError, match='9010093.*particles.define.*SecondaryMasses'):
        base().SecondaryMasses([siren.dataclasses.ParticleType(9010093)])


@pytest.mark.parametrize('mass', [-0.1, float('nan'), float('inf')])
def test_invalid_mass_does_not_register(mass):
    with pytest.raises(ConfigurationError, match='finite and non-negative'):
        particles.define('InvalidMassTest', 9010094, mass)
    with pytest.raises(ValueError, match='Unknown particle'):
        particles.resolve('InvalidMassTest')


def test_builtin_aliases_cannot_register_conflicting_masses():
    # Restore metadata so this test cannot change a later model's electron mass.
    before = dict(particles._name_to_mass)
    try:
        electron_mass = particles.mass('Electron')
        particles.define('Electron', int(particles.Electron), electron_mass)
        with pytest.raises(ConfigurationError, match='already has mass'):
            particles.define('EMinus', int(particles.Electron), 1.0)
    finally:
        particles._name_to_mass.clear()
        particles._name_to_mass.update(before)


@pytest.mark.parametrize('name', ['mass', 'define', '_internal', 'not a name'])
def test_registration_cannot_replace_registry_functions(name):
    with pytest.raises(ConfigurationError, match='unreserved public identifier'):
        particles.define(name, 9010095, 1.0)
    assert callable(particles.mass) and callable(particles.define)


@pytest.mark.parametrize('base', [siren.DecayModel, models.decay_model_base(siren.interactions.DarkNewsDecay)])
def test_native_sampling_uses_registered_daughter_masses(base):
    import math
    import numpy as np
    from test_authoring_bases import _template_record

    chi = particles.define('SampleMassChi', 9010096, 0.125)

    class Decay(base):
        parent = 'N4'
        daughters = (chi, 'Electron')
        measure = siren.Measure.SolidAngleRest()

        def total_width(self):
            return 1.0

        def differential_width(self, record):
            return 1.0 / (4 * math.pi)

        def sample(self, record, random):
            self.sample_isotropic(record, random)

    model = Decay()
    source = _template_record(model.GetPossibleSignatures()[0], energy=2, primary_mass=1)
    source.secondary_masses = []
    record = siren.dataclasses.CrossSectionDistributionRecord(source)
    siren.interactions.Decay.SampleFinalState(model, record, siren.utilities.SIREN_random(12))
    record.finalize(source)
    assert source.secondary_masses == [0.125, particles.mass('Electron')]
    momenta = np.array(source.secondary_momenta)
    np.testing.assert_allclose(momenta.sum(axis=0), source.primary_momentum, atol=1e-14)
    np.testing.assert_allclose(momenta[:, 0] ** 2 - (momenta[:, 1:] ** 2).sum(axis=1),
                               np.square(source.secondary_masses), atol=1e-14)
