"""Equality and sampling contracts shared by native and Python entry points."""

import math
import types

import pytest
import siren
from siren import _validation, models
from siren.errors import ConfigurationError
from test_authoring_bases import _LegacyBrokenDecay, _template_csdr
from test_vertex_spec import _KeepAliveCrossSection


@pytest.mark.parametrize("base", [
    siren.DecayModel, siren.CrossSectionModel,
    models.decay_model_base(siren.interactions.DarkNewsDecay),
    models.cross_section_model_base(siren.interactions.DarkNewsCrossSection),
])
@pytest.mark.parametrize("measure", [
    siren.Measure.Recursive2Body(2, 0, 1),
    siren.Measure.DalitzPair(1, 2, 0),
])
def test_selected_signature_preserves_declared_measure(base, measure):
    class Model(base):
        parent = primary = "N4"
        target = "PPlus"
        daughters = finals = ("NuLight", "EMinus", "EPlus")

    model = Model()
    model.measure = measure
    signature = model.GetPossibleSignatures()[0]
    assert model.DensityVariables() == []
    assert model.MeasureForSignature(signature) == measure
    assert model.TopologyForSignature(signature) == model.Topology()
    # Unrelated signatures must not inherit this model's coordinate convention.
    signature.primary_type = siren.dataclasses.Particle.ParticleType.Gamma
    assert model.MeasureForSignature(signature) == siren.Measure.Unspecified()
    assert model.TopologyForSignature(signature) == siren.Topology.Unspecified


@pytest.mark.parametrize("base", [
    siren.DecayModel, siren.CrossSectionModel,
    models.decay_model_base(siren.interactions.DarkNewsDecay),
    models.cross_section_model_base(siren.interactions.DarkNewsCrossSection),
])
def test_advertised_additional_signatures_keep_declared_measure(base):
    # A subclass may advertise more channels than the convenience attributes
    # build; each keeps the declared measure and its own final-state arity.
    P = siren.dataclasses.Particle.ParticleType
    is_decay = issubclass(base, siren.interactions.Decay)

    class Model(base):
        parent = primary = "N4"
        target = "PPlus"
        daughters = finals = ("NuLight", "EMinus", "EPlus")
        measure = siren.Measure.Recursive2Body(2, 0, 1)

        def GetPossibleSignatures(self):
            first = self._signature()
            second = siren.dataclasses.InteractionSignature()
            second.primary_type = first.primary_type
            second.target_type = first.target_type
            if is_decay:
                second.secondary_types = [P.NuLight, P.Gamma]
            else:
                second.target_type = P.Neutron
                second.secondary_types = list(first.secondary_types)
            return [first, second]

    model = Model()
    first, second = model.GetPossibleSignatures()
    assert second != first
    for signature in (first, second):
        assert model.MeasureForSignature(signature) == model.measure
    assert model.TopologyForSignature(first) == model.Topology()
    if is_decay:
        assert model.TopologyForSignature(second) == siren.Topology.Decay2Body
    else:
        assert model.TopologyForSignature(second) == siren.Topology.Scatter2to3
    second.primary_type = P.Gamma
    assert model.MeasureForSignature(second) == siren.Measure.Unspecified()
    assert model.TopologyForSignature(second) == siren.Topology.Unspecified


@pytest.mark.parametrize("base,root", [
    (siren.DecayModel, siren.interactions.Decay),
    (siren.CrossSectionModel, siren.interactions.CrossSection),
    (models.decay_model_base(siren.interactions.DarkNewsDecay), siren.interactions.Decay),
    (models.cross_section_model_base(siren.interactions.DarkNewsCrossSection), siren.interactions.CrossSection),
])
def test_equality_defaults_and_explicit_overrides_through_cpp(base, root):
    a, b = base(), base()
    assert a.equal(a) and not a.equal(b)
    assert root.__eq__(a, a) and not root.__eq__(a, b)

    class ByValue(base):
        def equal(self, other):
            return isinstance(other, ByValue) and self.value == other.value

    a, b = ByValue(), ByValue()
    a.value = b.value = 7
    assert root.__eq__(a, b)
    b.value = 8
    assert not root.__eq__(a, b)


@pytest.mark.parametrize("base", [siren.DecayModel, siren.CrossSectionModel])
def test_equality_typo_is_rejected(base):
    with pytest.raises(ConfigurationError, match="equal"):
        class Typo(base):
            def equals(self, other):
                return self is other


def test_direct_decay_without_equal_is_rejected_before_native_comparison():
    methods = {name: value for name, value in vars(_LegacyBrokenDecay).items()
               if name != "equal" and (not name.startswith("__") or name == "__init__")}
    methods["DifferentialDecayWidth"] = lambda self, record: 1.0
    RawDecay = type("RawDecay", (siren.interactions.Decay,), methods)
    with pytest.raises(ConfigurationError, match="required method 'equal'"):
        _validation.audit_overrides([RawDecay()])


@pytest.mark.parametrize("skew", [0.0, 0.8])
@pytest.mark.parametrize("base", [siren.DecayModel, siren.CrossSectionModel])
def test_a_measure_does_not_supply_a_sampler(base, skew):
    class MissingSampler(base):
        parent = primary = "N4"
        target = "PPlus"
        daughters = finals = ("NuLight", "Gamma")
        measure = siren.Measure.SolidAngleRest()

        def total_width(self):
            return 1.0

        def differential_width(self, record):
            return (1 + skew) / (4 * math.pi)

        def total_xs(self, record):
            return 1.0

        def differential_xs(self, record):
            return (1 + skew) / (4 * math.pi)

    model = MissingSampler()
    record, source = _template_csdr(model.GetPossibleSignatures()[0])
    expected = model.FinalStateProbability(source)
    with pytest.raises(ConfigurationError, match="implement sample"):
        _validation.audit_overrides([model])
    with pytest.raises(ConfigurationError, match="implement sample"):
        model.SampleFinalState(record, siren.utilities.SIREN_random(0))
    assert model.FinalStateProbability(source) == expected


def test_explicit_native_entry_point_is_accepted_on_factory_base():
    Base = models.decay_model_base(siren.interactions.DarkNewsDecay)

    class NativeEntry(Base):
        parent = "N4"
        daughters = ("NuLight", "Gamma")
        measure = siren.Measure.SolidAngleRest()

        def total_width(self):
            return 1.0

        def differential_width(self, record):
            return 1.0 / (4 * math.pi)

        def SampleFinalState(self, record, random):
            self.sample_isotropic(record, random)

    model = NativeEntry()
    _validation.audit_overrides([model])
    record, source = _template_csdr(model.GetPossibleSignatures()[0])
    model.SampleFinalState(record, siren.utilities.SIREN_random(0))
    record.finalize(source)
    assert source.secondary_momenta[0][0] > 0


def test_direct_cross_section_without_equal_is_rejected():
    methods = {name: value for name, value in vars(_KeepAliveCrossSection).items()
               if name != "equal" and not name.startswith("__")}
    RawXS = type("RawXS", (siren.interactions.CrossSection,), methods)
    with pytest.raises(ConfigurationError, match="required method 'equal'"):
        _validation.audit_overrides([RawXS()])


@pytest.mark.parametrize('native,factory', [
    (siren.interactions.DarkNewsDecay, models.decay_model_base),
    (siren.interactions.DarkNewsCrossSection, models.cross_section_model_base),
])
def test_darknews_default_equality_is_symmetric_across_entry_points(native, factory):
    class Legacy(native):
        pass

    Author = factory(native)
    objects = [native(), native(), Legacy(), Legacy(), Author(), Author()]
    for left in objects:
        for right in objects:
            assert (left == right) == (left is right)
            assert native.__eq__(left, right) == (left is right)
            assert left.equal(right) == (left is right)
    # Python overrides on legacy subclasses must still reach C++ comparisons.
    class ByValue(Legacy):
        def equal(self, other):
            return isinstance(other, ByValue) and self.value == other.value

    left, right = ByValue(), ByValue()
    left.value = right.value = 7
    assert native.__eq__(left, right) and native.__eq__(right, left)
    right.value = 8
    assert not native.__eq__(left, right)


@pytest.mark.parametrize('base', [
    siren.DecayModel, siren.CrossSectionModel,
    models.decay_model_base(siren.interactions.DarkNewsDecay),
    models.cross_section_model_base(siren.interactions.DarkNewsCrossSection),
])
@pytest.mark.parametrize('entry', ['default', 'class', 'bound', 'lambda', 'native', 'noncallable'])
def test_sampler_audit_matches_runtime_method_resolution(base, entry):
    class Model(base):
        parent = primary = 'N4'
        target = 'PPlus'
        daughters = finals = ('NuLight', 'Gamma')
        measure = siren.Measure.SolidAngleRest()

        def total_width(self):
            return 1.0

        def differential_width(self, record):
            return 1 / (4*math.pi)

        def total_xs(self, record):
            return 1.0

        def differential_xs(self, record):
            return 1 / (4*math.pi)

    def sample(self, record, random):
        if isinstance(self, siren.interactions.Decay):
            self.sample_isotropic(record, random)
        else:
            for i, secondary in enumerate(record.get_secondary_particle_records()):
                secondary.mass = 0.0
                secondary.four_momentum = [0.01, 0, 0, 0.01 if i == 0 else -0.01]

    if entry == 'class':
        Model = type('ClassSampler', (Model,), {'sample': sample})
    model = Model()
    if entry == 'bound':
        model.sample = (model.sample_isotropic if isinstance(model, siren.interactions.Decay)
                        else types.MethodType(sample, model))
    elif entry == 'lambda':
        model.sample = lambda record, random: sample(model, record, random)
    elif entry == 'native':
        model.SampleFinalState = types.MethodType(sample, model)
    elif entry == 'noncallable':
        model.sample = None

    if entry in ('default', 'noncallable'):
        with pytest.raises(ConfigurationError, match='implement sample'):
            _validation.audit_overrides([model])
        return
    _validation.audit_overrides([model])
    record, source = _template_csdr(model.GetPossibleSignatures()[0])
    root = (siren.interactions.Decay if isinstance(model, siren.interactions.Decay)
            else siren.interactions.CrossSection)
    root.SampleFinalState(model, record, siren.utilities.SIREN_random(0))
    record.finalize(source)
    assert all(momentum[0] > 0 for momentum in source.secondary_momenta)


@pytest.mark.parametrize('masses,message', [
    ([0.0], 'returned 1 masses for 2 secondary types'),
    ([0.0, object()], 'must return numeric masses'),
])
def test_isotropic_mass_errors_name_the_model_and_hook(masses, message):
    class InvalidMasses(siren.DecayModel):
        parent = 'N4'
        daughters = ('NuLight', 'Gamma')
        measure = siren.Measure.SolidAngleRest()

        def SecondaryMasses(self, types):
            return masses

    model = InvalidMasses()
    _, source = _template_csdr(model.GetPossibleSignatures()[0])
    source.secondary_masses = []
    record = siren.dataclasses.CrossSectionDistributionRecord(source)
    with pytest.raises(ConfigurationError, match='InvalidMasses.SecondaryMasses.*' + message):
        model.sample_isotropic(record, siren.utilities.SIREN_random(0))


def test_isotropic_mass_hook_retains_float_conversion():
    class ConvertibleMasses(siren.DecayModel):
        parent = 'N4'
        daughters = ('NuLight', 'Gamma')
        measure = siren.Measure.SolidAngleRest()

        def SecondaryMasses(self, types):
            return ['0.003', '0.004']

    model = ConvertibleMasses()
    _, source = _template_csdr(model.GetPossibleSignatures()[0])
    source.secondary_masses = []
    record = siren.dataclasses.CrossSectionDistributionRecord(source)
    model.sample_isotropic(record, siren.utilities.SIREN_random(0))
    record.finalize(source)
    assert list(source.secondary_masses) == [0.003, 0.004]
