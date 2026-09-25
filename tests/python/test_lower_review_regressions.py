"""Control-state, constrained-measure and caller-metadata review regressions."""
import pickle

import pytest
import siren

from test_on_shell_cascade import record, target


@pytest.mark.parametrize("mass", [float('nan'), -0.1, float('inf')])
def test_pair_mass_cannot_bypass_validated_factory(mass):
    measure = siren.Measure.OnShellCascade(.03)
    with pytest.raises(AttributeError):
        measure.pair_mass = mass
    with pytest.raises(ValueError):
        siren.Measure.OnShellCascade(mass)
    assert measure == measure


@pytest.mark.parametrize("measure", [siren.Measure.OnShellCascade(.03),
    siren.Measure.OnShellCascade(.06, 2, 0, 1), siren.Measure.SolidAngleRest(),
    siren.Measure.Recursive2Body(1, 0, 2), siren.Measure.SolidAngleLab(1),
    siren.Measure.Unspecified()])
def test_measure_pickle_preserves_constraint_and_dictionary_identity(measure):
    copy = pickle.loads(pickle.dumps(measure))
    assert copy == measure
    assert hash(copy) == hash(measure)
    assert copy.pair_mass == measure.pair_mass
    assert {measure: 1}[copy] == 1


class _CascadeModel(siren.DecayModel):
    parent = 'Pi0'
    daughters = ('Gamma', 'N4', 'N5')
    def __init__(self, measure, masses):
        super().__init__()
        self.measure = measure
        self.masses = masses
    def SecondaryMasses(self, types): return self.masses
    def total_width(self): return 1.
    def differential_width(self, record): return 1.
    def sample(self, record, random): pass


@pytest.mark.parametrize("measure,masses", [
    (siren.Measure.OnShellCascade(.03), [.01, .01, 0.]),
    (siren.Measure.OnShellCascade(.03, 2, 0, 1), [.01]*3),
    (siren.Measure.OnShellCascade(.06), [0., .01, .01]),
    (siren.Measure.Recursive2Body(), [0., .01, .01]),
])
def test_facade_rejects_incompatible_model_order_mass_or_measure(measure, masses):
    model = _CascadeModel(measure, masses)
    proposal = 1 * siren.channels.on_shell_cascade(target(), .03)
    with pytest.raises(siren.errors.ConfigurationError, match="on_shell_cascade"):
        proposal._build(model.GetPossibleSignatures()[0], models=[model])


def test_named_selectors_check_order_even_when_all_masses_equal():
    r = record()
    channel = 1 * siren.channels.on_shell_cascade(
        target(), .03, spectator='N5', pair=('Gamma', 'N4'))
    with pytest.raises(siren.errors.ConfigurationError, match="reorder"):
        channel._build(r.signature)
    # With all three masses equal, masses cannot reveal the order; only the
    # named selectors can, and an unsupported permutation must still fail.
    equal = _CascadeModel(siren.Measure.OnShellCascade(.03), [.01, .01, .01])
    with pytest.raises(siren.errors.ConfigurationError, match="reorder"):
        channel._build(r.signature, models=[equal])
    accepted = 1 * siren.channels.on_shell_cascade(
        target(), .03, spectator='Gamma', pair=('N4', 'N5'))
    assert accepted._build(r.signature, models=[equal]).Measure() == equal.Measure()
    model = _CascadeModel(siren.Measure.OnShellCascade(.03), [0., .01, .01])
    channel = 1 * siren.channels.on_shell_cascade(
        target(), .03, spectator='Gamma', pair=('N4', 'N5'))
    assert channel._build(r.signature, models=[model]).Measure() == model.Measure()
