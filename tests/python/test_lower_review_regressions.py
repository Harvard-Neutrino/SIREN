"""Control-state, constrained-measure and caller-metadata review regressions."""
import math
import pickle

import pytest
import siren

from test_injector_api import _chain_injector, _depth_ge_1
from test_on_shell_cascade import record, target
from test_phase_space_decay import decay, make_record


def _rules():
    pt = siren.particles.NuMu
    return siren.injection.SecondaryExpansion([[int(pt), int(pt), -1, 1, 0]])


def test_none_resets_built_injector_to_primary_only():
    inj = _chain_injector()
    assert all(len(t.tree) == 2 for t in inj.generate(3, on_shortfall="raise"))
    inj.stopping_condition = None
    assert inj.stopping_condition is None
    assert all(len(t.tree) == 1 for t in inj.generate(3, on_shortfall="raise"))


def test_native_rules_discard_previous_callback_and_wrapper_tracks_clear():
    inj = _chain_injector()
    inj._build()
    calls = []
    inj.stopping_condition = lambda *args: calls.append(args) or False
    inj.engine.SetSecondaryExpansion(_rules())
    assert inj.stopping_condition == _rules()
    assert all(len(t.tree) == 2 for t in inj.generate(3, on_shortfall="raise"))
    inj.engine.SetSecondaryExpansion(None)
    assert inj.stopping_condition is None
    assert all(len(t.tree) == 1 for t in inj.generate(3, on_shortfall="raise"))
    assert not calls
    assert pickle.loads(pickle.dumps(inj)).stopping_condition is None


@pytest.mark.parametrize("archived_rules", [False, True])
def test_raw_v3_load_replaces_control_policy(tmp_path, archived_rules):
    saved = _chain_injector()
    saved._build()
    saved.stopping_condition = _rules() if archived_rules else None
    path = str(tmp_path / "policy")
    saved.save(path)
    inj = _chain_injector()
    inj._build()
    calls = []
    inj.stopping_condition = lambda *args: calls.append(args) or False
    inj.engine.LoadInjector(path)
    assert inj.stopping_condition == (_rules() if archived_rules else None)
    assert all(len(t.tree) == (2 if archived_rules else 1)
               for t in inj.generate(3, on_shortfall="raise"))
    assert not calls


def test_raw_callback_replacement_is_visible_and_blocks_native_save(tmp_path):
    inj = _chain_injector()
    inj._build()
    inj.stopping_condition = _rules()
    inj.engine.SetStoppingCondition(_depth_ge_1)
    assert not isinstance(inj.stopping_condition, siren.injection.SecondaryExpansion)
    with pytest.raises(siren.errors.NotSerializableError):
        inj.save(str(tmp_path / "callback"))
    assert all(len(t.tree) == 2 for t in inj.generate(3, on_shortfall="raise"))


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


def test_native_sampler_preserves_pending_caller_parameters():
    model = decay()
    r = make_record(model)
    r.interaction_parameters = {"source": 8., "override": 1.}
    output = siren.dataclasses.CrossSectionDistributionRecord(r)
    output.interaction_parameters = {"source": 8., "override": 2., "caller": 17.}
    model.SampleFinalState(output, siren.utilities.SIREN_random(31))
    assert output.interaction_parameters == {"source": 8., "override": 2., "caller": 17.}


def test_unequal_masses_are_configuration_errors_and_classified_direct_failures():
    r = record()
    r.secondary_masses = [0., .01, math.nextafter(.01, math.inf)]
    channel = siren.injection.OnShellCascadeChannel(target(), .03)
    with pytest.raises(ValueError, match="exactly equal pair masses"):
        siren.injection.PhaseSpaceDecay(r.signature, r.secondary_masses, 1., 1., channel)
    assert channel.Density(None, r) == 0
    with pytest.raises(siren.utilities.InjectionFailure):
        channel.Sample(siren.utilities.SIREN_random(3), None, r)


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
