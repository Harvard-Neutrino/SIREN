"""Expansion control replacement and persistence contracts."""
import pickle
import pytest
import siren
from test_injector_api import _chain_injector, _depth_ge_1

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
def test_raw_v2_load_replaces_control_policy(tmp_path, archived_rules):
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
