"""Tables retain computed decay state and reload from the directory saved."""
import pickle
from types import SimpleNamespace

import numpy as np
import pytest
import siren

pytest.importorskip("DarkNews")


def resource():
    return siren._util.import_resource("processes", "DarkNewsTables")


def test_decay_table_restores_current_cache_without_recomputation(tmp_path):
    module = resource()
    from DarkNews.processes import FermionDileptonDecay
    from DarkNews import pdg
    case = object.__new__(FermionDileptonDecay)
    case.nu_parent = pdg.neutrino4
    case.nu_daughter = pdg.nulight
    case.secondaries = [pdg.electron, pdg.positron]
    case.vector_off_shell = case.scalar_off_shell = True
    case.total_width = lambda: pytest.fail("cached width was recomputed")
    decay = module.PyDarkNewsDecay(case)
    decay.total_width = 0.123
    decay.PS_samples = np.array([[0.1, 0.7]])
    decay.PS_weights = np.array([1., 2.])
    decay.save_to_table(tmp_path)
    restored = module.PyDarkNewsDecay(case)
    restored.load_from_table(tmp_path)
    assert restored.TotalDecayWidthAllFinalStates(siren.particles.N4) == 0.123
    record = siren.dataclasses.InteractionRecord()
    record.signature = restored.GetPossibleSignatures()[0]
    assert restored.TotalDecayWidth(record) == 0.123
    np.testing.assert_array_equal(restored.PS_samples, decay.PS_samples)
    np.testing.assert_array_equal(restored.PS_weights, decay.PS_weights)
    np.testing.assert_array_equal(restored.GetPSSample(SimpleNamespace(Uniform=lambda a, b: 0.)), [0.1])


def test_decay_table_reads_legacy_mapping_by_key(tmp_path):
    state = {"decay_norm": {"diff_decay_rate_0": 2.}, "decay_integrator": None}
    with (tmp_path / "decay.pkl").open("wb") as stream:
        pickle.dump(state, stream)
    decay = resource().PyDarkNewsDecay(SimpleNamespace())
    decay.load_from_table(tmp_path)
    assert decay.decay_norm == state["decay_norm"]
    assert decay.decay_integrator is None


@pytest.mark.parametrize("preference", ["table", "pickle"])
def test_decay_loader_uses_one_subdirectory(tmp_path, preference):
    module = resource()
    key = ("N4", "nu", "photon")
    directory = tmp_path / "Decay_N4_nu_photon"
    directory.mkdir()
    case = SimpleNamespace(marker=17)
    decay = module.PyDarkNewsDecay(case)
    decay.total_width = 0.25
    decay.save_to_table(directory)
    with (directory / "dec_object.pkl").open("wb") as stream:
        pickle.dump(decay, stream)
    restored = module.attempt_to_load_decay(
        SimpleNamespace(dec_cases={key: case}), key, str(tmp_path), [preference])
    assert restored.total_width == 0.25
    assert not (directory / directory.name).exists()
