from types import SimpleNamespace
import numpy as np
import pytest
import siren
pytest.importorskip('DarkNews')


def test_darknews_sampleps_returns_samples_without_return_flags(monkeypatch):
    mod=siren._util.import_resource('processes','DarkNewsTables')
    module=siren._util.load_module('siren.resources.processes.DarkNewsTables.DarkNewsDecay',
        str(__import__('pathlib').Path(mod.__file__).with_name('DarkNewsDecay.py')))
    calls=[]
    def sample(*,existing_integrator=None):
        calls.append(existing_integrator)
        return np.array([[0.1]]),{'diff_decay_rate_0':np.array([1.])}
    decay=module.PyDarkNewsDecay(SimpleNamespace(SamplePS=sample))
    monkeypatch.setattr(module,'get_decay_momenta_from_vegas_samples',lambda *args:{})
    record=SimpleNamespace(primary_momentum=[1.,0.,0.,1.],get_secondary_particle_records=lambda:[])
    random=SimpleNamespace(Uniform=lambda low,high:low)
    decay.SampleRecordFromDarkNews(record,random)
    assert calls==[None]
    decay.SampleRecordFromDarkNews(record,random)
    assert calls==[None]


def test_darknews_uncached_total_width_uses_scalar_api():
    mod=siren._util.import_resource('processes','DarkNewsTables')
    from DarkNews.processes import FermionDileptonDecay
    case=object.__new__(FermionDileptonDecay)
    case.nu_parent=5914
    case.vector_off_shell=case.scalar_off_shell=True
    case.total_width=lambda:0.123
    decay=mod.PyDarkNewsDecay(case)
    assert decay.TotalDecayWidthAllFinalStates(siren.particles.N4)==0.123
