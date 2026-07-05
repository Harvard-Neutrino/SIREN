"""Tests for HepMC3/NuHepMC event export (siren.io).

Exercises both the direct ``siren.io.SaveInteractionTreesAsHepMC3`` entry point
and the ``siren._util.SaveEvents(save_hepmc3=True)`` wiring. The tests skip
cleanly when SIREN was built without HepMC3 support (the bound functions raise a
RuntimeError in that configuration).
"""
import os

import pytest

siren = pytest.importorskip("siren")


@pytest.fixture(scope="module")
def dc():
    from siren import dataclasses
    return dataclasses


def _make_tree(dc):
    record = dc.InteractionRecord()
    sig = record.signature
    sig.primary_type = dc.ParticleType.NuMu
    sig.target_type = dc.ParticleType.EPlus
    sig.secondary_types = [dc.ParticleType.NuMu]
    record.signature = sig
    record.primary_momentum = [10.0, 0.0, 0.0, 10.0]
    record.target_mass = 0.000511
    record.interaction_vertex = [1.0, 2.0, 3.0]  # internal meters
    record.interaction_time = 5.0
    record.secondary_ids = [dc.ParticleID(1, 1)]
    record.secondary_momenta = [[10.0, 0.0, 0.0, 10.0]]
    record.secondary_masses = [0.0]
    record.secondary_helicities = [-1.0]
    record.secondary_times = [5.0]
    record.interaction_parameters = {"Q2": 0.7}

    tree = dc.InteractionTree()
    tree.add_entry(record, None)
    tree.header.event_number = 7
    tree.header.weights = [3.0]
    return tree


def _run_or_skip(callable_):
    try:
        callable_()
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise


def test_io_save_hepmc3(dc, tmp_path):
    from siren import hepmc3 as sio

    tree = _make_tree(dc)
    out = str(tmp_path / "events.hepmc3")
    _run_or_skip(lambda: sio.SaveInteractionTreesAsHepMC3([tree], out))

    with open(out) as f:
        text = f.read()
    assert "HepMC::Version" in text          # valid HepMC3 Ascii header
    assert "U GEV CM" in text                # units line
    assert text.count("\nE ") >= 1 or text.startswith("E ") or "\nE 7 " in text
    assert "siren.param.Q2" in text          # reweight-critical attribute
    # internal-meters -> CM conversion: vertex {1,2,3} m -> {100,200,300}
    assert "1.0000000000000000e+02" in text


def test_roundtrip_write_read(dc, tmp_path):
    from siren import hepmc3 as sio

    tree = _make_tree(dc)
    out = str(tmp_path / "roundtrip.hepmc3")
    _run_or_skip(lambda: sio.SaveInteractionTreesAsHepMC3([tree], out))

    loaded = sio.LoadInteractionTreesFromHepMC3(out)
    assert len(loaded) == 1
    lt = loaded[0]
    assert len(lt.tree) == 1
    assert lt.header.event_number == 7
    assert list(lt.header.weights) == [3.0]

    rec = lt.tree[0].record
    assert rec.signature.primary_type == dc.ParticleType.NuMu
    assert list(rec.primary_momentum) == pytest.approx([10.0, 0.0, 0.0, 10.0], abs=1e-6)
    assert list(rec.interaction_vertex) == pytest.approx([1.0, 2.0, 3.0], abs=1e-6)
    assert rec.interaction_parameters["Q2"] == pytest.approx(0.7, abs=1e-6)


def test_util_saveevents_hepmc3(dc, tmp_path):
    from siren import _util

    tree = _make_tree(dc)
    out = str(tmp_path / "run")
    _run_or_skip(lambda: _util.SaveEvents(
        [tree],
        weighter=None,
        gen_times=[0.0],
        save_hdf5=False,
        save_parquet=False,
        save_siren_events=False,
        save_hepmc3=True,
        fid_vol=None,
        output_filename=out))

    assert os.path.exists(out + ".hepmc3")
