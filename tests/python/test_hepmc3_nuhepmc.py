"""NuHepMC compliance + external-reader validation for the SIREN HepMC3 export.

Writes a two-level cascade (so a non-root interaction carries a real target) and
checks the NuHepMC metadata the writer now emits: the vertex/particle status
registries (G.R.9/G.R.10), the additional-particle-number registry (G.R.11), the
per-event lab position (E.R.5/E.C.5), and the one-beam/one-target invariant
(E.R.7). The file is re-read with an independent tool (pyhepmc) to confirm the
output is genuinely readable outside SIREN, not just by SIREN's own reader.

Skips cleanly when SIREN was built without HepMC3 support or when pyhepmc is not
installed.
"""
import pytest

siren = pytest.importorskip("siren")


@pytest.fixture(scope="module")
def dc():
    from siren import dataclasses
    return dataclasses


def _record(dc, primary_type, primary_id, target_type, target_mass,
            secondary_types, secondary_ids, vertex, time):
    r = dc.InteractionRecord()
    sig = r.signature
    sig.primary_type = primary_type
    sig.target_type = target_type
    sig.secondary_types = list(secondary_types)
    r.signature = sig
    r.primary_momentum = [10.0, 0.0, 0.0, 10.0]
    r.primary_id = primary_id
    r.target_mass = target_mass
    r.interaction_vertex = list(vertex)
    r.interaction_time = time
    r.secondary_ids = list(secondary_ids)
    r.secondary_momenta = [[5.0, 0.0, 0.0, 5.0] for _ in secondary_types]
    r.secondary_masses = [0.0 for _ in secondary_types]
    r.secondary_helicities = [0.0 for _ in secondary_types]
    r.secondary_times = [time for _ in secondary_types]
    return r


def _make_cascade(dc):
    """root (NuMu+EPlus -> NuMu) then child (NuMu+proton -> EMinus).

    The child's incoming primary is the root's outgoing secondary (shared id), so
    the child interaction has a genuine, non-root target -- exercising the
    status-22 encoding.
    """
    PT = dc.ParticleType
    root = _record(dc, PT.NuMu, dc.ParticleID(9, 9), PT.EPlus, 0.000511,
                   [PT.NuMu], [dc.ParticleID(1, 1)], [1.0, 2.0, 3.0], 5.0)
    child = _record(dc, PT.NuMu, dc.ParticleID(1, 1), PT.HNucleus, 0.938,
                    [PT.EMinus], [dc.ParticleID(1, 2)], [4.0, 5.0, 6.0], 6.0)
    tree = dc.InteractionTree()
    d0 = tree.add_entry(root, None)
    tree.add_entry(child, d0)
    tree.header.event_number = 7
    tree.header.weights = [3.0]
    return tree


def _write_or_skip(dc, path):
    from siren import hepmc3 as sio
    tree = _make_cascade(dc)
    try:
        sio.SaveInteractionTreesAsHepMC3([tree], path)
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise


def test_nuhepmc_registry_keys_in_ascii(dc, tmp_path):
    out = str(tmp_path / "cascade.hepmc3")
    _write_or_skip(dc, out)
    with open(out) as f:
        text = f.read()
    for key in (
        "NuHepMC.VertexStatusIDs",
        "NuHepMC.VertexStatusInfo[22].Name",
        "NuHepMC.ParticleStatusIDs",
        "NuHepMC.ParticleStatusInfo[22].Name",
        "NuHepMC.AdditionalParticleNumbers",
        "NuHepMC.AdditionalParticleNumber[5914].Name",   # writer stub
        "NuHepMC.AdditionalParticleInfo[5914].Name",     # reader stub (dual-emit)
        "NuHepMC.AdditionalParticleInfo[5914].Description",
        "NuHepMC.Conventions",
        "lab_pos",
    ):
        assert key in text, "missing NuHepMC key: " + key


def test_external_reader_pyhepmc(dc, tmp_path):
    pyhepmc = pytest.importorskip("pyhepmc")
    out = str(tmp_path / "cascade.hepmc3")
    _write_or_skip(dc, out)

    with pyhepmc.open(out) as f:
        events = [e for e in f]
    assert len(events) == 1
    evt = events[0]

    # Run-level NuHepMC registries are visible to the external reader.
    ri_attrs = evt.run_info.attributes
    vids = str(ri_attrs["NuHepMC.VertexStatusIDs"]).split()
    assert vids == ["1", "21", "22"]
    psids = str(ri_attrs["NuHepMC.ParticleStatusIDs"]).split()
    assert "22" in psids and "4" in psids and "20" in psids
    apn = str(ri_attrs["NuHepMC.AdditionalParticleNumbers"]).split()
    assert "5914" in apn                                  # the HNL N4 is declared
    assert "NuHepMC.AdditionalParticleInfo[5914].Name" in list(ri_attrs)

    # E.R.7: exactly one beam (status 4) and one target (status 20) per event;
    # the deeper-interaction target is marked 22.
    statuses = [p.status for p in evt.particles]
    assert statuses.count(4) == 1
    assert statuses.count(20) == 1
    assert statuses.count(22) >= 1

    # E.R.5: per-event lab position (three spatial slots; the lossless time lives
    # in the Minkowski vertex ct slot, not a truncation-prone lab_pos seconds slot).
    lab = str(evt.attributes["lab_pos"]).split()
    assert len(lab) == 3
    # vertex {1,2,3} internal meters -> {100,200,300} CM
    assert [round(float(v)) for v in lab[:3]] == [100, 200, 300]


def test_nonroot_target_roundtrips_through_siren(dc, tmp_path):
    from siren import hepmc3 as sio
    out = str(tmp_path / "cascade.hepmc3")
    _write_or_skip(dc, out)

    loaded = sio.LoadInteractionTreesFromHepMC3(out)
    assert len(loaded) == 1
    lt = loaded[0]
    assert len(lt.tree) == 2
    assert lt.tree[0].is_root()
    assert lt.tree[0].record.signature.target_type == dc.ParticleType.EPlus
    assert not lt.tree[1].is_root()
    assert lt.tree[1].parent_index == 0
    assert lt.tree[1].record.signature.target_type == dc.ParticleType.HNucleus
