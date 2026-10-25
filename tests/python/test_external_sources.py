"""Complete-file exposure, stable source IDs, raw records and weighted replay."""
import math
import pickle
import shutil

import numpy as np
import pytest
import siren
from siren import dk2nu
from siren.errors import ConfigurationError


def write_beam(path, *, times=True, pots=(3.0, 7.0), bad_ancestor=False):
    uproot = pytest.importorskip("uproot")
    ak = pytest.importorskip("awkward")
    columns = {name: np.zeros(3) for name in dk2nu._BRANCH_SETS[1].values()}
    columns.update({"decay.ptype": np.array([211, 321, 211]),
                    "decay.ntype": np.full(3, 14), "decay.ndecay": np.full(3, 13),
                    "decay.ppenergy": np.full(3, 99.0),
                    "decay.pdpz": np.array([1.0, 2.0, 3.0]),
                    "decay.nimpwt": np.array([0.2, 0.3, 0.5])})
    with uproot.recreate(path) as f:
        # TTrees exercise the standard dk2nu branch names and jagged ancestors.
        if times:
            columns["ancestor.pdg"] = ak.Array([[211, 14], [] if bad_ancestor else [321, 14], [211, 14]])
            columns["ancestor.startt"] = ak.Array([[0.0, 12.0], [] if bad_ancestor else [0.0, 13.0], [0.0, 14.0]])
        f.mktree("dk2nuTree", columns)
        f.mktree("dkmetaTree", {"dkmeta/pots": np.asarray(pots)})
    return path


def test_full_exposure_selection_and_order_independent_source_ids(tmp_path):
    a = write_beam(tmp_path / "a.root")
    b = write_beam(tmp_path / "b.root", pots=(20.0,))
    x = dk2nu.read_dk2nu([a, b], parent_pdg=211)
    y = dk2nu.read_dk2nu([b, a], parent_pdg=211)
    assert x["pot"] == y["pot"] == 30
    assert x["source_entry"].tolist() == [0, 2, 0, 2]
    for k in ("source_file_id", "source_entry", "t0", "E", "nimpwt"):
        assert np.array_equal(x[k], np.concatenate([y[k][2:], y[k][:2]]))
    assert all(len(s["sha256"]) == 64 and s["entries"] == 3 for s in x["source_files"])
    assert np.all(x["source_sampling_probability"] == 1)
    table = dk2nu.dk2nu_to_primary_distribution(x, None, frame="detector")
    assert table.normalization == pytest.approx(1.4 / 30)
    r = siren.dataclasses.PrimaryDistributionRecord(siren.particles.PiPlus)
    table.Sample(siren.utilities.SIREN_random(1), None, None, r)
    event = siren.dataclasses.InteractionRecord(); r.finalize(event)
    assert event.interaction_parameters["source_file_id"] in x["source_file_id"]
    assert event.interaction_parameters["source_entry"] in (0, 2)


def test_mixed_file_times_do_not_depend_on_order(tmp_path):
    a = write_beam(tmp_path / "timed.root")
    b = write_beam(tmp_path / "untimed.root", times=False)
    for files in ([a, b], [b, a]):
        data = dk2nu.read_dk2nu(files)
        assert data["time_available"].sum() == 3
        assert np.isnan(data["t0"]).sum() == 3
        with pytest.raises(ConfigurationError, match="t0"):
            dk2nu.dk2nu_to_primary_distribution(data, None, frame="detector")
        timeless = dk2nu.read_dk2nu(files, read_time=False)
        assert "t0" not in timeless
        assert len(timeless["ptype"]) == 6


def test_reader_refuses_slicing_duplicates_and_damaged_ancestors(tmp_path):
    a = write_beam(tmp_path / "a.root")
    b = tmp_path / "copy.root"; shutil.copyfile(a, b)
    with pytest.raises(ConfigurationError, match="Duplicate"):
        dk2nu.read_dk2nu([a, b])
    with pytest.raises(ConfigurationError, match="slicing"):
        dk2nu.read_dk2nu(a, entry_stop=1)
    bad = write_beam(tmp_path / "bad.root", bad_ancestor=True)
    with pytest.raises(ConfigurationError, match="ancestor"):
        dk2nu.read_dk2nu(bad)


def test_spectrum_uses_decay_energy_and_all_metadata_pot(tmp_path, capsys):
    data = dk2nu.read_dk2nu(write_beam(tmp_path / "beam.root"))
    centers, density = dk2nu.meson_energy_spectrum(data, 211, n_bins=5, E_range=(0, 5))
    assert centers[np.nonzero(density)].tolist() == [1.5, 3.5]
    assert density.sum() == pytest.approx(0.07)
    dk2nu.print_summary(data)
    assert "99.000" not in capsys.readouterr().out
    with pytest.raises(ConfigurationError, match="POT"):
        dk2nu.meson_energy_spectrum(dict(data, pot=0), 211)


@pytest.mark.parametrize("keys, rows", [(["E", "E"], [[1., 2.]]),
    (["E", "px"], [[1., 2.]]), (["E"], [[]]), (["E"], [[float("nan")]]),
    (["E"], []), (["PrimaryExternalDistribution_row"], [[0.]])])
def test_native_table_rejects_malformed_arrays(keys, rows):
    with pytest.raises(RuntimeError):
        siren.distributions.PrimaryExternalDistribution(keys, rows)


def table(weights=(0.2, 0.8), bias=None, reverse=False):
    keys = ["E", "px", "py", "pz", "m", "x", "y", "z", "weight", "source_entry"]
    rows = [[math.hypot(p, 1.), 0., 0., p, 1., 0., 0., float(i), weights[i], i]
            for i, p in enumerate((1., 3.))]
    if reverse: rows.reverse()
    args = [keys, rows] + ([] if bias is None else [list(bias)])
    return siren.distributions.PrimaryExternalDistribution(*args)


def setup(dist, physical=None, n=256, seed=717):
    from test_native_phase_space import fixed_decay_setup
    model = siren.interactions.HNLDipoleDecay(1., 1e-6, siren.interactions.HNLDipoleDecay.ChiralNature.Dirac)
    gen, unused = fixed_decay_setup(model, siren.particles.N4, 1., n)
    process = gen.GetPrimaryProcess()
    process.distributions = [dist, siren.distributions.PrimaryNeutrinoHelicityDistribution()]
    detector = siren.detector.DetectorModel()
    generator = siren.injection._Injector(n, detector, process, siren.utilities.SIREN_random(seed))
    phys = siren.injection.PhysicalProcess()
    phys.primary_type = siren.particles.N4
    phys.interactions = process.interactions
    phys.weighting_mode = siren.Fixed()
    phys.distributions = [physical or dist, process.distributions[1]]
    weighter = siren.injection._Weighter([generator], detector, phys)
    return generator, weighter, phys, detector


@pytest.mark.parametrize("bias", [None, (1., 4.)])
def test_absolute_row_weight_is_physical_over_proposal_and_survives_archive(bias, tmp_path):
    dist = table(bias=bias)
    gen, weighter, phys, detector = setup(dist)
    events = [gen.GenerateEvent() for _ in range(256)]
    gen.SaveInjector(str(tmp_path / "source"))
    restored = siren.injection._Injector(256, str(tmp_path / "source"), siren.utilities.SIREN_random(1))
    replay = siren.injection._Weighter([restored], detector, phys)
    weights = []
    for event in events:
        i = int(event.tree[0].record.interaction_parameters["source_entry"])
        q = 0.5 if bias is None else bias[i] / sum(bias)
        expected = (0.2, 0.8)[i] / q / 256
        assert weighter.EventWeight(event) == pytest.approx(expected, rel=2e-12)
        assert replay.EventWeight(event) == pytest.approx(expected, rel=2e-12)
        weights.append(expected)
    assert sum(weights) == pytest.approx(1, abs=0.1)


def test_reordered_physical_table_rejected():
    with pytest.raises(ConfigurationError, match="row layout"):
        setup(table(), table(reverse=True))


def test_raw_four_components_and_mass_are_not_projected():
    dist = siren.distributions.PrimaryExternalDistribution(
        ["E", "px", "py", "pz", "m", "x", "y", "z"], [[2., 1., 2., 3., 0.7, 0., 0., 0.]])
    r = siren.dataclasses.PrimaryDistributionRecord(siren.particles.N4)
    dist.Sample(siren.utilities.SIREN_random(7), None, None, r)
    event = siren.dataclasses.InteractionRecord(); r.finalize(event)
    assert list(event.primary_momentum) == [2., 1., 2., 3.]
    assert event.primary_mass == 0.7


@pytest.mark.parametrize("bias", [[1., float("nan")], [1., float("inf")], [1., 0.]])
def test_invalid_or_incomplete_native_proposal_support_rejected(bias):
    with pytest.raises(RuntimeError): table(bias=bias)


def test_dif_dar_split_uses_decay_momentum_and_preserves_columns(tmp_path):
    dif = dk2nu.read_dk2nu(write_beam(tmp_path / "dif.root"))
    dar = dk2nu.read_dk2nu(write_beam(tmp_path / "dar.root", pots=(20.,)))
    dar["pz"][:] = 0.
    merged = dk2nu.combine_dif_dar(dif, dar)
    assert merged["dar"].tolist() == [False] * 3 + [True] * 3
    assert merged["pot"] == 10
    assert merged["nimpwt"] / merged["pot"] == pytest.approx([.02, .03, .05, .01, .015, .025])
    assert len(merged["source_files"]) == 2
    assert len(merged["source_entry"]) == 6
    with pytest.raises(ConfigurationError, match="columns differ"):
        dk2nu.combine_dif_dar(dif, {k: v for k, v in dar.items() if k != "t0"})
