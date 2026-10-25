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


def boundary_data():
    return {'ptype':np.array([211,211]),'E':np.array([99.,99.]),'px':np.zeros(2),'py':np.zeros(2),'pz':np.array([0.,1.]),'vx':np.zeros(2),'vy':np.zeros(2),'vz':np.zeros(2),'nimpwt':np.array([.2,.8]),'pot':10.}


def test_masked_species_must_not_be_silently_unmasked():
    d=boundary_data();d['ptype']=np.ma.array([211,211],mask=[True,False])
    with pytest.raises(ConfigurationError):
        dk2nu.dk2nu_to_primary_distribution(d,None,frame='detector')


def test_dif_dar_invalid_momentum_must_not_disappear():
    d=boundary_data();d['pz'][0]=np.nan
    with pytest.raises(ConfigurationError):
        dk2nu.combine_dif_dar(d,boundary_data())


def test_spectrum_invalid_physical_weights_must_not_be_negative():
    d=boundary_data();d['nimpwt'][:]=-1
    with pytest.raises(ConfigurationError):
        dk2nu.meson_energy_spectrum(d,211)


def test_unrepresentable_proposal_support_must_not_be_claimed():
    with pytest.raises(RuntimeError):
        siren.distributions.PrimaryExternalDistribution(['E','weight'],[[1.,1.],[2.,1.]],[1.,1e-300])


def test_unweighted_table_removes_proposal_bias_per_event():
    dist = siren.distributions.PrimaryExternalDistribution(
        ["E", "px", "py", "pz", "m", "x", "y", "z", "source_entry"],
        [[math.sqrt(2), 0, 0, 1, 1, 0, 0, 0, 0],
         [math.sqrt(10), 0, 0, 3, 1, 0, 0, 0, 1]], [1., 3.])
    gen, weighter, _, _ = setup(dist, n=128)
    assert dist.PhysicalDensityDiffers()
    seen = set()
    events = [gen.GenerateEvent() for _ in range(128)]
    for event in events:
        row = int(event.tree[0].record.interaction_parameters["source_entry"])
        seen.add(row)
        assert weighter.EventWeight(event) == pytest.approx(.5 / (.25, .75)[row] / 128)
    assert seen == {0, 1}


def test_filter_rebuilds_physical_support_before_proposal():
    dist = siren.distributions.PrimaryExternalDistribution(
        ["E", "weight"], [[.5, 1.], [2., 0.], [3., 1.]], [1., 0., 1.], 1.)
    assert dist.GetPhysicalNumEvents() == 2
    assert dist.normalization == 1.


def test_pot_sum_is_stable_for_unequal_file_exposures(tmp_path):
    paths = [write_beam(tmp_path / f"beam-{i}.root", pots=(pot,))
             for i, pot in enumerate([1e16, 1., 1.])]
    assert dk2nu.read_dk2nu(paths)["pot"] == dk2nu.read_dk2nu(paths[::-1])["pot"] == 1e16 + 2


def test_nonzero_threshold_requires_energy_column():
    with pytest.raises(RuntimeError, match="E column"):
        siren.distributions.PrimaryExternalDistribution(["px", "py", "pz"], [[0., 0., .1]], .2)


@pytest.mark.parametrize("density", [-1., float("nan"), float("inf")])
def test_legacy_cached_proposal_density_must_be_valid(density):
    dist = siren.distributions.PrimaryExternalDistribution(["E"], [[1.], [2.]], [1., 3.])
    record = siren.dataclasses.InteractionRecord()
    record.primary_momentum = [1., 0., 0., 1.]
    record.interaction_parameters = {"PrimaryExternalDistribution_gen_prob": density}
    with pytest.raises(siren.errors.WeightCalculationError, match="legacy"):
        dist.GenerationProbability(None, None, record)


@pytest.mark.parametrize("label", ["positive", "negative"])
def test_version_zero_weight_column_stays_metadata_on_resave(label, tmp_path):
    import sys
    from pathlib import Path
    if sys.byteorder != "little":
        pytest.skip("Native archive fixture is little-endian")
    fixture = Path(__file__).parents[1] / "fixtures" / "external_v0" / f"{label}.injector"
    random = siren.utilities.SIREN_random(17)
    gen = siren.injection._Injector(1, str(fixture), random)
    assert gen.GetRandom() is random
    for _ in range(2):
        dist = gen.GetPrimaryProcess().distributions[0]
        assert not dist.IsNormalizationSet()
        assert dist.normalization == 1
        assert not dist.PhysicalDensityDiffers()
        gen.SaveInjector(str(tmp_path / "roundtrip"))
        gen = siren.injection._Injector(1, str(tmp_path / "roundtrip"), siren.utilities.SIREN_random(1))


def test_explicit_normalization_survives_native_archive(tmp_path):
    dist = table()
    dist.normalization = 9.
    gen, _, _, _ = setup(dist)
    gen.SaveInjector(str(tmp_path / "normalized"))
    restored = siren.injection._Injector(256, str(tmp_path / "normalized"), siren.utilities.SIREN_random(1))
    assert restored.GetPrimaryProcess().distributions[0].normalization == 9.


def test_initial_only_table_preserves_separately_sampled_vertex():
    dist = siren.distributions.PrimaryExternalDistribution(["x0", "y0", "z0"], [[1., 2., 3.]])
    record = siren.dataclasses.PrimaryDistributionRecord(siren.particles.N4)
    record.interaction_vertex = [4., 5., 6.]
    dist.Sample(siren.utilities.SIREN_random(1), None, None, record)
    assert list(record.interaction_vertex) == [4., 5., 6.]
    assert list(record.initial_position) == [1., 2., 3.]


def test_point_vertex_rejects_propagation_weighting_at_configuration():
    gen, _, phys, detector = setup(table())
    phys.weighting_mode = siren.Propagated()
    with pytest.raises(ConfigurationError, match="Fixed weighting"):
        siren.injection._Weighter([gen], detector, phys)


def test_reader_can_explicitly_disable_missing_polarization_fields(tmp_path):
    uproot = pytest.importorskip("uproot")
    path = write_beam(tmp_path / "full.root", times=False)
    with uproot.open(path) as f:
        columns = f["dk2nuTree"].arrays(library="np")
    for name in ("ppdxdz", "ppdydz", "pppz", "mupare", "muparpx", "muparpy", "muparpz"):
        del columns[f"decay.{name}"]
    slim = tmp_path / "slim.root"
    with uproot.recreate(slim) as f:
        f.mktree("dk2nuTree", columns)
        f.mktree("dkmetaTree", {"dkmeta/pots": np.array([10.])})
    data = dk2nu.read_dk2nu(slim, read_time=False, read_polarization=False)
    assert len(data["ptype"]) == 3
    with pytest.raises(ValueError, match="branches"):
        dk2nu.read_dk2nu(slim, read_time=False)


def test_dif_dar_at_rest_parents_run_through_public_simulation():
    dif = boundary_data(); dif["ptype"][:] = 111; dif["pz"][:] = [1., 2.]
    dar = boundary_data(); dar["ptype"][:] = 111; dar["pz"][:] = 0.; dar["pot"] = 20.
    merged = dk2nu.combine_dif_dar(dif, dar)
    dist = dk2nu.dk2nu_to_primary_distribution(merged, None, frame="detector", sampling_bias=merged["nimpwt"])
    signature = siren.dataclasses.InteractionSignature()
    signature.primary_type = siren.particles.Pi0
    signature.target_type = siren.particles.Decay
    signature.secondary_types = [siren.particles.Gamma] * 2
    decay = siren.injection.PhaseSpaceDecay(signature, [0., 0.], 1e-16, 1e-16, siren.injection.Isotropic2BodyChannel())
    vertex = siren.Vertex("Pi0", decay, distributions=[dist], physical=[dist],
                          kinematics=siren.channels.physical(), weighting=siren.Fixed())
    result = siren.Simulation(detector=siren.detector.DetectorModel(), primary=vertex, events=128, seed=37).run(
        on_failure="raise", on_shortfall="raise")
    assert float(result.weights.sum()) == pytest.approx(.15, rel=2e-12)


@pytest.mark.parametrize("weighted", [False, True])
def test_proposal_support_only_required_above_threshold(weighted):
    keys = ["E", "weight"] if weighted else ["E"]
    rows = [[.5, 1.], [2., 1.]] if weighted else [[.5], [2.]]
    dist = siren.distributions.PrimaryExternalDistribution(keys, rows, [0., 1.], 1.)
    assert dist.GetPhysicalNumEvents() == 1


def test_initial_only_source_requires_a_vertex_sampler_at_configuration():
    dist = siren.distributions.PrimaryExternalDistribution(
        ["E", "px", "py", "pz", "m", "x0", "y0", "z0"],
        [[2., 0., 0., math.sqrt(3), 1., 0., 0., 0.]])
    with pytest.raises(ConfigurationError, match="requires a separate vertex sampler"):
        setup(dist)


@pytest.mark.parametrize("bias", [None, []])
def test_physical_weight_validation_uses_threshold_survivors(bias, tmp_path):
    path = tmp_path / "threshold.csv"
    path.write_text("E,weight\n0.5,-1\n2,1\n")
    csv = siren.distributions.PrimaryExternalDistribution(str(path), 1.)
    args = [["E", "weight"], [[.5, -1.], [2., 1.]]]
    if bias is not None:
        args.append(bias)
    array = siren.distributions.PrimaryExternalDistribution(*args, 1.)
    assert array.GetPhysicalNumEvents() == csv.GetPhysicalNumEvents() == 1
    assert array.normalization == csv.normalization == 1.
