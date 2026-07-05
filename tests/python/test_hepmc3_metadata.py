"""Conformance checks for the NuHepMC run-level metadata SIREN now emits:
process-ID registry (G.R.8) + per-event signal_process_id (E.R.3), generation
counts, cross-section units (G.R.6), and the flux-averaged cross-section
diagnostics/estimator (E.C.4). Also a small write-throughput benchmark.

Skips cleanly when SIREN was built without HepMC3 support.
"""
import os
import tempfile
import time

import pytest

siren = pytest.importorskip("siren")


@pytest.fixture(scope="module")
def dc():
    from siren import dataclasses
    return dataclasses


def _tree(dc, primary_type, weight):
    r = dc.InteractionRecord()
    sig = r.signature
    sig.primary_type = primary_type
    sig.target_type = dc.ParticleType.EPlus
    sig.secondary_types = [dc.ParticleType.NuMu]
    r.signature = sig
    r.primary_momentum = [10.0, 0.0, 0.0, 10.0]
    r.target_mass = 0.000511
    r.interaction_vertex = [1.0, 2.0, 3.0]
    r.interaction_time = 5.0
    r.secondary_ids = [dc.ParticleID(1, 1)]
    r.secondary_momenta = [[10.0, 0.0, 0.0, 10.0]]
    r.secondary_masses = [0.0]
    r.secondary_helicities = [-1.0]
    r.secondary_times = [5.0]
    r.interaction_parameters = {"Q2": 0.7}
    t = dc.InteractionTree()
    t.add_entry(r, None)
    t.header.event_number = 1
    t.header.weights = [weight]
    return t


def _write_or_skip(dc, trees, path, options=None):
    from siren import hepmc3 as sio
    try:
        if options is not None:
            sio.SaveInteractionTreesAsHepMC3(trees, path, options)
        else:
            sio.SaveInteractionTreesAsHepMC3(trees, path)
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise


def test_process_registry_and_signal_process_id(dc, tmp_path):
    from siren import hepmc3
    out = str(tmp_path / "reg.hepmc3")
    # Two distinct root signatures (different primaries) -> two process ids.
    trees = [_tree(dc, dc.ParticleType.NuMu, 2.0), _tree(dc, dc.ParticleType.NuEBar, 4.0)]
    _write_or_skip(dc, trees, out, hepmc3.HepMC3WriterOptions())
    text = open(out).read()
    # G.R.8 registry declared, with two entries starting at the "Other" band (700).
    assert "NuHepMC.ProcessIDs" in text
    assert "NuHepMC.ProcessInfo[700].Name" in text
    assert "NuHepMC.ProcessInfo[701].Name" in text
    # E.R.3 per-event signal_process_id present on events.
    assert "signal_process_id" in text


def test_counts_and_fatx_metadata(dc, tmp_path):
    pyhepmc = pytest.importorskip("pyhepmc")
    from siren import hepmc3
    out = str(tmp_path / "fatx.hepmc3")
    opts = hepmc3.HepMC3WriterOptions()
    opts.attempted_events = 1000
    weights = (2.0, 4.0, 6.0)
    trees = [_tree(dc, dc.ParticleType.NuMu, w) for w in weights]
    # Default weights_state ("header") is a NuHepMC mode: the trees carry header CVs.
    _write_or_skip(dc, trees, out, opts)
    with pyhepmc.open(out) as f:
        attrs = f.read().run_info.attributes
    # Counts stored as metadata; accepted auto-fills to the number of trees.
    assert float(str(attrs["siren.attempted_events"])) == 1000
    assert float(str(attrs["siren.accepted_events"])) == len(trees)
    # FATX raw ingredient (siren.fatx.weight_sum) == sum of the CV weights.
    assert abs(float(str(attrs["siren.fatx.weight_sum"])) - sum(weights)) < 1e-9
    assert "siren.fatx.value" in attrs
    # Invariant: a declared NuHepMC.Version implies the mandatory FATX key + G.R.6 units.
    assert "NuHepMC.Version.Major" in attrs
    assert "NuHepMC.FluxAveragedTotalCrossSection" in attrs
    assert "NuHepMC.Units.CrossSection.Unit" in attrs
    assert "NuHepMC.Units.CrossSection.TargetScale" in attrs


def test_fatx_target_scale_labels_per_atom_optin(dc, tmp_path):
    """In a NuHepMC mode the reserved FATX key + units are always emitted; the
    fatx_per_atom opt-in only selects the TargetScale label (the claim that the
    value is a genuine per-atom sigma vs an unnormalized rate)."""
    from siren import hepmc3
    # With the opt-in the caller-provided TargetScale is used.
    out = str(tmp_path / "peratom.hepmc3")
    opts = hepmc3.HepMC3WriterOptions()
    opts.attempted_events = 1000
    opts.fatx_per_atom = True
    opts.target_scale = "PerAtom"
    _write_or_skip(dc, [_tree(dc, dc.ParticleType.NuMu, 2.0)], out, opts)
    text = open(out).read()
    assert "NuHepMC.FluxAveragedTotalCrossSection" in text
    assert "NuHepMC.Units.CrossSection.Unit" in text
    assert "PerAtom" in text
    # Without the opt-in the value is labelled Unnormalized (rate-weight caveat).
    out2 = str(tmp_path / "unnorm.hepmc3")
    opts2 = hepmc3.HepMC3WriterOptions()
    opts2.attempted_events = 1000
    opts2.fatx_per_atom = False
    _write_or_skip(dc, [_tree(dc, dc.ParticleType.NuMu, 2.0)], out2, opts2)
    text2 = open(out2).read()
    assert "NuHepMC.FluxAveragedTotalCrossSection" in text2
    assert "Unnormalized" in text2


def test_unweighted_mode_omits_nuhepmc_and_fatx(dc, tmp_path):
    """weights_state == "unweighted" produces a plain HepMC3 file: no NuHepMC.*
    keys and no siren.fatx.* keys, but the generation counts survive for pooling."""
    from siren import hepmc3
    out = str(tmp_path / "unweighted.hepmc3")
    opts = hepmc3.HepMC3WriterOptions()
    opts.weights_state = "unweighted"
    opts.attempted_events = 1000
    _write_or_skip(dc, [_tree(dc, dc.ParticleType.NuMu, 2.0)], out, opts)
    text = open(out).read()
    assert "siren.weights_state" in text and "unweighted" in text
    assert "NuHepMC." not in text
    assert "siren.fatx.value" not in text
    assert "siren.fatx.weight_sum" not in text
    # Counts still present (later pooling needs them).
    assert "siren.attempted_events" in text
    assert "siren.accepted_events" in text


def test_fatx_value_matches_estimator(dc, tmp_path):
    """FATX estimate == 3.894e8 pb/GeV^-2 * sum(weights)/attempted (via pyhepmc)."""
    pyhepmc = pytest.importorskip("pyhepmc")
    out = str(tmp_path / "est.hepmc3")
    from siren import hepmc3
    opts = hepmc3.HepMC3WriterOptions()
    opts.attempted_events = 10
    opts.fatx_per_atom = True
    weights = (2.0, 4.0, 6.0)
    _write_or_skip(dc, [_tree(dc, dc.ParticleType.NuMu, w) for w in weights], out, opts)
    with pyhepmc.open(out) as f:
        evt = f.read()
        ri = evt.run_info
    val = float(str(ri.attributes["NuHepMC.FluxAveragedTotalCrossSection"]))
    expected = 3.894e8 * sum(weights) / 10.0
    assert abs(val - expected) / expected < 1e-6


def test_norm_basis_records_denominator(dc, tmp_path):
    """siren.fatx.norm_basis records which count normalized the FATX: "attempted"
    when attempted_events was provided, "accepted" when the code fell back to the
    accepted count. A pooler must not mix the two bases."""
    pyhepmc = pytest.importorskip("pyhepmc")
    from siren import hepmc3
    # attempted provided -> attempted basis.
    out = str(tmp_path / "basis_att.hepmc3")
    opts = hepmc3.HepMC3WriterOptions()
    opts.attempted_events = 1000
    _write_or_skip(dc, [_tree(dc, dc.ParticleType.NuMu, 2.0)], out, opts)
    with pyhepmc.open(out) as f:
        attrs = f.read().run_info.attributes
    assert str(attrs["siren.fatx.norm_basis"]) == "attempted"
    # attempted absent -> silent fallback to the accepted (auto-filled) count.
    out2 = str(tmp_path / "basis_acc.hepmc3")
    opts2 = hepmc3.HepMC3WriterOptions()  # no attempted_events
    _write_or_skip(dc, [_tree(dc, dc.ParticleType.NuMu, 2.0)], out2, opts2)
    with pyhepmc.open(out2) as f:
        attrs2 = f.read().run_info.attributes
    assert str(attrs2["siren.fatx.norm_basis"]) == "accepted"


def test_events_to_inject_persisted(dc, tmp_path):
    """The pooled-weighting seed N_i (Injector EventsToInject) is persisted as
    siren.events_to_inject when set, and absent when left unset (< 0)."""
    from siren import hepmc3
    out = str(tmp_path / "n_i.hepmc3")
    opts = hepmc3.HepMC3WriterOptions()
    opts.attempted_events = 1000
    opts.events_to_inject = 500
    _write_or_skip(dc, [_tree(dc, dc.ParticleType.NuMu, 2.0)], out, opts)
    text = open(out).read()
    assert "siren.events_to_inject" in text
    # Left unset -> not emitted.
    out2 = str(tmp_path / "no_n_i.hepmc3")
    opts2 = hepmc3.HepMC3WriterOptions()
    opts2.attempted_events = 1000
    _write_or_skip(dc, [_tree(dc, dc.ParticleType.NuMu, 2.0)], out2, opts2)
    assert "siren.events_to_inject" not in open(out2).read()


def test_partitioned_fatx_carries_raw_sums(dc, tmp_path):
    """A partitioned (mixed-primary) file must be losslessly poolable: each primary
    carries the normalized siren.fatx.<pdg> value plus the raw ingredients
    (.weight_sum, .accepted) that formed it, so summing raw across files and
    renormalizing reproduces the combined estimate."""
    pyhepmc = pytest.importorskip("pyhepmc")
    from siren import hepmc3
    out = str(tmp_path / "part.hepmc3")
    opts = hepmc3.HepMC3WriterOptions()
    opts.attempted_events = 1000
    opts.fatx_partition_by_primary = True
    # Two primaries: NuMu with weights 2 + 6, NuEBar with weight 4.
    trees = [
        _tree(dc, dc.ParticleType.NuMu, 2.0),
        _tree(dc, dc.ParticleType.NuEBar, 4.0),
        _tree(dc, dc.ParticleType.NuMu, 6.0),
    ]
    _write_or_skip(dc, trees, out, opts)
    with pyhepmc.open(out) as f:
        attrs = f.read().run_info.attributes
    numu = int(dc.ParticleType.NuMu)
    nuebar = int(dc.ParticleType.NuEBar)
    # The partition marker flags that per-primary keys are present (a consumer must
    # read the per-pdg breakdown, not treat siren.fatx.value as the whole story).
    assert int(str(attrs["siren.fatx_partitioned"])) == 1
    # Raw per-primary CV weight sums.
    assert abs(float(str(attrs["siren.fatx.%d.weight_sum" % numu])) - 8.0) < 1e-9
    assert abs(float(str(attrs["siren.fatx.%d.weight_sum" % nuebar])) - 4.0) < 1e-9
    # Raw per-primary accepted counts.
    assert float(str(attrs["siren.fatx.%d.accepted" % numu])) == 2
    assert float(str(attrs["siren.fatx.%d.accepted" % nuebar])) == 1
    # Normalized values are still present, and equal 3.894e8 * weight_sum / norm.
    norm = 1000.0
    exp_numu = 3.894e8 * 8.0 / norm
    assert abs(float(str(attrs["siren.fatx.%d" % numu])) - exp_numu) / exp_numu < 1e-6
    # Poolability: the normalized value can be reconstructed from the raw ingredients.
    ws = float(str(attrs["siren.fatx.%d.weight_sum" % numu]))
    assert abs(3.894e8 * ws / norm - float(str(attrs["siren.fatx.%d" % numu]))) / exp_numu < 1e-6


def test_partitioned_fatx_suppressed_when_unweighted(dc, tmp_path):
    """fatx_partition_by_primary must be gated by the same NuHepMC-mode guard as
    the rest of the FATX block: with weights_state == "unweighted" a mixed-primary
    file emits NO per-primary siren.fatx.<pdg> keys and NO siren.fatx_partitioned
    marker, and it must not fall back to presenting a single conflated global FATX
    value. (Under "header" the same inputs DO emit the partition; this is purely
    the suppression axis.)"""
    from siren import hepmc3
    numu = int(dc.ParticleType.NuMu)
    nuebar = int(dc.ParticleType.NuEBar)
    trees = [
        _tree(dc, dc.ParticleType.NuMu, 2.0),
        _tree(dc, dc.ParticleType.NuEBar, 4.0),
        _tree(dc, dc.ParticleType.NuMu, 6.0),
    ]

    # Control: default (header) mode with the same inputs DOES partition.
    ctrl = str(tmp_path / "part_header.hepmc3")
    ctrl_opts = hepmc3.HepMC3WriterOptions()
    ctrl_opts.attempted_events = 1000
    ctrl_opts.fatx_partition_by_primary = True
    _write_or_skip(dc, trees, ctrl, ctrl_opts)
    ctrl_text = open(ctrl).read()
    assert "siren.fatx_partitioned" in ctrl_text
    assert ("siren.fatx.%d" % numu) in ctrl_text

    # Unweighted mode with partitioning requested: everything FATX is suppressed.
    out = str(tmp_path / "part_unweighted.hepmc3")
    opts = hepmc3.HepMC3WriterOptions()
    opts.weights_state = "unweighted"
    opts.attempted_events = 1000
    opts.fatx_partition_by_primary = True
    _write_or_skip(dc, trees, out, opts)
    text = open(out).read()
    assert "siren.weights_state" in text and "unweighted" in text
    # No partition marker and no per-primary keys.
    assert "siren.fatx_partitioned" not in text
    assert ("siren.fatx.%d" % numu) not in text
    assert ("siren.fatx.%d" % nuebar) not in text
    # And no conflated single global FATX value presented as the answer.
    assert "siren.fatx.value" not in text
    assert "NuHepMC.FluxAveragedTotalCrossSection" not in text
    # Counts still survive for later pooling.
    assert "siren.accepted_events" in text


def test_fatx_precision_sweep(dc, tmp_path):
    """The ascii writer must preserve the flux-averaged cross-section value across
    a wide magnitude range. Driving synthetic header CV weights so the writer's
    normal path emits FATX values from 1e-2 pb down to 1e-30 pb, the read-back
    value must match the writer's own estimate to a relative error < 1e-12 at
    every magnitude. A larger error would mean the ascii float formatting (default
    16 significant digits) is truncating precision on small cross sections."""
    pyhepmc = pytest.importorskip("pyhepmc")
    from siren import hepmc3
    # Same constant the writer uses: 1 GeV^-2 = 3.894e8 pb (PDG).
    kGeVm2_to_pb = 3.894e8
    norm = 1000  # attempted_events; exact integer denominator
    # Targets span 28 decades. For each we pick the CV weight that makes the
    # writer's estimate 3.894e8 * weight_sum / norm land on the target magnitude.
    for exponent in range(-2, -31, -1):
        target = 10.0 ** exponent
        weight = target * norm / kGeVm2_to_pb
        # The value the writer will actually compute and emit for this weight.
        expected = kGeVm2_to_pb * weight / norm

        out = str(tmp_path / ("sweep_%d.hepmc3" % exponent))
        opts = hepmc3.HepMC3WriterOptions()
        opts.attempted_events = norm
        opts.fatx_per_atom = True
        _write_or_skip(dc, [_tree(dc, dc.ParticleType.NuMu, weight)], out, opts)

        with pyhepmc.open(out) as f:
            attrs = f.read().run_info.attributes
        # Both the raw ingredient and the normalized estimate must survive.
        got_weight_sum = float(str(attrs["siren.fatx.weight_sum"]))
        got_fatx = float(str(attrs["NuHepMC.FluxAveragedTotalCrossSection"]))
        assert abs(got_weight_sum - weight) <= abs(weight) * 1e-12, \
            (exponent, got_weight_sum, weight)
        assert abs(got_fatx - expected) <= abs(expected) * 1e-12, \
            (exponent, got_fatx, expected)


def test_write_throughput_benchmark(dc, tmp_path, capsys):
    """Not a hard gate: report events/s and gzip-vs-ascii size where available."""
    out = str(tmp_path / "perf.hepmc3")
    n = 2000
    trees = [_tree(dc, dc.ParticleType.NuMu, 1.5) for _ in range(n)]
    t0 = time.perf_counter()
    _write_or_skip(dc, trees, out)
    dt = time.perf_counter() - t0
    size = os.path.getsize(out)
    rate = n / dt if dt > 0 else float("inf")
    with capsys.disabled():
        print(f"\n[hepmc3 perf] wrote {n} events in {dt*1e3:.1f} ms "
              f"({rate:.0f} ev/s), {size/1024:.1f} KiB ascii")
    assert rate > 100  # generous floor; the point is the printed number


def test_gzip_roundtrip(dc, tmp_path):
    """gzip write + auto-detected gzip read. Skips on a no-compression HepMC3 build."""
    from siren import hepmc3, _util
    out = str(tmp_path / "g.hepmc3.gz")
    opts = hepmc3.HepMC3WriterOptions()
    opts.gzip = True
    trees = [_tree(dc, dc.ParticleType.NuMu, 2.0), _tree(dc, dc.ParticleType.NuEBar, 3.0)]
    try:
        hepmc3.SaveInteractionTreesAsHepMC3(trees, out, opts)
    except RuntimeError as exc:
        msg = str(exc)
        if "compression" in msg:
            pytest.skip("HepMC3 build has no gzip/compression support")
        if "without HepMC3" in msg:
            pytest.skip("SIREN built without HepMC3 support")
        raise
    # Genuine gzip stream (magic bytes), not plaintext.
    with open(out, "rb") as f:
        assert f.read(2) == b"\x1f\x8b"
    back = _util.LoadEventsFromHepMC3(out)  # reader auto-detects the .gz suffix
    assert len(back) == len(trees)


def test_strict_reader_rejects_missing_helicity(dc, tmp_path):
    """siren.helicity is emitted unconditionally on every particle, so a written
    file with one helicity line removed is corrupt: a strict load (the default)
    raises naming the attribute, while strict=False tolerates it (0.0 fallback)."""
    from siren import hepmc3, _util
    out = str(tmp_path / "corrupt.hepmc3")
    _write_or_skip(dc, [_tree(dc, dc.ParticleType.NuMu, 2.0)], out)

    # Textually drop the first siren.helicity attribute line from the ascii file.
    with open(out) as f:
        lines = f.readlines()
    removed = 0
    kept = []
    for ln in lines:
        if "siren.helicity" in ln and removed == 0:
            removed += 1
            continue
        kept.append(ln)
    assert removed == 1, "fixture did not contain a siren.helicity line to remove"
    with open(out, "w") as f:
        f.writelines(kept)

    # Strict (default) raises, naming the missing attribute.
    with pytest.raises(RuntimeError) as exc:
        _util.LoadEventsFromHepMC3(out)  # strict defaults to True
    assert "siren.helicity" in str(exc.value)

    # Lax tolerates the corruption and still returns the tree.
    back = _util.LoadEventsFromHepMC3(out, strict=False)
    assert len(back) == 1
