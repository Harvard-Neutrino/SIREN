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
    _write_or_skip(dc, trees, out, opts)
    with pyhepmc.open(out) as f:
        attrs = f.read().run_info.attributes
    # Counts stored as metadata; accepted auto-fills to the number of trees.
    assert float(str(attrs["siren.attempted_events"])) == 1000
    assert float(str(attrs["siren.accepted_events"])) == len(trees)
    # FATX raw ingredient (siren.fatx.weight_sum) == sum of the CV weights.
    assert abs(float(str(attrs["siren.fatx.weight_sum"])) - sum(weights)) < 1e-9
    assert "siren.fatx.value" in attrs
    # Reserved per-atom key AND its G.R.6 units are gated OFF by default (rate-weight caveat).
    assert "NuHepMC.FluxAveragedTotalCrossSection" not in attrs
    assert "NuHepMC.Units.CrossSection.Unit" not in attrs


def test_fatx_reserved_key_requires_opt_in(dc, tmp_path):
    from siren import hepmc3
    out = str(tmp_path / "peratom.hepmc3")
    opts = hepmc3.HepMC3WriterOptions()
    opts.attempted_events = 1000
    opts.fatx_per_atom = True  # caller asserts the per-atom normalization holds
    _write_or_skip(dc, [_tree(dc, dc.ParticleType.NuMu, 2.0)], out, opts)
    text = open(out).read()
    # The reserved key AND its units are emitted only under the opt-in.
    assert "NuHepMC.FluxAveragedTotalCrossSection" in text
    assert "NuHepMC.Units.CrossSection.Unit" in text
    assert "NuHepMC.Units.CrossSection.TargetScale" in text


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
