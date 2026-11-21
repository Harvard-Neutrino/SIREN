"""Weight-policy + provenance tests for the HepMC3 export path (R2 stage P6).

The production defect these guard against: SaveEvents(..., hepmc3=True) used to
write a CV of 1.0 for every event because nothing populated tree.header.weights
before the file was written. The policy now computes the event weight once and
populates the headers before any file is saved, and stamps a siren.weights_state
provenance marker so a NuHepMC reader is never told an unweighted file is a
conforming, rate-normalized one.

Uses the same data-free DummyCrossSection injector+weighter fixture as
test_hepmc3_closure.py. Skips cleanly without HepMC3 support or when the CCM
detector material files are unavailable.
"""
import pytest
import numpy as np

siren = pytest.importorskip("siren")


def _build():
    """Build (weighter, events) for a data-free DummyCrossSection setup whose
    weight depends on energy and helicity (so a CV of 1.0 would be obviously
    wrong)."""
    import os
    from siren import dataclasses as dc
    from siren import injection
    from siren import interactions
    from siren import distributions
    from siren import detector
    from siren import math as smath
    from siren import utilities
    from siren import _util

    NuMu = dc.Particle.ParticleType.NuMu

    try:
        dm = detector.DetectorModel()
        det_dir = _util.get_detector_model_path("CCM")
        dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
        dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    except Exception as e:  # detector model files not available in this env
        pytest.skip(f"CCM detector model unavailable: {e}")

    xs = interactions.DummyCrossSection()
    int_col = interactions.InteractionCollection(NuMu, [xs])

    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = NuMu
    primary_inj.interactions = int_col
    primary_inj.distributions = [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]

    primary_phys = injection.PhysicalProcess()
    primary_phys.primary_type = NuMu
    primary_phys.interactions = int_col
    primary_phys.distributions = [
        distributions.PrimaryMass(0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
    ]

    n_inject = 40
    rand = utilities.SIREN_random(1234)
    inj = injection._Injector(n_inject, dm, primary_inj, rand)
    weighter = injection._Weighter([inj], dm, primary_phys)

    events = []
    for _ in range(n_inject):
        try:
            ev = inj.GenerateEvent()
        except RuntimeError:
            break
        if len(ev.tree) > 0:
            events.append(ev)
    if not events:
        pytest.skip("no non-empty events generated")

    return weighter, events


def _skip_if_no_hepmc3(callable_):
    try:
        return callable_()
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise


def test_production_path_cv_regression(tmp_path):
    """SaveEvents(hepmc3=True) under "auto" writes the weighter CV per event, not
    1.0. The restored header CV must equal the weighter's EventWeight of the
    original tree, and the ascii marker must read "computed"."""
    from siren import _util, hepmc3 as sio

    weighter, events = _build()
    expected = [weighter.EventWeight(ev) for ev in events]
    for w in expected:
        assert np.isfinite(w) and w > 0
    # A meaningful test needs weights that are not the 1.0 placeholder.
    assert any(abs(w - 1.0) > 1e-6 for w in expected)

    out = str(tmp_path / "prod")
    _skip_if_no_hepmc3(lambda: _util.SaveEvents(
        events,
        weighter=weighter,
        gen_times=[0.0] * len(events),
        save_hdf5=False,
        save_parquet=False,
        save_siren_events=False,
        save_hepmc3=True,
        hepmc3_weights="auto",
        output_filename=out))

    restored = sio.LoadInteractionTreesFromHepMC3(out + ".hepmc3")
    assert len(restored) == len(events)
    for tree, want in zip(restored, expected):
        assert len(tree.header.weights) >= 1
        assert tree.header.weights[0] == pytest.approx(want, rel=1e-9), \
            (tree.header.weights[0], want)

    text = open(out + ".hepmc3").read()
    assert "siren.weights_state" in text
    assert "computed" in text


def test_unweighted_closure(tmp_path):
    """hepmc3_weights="none" writes a plain HepMC3 file: state "unweighted", no
    NuHepMC.Version, no siren.fatx.value/weight_sum, but the accepted count
    survives (later pooling needs it). Reweighting the reconstructed trees with
    the live weighter must reproduce the original weights (comparing recomputed
    values, not the padded placeholder CVs that round-trip in the header).

    (The attempted-events count is only known to the controller/injector, not to
    _util.SaveEvents; that both counts survive an unweighted write is covered at
    the writer level by test_hepmc3_metadata.test_unweighted_mode_omits_...)."""
    from siren import _util, hepmc3 as sio

    weighter, events = _build()
    weights_before = [weighter.EventWeight(ev) for ev in events]

    out = str(tmp_path / "unw")
    _skip_if_no_hepmc3(lambda: _util.SaveEvents(
        events,
        weighter=weighter,
        gen_times=[0.0] * len(events),
        save_hdf5=False,
        save_parquet=False,
        save_siren_events=False,
        save_hepmc3=True,
        hepmc3_weights="none",
        output_filename=out))

    text = open(out + ".hepmc3").read()
    assert "siren.weights_state" in text and "unweighted" in text
    assert "NuHepMC.Version" not in text
    assert "siren.fatx.value" not in text
    assert "siren.fatx.weight_sum" not in text
    assert "siren.accepted_events" in text

    loaded = sio.LoadInteractionTreesFromHepMC3(out + ".hepmc3")
    assert len(loaded) == len(events)
    weights_after = [weighter.EventWeight(ev) for ev in loaded]
    for before, after in zip(weights_before, weights_after):
        assert after == pytest.approx(before, rel=1e-9), (before, after)
