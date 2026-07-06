"""Set-combination (pooling) closure + FATX-merge checks for combine_and_export_hepmc3.

Two simulation sets are generated with genuinely distinct injectors that have
OVERLAPPING support (two different power-law energy spectra over overlapping
ranges, same physical process + detector). Each set is saved as the native
.siren_events + .siren_weighter pair. combine_and_export_hepmc3 pools them into
one HepMC3 file.

The essential guarantee is pooled-weight closure: the C++ Weighter is natively
pooled (EventWeight = 1 / sum_i(N_i G_i / P_i) over its injector list), so
combining sets means building ONE Weighter whose injector list is the union and
re-running EventWeight on every event. The distinct-injector requirement is
load-bearing: with a single shared injector the cross-terms are trivial and the
test could not catch an implementation that wrongly reuses per-set-only weights.

Skips cleanly when SIREN was built without HepMC3, or when the CCM detector
material files are unavailable.
"""
import os

import pytest
import numpy as np

siren = pytest.importorskip("siren")

# 1 GeV^-2 -> pb, matching HepMC3Writer.cxx (kGeVm2_to_pb).
K_GEVM2_TO_PB = 3.894e8


def _detector():
    from siren import detector
    from siren import _util
    try:
        dm = detector.DetectorModel()
        det_dir = _util.get_detector_model_path("CCM")
        dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
        dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    except Exception as e:  # detector model files not available in this env
        pytest.skip(f"CCM detector model unavailable: {e}")
    return dm


def _make_injector(dm, power_index, e_min, e_max, seed, n_inject, primary_type=None):
    """A single-vertex injector whose only weight-relevant knob is the power law."""
    from siren import dataclasses as dc
    from siren import injection
    from siren import interactions
    from siren import distributions
    from siren import math as smath
    from siren import utilities

    ptype = dc.Particle.ParticleType.NuMu if primary_type is None else primary_type

    xs = interactions.DummyCrossSection()
    int_col = interactions.InteractionCollection(ptype, [xs])

    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = ptype
    primary_inj.interactions = int_col
    primary_inj.distributions = [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(power_index, e_min, e_max),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]

    rand = utilities.SIREN_random(seed)
    inj = injection._Injector(n_inject, dm, primary_inj, rand)
    return inj


def _shared_physical_process(dm, primary_type=None):
    from siren import dataclasses as dc
    from siren import injection
    from siren import interactions
    from siren import distributions

    ptype = dc.Particle.ParticleType.NuMu if primary_type is None else primary_type
    xs = interactions.DummyCrossSection()
    int_col = interactions.InteractionCollection(ptype, [xs])
    primary_phys = injection.PhysicalProcess()
    primary_phys.primary_type = ptype
    primary_phys.interactions = int_col
    primary_phys.distributions = [
        distributions.PrimaryMass(0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
    ]
    return primary_phys


def _generate(inj, n):
    events = []
    for _ in range(n):
        try:
            ev = inj.GenerateEvent()
        except RuntimeError:
            break
        if len(ev.tree) > 0:
            events.append(ev)
    return events


def _save_set(events, base):
    """Write <base>.siren_events + <base>.siren_weighter for a set."""
    from siren import dataclasses as dc
    dc.SaveInteractionTrees(events, base)


def _build_two_sets(tmp_path):
    """Two overlapping-support injectors + the shared physical process.

    Returns (injA, injB, primary_phys, dm, baseA, baseB, nA, nB).
    Each set is saved to disk as .siren_events + .siren_weighter.
    """
    from siren import injection

    dm = _detector()
    n_inject = 60

    # Two distinct power laws over OVERLAPPING energy ranges.
    injA = _make_injector(dm, power_index=1.5, e_min=0.5, e_max=4.0, seed=101, n_inject=n_inject)
    injB = _make_injector(dm, power_index=3.0, e_min=1.0, e_max=5.0, seed=202, n_inject=n_inject)

    physA = _shared_physical_process(dm)
    physB = _shared_physical_process(dm)

    eventsA = _generate(injA, n_inject)
    eventsB = _generate(injB, n_inject)
    if not eventsA or not eventsB:
        pytest.skip("no non-empty events generated for one of the sets")

    baseA = str(tmp_path / "setA")
    baseB = str(tmp_path / "setB")

    # Save each set's events and its single-injector weighter (SaveEvents naming
    # convention: <base>.siren_events alongside <base>.siren_weighter).
    _save_set(eventsA, baseA)
    weighterA = injection._Weighter([injA], dm, physA)
    weighterA.SaveWeighter(baseA)

    _save_set(eventsB, baseB)
    weighterB = injection._Weighter([injB], dm, physB)
    weighterB.SaveWeighter(baseB)

    return injA, injB, physA, dm, baseA, baseB, len(eventsA), len(eventsB)


def test_two_set_pooled_closure(tmp_path):
    """The essential test: combined CVs equal a fresh pooled Weighter's weights.

    Ground truth is a Weighter built directly over [injA, injB] with the shared
    physical config; combine_and_export_hepmc3 must reproduce it event-for-event.
    """
    from siren import injection, dataclasses as dc, _util

    injA, injB, physA, dm, baseA, baseB, nA, nB = _build_two_sets(tmp_path)

    out = str(tmp_path / "combined.hepmc3")
    try:
        out_path, pooled_weights = _util.combine_and_export_hepmc3([baseA, baseB], out)
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise

    # Ground truth: a fresh pooled Weighter over the two live injectors.
    ground = injection._Weighter([injA, injB], dm, physA)

    # Reload the flat tree order combine used (set A then set B) and compare.
    treesA = dc.LoadInteractionTrees(baseA)
    treesB = dc.LoadInteractionTrees(baseB)
    flat = list(treesA) + list(treesB)
    assert len(flat) == len(pooled_weights)
    assert len(pooled_weights) == nA + nB

    for tree, w_comb in zip(flat, pooled_weights):
        w_truth = ground.EventWeight(tree)
        assert np.isfinite(w_truth) and w_truth > 0
        assert w_comb == pytest.approx(w_truth, rel=1e-9), (w_comb, w_truth)


def test_two_set_pooling_is_not_single_set(tmp_path):
    """Sanity: the pooled weight actually differs from either single-set weight
    (overlapping support -> live cross-terms). Otherwise the closure above could
    pass for a broken combiner that reused per-set weights."""
    from siren import injection, dataclasses as dc, _util

    injA, injB, physA, dm, baseA, baseB, nA, nB = _build_two_sets(tmp_path)

    out = str(tmp_path / "combined.hepmc3")
    try:
        _out_path, pooled_weights = _util.combine_and_export_hepmc3([baseA, baseB], out)
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise

    single_A = injection._Weighter([injA], dm, physA)
    treesA = dc.LoadInteractionTrees(baseA)
    # For a set-A event, the pooled weight must be strictly smaller than the
    # single-injector-A weight (the injector-B term adds to the denominator).
    differed = False
    for tree, w_comb in zip(treesA, pooled_weights[:len(treesA)]):
        w_single = single_A.EventWeight(tree)
        if abs(w_comb - w_single) > 1e-12 * max(1.0, abs(w_single)):
            differed = True
            assert w_comb < w_single  # pooling adds a positive term to inv_weight
    assert differed, "pooled weights never differed from single-set weights"


def _fatx_ingredients(hepmc3_path):
    """Read (weight_sum, norm, value) from a written HepMC3 ascii file."""
    ws = norm = value = None
    attempted = accepted = None
    with open(hepmc3_path) as f:
        for line in f:
            if "siren.fatx.weight_sum" in line:
                ws = _last_float(line)
            elif "siren.fatx.value" in line:
                value = _last_float(line)
            elif "siren.attempted_events" in line:
                attempted = _last_float(line)
            elif "siren.accepted_events" in line:
                accepted = _last_float(line)
    norm = attempted if (attempted is not None and attempted > 0) else accepted
    return ws, norm, value


def _last_float(line):
    """Extract the trailing float token from a HepMC3 attribute ascii line."""
    tok = line.strip().split()[-1]
    return float(tok)


def test_fatx_merge_single_flux(tmp_path):
    """FATX poolability: the two per-set files' weight_sum ingredients pool to the
    combined file's siren.fatx.value via k*(ws1+ws2).

    Both per-set files and the combined file carry the SAME pooled central
    values, so the raw ingredient siren.fatx.weight_sum is additive across files
    and rescaling their sum by k reproduces the combined estimate exactly (the CV
    weight already carries 1/EventsToInject, so there is no count denominator).
    """
    from siren import hepmc3, dataclasses as dc, _util

    _injA, _injB, _phys, _dm, baseA, baseB, nA, nB = _build_two_sets(tmp_path)

    att_A, att_B = 1000, 3000

    out = str(tmp_path / "combined.hepmc3")
    sets = [
        {"events": baseA, "weighter": baseA, "attempted": att_A},
        {"events": baseB, "weighter": baseB, "attempted": att_B},
    ]
    try:
        _out_path, pooled_weights = _util.combine_and_export_hepmc3(sets, out)
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise

    # Combined file's FATX value + weight_sum.
    ws_comb, norm_comb, value_comb = _fatx_ingredients(out)
    assert value_comb is not None
    assert norm_comb == pytest.approx(att_A + att_B)
    assert ws_comb == pytest.approx(sum(pooled_weights), rel=1e-9)

    # Write each per-set file carrying its share of the SAME pooled CVs, with its
    # own attempted count, then read the raw ingredients back.
    treesA = dc.LoadInteractionTrees(baseA)
    treesB = dc.LoadInteractionTrees(baseB)
    wA = pooled_weights[:len(treesA)]
    wB = pooled_weights[len(treesA):len(treesA) + len(treesB)]
    for tree, w in zip(treesA, wA):
        _util._set_header_cv(tree, w)
    for tree, w in zip(treesB, wB):
        _util._set_header_cv(tree, w)

    outA = str(tmp_path / "onlyA.hepmc3")
    outB = str(tmp_path / "onlyB.hepmc3")
    optsA = hepmc3.HepMC3WriterOptions(); optsA.weights_state = "computed"; optsA.attempted_events = att_A
    optsB = hepmc3.HepMC3WriterOptions(); optsB.weights_state = "computed"; optsB.attempted_events = att_B
    hepmc3.SaveInteractionTreesAsHepMC3(treesA, outA, optsA)
    hepmc3.SaveInteractionTreesAsHepMC3(treesB, outB, optsB)

    ws1, _norm1, _v1 = _fatx_ingredients(outA)
    ws2, _norm2, _v2 = _fatx_ingredients(outB)

    merged = K_GEVM2_TO_PB * (ws1 + ws2)
    assert merged == pytest.approx(value_comb, rel=1e-9), (merged, value_comb)


def test_combine_refuses_unweighted_set(tmp_path):
    """A set declaring weights_state 'unweighted' must not be silently pooled with
    a computed set; the combiner raises."""
    from siren import _util

    _injA, _injB, _phys, _dm, baseA, baseB, nA, nB = _build_two_sets(tmp_path)

    out = str(tmp_path / "combined_bad.hepmc3")
    sets = [
        {"events": baseA, "weighter": baseA},
        {"events": baseB, "weighter": baseB, "weights_state": "unweighted"},
    ]
    with pytest.raises(ValueError) as exc:
        _util.combine_and_export_hepmc3(sets, out)
    assert "unweighted" in str(exc.value)


def test_combine_requires_weighter_file(tmp_path):
    """A set without a .siren_weighter companion cannot be pooled; the combiner
    raises rather than emitting garbage weights."""
    from siren import _util

    _injA, _injB, _phys, _dm, baseA, baseB, nA, nB = _build_two_sets(tmp_path)

    # Point set B's weighter at a base with no .siren_weighter file.
    missing = str(tmp_path / "does_not_exist")
    out = str(tmp_path / "combined_missing.hepmc3")
    sets = [baseA, {"events": baseB, "weighter": missing}]
    with pytest.raises(ValueError) as exc:
        _util.combine_and_export_hepmc3(sets, out)
    assert "weighter" in str(exc.value)


def test_weighter_exposes_physical_process(tmp_path):
    """The new C++ accessor surfaces the primary/secondary physical processes so
    the combiner can build a REAL physical-side signature, not just the
    injection-side proxy."""
    from siren import injection, dataclasses as dc, _util

    _injA, _injB, _physA, _dm, baseA, _baseB, _nA, _nB = _build_two_sets(tmp_path)

    w = injection._Weighter([], baseA)
    prim = w.GetPrimaryPhysicalProcess()
    assert prim is not None
    assert int(prim.primary_type) == int(dc.Particle.ParticleType.NuMu)
    # Single-vertex config: no secondary physical processes.
    assert list(w.GetSecondaryPhysicalProcesses()) == []
    # A detector model must come back too (used for parity with the injectors).
    assert w.GetDetectorModel() is not None
    # The physical-side signature must be readable (non-None), so the guard is a
    # real check rather than degrading to the injection-side warning fallback.
    assert _util._physical_config_signature(w) is not None


def _build_mismatched_physics_sets(tmp_path):
    """Two sets with genuinely different PHYSICAL processes (NuMu vs NuE), each
    internally consistent. Pooling them is physically meaningless and must be
    refused."""
    from siren import dataclasses as dc, injection

    dm = _detector()
    n_inject = 40
    NuMu = dc.Particle.ParticleType.NuMu
    NuE = dc.Particle.ParticleType.NuE

    injA = _make_injector(dm, 1.5, 0.5, 4.0, seed=303, n_inject=n_inject, primary_type=NuMu)
    injB = _make_injector(dm, 1.5, 0.5, 4.0, seed=404, n_inject=n_inject, primary_type=NuE)
    physA = _shared_physical_process(dm, primary_type=NuMu)
    physB = _shared_physical_process(dm, primary_type=NuE)

    eventsA = _generate(injA, n_inject)
    eventsB = _generate(injB, n_inject)
    if not eventsA or not eventsB:
        pytest.skip("no non-empty events generated for one of the mismatched sets")

    baseA = str(tmp_path / "physNuMu")
    baseB = str(tmp_path / "physNuE")
    _save_set(eventsA, baseA)
    injection._Weighter([injA], dm, physA).SaveWeighter(baseA)
    _save_set(eventsB, baseB)
    injection._Weighter([injB], dm, physB).SaveWeighter(baseB)
    return baseA, baseB


def test_combine_refuses_mismatched_physical_config(tmp_path):
    """Pooling sets with mismatched physical processes raises before any HepMC3 write."""
    from siren import _util

    baseA, baseB = _build_mismatched_physics_sets(tmp_path)
    out = str(tmp_path / "combined_bad_physics.hepmc3")
    with pytest.raises(ValueError) as exc:
        _util.combine_and_export_hepmc3([baseA, baseB], out)
    msg = str(exc.value).lower()
    assert "physical" in msg
