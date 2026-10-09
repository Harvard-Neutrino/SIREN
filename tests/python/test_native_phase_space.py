"""Native proposals, independent normalization, and archive continuation."""
import math
import pickle

import pytest

import siren
from siren import dataclasses as dc, distributions as dist, injection as inj
from siren import interactions as xs, utilities
from siren.math import Vector3D

P = dc.Particle.ParticleType


def decay_setup(n=256, mixture=True):
    detector = siren.detector.DetectorModel()
    model = xs.HNLDipoleDecay(1.0, 1e-6, xs.HNLDipoleDecay.ChiralNature.Dirac)
    collection = xs.InteractionCollection(P.N4, [model])
    process = inj.PrimaryInjectionProcess()
    process.primary_type = P.N4
    process.interactions = collection
    process.weighting_mode = inj.VertexWeightingMode.Fixed()
    process.distributions = [
        dist.PrimaryMass(1.0), dist.Monoenergetic(2.0),
        dist.PrimaryNeutrinoHelicityDistribution(),
        dist.FixedDirection(Vector3D(0, 0, 1)),
        dist.SphereVolumePositionDistribution(siren.geometry.Sphere(1.0, 0.0)),
    ]
    if mixture:
        for signature in model.GetPossibleSignaturesFromParent(P.N4):
            proposal = inj.MultiChannelPhaseSpace([
                inj.PhysicalDecayChannel(model, signature),
                inj.Isotropic2BodyChannel(),
            ], [0.3, 0.7])
            process.SetPhaseSpace(signature, proposal)
    generator = inj._Injector(n, detector, process, utilities.SIREN_random(6217))
    physical = inj.PhysicalProcess()
    physical.primary_type = P.N4
    physical.interactions = collection
    physical.weighting_mode = inj.VertexWeightingMode.Fixed()
    physical.distributions = process.distributions
    weighter = inj._Weighter([generator], detector, physical)
    return generator, weighter


def record_state(event):
    return [(r.signature, list(r.primary_momentum), list(r.secondary_momenta),
             list(r.secondary_masses), list(r.secondary_helicities),
             list(r.interaction_vertex), r.interaction_time)
            for r in (datum.record for datum in event.tree)]


def rest_cosine(record):
    gamma_index = list(record.signature.secondary_types).index(P.Gamma)
    energy, _, _, pz = record.secondary_momenta[gamma_index]
    # Parent is (2, 0, 0, sqrt(3)); inverse boost is independent of SIREN helpers.
    energy_rest = 2.0 * energy - math.sqrt(3.0) * pz
    pz_rest = 2.0 * pz - math.sqrt(3.0) * energy
    return pz_rest / energy_rest


def test_decay_mixture_matches_analytic_angular_density():
    generator, weighter = decay_setup(n=4096)
    events = [generator.GenerateEvent() for _ in range(4096)]
    assert all(len(event.tree) == 1 for event in events)
    assert generator.InjectionAttempts() == 4096
    weights = []
    cosines = []
    for event in events:
        record = event.tree[0].record
        cosine = rest_cosine(record)
        alpha = -math.copysign(1.0, record.primary_helicity)
        physical = (1.0 + alpha * cosine) / (4.0 * math.pi)
        proposal = 0.3 * physical + 0.7 / (4.0 * math.pi)
        weight = weighter.EventWeight(event)
        assert weight == pytest.approx(physical / proposal / 4096, rel=2e-12, abs=0)
        breakdown = weighter.EventWeightWithBreakdown(event)
        assert breakdown.total == pytest.approx(weight, rel=2e-12, abs=0)
        for component in range(4):
            assert sum(p[component] for p in record.secondary_momenta) == pytest.approx(
                record.primary_momentum[component], rel=2e-12, abs=2e-14)
        weights.append(weight)
        cosines.append(alpha * cosine)
    # Integral of f is one; its first angular moment is 1/3.
    assert sum(weights) == pytest.approx(1.0, abs=0.035)
    assert sum(w * c for w, c in zip(weights, cosines)) == pytest.approx(1 / 3, abs=0.035)


def test_decay_archive_and_pickle_continue_rng_and_weights(tmp_path):
    generator, weighter = decay_setup()
    for _ in range(7):
        generator.GenerateEvent()
    archive = str(tmp_path / 'decay.injector')
    generator.SaveInjector(archive)
    restored = inj._Injector(256, archive, utilities.SIREN_random(999))
    pickled = pickle.loads(pickle.dumps(generator))
    weight_archive = str(tmp_path / 'decay')
    weighter.SaveWeighter(weight_archive)
    restored_weighter = inj._Weighter([], weight_archive)
    for _ in range(12):
        original = generator.GenerateEvent()
        loaded = restored.GenerateEvent()
        copied = pickled.GenerateEvent()
        assert record_state(original) == record_state(loaded) == record_state(copied)
    assert restored.InjectionAttempts() == generator.InjectionAttempts() == 19
    # The archived weighter owns its own injector snapshot (seven attempts).
    event = restored.GenerateEvent()
    archived_weight = restored_weighter.EventWeight(event)
    expected_snapshot_weight = weighter.EventWeight(event) * 19 / 7
    assert archived_weight == pytest.approx(expected_snapshot_weight, rel=3e-12, abs=0)


def test_failed_load_preserves_live_configuration(tmp_path):
    generator, _ = decay_setup()
    twin = pickle.loads(pickle.dumps(generator))
    broken = tmp_path / 'broken.injector'
    broken.write_bytes(b'not an archive')
    with pytest.raises(RuntimeError, match='archive'):
        generator.LoadInjector(str(broken))
    assert record_state(generator.GenerateEvent()) == record_state(twin.GenerateEvent())


def test_native_archive_rejects_live_stopping_callback(tmp_path):
    generator, _ = decay_setup()
    generator.SetStoppingCondition(lambda *_: True)
    archive = tmp_path / 'existing.injector'
    archive.write_bytes(b'preserve me')
    with pytest.raises(RuntimeError, match='stopping-condition'):
        generator.SaveInjector(str(archive))
    assert archive.read_bytes() == b'preserve me'
    with pytest.raises(RuntimeError, match='stopping-condition'):
        pickle.dumps(generator)


def test_scattering_proposal_roundtrip_keeps_total_rate(tmp_path):
    from test_trivial_cross_section import (
        detector_model, _constant_cross_section, _make_injector, _physical_process)
    detector = detector_model.__wrapped__()
    model = _constant_cross_section()
    generator, distributions = _make_injector(detector, model, n_inject=64)
    process = generator.GetPrimaryProcess()
    for signature in model.GetPossibleSignaturesFromParents(P.NuMu, P.Nucleon):
        process.SetPhaseSpace(signature, inj.MultiChannelPhaseSpace([
            inj.PhysicalCrossSectionChannel(model, signature),
            inj.PhysicalCrossSectionChannel(model, signature)], [0.2, 0.8]))
    weighter = inj._Weighter([generator], detector, _physical_process(model, distributions))
    archive = str(tmp_path / 'scatter.injector')
    generator.SaveInjector(archive)
    loaded = inj._Injector(64, archive, utilities.SIREN_random(0))
    events = [generator.GenerateEvent() for _ in range(64)]
    for event in events:
        assert record_state(event) == record_state(loaded.GenerateEvent())
    position = siren.detector.DetectorPosition(Vector3D(*events[0].tree[0].record.interaction_vertex))
    density = detector.GetParticleDensity(position, P.Nucleon)
    # A 25 m uniform thin target: convert path length to cm explicitly.
    probability = -math.expm1(-density * 25.0 * 100.0 * 1e-38)
    assert sum(weighter.EventWeight(event) for event in events) == pytest.approx(
        probability, rel=2e-10, abs=0)


def test_isotropic_proposal_rejects_nonphysical_parent_and_preserves_record():
    record = dc.InteractionRecord()
    record.signature.secondary_types = [P.Gamma, P.NuMu]
    record.primary_mass = 1.0
    record.primary_momentum = [0.5, 0, 0, 0]
    record.secondary_masses = [0.0, 0.0]
    record.secondary_momenta = [[11, 12, 13, 14], [21, 22, 23, 24]]
    before = list(record.secondary_momenta)
    channel = inj.Isotropic2BodyChannel()
    assert channel.Density(None, record) == 0
    with pytest.raises(RuntimeError):
        channel.Sample(utilities.SIREN_random(0), None, record)
    assert list(record.secondary_momenta) == before
