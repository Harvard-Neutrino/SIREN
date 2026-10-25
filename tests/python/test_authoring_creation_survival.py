"""Public compiler witnesses for U03's physical survival contract."""

import math
import pickle

import pytest
import siren
from siren import distributions as dist, injection as inj
from test_creation_survival import setup


def facade(tmp_path, *, volume=False, selected=True):
    raw, _ = setup(tmp_path, 5.0, True, branch=0.2, selected=selected, volume=volume)
    process = raw.GetPrimaryProcess()
    vertex = siren.Vertex(
        process.primary_type,
        list(process.interactions.GetDecays()),
        decay_channels=process.interactions.GetDecayChannels(),
        distributions=list(process.distributions),
        weighting=siren.PropagatedFromCreation(),
    )
    return inj.Injector(
        events=8, detector=raw.GetDetectorModel(), primary=vertex, seed=31
    )


def test_selected_external_source_survives_archive_and_pickle(tmp_path):
    generator = facade(tmp_path)
    events = generator.generate(on_failure="raise")
    weights = inj.Weighter(generator)
    expected = 0.2 * -math.expm1(-1.0) * math.exp(-5.0) / 8
    for event in events:
        assert list(event.tree[0].record.primary_initial_position) == [0, 0, 0]
        assert weights(event) == pytest.approx(expected, rel=2e-12, abs=0)

    path = str(tmp_path / "source-origin")
    generator.save(path)
    loaded = inj.Injector.load(path)
    copied = pickle.loads(pickle.dumps(generator))
    for restored in (loaded, copied):
        assert restored.engine.GetPrimaryProcess().weighting_mode.survival_from_creation
        weighter = inj.Weighter(restored)
        for event in events:
            assert weighter(event) == pytest.approx(expected, rel=2e-12, abs=0)
    next_original = generator.generate(3, on_failure="raise")
    for restored in (loaded, copied):
        next_restored = restored.generate(3, on_failure="raise")
        assert [list(e.tree[0].record.interaction_vertex) for e in next_restored] == [
            list(e.tree[0].record.interaction_vertex) for e in next_original
        ]


@pytest.mark.parametrize("surface", ["vertex", "legacy", "raw_weighter"])
def test_lost_external_origin_is_rejected_on_every_compiler_path(tmp_path, surface):
    raw, _ = setup(tmp_path, 5.0, True, volume=True)
    process = raw.GetPrimaryProcess()
    with pytest.raises(siren.errors.ConfigurationError, match="creation point"):
        if surface == "vertex":
            facade(tmp_path, volume=True).engine
        elif surface == "legacy":
            inj.Injector(
                events=8,
                detector=raw.GetDetectorModel(),
                primary_type=process.primary_type,
                primary_interactions=list(process.interactions.GetDecays()),
                primary_injection_distributions=list(process.distributions),
                primary_weighting_mode=siren.PropagatedFromCreation(),
            ).engine
        else:
            inj.Weighter(raw).engine


def test_clipped_point_source_is_rejected(tmp_path):
    generator = facade(tmp_path)
    generator.primary.distributions = [
        dist.PrimaryMass(1.0),
        dist.Monoenergetic(math.sqrt(10)),
        dist.FixedDirection([0, 0, 1]),
        dist.PointSourcePositionDistribution([0, 0, -1e9], 1e10),
    ]
    with pytest.raises(siren.errors.ConfigurationError, match="creation point"):
        generator.engine


def test_secondary_replacement_is_transactional(tmp_path):
    generator = facade(tmp_path).engine
    primary = generator.GetPrimaryProcess()
    valid = inj.SecondaryInjectionProcess(primary.primary_type, primary.interactions)
    valid.distributions = [dist.SecondaryPhysicalVertexDistribution()]
    generator.SetSecondaryProcesses([valid])
    bad = inj.SecondaryInjectionProcess(siren.particles.NuMu, primary.interactions)
    with pytest.raises(RuntimeError):
        generator.SetSecondaryProcesses([valid, bad])
    assert list(generator.GetSecondaryProcessMap()) == [primary.primary_type]
    with pytest.raises(RuntimeError, match="Duplicate"):
        generator.SetSecondaryProcesses([valid, valid])
    assert list(generator.GetSecondaryProcessMap()) == [primary.primary_type]
