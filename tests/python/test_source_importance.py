import numpy as np
import pytest
from siren.source_importance import SourceImportanceTable


def test_exact_expectation_support_and_invalidation(tmp_path):
    definition = {
        "source": {"hash": "rows-v1"},
        "physics": {"mass": 0.03},
        "geometry": "cylinder-v1",
        "scoring": "both-chords",
    }
    table = SourceImportanceTable([2, 3, 0, 5], [0, 4, 2, 1], definition)
    truth = np.array([7, 4, 2, 1])
    positive = table.physical > 0
    assert np.all(table.proposal[positive] > 0)
    assert table.proposal.sum() == pytest.approx(1)
    correction = np.divide(
        table.physical, table.proposal, out=np.zeros(4), where=positive
    )
    assert table.proposal @ (truth * correction) == pytest.approx(
        table.physical @ truth
    )
    assert correction.max() <= 20
    path = tmp_path / "table.npz"
    table.save(path)
    np.testing.assert_array_equal(
        SourceImportanceTable.load(path, definition).proposal, table.proposal
    )
    with pytest.raises(ValueError, match="definition changed"):
        SourceImportanceTable.load(path, {**definition, "physics": {"mass": 0.06}})


def test_empty_pilot_is_baseline_and_invalid_scores_fail():
    metadata = {"source": 1, "physics": 2, "geometry": 3, "scoring": 4}
    t = SourceImportanceTable([1, 9], [0, 0], metadata)
    np.testing.assert_array_equal(t.physical, t.proposal)
    for score in [[-1, 0], [np.nan, 0], [np.inf, 1]]:
        with pytest.raises(ValueError):
            SourceImportanceTable([1, 9], score, metadata)


def test_native_distribution_owns_row_sampling_and_correction():
    import siren
    metadata = {"source": 1, "physics": 2, "geometry": 3, "scoring": 4}
    table = SourceImportanceTable([1, 2, 3], [5, 0, 2], metadata)
    keys = ['E', 'weight', 'x', 'y', 'z', 'px', 'py', 'pz', 'm']
    rows = [[e, e*.001, 0, 0, 0, 0, 0, e, 0] for e in [1., 2., 3.]]
    native = table.to_distribution(keys, rows, metadata=metadata)
    assert native.normalization == pytest.approx(.006)
    rng = siren.utilities.SIREN_random(198)
    counts = np.zeros(3)
    for _ in range(4000):
        draw = siren.dataclasses.PrimaryDistributionRecord(siren.particles.NuMu)
        native.Sample(rng, None, None, draw)
        record = siren.dataclasses.InteractionRecord()
        draw.finalize(record)
        i = int(record.interaction_parameters['PrimaryExternalDistribution_row'])
        counts[i] += 1
        assert native.GenerationProbability(None, None, record) == pytest.approx(3 * table.proposal[i])
        assert native.PhysicalDensity(None, None, record) == pytest.approx(3 * table.physical[i])
    expected = 4000 * table.proposal
    assert np.all(abs(counts - expected) < 6 * np.sqrt(expected))
    for keys, data in [(['E'], [[1], [2], [3]]),
                       (keys, rows[:2]),
                       (keys, rows[::-1])]:
        with pytest.raises(ValueError):
            table.to_distribution(keys, data, metadata=metadata)
    with pytest.raises(ValueError, match='definition changed'):
        table.to_distribution(keys, rows, metadata={**metadata, 'source': 2})


@pytest.mark.parametrize('filename', ['pilot_table', 'pilot.npz'])
def test_canonical_metadata_and_exact_path(tmp_path, filename):
    definition = {'source': {'files': {np.int64(2): 'b', 10: 'k'}},
                  'physics': np.float32(.03), 'geometry': 1, 'scoring': np.int64(4)}
    table = SourceImportanceTable([1, 2], [3, 1], definition)
    table.validate_definition(definition)
    table.to_distribution(['E', 'weight'], [[1, .01], [2, .02]], metadata=definition)
    path = tmp_path / filename
    table.save(path)
    assert path.exists()
    copy = SourceImportanceTable.load(path, definition)
    assert copy.definition_hash == table.definition_hash
    np.testing.assert_array_equal(copy.proposal, table.proposal)
    with pytest.raises(ValueError, match='collide'):
        SourceImportanceTable([1], [1], {**definition, 'source': {1: 1, '1': 2}})
