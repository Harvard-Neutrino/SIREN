"""Source identity and small physical weights survive public table conversion."""

import io

import numpy as np
import pytest

import siren
from siren import dk2nu
from siren.errors import ConfigurationError


def _data(n=2):
    return {
        "ptype": np.full(n, 211),
        "E": np.full(n, 2.0),
        "px": np.full(n, 0.1),
        "py": np.full(n, 0.2),
        "pz": np.full(n, 1.5),
        "vx": np.full(n, 100.0),
        "vy": np.full(n, -50.0),
        "vz": np.full(n, 3000.0),
        "nimpwt": np.full(n, 0.5),
        "pot": 2.0,
    }


def _export(kind, data, tmp_path, **kwargs):
    if kind == "direct":
        return dk2nu.dk2nu_to_primary_distribution(
            data, None, frame="detector", **kwargs
        )
    path = tmp_path / "parents.csv"
    dk2nu.dk2nu_to_csv(data, path, units_cm=False, **kwargs)
    return siren.distributions.PrimaryExternalDistribution(str(path))


def _samples(dist):
    random = siren.utilities.SIREN_random(7)
    seen = {}
    for _ in range(32):
        record = siren.dataclasses.PrimaryDistributionRecord(siren.particles.PiPlus)
        dist.Sample(random, None, None, record)
        event = siren.dataclasses.InteractionRecord()
        record.finalize(event)
        seen[int(event.interaction_parameters["PrimaryExternalDistribution_row"])] = (
            event
        )
    return seen


@pytest.mark.parametrize("kind", ["direct", "csv"])
@pytest.mark.parametrize(
    "ids",
    [
        [np.uint64(2**53 + 1), np.int64(2**53)],
        [2**53 + 1, 1.0],
        [-(2**53) - 1, 1.0],
        [np.array(2**53 + 1), np.array(1)],
        [np.array(2**53 + 1), np.array(1.0)],
        [np.array(-(2**53) - 1), np.array(1.0)],
        np.array([np.iinfo(np.int64).min, 1], dtype=np.int64),
        np.array([np.iinfo(np.uint64).max, 1], dtype=np.uint64),
    ],
)
def test_integer_ids_are_checked_before_lossy_array_promotion(kind, ids, tmp_path):
    with pytest.raises(ConfigurationError, match="Integer extra column"):
        _export(kind, _data(), tmp_path, extra_columns={"source_id": ids})


@pytest.mark.parametrize("kind", ["direct", "csv"])
@pytest.mark.parametrize(
    "ids",
    [
        [np.uint64(2**53 + 1), np.int64(2**53)],
        [np.array(2**53 + 1), np.array(float(2**53))],
    ],
)
def test_filtered_integer_ids_do_not_invalidate_retained_metadata(kind, ids, tmp_path):
    data = _data()
    data["ptype"][0] = 111
    dist = _export(
        kind,
        data,
        tmp_path,
        parent_pdg=211,
        extra_columns={"source_id": ids},
    )
    assert _samples(dist)[0].interaction_parameters["source_id"] == 2**53


@pytest.mark.parametrize("kind", ["direct", "csv"])
def test_scalar_array_integer_ids_preserve_both_exact_boundaries(kind, tmp_path):
    ids = [np.array(2**53), np.array(-(2**53))]
    dist = _export(kind, _data(), tmp_path, extra_columns={"source_id": ids})
    assert {
        row: event.interaction_parameters["source_id"]
        for row, event in _samples(dist).items()
    } == {0: 2**53, 1: -(2**53)}


@pytest.mark.parametrize("kind", ["direct", "csv"])
def test_selected_masked_metadata_is_rejected(kind, tmp_path):
    ids = np.ma.array([101, 202], mask=[False, True])
    with pytest.raises(ConfigurationError, match="masked"):
        _export(kind, _data(), tmp_path, extra_columns={"source_id": ids})


@pytest.mark.parametrize("kind", ["direct", "csv"])
@pytest.mark.parametrize("filter_by", ["species", "importance"])
def test_discarded_masked_metadata_does_not_invalidate_retained_rows(
    kind, filter_by, tmp_path
):
    data = _data()
    ids = np.ma.array([101, 202], mask=[False, True])
    if filter_by == "species":
        data["ptype"][1] = 111
    else:
        data["nimpwt"][1] = -1
    dist = _export(
        kind, data, tmp_path, parent_pdg=211, extra_columns={"source_id": ids}
    )
    assert _samples(dist)[0].interaction_parameters["source_id"] == 101


@pytest.mark.parametrize("kind", ["direct", "csv"])
def test_unmasked_array_preserves_integer_identity(kind, tmp_path):
    ids = np.ma.array([2**53 - 1, 2**53], mask=False)
    dist = _export(kind, _data(), tmp_path, extra_columns={"source_id": ids})
    assert {
        row: event.interaction_parameters["source_id"]
        for row, event in _samples(dist).items()
    } == {0: 2**53 - 1, 1: 2**53}


@pytest.mark.parametrize("kind", ["direct", "csv"])
def test_floating_metadata_is_not_restricted_to_integer_id_range(kind, tmp_path):
    dist = _export(
        kind, _data(1), tmp_path, extra_columns={"exposure": np.array([2.0**60])}
    )
    assert _samples(dist)[0].interaction_parameters["exposure"] == 2.0**60


@pytest.mark.parametrize("kind", ["direct", "csv"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("rows", [1, 2])
def test_small_positive_physical_weights_survive_normalization(
    kind, dtype, rows, tmp_path
):
    data = _data(rows)
    data["nimpwt"] = np.array([1e-30, 1.0][:rows], dtype=dtype)
    data["pot"] = 1e20
    expected = [float(w) / float(data["pot"]) for w in data["nimpwt"]]
    dist = _export(kind, data, tmp_path)
    seen = _samples(dist)
    assert set(seen) == set(range(rows))
    for row, event in seen.items():
        assert expected[row] > 0
        assert event.interaction_parameters["weight"] == expected[row]
        assert dist.PhysicalDensity(None, None, event) == pytest.approx(
            rows * expected[row] / sum(expected), rel=2e-15, abs=0
        )
    assert dist.normalization == pytest.approx(sum(expected), rel=2e-15, abs=0)


@pytest.mark.parametrize("suffix", [".csv", ".csv.gz", ".csv.bz2"])
@pytest.mark.parametrize("string_path", [False, True])
def test_csv_suffix_does_not_enable_compression(suffix, string_path, tmp_path):
    path = tmp_path / ("parents" + suffix)
    dk2nu.dk2nu_to_csv(_data(1), str(path) if string_path else path, units_cm=False)
    dist = siren.distributions.PrimaryExternalDistribution(str(path))
    assert dist.GetPhysicalNumEvents() == 1
    assert path.read_text().startswith("E,px,py,pz,")
    assert b"\r" not in path.read_bytes()
    assert _samples(dist)[0].interaction_parameters["weight"] == 0.25


def test_invalid_csv_columns_do_not_truncate_an_existing_file(tmp_path):
    path = tmp_path / "parents.csv"
    path.write_text("preserve existing table\n")

    def malformed_positions(x, y, z):
        return np.repeat(x, 2), y, z

    with pytest.raises(ValueError):
        dk2nu.dk2nu_to_csv(_data(1), path, position_transform=malformed_positions)
    assert path.read_text() == "preserve existing table\n"


def test_csv_output_requires_a_filesystem_path():
    output = io.StringIO("preserve existing content")
    with pytest.raises(TypeError):
        dk2nu.dk2nu_to_csv(_data(1), output)
    assert output.getvalue() == "preserve existing content"
