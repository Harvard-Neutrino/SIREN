"""Public beam tables retain mass shells and source identities through selection."""

import copy
import math
import pickle

import numpy as np
import pytest

import siren
from siren import dk2nu
from siren.errors import ConfigurationError
from test_dk2nu import ROT_Y90, _one_row_data


def _rows():
    data = {
        k: np.repeat(v, 5) if isinstance(v, np.ndarray) else v
        for k, v in _one_row_data().items()
    }
    data["ptype"] = np.array([211, 321, 211, 111, 211])
    data["nimpwt"] = np.array([1.0, 50.0, -1.0, 8.0, 3.0])
    data["pz"] = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    data["t0"] = np.arange(5.0) + 0.125
    return data


def _sample(dist, random):
    record = siren.dataclasses.PrimaryDistributionRecord(siren.particles.PiPlus)
    dist.Sample(random, None, None, record)
    event = siren.dataclasses.InteractionRecord()
    record.finalize(event)
    return event


def _export(data, kind, tmp_path, **kwargs):
    if kind == "direct":
        return dk2nu.dk2nu_to_primary_distribution(
            data, None, frame="detector", **kwargs
        )
    path = tmp_path / "parents.csv"
    dk2nu.dk2nu_to_csv(data, path, units_cm=False, **kwargs)
    return np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True))


@pytest.mark.parametrize(
    "pdg,mass",
    [
        (211, 0.13957039),
        (-211, 0.13957039),
        (321, 0.49368),
        (-321, 0.49368),
        (130, 0.49761),
        (111, 0.1349768),
        (13, 0.10566),
        (-13, 0.10566),
    ],
)
@pytest.mark.parametrize("kind", ["direct", "csv"])
def test_parent_mass_and_decay_energy(pdg, mass, kind, tmp_path):
    data = _one_row_data()
    data["ptype"][:] = pdg
    data["E"][:] = 90.0  # production energy must not be exported at the decay point
    table = _export(data, kind, tmp_path)
    expected = math.sqrt(0.1**2 + 0.2**2 + 1.5**2 + mass**2)
    if kind == "csv":
        assert table["m"][0] == mass
        assert table["E"][0] == pytest.approx(expected, rel=2e-15)
        assert table["weight"][0] == 0.25
    else:
        event = _sample(table, siren.utilities.SIREN_random(7))
        assert event.primary_mass == mass
        assert event.primary_momentum[0] == pytest.approx(expected, rel=2e-15)


@pytest.mark.parametrize("kind", ["direct", "csv"])
def test_unknown_mass_fails_before_export(kind, tmp_path):
    data = _one_row_data()
    data["ptype"][:] = 8912345
    path = tmp_path / "parents.csv"
    path.write_text("keep existing file\n")
    with pytest.raises(ConfigurationError, match="8912345"):
        _export(data, kind, tmp_path)
    assert path.read_text() == "keep existing file\n"


@pytest.mark.parametrize("kind", ["direct", "csv"])
@pytest.mark.parametrize("registered", [False, True])
def test_explicit_and_registered_mass(kind, registered, tmp_path):
    data = _one_row_data()
    data["ptype"][:] = 8912346
    if registered:
        siren.particles.define("SIR11Parent", 8912346, 0.75)
    overrides = {} if registered else {8912346: 0.625}
    table = _export(data, kind, tmp_path, masses=overrides)
    expected = 0.75 if registered else 0.625
    mass = (
        table["m"][0]
        if kind == "csv"
        else _sample(table, siren.utilities.SIREN_random(7)).primary_mass
    )
    assert mass == expected


@pytest.mark.parametrize("kind", ["direct", "csv"])
def test_mass_override_takes_precedence(kind, tmp_path):
    table = _export(_one_row_data(), kind, tmp_path, masses={211: np.float32(0.8)})
    mass = (
        table["m"][0]
        if kind == "csv"
        else _sample(table, siren.utilities.SIREN_random(7)).primary_mass
    )
    assert mass == float(np.float32(0.8))


@pytest.mark.parametrize("kind", ["direct", "csv"])
@pytest.mark.parametrize(
    "masses",
    [
        {211: -1},
        {211: np.nan},
        {211: np.inf},
        {211: "invalid"},
        {"211": 0.2},
        {True: 0.2},
        [0.2],
    ],
)
def test_invalid_mass_override(kind, masses, tmp_path):
    with pytest.raises(ConfigurationError, match="[Mm]ass|PDG"):
        _export(_one_row_data(), kind, tmp_path, masses=masses)


@pytest.mark.parametrize("kind", ["direct", "csv"])
@pytest.mark.parametrize("pot", [0.0, -1.0, np.inf, np.nan])
def test_invalid_exposure(kind, pot, tmp_path):
    data = _one_row_data()
    data["pot"] = pot
    with pytest.raises(ConfigurationError, match="POT"):
        _export(data, kind, tmp_path)


@pytest.mark.parametrize("kind", ["direct", "csv"])
@pytest.mark.parametrize(
    "column,values",
    [
        ("weight", [1]),
        ("x0", [1]),
        ("pol_x", [1]),
        ("PrimaryExternalDistribution_row", [1]),
        ("PrimaryExternalDistribution_gen_prob", [1]),
        ("bad,name", [1]),
        ("bad\nname", [1]),
        ("id", [2**53 + 1]),
        ("id", np.array([2**63 + 1], dtype=np.uint64)),
        ("id", [np.inf]),
        ("id", [np.nan]),
        ("id", [1, 2]),
        ("id", [[1]]),
        ("id", 1),
        ("id", ["text"]),
        ("id", [1j]),
    ],
)
def test_invalid_metadata_is_not_silently_coerced(kind, column, values, tmp_path):
    with pytest.raises(ConfigurationError, match="column|columns"):
        _export(_one_row_data(), kind, tmp_path, extra_columns={column: values})


@pytest.mark.parametrize("kind", ["direct", "csv"])
def test_metadata_filters_and_large_ids_round_trip(kind, tmp_path):
    data = _rows()
    ids = np.array([2**53 - 1, 2**53 - 2, 2**53 - 3, 2**53 - 4, 2**53], dtype=np.int64)
    extra = {
        "source_id": ids,
        "inclusion_probability": np.array([0.25, 0.5, np.nan, 0.75, 1.0]),
    }
    table = _export(data, kind, tmp_path, parent_pdg=211, extra_columns=extra)
    if kind == "csv":
        np.testing.assert_array_equal(table["source_id"], ids[[0, 4]])
        np.testing.assert_array_equal(table["weight"], [0.5, 1.5])
        np.testing.assert_array_equal(table["t0"], data["t0"][[0, 4]])
        return
    seen = set()
    random = siren.utilities.SIREN_random(18)
    for _ in range(40):
        event = _sample(table, random)
        params = dict(event.interaction_parameters)
        original = {int(ids[0]): 0, int(ids[4]): 4}[int(params["source_id"])]
        seen.add(original)
        assert params["weight"] == data["nimpwt"][original] / data["pot"]
        assert (
            params["inclusion_probability"] == extra["inclusion_probability"][original]
        )
        assert event.primary_initial_time == data["t0"][original]
    assert seen == {0, 4}
    assert (
        table.normalization == 2.0
    )  # metadata must not apply another 1/inclusion factor


@pytest.mark.parametrize("copy_kind", ["original", "archive", "pickle_rejected"])
def test_rotated_table_bias_metadata_and_persistence(copy_kind, tmp_path):
    data = _rows()
    ids = np.arange(5) + 100
    bias = np.array([1.0, 500.0, 9.0, 2.0, 4.0])
    extra = {
        "beam_event_id": ids,
        "beam_file_id": np.array([0, 0, 1, 1, 1]),
        "source_sampling_probability": np.array([0.2, 0.3, 0.4, 0.5, 0.6]),
    }
    dist = dk2nu.dk2nu_to_primary_distribution(
        data,
        None,
        parent_pdg=[211, 111],
        frame=dk2nu.FrameTransform(ROT_Y90, [3.0, 4.0, 5.0], target="detector"),
        sampling_bias=bias,
        extra_columns=extra,
    )
    if copy_kind == "archive":
        ptype = siren.particles.NuMu
        collection = siren.interactions.InteractionCollection(
            ptype, [siren.interactions.DummyCrossSection()]
        )
        primary = siren.injection.PrimaryInjectionProcess(ptype, collection)
        primary.distributions = [dist]
        injector = siren.injection._Injector(
            1, siren.detector.DetectorModel(), primary, siren.utilities.SIREN_random(0)
        )
        path = str(tmp_path / "table-injector")
        injector.SaveInjector(path)
        restored = siren.injection._Injector(1, path, siren.utilities.SIREN_random(0))
        dist = restored.GetPrimaryProcess().distributions[0]
    if copy_kind == "pickle_rejected":
        # Native distributions persist through injector archives; direct pickle
        # must reject explicitly, rather than losing any table state.
        with pytest.raises(
            RuntimeError, match=r"Cannot pickle a C\+\+-defined instance"
        ):
            pickle.dumps(dist)
    selected = np.array([0, 3, 4])
    sw = bias[selected]
    weights = data["nimpwt"][selected] / data["pot"]
    random = siren.utilities.SIREN_random(33)
    reference = siren.utilities.SIREN_random(33)
    seen = set()
    for _ in range(120):
        j = np.searchsorted(np.cumsum(sw) / sum(sw), reference.Uniform(0.0, 1.0))
        original = selected[j]
        seen.add(original)
        event = _sample(dist, random)
        params = dict(event.interaction_parameters)
        assert params["beam_event_id"] == ids[original]
        assert params["beam_file_id"] == extra["beam_file_id"][original]
        assert (
            params["source_sampling_probability"]
            == extra["source_sampling_probability"][original]
        )
        assert params["weight"] == weights[j]
        assert event.primary_initial_time == data["t0"][original]
        assert event.interaction_vertex == pytest.approx(
            ROT_Y90 @ [1.0, -0.5, 30.0] + [3, 4, 5]
        )
        assert list(event.primary_momentum)[1:] == pytest.approx(
            ROT_Y90 @ [0.1, 0.2, data["pz"][original]]
        )
        assert [params[k] for k in ("pol_x", "pol_y", "pol_z")] == pytest.approx(
            ROT_Y90 @ [0.6, 0, 0.8]
        )
        assert dist.GenerationProbability(None, None, event) == pytest.approx(
            3 * sw[j] / sum(sw)
        )
        assert dist.PhysicalDensity(None, None, event) == pytest.approx(
            3 * weights[j] / sum(weights)
        )
    assert seen == {0, 3, 4}
    assert dist.normalization == sum(weights)


@pytest.mark.parametrize(
    "bias",
    [
        [0, 1],
        [-1, 1],
        [np.nan, 1],
        [np.inf, 1],
        [0, 0],
        [1e308, 1e308],
        1.0,
        [[1], [1]],
    ],
)
def test_invalid_sampling_bias_is_rejected(bias):
    data = {
        k: np.repeat(v, 2) if isinstance(v, np.ndarray) else v
        for k, v in _one_row_data().items()
    }
    with pytest.raises(ConfigurationError, match="sampling_bias"):
        dk2nu.dk2nu_to_primary_distribution(
            data, None, frame="detector", sampling_bias=bias
        )


def test_zero_bias_is_allowed_only_for_zero_physical_weight():
    data = {
        k: np.repeat(v, 2) if isinstance(v, np.ndarray) else v
        for k, v in _one_row_data().items()
    }
    data["nimpwt"][0] = 0.0
    dist = dk2nu.dk2nu_to_primary_distribution(
        data, None, frame="detector", sampling_bias=[0, 1]
    )
    event = _sample(dist, siren.utilities.SIREN_random(7))
    assert event.interaction_parameters["PrimaryExternalDistribution_row"] == 1


def test_callable_bias_receives_selected_original_kinematics():
    data = _rows()
    received = []

    def bias(*args):
        received.append(args)
        return 1.0

    dk2nu.dk2nu_to_primary_distribution(
        data, None, parent_pdg=211, frame="detector", sampling_bias=bias
    )
    for actual, name in zip(received[0], ("E", "px", "py", "pz", "vx", "vy", "vz")):
        np.testing.assert_array_equal(actual, data[name][[0, 4]])


def test_csv_matches_direct_native_table(tmp_path):
    data = _one_row_data()
    extra = {"event_id": np.array([2**53 - 1]), "source_probability": np.array([0.125])}
    # CSV's legacy callback transforms positions only. Express all vectors in
    # the output basis before export; the direct builder performs this rotation.
    rotated = copy.deepcopy(data)
    for group in [("px", "py", "pz"), ("pol_x", "pol_y", "pol_z")]:
        values = np.column_stack([data[k] for k in group]) @ ROT_Y90.T
        for j, k in enumerate(group):
            rotated[k] = values[:, j]

    def positions(x, y, z):
        return (np.column_stack([x, y, z]) @ ROT_Y90.T + [300, 400, 500]).T

    path = tmp_path / "parents.csv"
    dk2nu.dk2nu_to_csv(
        rotated, path, position_transform=positions, units_cm=False, extra_columns=extra
    )
    csv = siren.distributions.PrimaryExternalDistribution(str(path))
    direct = dk2nu.dk2nu_to_primary_distribution(
        data,
        None,
        frame=dk2nu.FrameTransform(ROT_Y90, [3, 4, 5], target="detector"),
        extra_columns=extra,
    )
    a = _sample(csv, siren.utilities.SIREN_random(7))
    b = _sample(direct, siren.utilities.SIREN_random(7))
    assert a.primary_mass == b.primary_mass
    np.testing.assert_array_equal(a.primary_momentum, b.primary_momentum)
    np.testing.assert_array_equal(a.interaction_vertex, b.interaction_vertex)
    assert a.primary_initial_time == b.primary_initial_time
    assert {k: v for k, v in a.interaction_parameters.items() if k != "nimpwt"} == dict(
        b.interaction_parameters
    )
    assert csv.normalization == direct.normalization
    assert csv.PhysicalDensity(None, None, a) == direct.PhysicalDensity(None, None, b)
    raw = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True))
    assert raw["E"][0] == pytest.approx(b.primary_momentum[0], rel=2e-15)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_csv_promotes_momentum_before_computing_energy(dtype, tmp_path):
    data = _one_row_data()
    for name in ("px", "py", "pz"):
        data[name] = data[name].astype(dtype)
    table = _export(data, "csv", tmp_path)
    expected = math.sqrt(
        sum(float(data[k][0]) ** 2 for k in ("px", "py", "pz")) + 0.13957039**2
    )
    assert table["E"][0] == pytest.approx(expected, rel=2e-15)
