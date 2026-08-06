"""Regression test for the DarkNews decay phase-space sample shape."""

import os
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("siren")
pytest.importorskip("DarkNews")

from siren import _util


DARKNEWS_DIR = os.path.join(
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))
        )
    ),
    "resources",
    "processes",
    "DarkNewsTables",
)


def _load_darknews_module(name):
    path = os.path.join(DARKNEWS_DIR, f"{name}.py")

    if not os.path.isfile(path):
        pytest.skip(f"{name}.py not present")

    try:
        return _util.load_module(name, path)
    except (ImportError, OSError, RuntimeError) as exc:
        pytest.skip(f"could not load {name}: {exc}")


def test_phase_space_sample_shape(monkeypatch):
    """A four-variable phase-space sample must have shape (4, 1)."""
    mod = _load_darknews_module("DarkNewsDecay")

    decay = mod.PyDarkNewsDecay(SimpleNamespace())

    # Four phase-space variables for one sampled event example.
    decay.PS_samples = np.array(
        [
            [0.1],
            [0.2],
            [0.3],
            [0.4],
        ]
    )
    decay.PS_weights = np.array([1.0])

    def check_shape(vsamples, *_args):
        assert vsamples.shape == (4, 1)
        return {}

    monkeypatch.setattr(
        mod,
        "get_decay_momenta_from_vegas_samples",
        check_shape,
    )

    record = SimpleNamespace(
        primary_momentum=[1.0, 0.0, 0.0, 1.0],
        get_secondary_particle_records=lambda: [],
    )

    random = SimpleNamespace(
        Uniform=lambda minimum, maximum: minimum
    )

    decay.SampleRecordFromDarkNews(record, random)