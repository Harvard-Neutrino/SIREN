"""Example4 decay-in-flight / decay-at-rest merging helper.

The helper is analysis-level code living beside the beam-timing example
scripts; these tests load it by path, the same way the examples do.
"""

import os

import numpy as np
import pytest

siren = pytest.importorskip("siren")
from siren import _util

beam_samples = _util.load_module(
    "test_example4_beam_samples_helper",
    os.path.join(_util.resource_package_dir(),
                 "examples", "example4", "_beam_samples.py"))


def _ke_sample(ptypes, kinetic_energies, weights, pot, extra=None):
    """Synthetic read_dk2nu-style dict with parents at given kinetic energies."""
    ptypes = np.asarray(ptypes)
    masses = np.array(
        [beam_samples.PARENT_MASSES[abs(int(p))] for p in ptypes])
    data = {
        "ptype": ptypes,
        "E": masses + np.asarray(kinetic_energies, dtype=float),
        "nimpwt": np.asarray(weights, dtype=float),
        "pot": float(pot),
    }
    if extra:
        data.update(extra)
    return data


def test_combine_dif_dar_classifies_by_kinetic_energy():
    dif = _ke_sample([211, 211, 321], [0.2, 0.01, 1.0], [1.0, 1.0, 1.0], 4.0)
    dar = _ke_sample([211, 211], [0.2, 0.01], [1.0, 1.0], 2.0)
    merged = beam_samples.combine_dif_dar(dif, dar, kinetic_energy_cut=0.05)
    # The decay-in-flight sample keeps its two rows at or above the cut,
    # the decay-at-rest sample keeps its one row below it.
    assert len(merged["E"]) == 3
    assert merged["dar"].tolist() == [False, False, True]
    assert merged["pot"] == 4.0
    assert merged["dif_pot"] == 4.0
    assert merged["dar_pot"] == 2.0


def test_combine_dif_dar_weights_are_per_sample_pot():
    dif = _ke_sample([211], [0.2], [3.0], 4.0)
    dar = _ke_sample([211], [0.01], [5.0], 2.0)
    merged = beam_samples.combine_dif_dar(dif, dar)
    # Dividing by the merged POT must recover each row's weight per POT of
    # its own sample; this is what dk2nu_to_primary_distribution computes.
    per_pot = merged["nimpwt"] / merged["pot"]
    assert per_pot[0] == pytest.approx(3.0 / 4.0)
    assert per_pot[1] == pytest.approx(5.0 / 2.0)


def test_combine_dif_dar_drops_one_sided_columns():
    dif = _ke_sample([211], [0.2], [1.0], 4.0,
                     extra={"t0": np.array([1.0])})
    dar = _ke_sample([211], [0.01], [1.0], 2.0)
    merged = beam_samples.combine_dif_dar(dif, dar)
    assert "t0" not in merged


def test_combine_dif_dar_requires_positive_pot():
    dif = _ke_sample([211], [0.2], [1.0], 4.0)
    dar = _ke_sample([211], [0.01], [1.0], 0.0)
    with pytest.raises(ValueError, match="POT"):
        beam_samples.combine_dif_dar(dif, dar)


def test_parent_kinetic_energy_rejects_unknown_species():
    data = {"ptype": np.array([9999]), "E": np.array([1.0])}
    with pytest.raises(KeyError, match="9999"):
        beam_samples.parent_kinetic_energy(data)
