"""Golden regression pin for the assembled on-shell Dutta-Kim chain.

The chain tests in test_dutta_kim_chain.py check the assembled
inject-and-weight product (dk2nu-free, DarkNews-free: the self-contained
VectorPortal and MesonProduction models over an SBND geometry) only for
weight finiteness, non-negativity, and internal self-consistency. This
module adds an absolute numeric pin: the full on-shell chain
(pi -> V1 -> chi -> chi' -> V1_sig -> e+ e-) run at a fixed seed must
reproduce a committed per-event weight vector, so a silent regression in
the production, decay, upscattering, or flux physics -- or in the channel
mixture that biases it -- is caught rather than passing the weaker
finiteness checks.

The chain uses no data files and no external generator, so the pin is
deterministic and machine-independent (verified run-to-run identical).
Re-bless GOLDEN_WEIGHTS only for an intentional, changelog-labeled change
to the on-shell chain physics, exactly as test_golden_regression.py's
archive is re-blessed.
"""
import math
import os

import pytest

siren = pytest.importorskip("siren")
from siren import _util
from siren.Injector import Injector
from siren.Weighter import Weighter

_SEED = 99
_N_EVENTS = 12

# Per-event central-value weights of the on-shell chain at seed 99 (see the
# module docstring for the re-bless policy).
GOLDEN_WEIGHTS = [
    4.52067640585333574e-20,
    3.23170162923930390e-20,
    3.45834698021994380e-20,
    3.75005453783721450e-20,
    4.44326497234204670e-20,
    3.88145861646823405e-20,
    3.53332366390994690e-20,
    3.74786803580858473e-20,
    4.48869762751869907e-20,
    3.46571023670208949e-20,
    4.63210154125096197e-20,
    4.84328010146990786e-20,
]


@pytest.fixture(scope="module")
def chain_module():
    return _util.load_module(
        "golden_DuttaKim_chain",
        os.path.join(_util.resource_package_dir(), "examples", "example4",
                     "DuttaKim_SBND_full_chain.py"))


def test_onshell_chain_weights_match_golden(chain_module):
    """The full on-shell chain reproduces the committed weight vector."""
    detector_model = siren.utilities.load_detector("SBN", detector="SBND")
    fiducial = siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0, 0, 0)), 4.0, 4.0, 4.0)

    models = chain_module.build_onshell_models()
    primary, secondaries = chain_module.build_vertices(models, fiducial)

    injector = Injector(detector=detector_model, primary=primary,
                        secondaries=secondaries, events=_N_EVENTS, seed=_SEED)
    events = injector.generate(_N_EVENTS, on_shortfall="warn")
    weighter = Weighter(injector, primary_physical=primary.physical)
    weights = [weighter(ev) for ev in events]

    assert len(weights) == len(GOLDEN_WEIGHTS), (
        "event count drifted: got %d, pinned %d" % (len(weights), len(GOLDEN_WEIGHTS)))
    for i, (got, want) in enumerate(zip(weights, GOLDEN_WEIGHTS)):
        assert math.isfinite(got) and got >= 0.0, "event %d: %r" % (i, got)
        assert got == pytest.approx(want, rel=1e-6), (
            "event %d weight %.17e != golden %.17e (re-bless only for a "
            "labeled physics change)" % (i, got, want))
