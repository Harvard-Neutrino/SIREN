"""Golden-physics regression pin for the assembled SIREN pipeline.

This is a fixed-seed numeric anchor. It archives weighted distributions from a
real assembled injection + weighting pipeline, and the converged per-channel
alphas of the Kleiss-Pittau optimizer, into a checked-in .npz. The tests
regenerate those arrays from scratch (same seeds) and assert they match the
archive to rtol=1e-12. The point is to make any future engine change PROVE it
left the numerics untouched: if a refactor is meant to be behavior-preserving
and this test moves, the refactor changed physics.

RE-BLESS POLICY (read before touching the archive)
--------------------------------------------------
The archive tests/python/data/golden_regression.npz is a frozen physics anchor.
It may be regenerated ONLY as part of a deliberate, changelog-labeled physics
fix -- a change that is INTENDED to alter the sampled/weighted numbers. Every
regeneration MUST ship with its own changelog entry naming the physics change
and why the golden values move. A regeneration for any other reason (a
"the test is failing, just re-bless it" reflex during an unrelated refactor)
defeats the entire purpose of this pin and is not allowed. Refactors that are
meant to preserve behavior must leave this archive byte-for-byte reproducible.

To regenerate after an approved, changelog-labeled physics change:

    python tests/python/test_golden_regression.py --regenerate

The archive is NEVER auto-regenerated from inside pytest. If it is missing the
tests FAIL loudly with the regeneration command, rather than silently minting a
new "golden" from whatever the current (possibly broken) code produces.

What is pinned
--------------
1. test_golden_kp_alphas -- a fixed-seed Kleiss-Pittau optimization of a toy
   3-channel mixture (physical/isotropic + two detector-directed channels) via
   tune.optimize_multichannel_weights, which drives the bare
   W_i = mean(w^2 g_i/g) statistic through the C++ Accumulate/UpdateWeights
   path. Archived: the converged per-channel alphas. This pins that bare
   statistic (it must never be silently replaced by an alpha-weighted variant).

ASCII only. No network. No external tables (the DarkNews path is deliberately
avoided so this pin is fully deterministic and machine-independent).
"""
from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
import pytest

import siren
from siren import dataclasses as dc
from siren import injection
from siren import interactions
from siren import distributions
from siren import detector
from siren import math as smath
from siren import utilities
from siren import _util
from siren.tune import optimize_multichannel_weights


# --------------------------------------------------------------------------- #
# Fixed configuration (seeds and sizes are part of the pin -- do not change    #
# them without a changelog-labeled physics-fix regeneration).                  #
# --------------------------------------------------------------------------- #

ARCHIVE = Path(__file__).resolve().parent / "data" / "golden_regression.npz"

_REGEN_MSG = (
    "Golden regression archive is missing:\n    {path}\n"
    "This archive is a checked-in physics anchor; pytest never creates it.\n"
    "If (and only if) you are making a deliberate, changelog-labeled physics\n"
    "change, regenerate it with:\n"
    "    python tests/python/test_golden_regression.py --regenerate\n"
    "Otherwise the archive should have been present -- check your checkout."
)

RTOL = 1e-12
ATOL = 0.0

# Part 1: assembled chain
CHAIN_SEED = 1234
N_EVENTS = 500          # accepted (non-empty) events to pin
MAX_ATTEMPTS = 2000     # attempt budget; the chain accepts ~everything, so
                        # this is a comfortable ceiling for N_EVENTS successes
E_BINS = np.linspace(0.5, 5.0, 13)     # PowerLaw support [0.5, 5.0]
R_BINS = np.linspace(0.0, 25.0, 11)    # PointSource max_dist = 25

# Part 2: Kleiss-Pittau optimizer toy
KP_SEED = 4242
KP_ITERS = 8
KP_BATCH = 3000
KP_DAMPING = 0.5
KP_MIN_WEIGHT = 0.01
KP_RULE = "alpha_sqrt_W"

_NuMu = dc.Particle.ParticleType.NuMu
_PT = siren.dataclasses.ParticleType
_M_V1 = 0.017


# --------------------------------------------------------------------------- #
# Part 2: Kleiss-Pittau optimizer toy (pins the bare W_i statistic)            #
# --------------------------------------------------------------------------- #

def _kp_make_record(m_chi, E_V1):
    pz = math.sqrt(max(E_V1 * E_V1 - _M_V1 * _M_V1, 0.0))
    r = siren.dataclasses.InteractionRecord()
    r.signature.primary_type = _PT.N4
    r.signature.target_type = _PT.Decay
    r.signature.secondary_types = [_PT.NuLight, _PT.Gamma]
    r.primary_mass = _M_V1
    r.primary_momentum = [E_V1, 0.0, 0.0, pz]
    r.secondary_masses = [m_chi, m_chi]
    r.secondary_momenta = [[0, 0, 0, 0], [0, 0, 0, 0]]
    r.secondary_helicities = [0, 0]
    r.interaction_vertex = [0.0, 0.0, 0.0]
    r.primary_initial_position = [0.0, 0.0, 0.0]
    return r


def _kp_box_at(x, y, z, full_width):
    placement = siren.geometry.Placement(siren.math.Vector3D(x, y, z))
    return siren.geometry.Box(placement, full_width, full_width, full_width)


def _generate_kp_alphas():
    """Run the fixed-seed KP optimization and return the converged alphas.

    A 3-channel mixture: an isotropic/physical channel plus two symmetric
    detector-directed channels aimed at boxes at +-30 deg. In this (collimated)
    regime the directed channels are pure isotropic fallback, so the canonical
    alpha_sqrt_W rule drives them to the floor and the physical channel takes
    almost all the weight -- a stable fixed point the bare W_i statistic must
    keep reproducing exactly.
    """
    th = math.radians(30.0)
    D = 20.0
    A = _kp_box_at(D * math.sin(th), 0.0, D * math.cos(th), 2.0)
    B = _kp_box_at(-D * math.sin(th), 0.0, D * math.cos(th), 2.0)

    iso = siren.injection.Isotropic2BodyChannel(0)
    dA = siren.injection.DetectorDirected2BodyChannel(A, 0)
    dB = siren.injection.DetectorDirected2BodyChannel(B, 0)
    mc = siren.injection.MultiChannelPhaseSpace()
    mc.channels = [iso, dA, dB]
    mc.weights = [0.6, 0.3, 0.1]

    optimize_multichannel_weights(
        mc,
        lambda r: iso.Density(None, r),
        siren.utilities.SIREN_random(KP_SEED),
        None,
        _kp_make_record(0.008, 0.018),
        n_iterations=KP_ITERS,
        batch_size=KP_BATCH,
        damping=KP_DAMPING,
        min_weight=KP_MIN_WEIGHT,
        update_rule=KP_RULE,
    )
    return {"kp_alphas": np.asarray(list(mc.weights), dtype=np.float64)}


# --------------------------------------------------------------------------- #
# Archive helpers                                                             #
# --------------------------------------------------------------------------- #

def _load_archive():
    if not ARCHIVE.exists():
        pytest.fail(_REGEN_MSG.format(path=ARCHIVE))
    with np.load(ARCHIVE) as npz:
        return {k: npz[k] for k in npz.files}


# --------------------------------------------------------------------------- #
# Tests                                                                       #
# --------------------------------------------------------------------------- #

def test_golden_kp_alphas():
    """Fixed-seed Kleiss-Pittau optimizer reproduces the archived converged
    per-channel alphas to rtol=1e-12 -- pinning the bare W_i = mean(w^2 g_i/g)
    statistic."""
    archive = _load_archive()
    fresh = _generate_kp_alphas()

    assert "kp_alphas" in archive, "archive missing 'kp_alphas'; regenerate it"
    np.testing.assert_allclose(
        fresh["kp_alphas"], archive["kp_alphas"], rtol=RTOL, atol=ATOL,
        err_msg=("converged KP alphas changed. If this is an intended, "
                 "changelog-labeled physics change, regenerate the archive; "
                 "otherwise the bare W_i statistic or the sampler moved."))

    alphas = fresh["kp_alphas"]
    assert alphas.shape == (3,)
    assert abs(float(alphas.sum()) - 1.0) < 1e-9
    assert np.all(alphas >= 0.0)
