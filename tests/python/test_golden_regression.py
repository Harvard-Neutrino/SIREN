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
1. test_golden_weighted_distributions -- a fixed-seed multi-vertex chain built
   from a real siren injection.Injector + injection.Weighter (the same
   data-free DummyCrossSection + CCM-detector assembly used by the HepMC3
   closure tests). N=500 accepted events. Archived: the per-event weights, the
   primary energies, the interaction-vertex coordinates, the per-event tree
   depth (a cheap kinematic-topology observable available from the chain), and
   weighted histograms of energy and of vertex radius.
2. test_golden_kp_alphas -- a fixed-seed Kleiss-Pittau optimization of a toy
   3-channel mixture (physical/isotropic + two detector-directed channels) via
   optimize.optimize_multichannel_weights, which drives the bare
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
from siren.optimize import optimize_multichannel_weights


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
# Part 1: a real assembled injector + weighter (multi-vertex, data-free)       #
# --------------------------------------------------------------------------- #

def _skip_unless_ccm_data():
    """Skip only when the CCM detector data files are absent.

    Missing files are an environment condition; any error raised while
    PARSING present files is a real regression and must propagate.
    """
    try:
        det_dir = _util.get_detector_model_path("CCM")
    except ValueError as e:
        pytest.skip(f"CCM detector model path unavailable: {e}")
    missing = [p for p in (os.path.join(det_dir, "materials.dat"),
                           os.path.join(det_dir, "densities.dat"))
               if not os.path.exists(p)]
    if missing:
        pytest.skip("CCM detector data missing: " + ", ".join(missing))


def _load_ccm_detector():
    """Load the CCM detector model; errors propagate to the caller."""
    dm = detector.DetectorModel()
    det_dir = _util.get_detector_model_path("CCM")
    dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
    dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    return dm


def _build_chain(detector_model, n_inject, seed):
    """Assemble a real multi-vertex Injector + Weighter over a data-free
    DummyCrossSection.

    Generation flux: a power-law energy spectrum, a neutrino-helicity
    distribution, isotropic direction, and a point-source position -- so the
    weight genuinely depends on energy, helicity, and geometry (a monoenergetic
    setup would not exercise the reweighting). One generation of secondaries is
    simulated so every event is a multi-vertex cascade. Returns
    (injector, weighter, keepalive) where keepalive holds Python objects that
    must outlive the C++ engines.
    """
    xs = interactions.DummyCrossSection()
    int_col = interactions.InteractionCollection(_NuMu, [xs])

    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = _NuMu
    primary_inj.interactions = int_col
    primary_inj.distributions = [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]

    primary_phys = injection.PhysicalProcess()
    primary_phys.primary_type = _NuMu
    primary_phys.interactions = int_col
    primary_phys.distributions = [
        distributions.PrimaryMass(0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
    ]

    sec_xs = interactions.DummyCrossSection()
    sec_col = interactions.InteractionCollection(_NuMu, [sec_xs])

    sec_inj = injection.SecondaryInjectionProcess()
    sec_inj.primary_type = _NuMu
    sec_inj.interactions = sec_col
    sec_inj.distributions = [distributions.SecondaryPhysicalVertexDistribution()]

    sec_phys = injection.PhysicalProcess()
    sec_phys.primary_type = _NuMu
    sec_phys.interactions = sec_col
    sec_phys.distributions = [distributions.SecondaryPhysicalVertexDistribution()]

    rand = utilities.SIREN_random(seed)
    inj = injection._Injector(n_inject, detector_model, primary_inj, [sec_inj], rand)
    # Simulate exactly one generation of secondaries (stop at depth >= 1).
    inj.SetStoppingCondition(lambda tree, datum, i: datum.depth(tree) >= 1)
    weighter = injection._Weighter([inj], detector_model, primary_phys, [sec_phys])

    keepalive = (xs, int_col, sec_xs, sec_col, primary_inj, primary_phys,
                 sec_inj, sec_phys, rand)
    return inj, weighter, keepalive


def _generate_chain_arrays(detector_model):
    """Generate N_EVENTS accepted events and return the pinned per-event arrays
    plus the derived weighted histograms.

    Deterministic: identical seeds -> identical arrays (verified by the second
    pytest run in CI and by the determinism assertions here).
    """
    inj, weighter, _keepalive = _build_chain(detector_model, MAX_ATTEMPTS, CHAIN_SEED)

    events = []
    energy = []
    vx, vy, vz = [], [], []
    depth = []

    for _ in range(MAX_ATTEMPTS):
        if len(events) >= N_EVENTS:
            break
        try:
            ev = inj.GenerateEvent()
        except RuntimeError as err:
            # Only the max-attempts sentinel ends the loop; any other typed
            # engine RuntimeError re-raises so it cannot be masked here.
            if "maximum number of injection attempts" not in str(err):
                raise
            break
        if len(ev.tree) == 0:
            continue
        root = ev.tree[0].record
        p = root.primary_momentum
        v = root.interaction_vertex
        events.append(ev)
        energy.append(p[0])
        vx.append(v[0])
        vy.append(v[1])
        vz.append(v[2])
        depth.append(len(ev.tree))

    if len(events) < N_EVENTS:
        raise RuntimeError(
            f"golden chain produced only {len(events)} accepted events "
            f"(< N_EVENTS={N_EVENTS}) within {MAX_ATTEMPTS} attempts")

    # Weight after generation completes: EventWeight normalizes by the realized
    # injected count, so all events must share the same (final) normalization.
    # Weighting mid-generation would bake the running count into each weight.
    weights = [weighter.EventWeight(ev) for ev in events]

    weights = np.asarray(weights[:N_EVENTS], dtype=np.float64)
    energy = np.asarray(energy[:N_EVENTS], dtype=np.float64)
    vx = np.asarray(vx[:N_EVENTS], dtype=np.float64)
    vy = np.asarray(vy[:N_EVENTS], dtype=np.float64)
    vz = np.asarray(vz[:N_EVENTS], dtype=np.float64)
    depth = np.asarray(depth[:N_EVENTS], dtype=np.float64)

    radius = np.sqrt(vx * vx + vy * vy + vz * vz)
    energy_hist, _ = np.histogram(energy, bins=E_BINS, weights=weights)
    radius_hist, _ = np.histogram(radius, bins=R_BINS, weights=weights)

    return {
        "chain_weights": weights,
        "chain_energy": energy,
        "chain_vx": vx,
        "chain_vy": vy,
        "chain_vz": vz,
        "chain_depth": depth,
        "chain_energy_hist": energy_hist.astype(np.float64),
        "chain_radius_hist": radius_hist.astype(np.float64),
    }


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

def _build_archive():
    """Compute every pinned array. Used both by --regenerate and by the tests."""
    detector_model = _load_ccm_detector()
    data = {}
    data.update(_generate_chain_arrays(detector_model))
    data.update(_generate_kp_alphas())
    return data


def _load_archive():
    if not ARCHIVE.exists():
        pytest.fail(_REGEN_MSG.format(path=ARCHIVE))
    with np.load(ARCHIVE) as npz:
        return {k: npz[k] for k in npz.files}


# --------------------------------------------------------------------------- #
# Tests                                                                       #
# --------------------------------------------------------------------------- #

def test_golden_weighted_distributions():
    """Fixed-seed assembled chain reproduces the archived weighted distributions
    to rtol=1e-12."""
    _skip_unless_ccm_data()
    detector_model = _load_ccm_detector()

    archive = _load_archive()
    fresh = _generate_chain_arrays(detector_model)

    for key in ("chain_weights", "chain_energy", "chain_vx", "chain_vy",
                "chain_vz", "chain_depth", "chain_energy_hist",
                "chain_radius_hist"):
        assert key in archive, f"archive missing key {key!r}; regenerate it"
        np.testing.assert_allclose(
            fresh[key], archive[key], rtol=RTOL, atol=ATOL,
            err_msg=(f"golden chain array {key!r} changed. If this is an "
                     "intended, changelog-labeled physics change, regenerate "
                     "the archive; otherwise the engine numerics moved."))

    # Sanity: the pinned sample really is the assembled multi-vertex pipeline.
    assert fresh["chain_weights"].shape == (N_EVENTS,)
    assert np.all(np.isfinite(fresh["chain_weights"]))
    assert np.all(fresh["chain_weights"] > 0.0)
    assert np.all(fresh["chain_depth"] >= 2.0)   # secondaries -> multi-vertex


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


# --------------------------------------------------------------------------- #
# Regeneration entry point (never invoked by pytest)                          #
# --------------------------------------------------------------------------- #

def _regenerate():
    data = _build_archive()
    ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
    np.savez(ARCHIVE, **data)
    print(f"wrote {ARCHIVE}")
    for key in sorted(data):
        arr = np.asarray(data[key])
        print(f"  {key}: shape={arr.shape} dtype={arr.dtype}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Regenerate the golden-physics regression archive. "
                    "Only run this as part of a deliberate, changelog-labeled "
                    "physics change.")
    parser.add_argument("--regenerate", action="store_true",
                        help="recompute and overwrite the archive .npz")
    args = parser.parse_args()

    if args.regenerate:
        _regenerate()
    else:
        parser.error("pass --regenerate to write the archive")
