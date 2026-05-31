"""Directed-2-body channel across kinematic regimes (closes coverage gap R5).

The Dutta-Kim chain runs in the *collimated* regime (a highly boosted V1
gives a sub-degree chi cone), where detector-directing is inert and the
example uses physical/isotropic channels.  The opposite, non-collimated
regime -- where the daughter's kinematic cone is wider than the target, so
the BoundInKin / full-sphere directed path actually does the work -- was
previously untested.

Collimation is governed by the parent boost vs the daughter CM speed
beta* = p_rest/E_rest; for a symmetric V1 -> chi chi decay,
beta* = sqrt(1 - 4 m_chi^2/m_V1^2).  With m_V1 = 17 MeV and m_chi = 8 MeV
(beta* = 0.34), a slow V1 (E ~ m_V1) gives the full kinematic sphere, while
E_V1 = 30 MeV gives a ~14 deg cone -- enough to place a target inside,
straddling, or outside it.

These tests assert that across all geometric regimes the directed channel
stays consistent (closure E_g[f/g] -> 1, the weighted hit fraction matches
the isotropic baseline, no NaN), and that directing concentrates exactly on
the kinematically reachable part of the target.
"""
import math

import numpy as np
import pytest

import siren

PT = siren.dataclasses.ParticleType
M_V1 = 0.017


def _make_record(m_chi, E_V1, vertex=(0.0, 0.0, 0.0)):
    pz = math.sqrt(max(E_V1 * E_V1 - M_V1 * M_V1, 0.0))
    rec = siren.dataclasses.InteractionRecord()
    rec.signature.primary_type = PT.N4
    rec.signature.target_type = PT.Decay
    rec.signature.secondary_types = [PT.NuLight, PT.Gamma]
    rec.primary_mass = M_V1
    rec.primary_momentum = [E_V1, 0.0, 0.0, pz]
    rec.secondary_masses = [m_chi, m_chi]
    rec.secondary_momenta = [[0, 0, 0, 0], [0, 0, 0, 0]]
    rec.secondary_helicities = [0, 0]
    rec.interaction_vertex = list(vertex)
    rec.primary_initial_position = list(vertex)
    return rec


def _box_at(x, y, z, full_width):
    placement = siren.geometry.Placement(siren.math.Vector3D(x, y, z))
    return siren.geometry.Box(placement, full_width, full_width, full_width)


def _hits(rec, idx, box):
    m = rec.secondary_momenta[idx]
    p = math.sqrt(m[1] ** 2 + m[2] ** 2 + m[3] ** 2)
    if p < 1e-30:
        return False
    pos = siren.math.Vector3D(*rec.interaction_vertex)
    d = siren.math.Vector3D(m[1] / p, m[2] / p, m[3] / p)
    return len(box.Intersections(pos, d)) > 0


def _run(target, m_chi, E_V1, N=15000, seed=21):
    """Isotropic baseline + 50/50 (isotropic, directed) mixture metrics."""
    iso = siren.injection.Isotropic2BodyChannel(0)
    directed = siren.injection.DetectorDirected2BodyChannel(target, 0)
    mc = siren.injection.MultiChannelPhaseSpace()
    mc.channels = [iso, directed]
    mc.weights = [0.5, 0.5]

    rng_a = siren.utilities.SIREN_random(seed)
    hit_iso = 0
    for _ in range(N):
        r = _make_record(m_chi, E_V1)
        iso.Sample(rng_a, None, r)
        if _hits(r, 0, target):
            hit_iso += 1

    rng_b = siren.utilities.SIREN_random(seed + 1)
    ws, whit, sumw, n_dir, n_dir_hit, n_bad = [], 0.0, 0.0, 0, 0, 0
    for _ in range(N):
        r = _make_record(m_chi, E_V1)
        ch = mc.Sample(rng_b, None, r)
        f = iso.Density(None, r)
        g = mc.Density(None, r)
        if not (math.isfinite(g) and math.isfinite(f)):
            n_bad += 1
            continue
        if g <= 0:
            continue
        w = f / g
        ws.append(w)
        sumw += w
        if _hits(r, 0, target):
            whit += w
        if ch == 1:
            n_dir += 1
            n_dir_hit += 1 if _hits(r, 0, target) else 0
    ws = np.array(ws)
    return dict(
        iso=hit_iso / N,
        closure=(float(ws.mean()) if len(ws) else float("nan")),
        wt_hit=(whit / sumw if sumw > 0 else 0.0),
        dir_eff=(n_dir_hit / n_dir if n_dir else 0.0),
        n_bad=n_bad,
    )


def test_noncollimated_full_sphere_closure_and_efficiency():
    """Full kinematic sphere (slow V1): directing is the obvious win and
    stays unbiased -- the BoundInKin / full-sphere path (R5)."""
    target = _box_at(0.0, 0.0, 20.0, 2.0)  # small distant on-axis target
    m = _run(target, m_chi=0.008, E_V1=0.0180)

    assert m["n_bad"] == 0
    assert 0.9 <= m["closure"] <= 1.1, f"closure {m['closure']}"
    # Volume-mode directing hits the target every time; isotropic rarely.
    assert m["dir_eff"] > 0.95
    assert m["iso"] < 0.05
    assert m["dir_eff"] > 20 * m["iso"]  # large efficiency gain
    # Unbiased: reweighted hit fraction reproduces the isotropic baseline.
    assert abs(m["wt_hit"] - m["iso"]) <= 0.5 * m["iso"] + 0.002


def test_offaxis_target_inside_cone():
    """Off-axis target fully inside a ~14 deg cone: directing steers off the
    beam axis, hits every time, and stays unbiased."""
    D, th = 20.0, math.radians(6.0)
    target = _box_at(D * math.sin(th), 0.0, D * math.cos(th), 2.0)
    m = _run(target, m_chi=0.008, E_V1=0.030)

    assert m["n_bad"] == 0
    assert 0.9 <= m["closure"] <= 1.1, f"closure {m['closure']}"
    assert m["dir_eff"] > 0.9
    assert 0.0 < m["iso"] < 0.15
    assert abs(m["wt_hit"] - m["iso"]) <= 0.5 * m["iso"] + 0.005


def test_offaxis_partial_overlap_closure():
    """Off-axis target straddling the kinematic-cone edge (the lens-overlap
    BoundInKin case): directing reaches only the kinematically allowed part
    of the target (strictly between 0 and 1), and closure still holds."""
    D, th = 20.0, math.radians(14.0)
    target = _box_at(D * math.sin(th), 0.0, D * math.cos(th), 2.0)
    m = _run(target, m_chi=0.008, E_V1=0.030)

    assert m["n_bad"] == 0
    assert 0.9 <= m["closure"] <= 1.1, f"closure {m['closure']}"
    # Partial: not inert (>iso) but not the full target (<1) -- only the
    # kinematically reachable overlap is hit.
    assert 0.1 < m["dir_eff"] < 0.95
    assert abs(m["wt_hit"] - m["iso"]) <= 0.5 * m["iso"] + 0.005


def test_disjoint_unreachable_target_graceful():
    """Off-axis target fully beyond the kinematic cone (unreachable): both
    channels give zero hits, the directed channel falls back gracefully
    (no NaN, closure -> 1), and the estimate is an unbiased zero."""
    D, th = 20.0, math.radians(25.0)
    target = _box_at(D * math.sin(th), 0.0, D * math.cos(th), 2.0)
    m = _run(target, m_chi=0.008, E_V1=0.030)

    assert m["n_bad"] == 0
    assert 0.9 <= m["closure"] <= 1.1, f"closure {m['closure']}"
    assert m["iso"] < 0.005       # target is kinematically unreachable
    assert m["dir_eff"] < 0.05    # directing cannot reach it either
    assert m["wt_hit"] < 0.01     # unbiased zero
