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


# ------------------------------------------------------------------ #
#  Nestable sub-mixture channel (prototype A)                         #
# ------------------------------------------------------------------ #

def _separated_AB(D=20.0, deg=30.0):
    t = math.radians(deg)
    A = _box_at(D * math.sin(t), 0.0, D * math.cos(t), 2.0)
    B = _box_at(-D * math.sin(t), 0.0, D * math.cos(t), 2.0)
    return A, B


def test_nested_mixture_channel_closure_and_composition():
    """A directed sub-mixture wrapped as ONE channel preserves closure, and
    the outer density composes transparently as a weighted sum."""
    A, B = _separated_AB()
    iso = siren.injection.Isotropic2BodyChannel(0)
    dirA = siren.injection.DetectorDirected2BodyChannel(A, 0)
    dirB = siren.injection.DetectorDirected2BodyChannel(B, 0)
    inner = siren.injection.MultiChannelPhaseSpace()
    inner.channels = [dirA, dirB]
    inner.weights = [0.5, 0.5]
    nested = siren.injection.NestedMixtureChannel(inner)
    outer = siren.injection.MultiChannelPhaseSpace()
    outer.channels = [iso, nested]
    outer.weights = [0.5, 0.5]

    # The nested channel presents the inner common topology/measure, so the
    # outer mixture treats it like any other channel.
    assert nested.Topology() == iso.Topology()
    assert nested.Measure() == iso.Measure()
    assert outer.CommonTopology() == iso.Topology()

    rng = siren.utilities.SIREN_random(3)
    ws = []
    for n in range(15000):
        r = _make_record(0.008, 0.018)  # full kinematic sphere
        outer.Sample(rng, None, r)
        f = iso.Density(None, r)
        g = outer.Density(None, r)
        if g <= 0 or not math.isfinite(g):
            continue
        ws.append(f / g)
        if n < 500:
            # all channels share the measure, so the density composes exactly
            g_nested = nested.Density(None, r)
            assert g_nested == pytest.approx(
                0.5 * dirA.Density(None, r) + 0.5 * dirB.Density(None, r),
                rel=1e-9, abs=1e-30)
            assert g == pytest.approx(0.5 * f + 0.5 * g_nested, rel=1e-9, abs=1e-30)
    assert 0.9 <= np.mean(ws) <= 1.1   # closure E_g[f/g] -> 1


def test_nested_inner_weights_are_optimized():
    """optimize_multichannel_weights tunes a nested channel's INNER weights.
    Starting from an imbalanced inner mixture of two symmetric reachable
    targets, the optimizer moves them back toward balance (the symmetric
    optimum) and the integral stays unbiased."""
    from siren.optimize import optimize_multichannel_weights

    A, B = _separated_AB()
    iso = siren.injection.Isotropic2BodyChannel(0)
    inner = siren.injection.MultiChannelPhaseSpace()
    inner.channels = [siren.injection.DetectorDirected2BodyChannel(A, 0),
                      siren.injection.DetectorDirected2BodyChannel(B, 0)]
    inner.weights = [0.9, 0.1]   # deliberately imbalanced
    nested = siren.injection.NestedMixtureChannel(inner)
    outer = siren.injection.MultiChannelPhaseSpace()
    outer.channels = [iso, nested]
    outer.weights = [0.5, 0.5]

    template = _make_record(0.008, 0.018)
    before = list(inner.weights)
    optimize_multichannel_weights(
        outer, lambda r: iso.Density(None, r),
        siren.utilities.SIREN_random(13), None, template,
        n_iterations=6, batch_size=2000, damping=0.5, min_weight=0.01,
        update_rule="sqrt_W")
    after = list(inner.weights)

    # The recursion ran and produced valid inner weights.
    assert abs(sum(after) - 1.0) < 1e-9
    assert all(math.isfinite(w) and w >= 0.0 for w in after)
    # Symmetric reachable targets -> balanced optimum: the under-weighted
    # channel gains, the over-weighted one drops.
    assert after[1] > before[1] + 0.05
    assert after[0] < before[0] - 0.05

    # Integral remains unbiased after optimization.
    rng = siren.utilities.SIREN_random(14)
    ws = []
    for _ in range(8000):
        r = _make_record(0.008, 0.018)
        outer.Sample(rng, None, r)
        f = iso.Density(None, r)
        g = outer.Density(None, r)
        if g > 0 and math.isfinite(g):
            ws.append(f / g)
    assert 0.9 <= np.mean(ws) <= 1.1


# ------------------------------------------------------------------ #
#  Disjoint tiling of the reachable region (Option B)                 #
# ------------------------------------------------------------------ #
#
# The reachable region is presented as a DISJOINT PARTITION of directed
# sub-channels instead of a set of overlapping ones, so the channels are a
# non-degenerate basis.  The tiles are built from existing pieces
# (Box sub-cells and BooleanGeometry union/subtraction) and wrapped in the
# existing NestedMixtureChannel, so C1 holds automatically.  See
# docs/volume_aware_directed_channel_plan.md.

from siren import directed_tiling as dt  # noqa: E402

G = siren.geometry


def _run_channel(directed, hit_geo, m_chi, E_V1, N=15000, seed=21):
    """Like _run, but for a pre-built directed channel (e.g. a tiling)."""
    iso = siren.injection.Isotropic2BodyChannel(0)
    mc = siren.injection.MultiChannelPhaseSpace()
    mc.channels = [iso, directed]
    mc.weights = [0.5, 0.5]

    rng_a = siren.utilities.SIREN_random(seed)
    hit_iso = 0
    for _ in range(N):
        r = _make_record(m_chi, E_V1)
        iso.Sample(rng_a, None, r)
        if _hits(r, 0, hit_geo):
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
        if _hits(r, 0, hit_geo):
            whit += w
        if ch == 1:
            n_dir += 1
            n_dir_hit += 1 if _hits(r, 0, hit_geo) else 0
    ws = np.array(ws)
    return dict(
        iso=hit_iso / N,
        closure=(float(ws.mean()) if len(ws) else float("nan")),
        wt_hit=(whit / sumw if sumw > 0 else 0.0),
        dir_eff=(n_dir_hit / n_dir if n_dir else 0.0),
        n_bad=n_bad,
    )


def test_geometry_volume_helpers():
    """Exact volumes for primitives; Monte-Carlo for composites."""
    box = _box_at(0, 0, 20, 4.0)
    assert dt.geometry_volume(box) == pytest.approx(64.0, rel=1e-9)
    sph = G.Sphere(G.Placement(siren.math.Vector3D(0, 0, 0)), 5.0, 0.0)
    assert dt.geometry_volume(sph) == pytest.approx(4 / 3 * math.pi * 125, rel=1e-9)
    # Two 4-cubes offset by 2 in x: |union| = 64 + 64 - 32 = 96.
    A = _box_at(-1, 0, 20, 4.0)
    B = _box_at(1, 0, 20, 4.0)
    u = G.BooleanGeometry(G.BooleanOperation.UNION, A, B)
    assert dt.geometry_volume(u, n_mc=200000) == pytest.approx(96.0, rel=0.02)


def test_dedup_contained_collapses_nested_spheres():
    spheres = [G.Sphere(G.Placement(siren.math.Vector3D(0, 0, 0)), r, 0.0)
               for r in (2, 5, 10, 20)]
    kept = dt.dedup_contained(spheres)
    assert len(kept) == 1
    assert kept[0] is spheres[-1]   # the outermost survives


def test_cluster_by_cone_overlap_separates_and_groups():
    v = [siren.math.Vector3D(0, 0, 0)]
    sepA, sepB = _separated_AB(D=20.0, deg=30.0)   # far apart from the vertex
    assert len(dt.cluster_by_cone_overlap([sepA, sepB], vertices=v)) == 2
    nearA = _box_at(-1, 0, 20, 2.0)
    nearB = _box_at(1, 0, 20, 2.0)                  # overlapping cones
    assert len(dt.cluster_by_cone_overlap([nearA, nearB], vertices=v)) == 1


def test_grid_tiling_closure_full_sphere():
    """A grid of disjoint sub-box cells tiles a box; closure E_g[f/g] -> 1."""
    box = _box_at(0.0, 0.0, 20.0, 4.0)
    tiled = dt.build_grid_tiling(box, 2)             # 8 cells
    m = _run_channel(tiled, box, m_chi=0.008, E_V1=0.018)
    assert m["n_bad"] == 0
    assert 0.9 <= m["closure"] <= 1.1, f"closure {m['closure']}"
    assert m["dir_eff"] > 0.95                       # Volume mode hits every time
    assert abs(m["wt_hit"] - m["iso"]) <= 0.5 * m["iso"] + 0.003


def test_subtraction_tiling_closure_and_disjoint():
    """BooleanGeometry subtraction shells T_k = V_k \\ union(earlier) are
    disjoint and tile the union; closure holds."""
    A = _box_at(-1.0, 0.0, 20.0, 4.0)
    B = _box_at(1.0, 0.0, 20.0, 4.0)
    union = G.BooleanGeometry(G.BooleanOperation.UNION, A, B)
    tiles = dt.disjointify([A, B])                   # [A, B\A]
    # The shells are disjoint: a point in A is not in B\A.
    shell = tiles[1]
    rng = siren.utilities.SIREN_random(1)
    for _ in range(500):
        p = siren.math.Vector3D(
            -2.0 + 0.01 * rng.Uniform(0, 1), 0.0, 20.0)  # inside A
        if A.IsInside(p):
            assert not shell.IsInside(p)
    tiled = dt.build_subtraction_tiling([A, B])
    m = _run_channel(tiled, union, m_chi=0.008, E_V1=0.018)
    assert m["n_bad"] == 0
    assert 0.9 <= m["closure"] <= 1.1, f"closure {m['closure']}"


def test_union_tile_closure():
    """A single union tile (dedup + BooleanGeometry UNION) closes."""
    A = _box_at(-1.0, 0.0, 20.0, 4.0)
    B = _box_at(1.0, 0.0, 20.0, 4.0)
    union = G.BooleanGeometry(G.BooleanOperation.UNION, A, B)
    ch = dt.build_union_tile([A, B])
    m = _run_channel(ch, union, m_chi=0.008, E_V1=0.018)
    assert m["n_bad"] == 0
    assert 0.9 <= m["closure"] <= 1.1, f"closure {m['closure']}"


def test_angular_sector_full_cone_closure():
    """A single full-cone angular sector closes (the C++ FullConeClosure
    test, exercised through the Python binding)."""
    box = _box_at(0.0, 0.0, 20.0, 4.0)
    sector = siren.injection.DetectorDirectedAngularSectorChannel(
        box, 0.0, 1.0, 0.0, 2.0 * math.pi, 0)
    m = _run_channel(sector, box, m_chi=0.008, E_V1=0.018)
    assert m["n_bad"] == 0
    assert 0.9 <= m["closure"] <= 1.1, f"closure {m['closure']}"


def test_angular_sector_density_zero_outside_own_bin():
    """Sample == Density mechanism: a direction drawn by one sector lands in
    that sector's bin (positive own-density) and outside a disjoint sector's
    bin (zero cross-density)."""
    box = _box_at(0.0, 0.0, 20.0, 4.0)
    Sector = siren.injection.DetectorDirectedAngularSectorChannel
    s0 = Sector(box, 0.0, 1.0, 0.0, math.pi, 0)            # phi in [0, pi)
    s1 = Sector(box, 0.0, 1.0, math.pi, 2.0 * math.pi, 0)  # phi in [pi, 2pi)
    rng = siren.utilities.SIREN_random(5)
    n_self_pos, n_cross_zero, n = 0, 0, 1500
    for _ in range(n):
        r = _make_record(0.008, 0.018)
        s0.Sample(rng, None, r)
        if s0.Density(None, r) > 0.0:
            n_self_pos += 1
        if s1.Density(None, r) == 0.0:
            n_cross_zero += 1
    assert n_self_pos > 0.99 * n      # own sector claims its samples
    assert n_cross_zero > 0.99 * n    # the disjoint sector does not


def test_angular_tiling_closure():
    """build_angular_tiling tiles the bounding cone into sectors; closure holds."""
    box = _box_at(0.0, 0.0, 20.0, 4.0)
    tiled = dt.build_angular_tiling(box, n_theta=2, n_phi=4)
    m = _run_channel(tiled, box, m_chi=0.008, E_V1=0.018)
    assert m["n_bad"] == 0
    assert 0.9 <= m["closure"] <= 1.1, f"closure {m['closure']}"


def test_tiling_covers_separated_regions_single_misses():
    """Coverage mechanism: a tiled [A, B] channel directs at BOTH separated
    reachable regions, while a single channel aimed at A directs at A only --
    the multi-region coverage that single-fiducial gives up."""
    D = 20.0
    th = math.radians(6.0)   # both within the kinematic cone at E_V1=0.030
    A = _box_at(D * math.sin(th), 0.0, D * math.cos(th), 2.0)
    B = _box_at(-D * math.sin(th), 0.0, D * math.cos(th), 2.0)

    def directed_hits(channel, seed):
        rng = siren.utilities.SIREN_random(seed)
        hitA = hitB = 0
        for _ in range(6000):
            r = _make_record(0.008, 0.030)
            channel.Sample(rng, None, r)
            hitA += 1 if _hits(r, 0, A) else 0
            hitB += 1 if _hits(r, 0, B) else 0
        return hitA / 6000, hitB / 6000

    single_a = siren.injection.DetectorDirected2BodyChannel(A, 0)
    sa_A, sa_B = directed_hits(single_a, 31)
    assert sa_A > 0.8           # single channel directs at A
    assert sa_B < 0.02          # ... and never at B

    tiled = dt.build_subtraction_tiling([A, B])   # disjoint [A, B]
    ti_A, ti_B = directed_hits(tiled, 33)
    assert ti_A > 0.2 and ti_B > 0.2              # tiling directs at BOTH


def test_separated_region_tiling_beats_single_and_flat_ess():
    """Headline ESS A/B (2-body).  In the non-collimated, multi-region regime
    (full kinematic sphere, two widely separated reachable regions A,B at
    +-30 deg) the volume-aware Option B channels -- a disjoint tiling [A,B] and
    a union(A,B) -- beat BOTH a single fiducial box spanning the gap and a flat
    nested set, after Kleiss-Pittau optimization.  All configs stay unbiased
    (closure -> 1)."""
    from siren.optimize import optimize_multichannel_weights

    D = 20.0
    th = math.radians(30.0)
    x = D * math.sin(th)
    z = D * math.cos(th)
    A = _box_at(x, 0.0, z, 2.0)
    B = _box_at(-x, 0.0, z, 2.0)
    # one big box spanning both A and B (the gap is ~80% of its width)
    C = G.Box(G.Placement(siren.math.Vector3D(0.0, 0.0, z)), 2 * x + 2.0, 2.0, 2.0)
    C2 = G.Box(G.Placement(siren.math.Vector3D(0.0, 0.0, z)), 2 * x + 5.0, 5.0, 5.0)
    f2 = dt.default_2body_factory(0)

    def signal_hit(r):
        return _hits(r, 0, A) or _hits(r, 0, B)

    def ess_of(directed):
        iso = siren.injection.Isotropic2BodyChannel(0)
        outer = siren.injection.MultiChannelPhaseSpace()
        outer.channels = [iso, directed]
        outer.weights = [0.5, 0.5]
        optimize_multichannel_weights(
            outer, lambda r: iso.Density(None, r),
            siren.utilities.SIREN_random(5), None, _make_record(0.008, 0.018),
            n_iterations=5, batch_size=1500, damping=0.6, min_weight=0.01,
            update_rule="alpha_sqrt_W")
        rng = siren.utilities.SIREN_random(77)
        w_all, w_hit = [], []
        for _ in range(10000):
            r = _make_record(0.008, 0.018)
            outer.Sample(rng, None, r)
            f = iso.Density(None, r)
            g = outer.Density(None, r)
            if g <= 0 or not (math.isfinite(g) and math.isfinite(f)):
                continue
            w = f / g
            w_all.append(w)
            w_hit.append(w if signal_hit(r) else 0.0)
        w_all = np.array(w_all)
        w_hit = np.array(w_hit)
        closure = float(w_all.mean())
        denom = float((w_hit ** 2).sum())
        ess = (float(w_hit.sum()) ** 2 / denom / len(w_all)) if denom > 0 else 0.0
        return closure, ess

    c_single, e_single = ess_of(siren.injection.DetectorDirected2BodyChannel(C, 0))
    c_flat, e_flat = ess_of(dt._wrap([A, B, C, C2], f2, True, 100000))
    c_union, e_union = ess_of(dt.build_union_tile([A, B]))
    c_tiled, e_tiled = ess_of(dt.build_subtraction_tiling([A, B]))

    # Every config is unbiased.
    for c in (c_single, c_flat, c_union, c_tiled):
        assert 0.9 <= c <= 1.1, f"closure {c}"

    # Option B (union and disjoint tiling) beats the single fiducial (gap waste)
    # and the flat nested set (budget dilution).
    assert e_union > 1.25 * e_single
    assert e_tiled > 1.15 * e_single
    assert e_union > 1.10 * e_flat
    assert e_tiled > 0.98 * e_flat   # tiling at least matches the flat set


def test_group_directed_channels_preserves_density():
    """group_directed_channels wraps the directed channels into one
    NestedMixtureChannel WITHOUT changing the mixture density g(x) (the group's
    outer weight is the sum of the directed weights, inner weights renormalized)."""
    from siren.optimize import group_directed_channels

    box = _box_at(0.0, 0.0, 20.0, 2.0)
    iso = siren.injection.Isotropic2BodyChannel(0)
    flat = siren.injection.MultiChannelPhaseSpace()
    flat.channels = [iso] + [siren.injection.DetectorDirected2BodyChannel(box, 0)
                             for _ in range(4)]
    flat.weights = [0.2, 0.2, 0.2, 0.2, 0.2]

    grouped = group_directed_channels(flat)
    assert len(grouped.channels) == 2                       # physical + one group
    assert grouped.weights[1] == pytest.approx(0.8)         # sum of directed weights

    rng = siren.utilities.SIREN_random(2)
    for _ in range(500):
        r = _make_record(0.008, 0.018)
        flat.Sample(rng, None, r)
        gf = flat.Density(None, r)
        if gf > 0 and math.isfinite(gf):
            assert grouped.Density(None, r) == pytest.approx(gf, rel=1e-9, abs=1e-30)


def test_grouping_drives_fallback_down_further_than_flat():
    """With N pure-fallback directed channels, grouping floors the directed
    weight at one min_weight instead of N*min_weight, driving the redundant
    fallback substantially closer to zero (the many-channels convergence fix)."""
    from siren.optimize import (optimize_multichannel_weights,
                                group_directed_channels)

    box = _box_at(0.0, 0.0, 20.0, 2.0)
    N, mw, E_V1 = 8, 0.01, 0.5      # collimated -> all N directed are pure fallback

    def build_flat():
        iso = siren.injection.Isotropic2BodyChannel(0)
        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso] + [siren.injection.DetectorDirected2BodyChannel(box, 0)
                               for _ in range(N)]
        mc.weights = [1.0 / (N + 1)] * (N + 1)
        return mc, iso

    flat, iso_f = build_flat()
    grouped = group_directed_channels(build_flat()[0])
    iso_g = grouped.channels[0]

    for mc, isoc in ((flat, iso_f), (grouped, iso_g)):
        optimize_multichannel_weights(
            mc, lambda r, c=isoc: c.Density(None, r),
            siren.utilities.SIREN_random(1), None, _make_record(0.008, E_V1),
            n_iterations=8, batch_size=1500, damping=0.6, min_weight=mw,
            update_rule="alpha_sqrt_W", discount_fallback=True)

    flat_directed = sum(flat.weights[1:])     # N independent directed channels
    grouped_directed = grouped.weights[1]     # the single group
    assert grouped_directed < flat_directed   # grouping drives it lower
    assert grouped_directed < 3 * mw          # group lands near one min_weight
    assert flat_directed > 3 * mw             # flat stuck near N*min_weight


def test_discount_fallback_turns_off_pure_fallback_director():
    """discount_fallback (optimizer): a directed channel that is always in its
    isotropic 1/4pi fallback (collimated regime) contributes only the shared
    floor.  Crediting its variance only on directing-active events drives it to
    min_weight in favor of the physical/isotropic channel; without discounting
    it stays degenerate with the isotropic channel.  Both stay unbiased."""
    from siren.optimize import optimize_multichannel_weights

    box = _box_at(0.0, 0.0, 20.0, 2.0)
    # E_V1 = 0.5: highly boosted V1 -> sub-degree chi cone inside the target cone
    # -> the directed channel is always in its isotropic fallback.
    assert not siren.injection.DetectorDirected2BodyChannel(box, 0).DirectingActive(
        _make_record(0.008, 0.5))

    def optimized_directed_weight(discount):
        iso = siren.injection.Isotropic2BodyChannel(0)
        d = siren.injection.DetectorDirected2BodyChannel(box, 0)
        mc = siren.injection.MultiChannelPhaseSpace()
        mc.channels = [iso, d]
        mc.weights = [0.5, 0.5]
        optimize_multichannel_weights(
            mc, lambda r: iso.Density(None, r),
            siren.utilities.SIREN_random(7), None, _make_record(0.008, 0.5),
            n_iterations=6, batch_size=1500, damping=0.6, min_weight=0.01,
            update_rule="alpha_sqrt_W", discount_fallback=discount)
        # closure: the integral stays unbiased either way
        rng = siren.utilities.SIREN_random(8)
        ws = []
        for _ in range(4000):
            r = _make_record(0.008, 0.5)
            mc.Sample(rng, None, r)
            f = iso.Density(None, r)
            g = mc.Density(None, r)
            if g > 0 and math.isfinite(g):
                ws.append(f / g)
        return mc.weights[1], float(np.mean(ws))

    w_on, closure_on = optimized_directed_weight(True)
    w_off, closure_off = optimized_directed_weight(False)

    assert w_on < 0.1      # discounting drives the pure-fallback director down
    assert w_off > 0.3     # without it, it stays degenerate with the iso channel
    assert 0.9 <= closure_on <= 1.1
    assert 0.9 <= closure_off <= 1.1


def test_builder_tiles_optimize_and_stay_unbiased():
    """A builder-produced tiling integrates with the recursive Kleiss-Pittau
    optimizer: an imbalanced inner weight set is retuned and the integral
    stays unbiased."""
    from siren.optimize import optimize_multichannel_weights

    A, B = _separated_AB(D=20.0, deg=30.0)        # both reachable (full sphere)
    tiled = dt.build_subtraction_tiling([A, B])
    inner = tiled.mixture
    inner.weights = [0.9, 0.1]                     # imbalanced
    iso = siren.injection.Isotropic2BodyChannel(0)
    outer = siren.injection.MultiChannelPhaseSpace()
    outer.channels = [iso, tiled]
    outer.weights = [0.5, 0.5]

    before = list(inner.weights)
    optimize_multichannel_weights(
        outer, lambda r: iso.Density(None, r),
        siren.utilities.SIREN_random(13), None, _make_record(0.008, 0.018),
        n_iterations=6, batch_size=2000, damping=0.5, min_weight=0.01,
        update_rule="sqrt_W")
    after = list(inner.weights)

    assert abs(sum(after) - 1.0) < 1e-9
    assert all(math.isfinite(w) and w >= 0.0 for w in after)
    assert after[1] > before[1] + 0.05            # under-weighted tile gains
    assert after[0] < before[0] - 0.05

    rng = siren.utilities.SIREN_random(14)
    ws = []
    for _ in range(8000):
        r = _make_record(0.008, 0.018)
        outer.Sample(rng, None, r)
        f = iso.Density(None, r)
        g = outer.Density(None, r)
        if g > 0 and math.isfinite(g):
            ws.append(f / g)
    assert 0.9 <= np.mean(ws) <= 1.1
