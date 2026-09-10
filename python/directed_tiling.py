"""
Disjoint tiling of the reachable region for detector-directed phase-space
biasing.

The detector-directed channels point a daughter at a geometric target volume so
the downstream signal lands in the detector.  When several target volumes
overlap/nest, presenting them as separate channels makes the Kleiss-Pittau
optimizer degenerate (it cannot tell near-identical proposals apart) and dilutes
the sampling budget -- even though the mixture denominator already counts overlap
exactly once (so the estimate stays unbiased).

This module presents the reachable region as a DISJOINT PARTITION (tiling) of
directed sub-channels instead.  Each tile runs the SAME existing
kinematic-cone-intersect-geometry regime machinery at finer resolution and gets
its own optimizer weight; the recursive Kleiss-Pittau optimizer then allocates
sampling across tiles, concentrating on the tiles inside the per-event kinematic
cone and suppressing the rest.  Disjoint tiles are a non-degenerate channel basis
(the MadGraph/Kleiss-Pittau "one distinct channel per region" ideal).

The tiles are built from existing pieces -- `siren.geometry.Box` sub-cells and
`siren.geometry.BooleanGeometry` (UNION to dedup/collapse nested inputs,
SUBTRACTION to disjointify overlapping inputs) -- and wrapped in the existing
`NestedMixtureChannel`, so C1 (Sample == Density) holds automatically and the
existing `optimize_multichannel_weights` recursion tunes the tile weights.
"""

from __future__ import annotations

import inspect
import math
from functools import reduce
from typing import Callable, List, Optional, Sequence, Tuple

__all__ = [
    "geometry_volume",
    "bounding_sphere",
    "dedup_contained",
    "cluster_by_cone_overlap",
    "make_union",
    "default_2body_factory",
    "build_grid_tiling",
    "build_subtraction_tiling",
    "build_union_tile",
    "build_angular_tiling",
    "build_directed_tiling",
]


# ------------------------------------------------------------------ #
#  Geometry helpers                                                   #
# ------------------------------------------------------------------ #

def _siren():
    import siren
    return siren


def _center(geo) -> Tuple[float, float, float]:
    """Center of a geometry's placement (world coords)."""
    pos = geo.placement.Position
    pos = pos() if callable(pos) else pos
    return pos.GetX(), pos.GetY(), pos.GetZ()


def _aabb(geo):
    box = geo.GetWorldBoundingBox()
    mn, mx = box.min_corner, box.max_corner
    return (mn.GetX(), mn.GetY(), mn.GetZ()), (mx.GetX(), mx.GetY(), mx.GetZ())


def bounding_sphere(geo) -> Tuple[Tuple[float, float, float], float]:
    """(center, radius) of the world-AABB bounding sphere of a geometry."""
    (mnx, mny, mnz), (mxx, mxy, mxz) = _aabb(geo)
    cx, cy, cz = 0.5 * (mnx + mxx), 0.5 * (mny + mxy), 0.5 * (mnz + mxz)
    r = 0.5 * math.sqrt((mxx - mnx) ** 2 + (mxy - mny) ** 2 + (mxz - mnz) ** 2)
    return (cx, cy, cz), r


def geometry_volume(geo, n_mc: int = 100000, seed: int = 20240531) -> float:
    """Return the C++ directed-channel normalization volume.

    The implementation is exported by ``siren.injection`` so Python tiling and
    C++ channel construction cannot drift. Composite and unsupported
    geometries deliberately raise instead of using a Monte-Carlo estimate.
    ``n_mc`` and ``seed`` remain accepted for source compatibility.
    """
    siren = _siren()
    del n_mc, seed
    try:
        return float(siren.injection.geometry_volume(geo))
    except (TypeError, ValueError, RuntimeError) as error:
        raise ValueError(str(error)) from error


def _sample_inside(geo, n: int, seed: int) -> list:
    """Up to `n` points uniformly inside `geo` (AABB rejection)."""
    siren = _siren()
    V3 = siren.math.Vector3D
    (mnx, mny, mnz), (mxx, mxy, mxz) = _aabb(geo)
    dx, dy, dz = mxx - mnx, mxy - mny, mxz - mnz
    rng = siren.utilities.SIREN_random(seed)
    pts = []
    for _ in range(50 * n):
        p = V3(mnx + rng.Uniform(0, 1) * dx,
               mny + rng.Uniform(0, 1) * dy,
               mnz + rng.Uniform(0, 1) * dz)
        if geo.IsInside(p):
            pts.append(p)
            if len(pts) >= n:
                break
    return pts


# ------------------------------------------------------------------ #
#  Dedup + clustering                                                 #
# ------------------------------------------------------------------ #

def dedup_contained(volumes: Sequence, n_test: int = 200, seed: int = 7) -> List:
    """Drop volumes that are contained in another (conservative).

    `volumes[i]` is dropped if a sample of `n_test` points uniformly inside it
    are ALL inside some other `volumes[j]`.  Requiring every sampled point to be
    contained makes the test conservative -- a non-contained volume is very
    unlikely to be dropped.  Used to collapse nested input probes
    (e.g. concentric spheres) before tiling.
    """
    vols = list(volumes)
    n = len(vols)
    drop = [False] * n
    pts = [_sample_inside(vols[i], n_test, seed + i) for i in range(n)]
    for i in range(n):
        if drop[i] or not pts[i]:
            continue
        for j in range(n):
            if i == j or drop[j]:
                continue
            if all(vols[j].IsInside(p) for p in pts[i]):
                drop[i] = True
                break
    return [vols[i] for i in range(n) if not drop[i]]


def cluster_by_cone_overlap(
    volumes: Sequence,
    vertices: Optional[Sequence] = None,
) -> List[List]:
    """Group volumes whose bounding cones overlap, via union-find.

    From a representative vertex `v`, volume `V` subtends a cone with axis
    `center(V) - v` and half-angle `asin(R_bound / dist)`.  Two volumes are
    unioned into a cluster if their cones overlap (axis separation < sum of
    half-angles) from the representative vertex.  Keeping a union within one
    cluster keeps its bounding cone tight, so directing stays active (the
    over-inflated-cone failure is avoided).

    `vertices`: representative production vertices (Vector3D or (x,y,z)); their
    mean is used.  If None, falls back to bounding-sphere overlap (geometric
    proximity), which is vertex-independent.
    """
    vols = list(volumes)
    n = len(vols)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    spheres = [bounding_sphere(v) for v in vols]

    if vertices is None:
        # bounding-sphere overlap
        for i in range(n):
            (ci, ri) = spheres[i]
            for j in range(i + 1, n):
                (cj, rj) = spheres[j]
                d = math.dist(ci, cj)
                if d < ri + rj:
                    union(i, j)
    else:
        verts = []
        for v in vertices:
            if hasattr(v, "GetX"):
                verts.append((v.GetX(), v.GetY(), v.GetZ()))
            else:
                verts.append(tuple(v))
        vx = sum(p[0] for p in verts) / len(verts)
        vy = sum(p[1] for p in verts) / len(verts)
        vz = sum(p[2] for p in verts) / len(verts)
        vmean = (vx, vy, vz)

        def cone(i):
            (c, r) = spheres[i]
            ax = (c[0] - vmean[0], c[1] - vmean[1], c[2] - vmean[2])
            dist = math.sqrt(ax[0] ** 2 + ax[1] ** 2 + ax[2] ** 2)
            if dist < 1e-12:
                return (0.0, 0.0, 1.0), math.pi
            ax = (ax[0] / dist, ax[1] / dist, ax[2] / dist)
            half = math.asin(min(r / dist, 1.0))
            return ax, half

        cones = [cone(i) for i in range(n)]
        for i in range(n):
            ax_i, h_i = cones[i]
            for j in range(i + 1, n):
                ax_j, h_j = cones[j]
                dot = max(-1.0, min(1.0, ax_i[0] * ax_j[0] + ax_i[1] * ax_j[1] + ax_i[2] * ax_j[2]))
                sep = math.acos(dot)
                if sep < h_i + h_j:
                    union(i, j)

    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(vols[i])
    return list(groups.values())


# ------------------------------------------------------------------ #
#  Tile construction                                                  #
# ------------------------------------------------------------------ #

def make_union(volumes: Sequence):
    """A single geometry = the union of `volumes` (BooleanGeometry tree)."""
    siren = _siren()
    G = siren.geometry
    vols = list(volumes)
    if not vols:
        raise ValueError("make_union: empty volume list")
    if len(vols) == 1:
        return vols[0]
    return reduce(
        lambda a, b: G.BooleanGeometry(G.BooleanOperation.UNION, a, b), vols)


def default_2body_factory(daughter_index: int = 0, mode=None) -> Callable:
    """A ``factory(geometry, exact_volume=None)`` for two-body channels."""
    siren = _siren()
    inj = siren.injection
    if mode is None:
        mode = inj.DirectedMode.Volume

    def factory(geo, exact_volume=None):
        volume = -1.0 if exact_volume is None else exact_volume
        return inj.DetectorDirected2BodyChannel(
            geo, daughter_index, mode, volume)

    return factory


def _resolved_volume(geo, volume_fn, required: bool):
    try:
        volume = geometry_volume(geo)
    except ValueError:
        if volume_fn is None:
            if required:
                raise ValueError(
                    "Volume-mode composite tiles require an exact volume_fn")
            return None
        volume = volume_fn(geo)
    volume = float(volume)
    if not math.isfinite(volume) or volume < 0.0:
        raise ValueError("exact tile volume must be finite and non-negative")
    return volume


def _factory_accepts_exact_volume(factory) -> bool:
    """Whether ``factory`` can be called as ``factory(geometry, exact_volume)``.

    The arity is probed through ``inspect.signature`` so that a ``TypeError``
    raised inside a two-argument factory propagates to the caller unmasked.
    Callables whose signature cannot be introspected (extension types) are
    assumed to accept the two-argument form; a genuine arity error then
    surfaces from the call itself.
    """
    try:
        signature = inspect.signature(factory)
    except (TypeError, ValueError):
        return True
    try:
        signature.bind(None, None)
    except TypeError:
        return False
    return True


def _channel_for(geo, factory, set_volume: bool, volume_fn=None):
    volume = _resolved_volume(geo, volume_fn, required=set_volume)
    if volume == 0.0:
        return None
    if not set_volume:
        return factory(geo)

    # Passing the volume into construction is required because the C++ channel
    # rejects unsupported Volume-mode geometries at construction, before
    # SetVolume can be called. Custom factories should therefore accept
    # (geometry, volume).
    if not _factory_accepts_exact_volume(factory):
        siren = _siren()
        configuration_error = getattr(
            getattr(siren, "utilities", None), "ConfigurationError", TypeError)
        raise configuration_error(
            "set_volume=True requires factory(geometry, exact_volume)")
    return factory(geo, volume)


def _wrap(tiles: Sequence, factory, set_volume: bool, volume_fn=None):
    """Wrap tile geometries as directed channels in a NestedMixtureChannel.

    Returns the NestedMixtureChannel (for >1 tile) or the bare channel (1 tile).
    """
    siren = _siren()
    inj = siren.injection
    chans = [
        channel for channel in
        (_channel_for(g, factory, set_volume, volume_fn) for g in tiles)
        if channel is not None
    ]
    if not chans:
        raise ValueError("directed tiling contains no non-empty tiles")
    if len(chans) == 1:
        return chans[0]
    inner = inj.MultiChannelPhaseSpace()
    inner.channels = chans
    inner.weights = [1.0 / len(chans)] * len(chans)
    return inj.NestedMixtureChannel(inner)


def grid_cells(box, n) -> List:
    """Partition an axis-aligned Box into n_x*n_y*n_z disjoint sub-box cells."""
    siren = _siren()
    G = siren.geometry
    V3 = siren.math.Vector3D
    if not isinstance(box, G.Box):
        from siren.utilities import ConfigurationError
        raise ConfigurationError("grid_cells requires a siren.geometry.Box region")
    q = box.placement.Quaternion
    q = q() if callable(q) else q
    # Sub-cells are emitted axis-aligned in world coordinates, so a rotated box
    # would neither partition nor cover its true (rotated) volume. Require the
    # placement rotation to be the identity: (x, y, z) == 0 and |w| == 1.
    if not (abs(q.X) < 1e-9 and abs(q.Y) < 1e-9 and abs(q.Z) < 1e-9
            and abs(abs(q.W) - 1.0) < 1e-9):
        configuration_error = getattr(
            getattr(siren, "utilities", None), "ConfigurationError", ValueError)
        raise configuration_error(
            "grid_cells requires an axis-aligned Box; the region carries a "
            "non-identity placement rotation (quaternion x={:.3g} y={:.3g} "
            "z={:.3g} w={:.3g})".format(q.X, q.Y, q.Z, q.W))
    nx, ny, nz = (n, n, n) if isinstance(n, int) else n
    cx, cy, cz = _center(box)
    wx, wy, wz = box.X, box.Y, box.Z
    sx, sy, sz = wx / nx, wy / ny, wz / nz
    cells = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                ox = cx - wx / 2 + (i + 0.5) * sx
                oy = cy - wy / 2 + (j + 0.5) * sy
                oz = cz - wz / 2 + (k + 0.5) * sz
                cells.append(G.Box(G.Placement(V3(ox, oy, oz)), sx, sy, sz))
    return cells


def disjointify(volumes: Sequence) -> List:
    """Disjoint shells T_k = V_k \\ (V_1 union ... union V_{k-1}).

    Union of shells = union of inputs (full coverage); the shells are pairwise
    disjoint, so they are a non-degenerate channel basis.
    """
    siren = _siren()
    G = siren.geometry
    vols = list(volumes)
    tiles = [vols[0]]
    for k in range(1, len(vols)):
        earlier = make_union(vols[:k])
        tiles.append(
            G.BooleanGeometry(G.BooleanOperation.SUBTRACTION, vols[k], earlier))
    return tiles


# ------------------------------------------------------------------ #
#  Public builders                                                    #
# ------------------------------------------------------------------ #

def build_grid_tiling(box, n=2, factory=None, daughter_index: int = 0,
                      set_volume: bool = True, n_mc: int = 100000,
                      volume_fn=None):
    """Tile an axis-aligned Box region into n^3 disjoint sub-box directed
    channels, wrapped in a NestedMixtureChannel.  The primary workhorse."""
    if factory is None:
        factory = default_2body_factory(daughter_index)
    del n_mc
    return _wrap(grid_cells(box, n), factory, set_volume, volume_fn)


def build_subtraction_tiling(volumes: Sequence, factory=None,
                             daughter_index: int = 0, set_volume: bool = True,
                             n_mc: int = 100000, volume_fn=None):
    """Disjointify overlapping `volumes` into shells and wrap them as a
    NestedMixtureChannel of directed channels.

    With ``set_volume=True``, ``volume_fn`` must return each composite tile's
    exact volume. Zero-volume tiles are removed before mixture weights are
    assigned. ``n_mc`` is retained only for source compatibility.
    """
    if factory is None:
        factory = default_2body_factory(daughter_index)
    del n_mc
    return _wrap(disjointify(list(volumes)), factory, set_volume, volume_fn)


def build_union_tile(volumes: Sequence, factory=None, daughter_index: int = 0,
                     dedup: bool = True, set_volume: bool = True,
                     n_mc: int = 100000, volume_fn=None):
    """Collapse `volumes` (optionally deduped) into ONE union directed channel.

    Returns a single channel (the union tile).  Use for a tightly-nested cluster
    where partitioning adds no optimizer value. With ``set_volume=True``, an
    exact ``volume_fn`` is required when the result is composite.
    """
    if factory is None:
        factory = default_2body_factory(daughter_index)
    vols = dedup_contained(volumes) if dedup else list(volumes)
    del n_mc
    channel = _channel_for(
        make_union(vols), factory, set_volume, volume_fn)
    if channel is None:
        raise ValueError("union tile is empty")
    return channel


def build_angular_tiling(target, n_theta: int = 2, n_phi: int = 4, factory=None,
                         daughter_index: int = 0, mode=None):
    """Tile the cone subtended by `target` into n_theta*n_phi (cos theta, phi)
    angular-sector directed channels, wrapped in a NestedMixtureChannel.

    Requires `siren.injection.DetectorDirectedAngularSectorChannel` (Phase 3).
    """
    siren = _siren()
    inj = siren.injection
    if not hasattr(inj, "DetectorDirectedAngularSectorChannel"):
        raise NotImplementedError(
            "DetectorDirectedAngularSectorChannel is not available; build and "
            "install the injection module with the angular-sector channel "
            "(Phase 3).")
    Sector = inj.DetectorDirectedAngularSectorChannel
    chans = []
    for it in range(n_theta):
        u_lo = it / n_theta                  # fraction of the bounding cone
        u_hi = (it + 1) / n_theta
        for ip in range(n_phi):
            phi_lo = 2.0 * math.pi * ip / n_phi
            phi_hi = 2.0 * math.pi * (ip + 1) / n_phi
            chans.append(Sector(target, u_lo, u_hi, phi_lo, phi_hi, daughter_index))
    inner = inj.MultiChannelPhaseSpace()
    inner.channels = chans
    inner.weights = [1.0 / len(chans)] * len(chans)
    return inj.NestedMixtureChannel(inner)


def build_directed_tiling(volumes: Sequence, factory=None, daughter_index: int = 0,
                          strategy: str = "union", vertices=None,
                          dedup: bool = True, grid_n=2, set_volume: bool = True,
                          n_mc: int = 100000, volume_fn=None):
    """Top-level: dedup -> cluster by cone overlap -> tile each cluster -> wrap
    all tiles (flattened) in one NestedMixtureChannel.

    strategy per cluster:
      "union"       -- one union tile per cluster (collapse; default).
      "subtraction" -- disjointify the cluster's volumes into shells.
      "grid"        -- grid-tile each single-Box cluster (else union).

    With ``set_volume=True``, ``volume_fn`` must return the exact volume of
    every composite tile. Empty tiles are dropped before weights are assigned;
    ``n_mc`` is retained only for source compatibility and is not used.

    Flattening all clusters' tiles into one inner mixture lets the existing
    single-level optimizer recursion tune every tile weight.
    """
    if factory is None:
        factory = default_2body_factory(daughter_index)
    vols = dedup_contained(volumes) if dedup else list(volumes)
    clusters = cluster_by_cone_overlap(vols, vertices)

    siren = _siren()
    G = siren.geometry
    tiles: List = []
    for cluster in clusters:
        if strategy == "union":
            tiles.append(make_union(cluster))
        elif strategy == "subtraction":
            tiles.extend(disjointify(cluster))
        elif strategy == "grid":
            if len(cluster) == 1 and isinstance(cluster[0], G.Box):
                tiles.extend(grid_cells(cluster[0], grid_n))
            else:
                tiles.append(make_union(cluster))
        else:
            raise ValueError(f"unknown strategy {strategy!r}")
    del n_mc
    return _wrap(tiles, factory, set_volume, volume_fn)
