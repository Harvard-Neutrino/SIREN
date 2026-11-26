"""Filled plane sections of placed triangle meshes (numpy only, display only).

Every volume fills the part of its section that none of its daughters occupies,
so hiding a daughter leaves an opening instead of showing its mother's material.
"""
import numpy as np

_GRID_MM = 1e-5  # endpoints shared by two triangles must agree to this


def plane_axes(normal, up=None):
    """Right-handed in-plane axes (u, w): u x w = normal, w along *up* if given."""
    n = np.asarray(normal, dtype=float)
    length = np.linalg.norm(n)
    if not np.isfinite(length) or length == 0:
        raise ValueError("section normal must be a nonzero vector")
    n = n / length
    for hint in (up, (0., 0., 1.), (0., 1., 0.)):
        if hint is None:
            continue
        w = np.asarray(hint, dtype=float)
        w = w - n * (w @ n)
        if np.linalg.norm(w) > 1e-6:
            w = w / np.linalg.norm(w)
            return n, np.cross(w, n), w
    raise ValueError("section 'up' is parallel to its normal")


def mesh_segments(vertices, faces, normal, offset):
    """End points of the cut through a mesh: where normal.x == offset.

    A vertex on the plane counts as above it, which keeps the outline closed
    when the plane touches vertices, edges or whole faces; a plane holding a
    face shows the material below it, the side a section keeps.
    """
    distance = vertices @ normal - offset
    above = (distance >= 0)[faces]
    count = above.sum(axis=1)
    rows = np.flatnonzero((count == 1) | (count == 2))
    if not len(rows):
        return np.empty((0, 2, 3))
    above = above[rows]
    lone = np.where(count[rows] == 1, above.argmax(axis=1), (~above).argmax(axis=1))
    k = np.arange(len(rows))
    tri = faces[rows]
    a, b, c = (tri[k, (lone + i) % 3] for i in range(3))
    def crossing(p, q):
        t = distance[p] / (distance[p] - distance[q])
        return vertices[p] + t[:, None] * (vertices[q] - vertices[p])
    return np.stack((crossing(a, b), crossing(a, c)), axis=1)


def fill(own, holes=()):
    """Trapezoids covering the inside of *own* minus the inside of every hole.

    *own* and each of *holes* are unordered (N, 2, 2) segment arrays.
    """
    if not len(own):
        return np.empty((0, 4, 2))
    # What lies wholly beside the outline cannot change it.
    low, high = own.reshape(-1, 2).min(axis=0), own.reshape(-1, 2).max(axis=0)
    holes = [h for h in holes if len(h) and (h.reshape(-1, 2).min(axis=0) < high).all() and
             (h.reshape(-1, 2).max(axis=0) > low).all()]
    return _sweep([own, *holes], lambda who, cover: cover(who == 0) & ~cover(who != 0))


def combine(shape):
    """Trapezoids of a Boolean solid's cut, from the cuts of its operands.

    *shape* is a segment array or (operation, left, right) with "union",
    "intersection" or "subtraction". The regions are combined in the plane, so
    no Boolean of meshes is needed.
    """
    leaves = []

    def index(node):
        if isinstance(node, tuple):
            return node[0], index(node[1]), index(node[2])
        leaves.append(node)
        return len(leaves) - 1
    tree = index(shape)

    def rule(who, cover):
        def inside(node):
            if not isinstance(node, tuple):
                return cover(who == node)
            operation, left, right = node[0], inside(node[1]), inside(node[2])
            if operation == "union":
                return left | right
            return left & (right if operation == "intersection" else ~right)
        return inside(tree)
    return _sweep(leaves, rule)


def sides(quads):
    """The slanted sides of trapezoids: segments that fill() turns back into them."""
    return np.concatenate((quads[:, (0, 3)], quads[:, (1, 2)]))


def _sweep(volumes, rule):
    """Cover the part of the plane where rule(who, cover) holds.

    The even-odd rule is applied to each volume strip by strip, so outlines need
    not be chained into loops, a mesh with T-junctions still fills, and volumes
    may overlap. ``who`` is the volume crossed at each step along the strips and
    ``cover(mask)`` tells where a volume in the mask covers the interval that
    follows. Returns (K, 4, 2) corners, counter-clockwise.
    """
    empty = np.empty((0, 4, 2))
    segments = np.concatenate(volumes)
    owner = np.repeat(np.arange(len(volumes)), [len(v) for v in volumes])
    segments = np.round(segments / _GRID_MM) * _GRID_MM
    keep = segments[:, 0, 1] != segments[:, 1, 1]
    segments, owner = segments[keep], owner[keep]
    if not len(segments):
        return empty
    swap = segments[:, 0, 1] > segments[:, 1, 1]
    segments[swap] = segments[swap][:, ::-1]
    x0, y0 = segments[:, 0].T
    x1, y1 = segments[:, 1].T
    levels = np.unique(np.concatenate((y0, y1)))
    first, last = np.searchsorted(levels, y0), np.searchsorted(levels, y1)
    spans = last - first
    segment = np.repeat(np.arange(len(segments)), spans)
    strip = np.arange(spans.sum()) - np.repeat(np.cumsum(spans) - spans, spans) + first[segment]
    middle = .5 * (levels[strip] + levels[strip + 1])
    slope = (x1 - x0) / (y1 - y0)
    x = x0[segment] + slope[segment] * (middle - y0[segment])
    # Along each strip a volume is entered and left in turn.
    who = owner[segment]
    order = np.lexsort((x, who, strip))
    first = np.r_[True, (strip[order][1:] != strip[order][:-1]) | (who[order][1:] != who[order][:-1])]
    place = np.arange(len(order))
    step = np.empty(len(order), dtype=int)
    step[order] = np.where((place - np.maximum.accumulate(np.where(first, place, 0))) % 2 == 0, 1, -1)
    order = np.lexsort((x, strip))
    strip, segment, step, who = strip[order], segment[order], step[order], who[order]
    starts = np.r_[True, strip[1:] != strip[:-1]]
    begin = np.maximum.accumulate(np.where(starts, place, 0))
    def cover(chosen):
        # Counted from the start of each strip: a gap in one must not spread.
        total = np.r_[0, np.cumsum(step * chosen)]
        return total[1:] - total[begin] > 0
    inside = rule(who, cover)
    inside[-1] = False
    inside[:-1] &= strip[1:] == strip[:-1]
    left = np.flatnonzero(inside & np.r_[True, ~inside[:-1] | starts[1:]])
    right = np.flatnonzero(inside & np.r_[~inside[1:], True]) + 1
    if not len(left):
        return empty
    # Join the same pair of sides over consecutive strips into one trapezoid.
    ls, rs, k = segment[left], segment[right], strip[left]
    order = np.lexsort((k, rs, ls))
    ls, rs, k = ls[order], rs[order], k[order]
    new = np.r_[True, (ls[1:] != ls[:-1]) | (rs[1:] != rs[:-1]) | (k[1:] != k[:-1] + 1)]
    begin = np.flatnonzero(new)
    end = np.r_[begin[1:], len(k)] - 1
    ls, rs, low, high = ls[begin], rs[begin], levels[k[begin]], levels[k[end] + 1]
    def at(s, y):
        return x0[s] + slope[s] * (y - y0[s])
    quads = np.empty((len(low), 4, 2))
    quads[:, 0, 0], quads[:, 1, 0] = at(ls, low), at(rs, low)
    quads[:, 2, 0], quads[:, 3, 0] = at(rs, high), at(ls, high)
    quads[:, :2, 1], quads[:, 2:, 1] = low[:, None], high[:, None]
    width = quads[:, 1, 0] - quads[:, 0, 0] + quads[:, 2, 0] - quads[:, 3, 0]
    return quads[width > 4 * _GRID_MM]


def area(quads):
    """Total area of fill()'s trapezoids."""
    width = quads[:, 1, 0] - quads[:, 0, 0] + quads[:, 2, 0] - quads[:, 3, 0]
    return float((.5 * width * (quads[:, 2, 1] - quads[:, 1, 1])).sum())


class Section:
    """One plane cut through a scene, as 2D pieces in the (u, w) plane axes."""

    def __init__(self, origin, normal, up=None):
        self.normal, self.u, self.w = plane_axes(normal, up)
        self.origin = np.asarray(origin, dtype=float)
        self.offset = float(self.origin @ self.normal)
        self.pieces = {}    # instance index -> (K, 4, 2) trapezoids
        self.outlines = {}  # instance index -> (N, 2, 2) segments

    def to_plane(self, points):
        d = np.asarray(points) - self.origin
        return np.stack((d @ self.u, d @ self.w), axis=-1)

    def to_world(self, points):
        points = np.asarray(points)
        return self.origin + points[..., :1] * self.u + points[..., 1:] * self.w

    def find(self, point):
        """Index of the placed volume filling this point of the plane, or None."""
        x, y = self.to_plane(point)
        for index, q in self.pieces.items():
            row = q[(q[:, 0, 1] <= y) & (y <= q[:, 2, 1])]
            if not len(row):
                continue
            t = (y - row[:, 0, 1]) / (row[:, 2, 1] - row[:, 0, 1])
            lo = row[:, 0, 0] + t * (row[:, 3, 0] - row[:, 0, 0])
            hi = row[:, 1, 0] + t * (row[:, 2, 0] - row[:, 1, 0])
            if ((lo <= x) & (x <= hi)).any():
                return index
        return None


class SceneSections:
    """Cuts a scene's instances, holding only their placements and bounds."""

    def __init__(self, instances, meshes, exact=None):
        """*exact* replaces shapes: instance index -> function giving a mesh
        (vertices, faces) in scene coordinates, a tree of them as in combine(),
        or None to leave that instance out. Its mesh still bounds the instance."""
        self.instances, self.meshes = instances, meshes
        self.exact = {} if exact is None else exact
        self.daughters = {}
        for index, row in enumerate(instances):
            self.daughters.setdefault(row.get("parent", -1), []).append(index)
        # Sectors of a SIREN model overlap; the one with the higher level wins.
        self.levels = {i: row["level"] for i, row in enumerate(instances) if "level" in row}
        self.known = np.array([i for i, row in enumerate(instances)
                               if row["prototype"] in meshes], dtype=int)
        bounds = {name: (v.min(axis=0), v.max(axis=0)) for name, (v, _) in meshes.items()}
        rows = [instances[i] for i in self.known]
        self.matrices = np.array([r["matrix"] for r in rows], dtype=float).reshape(-1, 4, 4)
        self.low = np.array([bounds[r["prototype"]][0] for r in rows]).reshape(-1, 3)
        self.high = np.array([bounds[r["prototype"]][1] for r in rows]).reshape(-1, 3)

    def _cut_shape(self, shape, section):
        if len(shape) == 3:
            return shape[0], self._cut_shape(shape[1], section), self._cut_shape(shape[2], section)
        return section.to_plane(mesh_segments(*shape, section.normal, section.offset))

    def cut(self, origin, normal, up=None, drawn=None):
        """Cut every instance in *drawn* (default: all with meshes)."""
        section = Section(origin, normal, up)
        # The plane n.(Mx + t) = c in each instance's own coordinates.
        local = np.einsum("kji,j->ki", self.matrices[:, :3, :3], section.normal)
        offsets = section.offset - self.matrices[:, :3, 3] @ section.normal
        center = (local * (self.low + self.high)).sum(axis=1) / 2 - offsets
        reach = (np.abs(local) * (self.high - self.low)).sum(axis=1) / 2
        segments = {}
        for row in np.flatnonzero(np.abs(center) <= reach + 1e-9):
            index = int(self.known[row])
            if index in self.exact:
                shape = self.exact[index]()
                if shape is not None:
                    cut = sides(combine(self._cut_shape(shape, section)))
                    if len(cut):
                        segments[index] = cut
                continue
            vertices, faces = self.meshes[self.instances[index]["prototype"]]
            cut = mesh_segments(vertices, faces, local[row], offsets[row])
            if len(cut):
                world = cut @ self.matrices[row, :3, :3].T + self.matrices[row, :3, 3]
                segments[index] = section.to_plane(world)
        levelled = np.array([i for i in segments if i in self.levels], dtype=int)
        levels = np.array([self.levels[i] for i in levelled], dtype=float)
        low = np.array([segments[i].min(axis=(0, 1)) for i in levelled]).reshape(-1, 2)
        high = np.array([segments[i].max(axis=(0, 1)) for i in levelled]).reshape(-1, 2)
        for index, own in segments.items():
            if drawn is not None and self.instances[index]["prototype"] not in drawn:
                continue
            if index in self.levels:
                over = levelled[(levels > self.levels[index]) &
                                (low < own.max(axis=(0, 1))).all(axis=1) &
                                (high > own.min(axis=(0, 1))).all(axis=1)]
            else:
                over = self.daughters.get(index, ())
            quads = fill(own, [segments[i] for i in over if i in segments])
            if len(quads):
                section.pieces[index] = quads
            section.outlines[index] = own
        return section
