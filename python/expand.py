"""Expansion vocabulary for chain-vertex secondary recursion.

An expansion rule is an instruction attached to a chain vertex saying
which of its secondaries should recurse into their own vertex. Rules are
combined with the vertex's optional ``continue_if`` predicate and
compiled into one ``stopping_condition(tree, parent, i)`` callable with
the engine's polarity: the engine calls it for secondary index ``i`` of
``parent`` and, when it returns True, that secondary is pruned (not
expanded). False means the secondary is expanded. See
``Injector.SetStoppingCondition`` and the call site in
``projects/injection/private/Injector.cxx``.

Rule types
----------
``child(name, index=None)``
    Expand the secondary named ``name`` (a particle name or
    ParticleType, resolved against the parent signature at compile
    time). If ``index`` is given, only the occurrence of that particle
    at that position among the parent's secondaries matches.

``depth_below(n)``
    Expand every secondary while ``parent.depth(tree) < n``. Independent
    of particle identity.

Vertex description (duck-typed)
--------------------------------
``compile_expansion`` and ``validate_expansion_wiring`` (in
``_validation.py``) accept ``vertices`` as any iterable of objects with:

``particle``
    The ParticleType or particle name this vertex produces/handles
    (i.e. the parent particle whose secondaries this vertex's ``expand``
    list governs).
``expand``
    A tuple of rule objects from ``child()`` / ``depth_below()``.
``continue_if``
    ``None`` or a callable ``(tree, parent, i) -> bool``; True means the
    secondary may proceed (subject also to the expand-list match), False
    forces a prune regardless of the expand list.

A plain namedtuple (``VertexSpec`` below) satisfies this; any object
exposing the same three attributes works too. This module never imports
``python/vertex.py`` -- the real Vertex compiler builds these records
from Vertex objects itself.
"""

from __future__ import annotations

from collections import namedtuple

from . import particles as _particles

__all__ = [
    "child",
    "depth_below",
    "compile_expansion",
    "VertexSpec",
]


_VertexSpecBase = namedtuple(
    "VertexSpec", ["particle", "expand", "continue_if", "secondary_types"]
)
_VertexSpecBase.__new__.__defaults__ = (None,)


class VertexSpec(_VertexSpecBase):
    """Minimal duck-typed vertex description consumed by this module.

    `secondary_types` is optional and only consulted by
    `_validation.validate_expansion_wiring` for its "secondaries with no
    expansion declaration" check; `compile_expansion` never reads it.
    """


class _ChildRule:
    """Expansion rule matching a secondary by particle type (and index)."""

    __slots__ = ("name", "index")

    def __init__(self, name, index=None):
        self.name = name
        self.index = index

    def _matches(self, ptype, i):
        if self.index is not None and i != self.index:
            return False
        return ptype == _particles.resolve(self.name)

    def __repr__(self):
        if self.index is None:
            return f"child({self.name!r})"
        return f"child({self.name!r}, index={self.index!r})"


class _DepthBelowRule:
    """Expansion rule matching every secondary while parent depth < n."""

    __slots__ = ("n",)

    def __init__(self, n):
        self.n = n

    def _matches_depth(self, depth):
        return depth < self.n

    def __repr__(self):
        return f"depth_below({self.n!r})"


def child(name, index=None):
    """Return a rule expanding the secondary named `name`.

    Parameters
    ----------
    name : str or ParticleType
        Particle identifying which secondary to expand.
    index : int, optional
        Restrict the match to the secondary at this position among the
        parent's ``secondary_types``. Without it, every occurrence of
        `name` matches.
    """
    return _ChildRule(name, index=index)


def depth_below(n):
    """Return a rule expanding all secondaries while parent depth < n."""
    return _DepthBelowRule(n)


def _index_vertices_by_particle(vertices):
    """Map resolved ParticleType -> vertex spec, first match wins."""
    by_particle = {}
    for v in vertices:
        ptype = _particles.resolve(v.particle)
        if ptype not in by_particle:
            by_particle[ptype] = v
    return by_particle


def _rule_permits(rule, ptype, i, depth):
    if isinstance(rule, _ChildRule):
        return rule._matches(ptype, i)
    if isinstance(rule, _DepthBelowRule):
        return rule._matches_depth(depth)
    return False


def compile_expansion(vertices):
    """Compile per-vertex expand rules into one stopping_condition callable.

    Parameters
    ----------
    vertices : iterable
        Objects with ``.particle``, ``.expand`` (tuple of rules from
        ``child()``/``depth_below()``), and ``.continue_if``
        (callable or None). See module docstring for the exact contract.

    Returns
    -------
    callable
        ``stopping_condition(tree, parent, i) -> bool`` in the engine's
        native polarity: True prunes secondary `i` of `parent`, False
        expands it. Safe to pass directly to
        ``Injector.SetStoppingCondition`` / the ``stopping_condition``
        constructor argument.
    """
    by_particle = _index_vertices_by_particle(vertices)

    def stopping_condition(tree, parent, i):
        primary_type = parent.record.signature.primary_type
        vertex = by_particle.get(primary_type)
        if vertex is None:
            return True

        secondary_types = parent.record.signature.secondary_types
        ptype = secondary_types[i]
        depth = parent.depth(tree)

        matched = any(
            _rule_permits(rule, ptype, i, depth) for rule in vertex.expand
        )
        if not matched:
            return True

        if vertex.continue_if is not None:
            if not vertex.continue_if(tree, parent, i):
                return True

        return False

    return stopping_condition
