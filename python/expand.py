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
import math

from . import particles as _particles

__all__ = [
    "child",
    "depth_below",
    "compile_expansion",
    "VertexSpec",
]


_VertexSpecBase = namedtuple(
    "VertexSpec", ["particle", "expand", "continue_if", "secondary_types",
                   "secondary_signatures"]
)
_VertexSpecBase.__new__.__defaults__ = (None, None)


class VertexSpec(_VertexSpecBase):
    """Minimal duck-typed vertex description consumed by this module.

    Optional daughter metadata lets validation check actual reachability,
    including positional child rules. Compilation uses the expansion policy.
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
    """Separate a repeated primary type from its secondary declaration."""
    by_particle = {}
    for vertex in vertices[1:]:
        ptype = _particles.resolve(vertex.particle)
        if ptype in by_particle:
            from .errors import ConfigurationError
            raise ConfigurationError("Duplicate secondary vertex for {!r}".format(ptype))
        by_particle[ptype] = vertex
    primary = None
    if vertices:
        ptype = _particles.resolve(vertices[0].particle)
        if ptype in by_particle:
            primary = vertices[0]
        else:
            by_particle[ptype] = vertices[0]
    return by_particle, primary


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

    The first vertex is the primary. A secondary with the same particle type
    keeps its own policy; otherwise policies retain their lookup by type.
    """
    vertices = list(vertices)
    from .errors import ConfigurationError
    for vertex in vertices:
        if any(not isinstance(rule, (_ChildRule, _DepthBelowRule)) for rule in vertex.expand):
            raise ConfigurationError("expand entries must be child() or depth_below() rules")
        if vertex.continue_if is not None and not callable(vertex.continue_if):
            raise ConfigurationError("continue_if must be callable or None")
    by_particle, primary = _index_vertices_by_particle(vertices)

    if all(v.continue_if is None for v in vertices):
        from . import injection
        rules = _native_rules(by_particle)
        if primary is not None:
            return injection.SecondaryExpansion(
                rules, _native_rules({_particles.resolve(primary.particle): primary}))
        return injection.SecondaryExpansion(rules)
    return _ExpansionCallback(by_particle, primary)


def _native_rules(by_particle):
    rules = []
    for parent, vertex in by_particle.items():
        for rule in vertex.expand:
            if isinstance(rule, _ChildRule):
                index = rule.index
                if index is not None:
                    if not 0 <= index < 2**63 or index != int(index):
                        continue  # no possible child index matches this rule
                    index = int(index)
                rules.append((int(parent), int(_particles.resolve(rule.name)),
                              -1 if index is None else index, -1, 0))
            elif isinstance(rule, _DepthBelowRule):
                n = rule.n
                if not n > 0:  # includes NaN, matching depth < NaN == False
                    continue
                limit = -1 if n >= 2**63 else math.ceil(n)
                rules.append((int(parent), 0, -1, limit, 1))
    return rules

class _ExpansionCallback:
    """Picklable policy when its user-supplied predicates are picklable."""
    def __init__(self, by_particle, primary):
        self.by_particle = by_particle
        self.primary = primary

    def __call__(self, tree, parent, i):
        primary_type = parent.record.signature.primary_type
        depth = parent.depth(tree)
        vertex = (self.primary if self.primary is not None and depth == 0
                  else self.by_particle.get(primary_type))
        if vertex is None or _particles.resolve(vertex.particle) != primary_type:
            return True

        secondary_types = parent.record.signature.secondary_types
        ptype = secondary_types[i]

        matched = any(
            _rule_permits(rule, ptype, i, depth) for rule in vertex.expand
        )
        if not matched:
            return True

        if vertex.continue_if is not None:
            if not vertex.continue_if(tree, parent, i):
                return True

        return False
