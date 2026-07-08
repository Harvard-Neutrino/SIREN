"""Selection predicates for use as ``tune()`` metrics or standalone cuts.

Each function here returns a predicate ``callable(tree) -> bool``, not a
``tune()`` metric directly (``tune()``'s ``metric`` is
``callable(event, weight) -> float``).  Wrap a predicate to use it as a
tune metric, e.g.::

    is_in = vertex_in(fiducial_geometry)
    report = tune(injector, weighter,
                  metric=lambda ev, w: w if is_in(ev) else 0.0)
"""

from __future__ import annotations

from typing import Callable

from . import particles as _particles


def vertex_in(geometry, particle=None) -> Callable[[object], bool]:
    """Build a predicate: does any vertex of a tree lie inside ``geometry``.

    Parameters
    ----------
    geometry : siren geometry object
        Anything exposing ``IsInside(position)`` (e.g. a ``siren.geometry``
        shape or a ``BooleanGeometry`` combination).
    particle : str or ParticleType, optional
        If given (resolved via ``siren.particles.resolve``), only vertices
        whose record involves this particle (as primary, target, or a
        secondary) are checked.  If ``None`` (default), every vertex is
        checked regardless of which particles it involves.

    Returns
    -------
    callable(tree) -> bool
        True if at least one qualifying interaction vertex in ``tree``
        (an ``InteractionTree``, i.e. a list of ``InteractionTreeDatum``)
        lies inside ``geometry``.  A tree with no qualifying vertex, or an
        empty tree, returns False.
    """
    from . import math as _math

    ptype = _particles.resolve(particle) if particle is not None else None

    def _involves(record) -> bool:
        if ptype is None:
            return True
        sig = record.signature
        if sig.primary_type == ptype or sig.target_type == ptype:
            return True
        return any(t == ptype for t in sig.secondary_types)

    def _predicate(tree) -> bool:
        # An InteractionTree exposes its datum list as .tree; accept a bare
        # datum list too.
        data = getattr(tree, "tree", tree)
        for datum in data:
            record = datum.record
            if not _involves(record):
                continue
            pos = _math.Vector3D(record.interaction_vertex)
            if geometry.IsInside(pos):
                return True
        return False

    return _predicate
