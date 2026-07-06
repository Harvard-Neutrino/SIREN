"""Vertex: one chain node's interactions, distributions, and phase space.

A Vertex describes everything the engine needs to compile one
PrimaryInjectionProcess (root) or SecondaryInjectionProcess (non-root)
node: which particle it is, which Decay/CrossSection models drive it,
which injection distributions sample its free parameters, and which
channels.Mixture/_directed.Directed biases its phase-space sampling.
Vertex.compile() builds the engine process object; Vertex.as_vertex_spec()
builds the expand.VertexSpec consumed by expand.compile_expansion().

A Vertex holds strong references to every Python object it is handed
(models, distributions, mixtures) for the lifetime of the Vertex, so a
Python-subclassed trampoline (a DarkNews model, a custom distribution)
is not collected out from under the compiled C++ process, which only
holds it via a pybind trampoline pointer.
"""

from __future__ import annotations

from . import particles as _particles
from . import expand as _expand
from .errors import ConfigurationError

__all__ = ["Vertex"]


def _siren():
    import siren
    return siren


def _resolve_weighting_mode(weighting):
    """Resolve `weighting` to an engine VertexWeightingMode.

    Accepts an already-built VertexWeightingMode (from either
    siren.injection.VertexWeightingMode -- where it actually lives -- or
    a future siren.dataclasses.VertexWeightingMode / siren.Propagated /
    siren.Fixed re-export), or None (defaults to Propagated()).
    """
    siren = _siren()
    if weighting is None:
        return siren.injection.VertexWeightingMode.Propagated()
    return weighting


def _as_list(interactions):
    if interactions is None:
        return []
    if isinstance(interactions, (list, tuple)):
        return list(interactions)
    return [interactions]


def _model_signatures(model):
    """GetPossibleSignatures() for a Decay or CrossSection model.

    Both bases expose the same method name; this just centralizes the
    call so Vertex.compile() does not care which base `model` derives
    from.
    """
    return list(model.GetPossibleSignatures())


def _as_compilable(kinematics):
    """Normalize `kinematics` to something exposing compile(sig, ...).

    channels.Mixture, _directed.Directed, and _directed's
    by_signature() result all already expose compile(signature, *,
    detector, models); a bare channels.Channel/Tiling (the un-mixed
    single-channel case) does not, so it is wrapped in a
    one-entry Mixture.
    """
    if hasattr(kinematics, "compile"):
        return kinematics
    from . import channels as _channels
    if isinstance(kinematics, _channels._MixtureNode):
        return _channels.Mixture([(kinematics.weight, kinematics)])
    raise ConfigurationError(
        "Vertex(kinematics=...): expected a channels.Mixture, "
        "channels.Channel/Tiling, or _directed.Directed, got {!r}".format(
            type(kinematics).__name__))


class Vertex:
    """One chain node: particle, interactions, distributions, phase space.

    Parameters
    ----------
    particle : str or ParticleType
        The particle this vertex is compiled for (the primary of the
        root vertex, or the secondary type of a non-root vertex).
    interactions : Decay/CrossSection or list thereof
        The interaction model(s) driving this vertex.
    distributions : list, optional
        PrimaryInjectionDistribution / SecondaryInjectionDistribution
        objects sampling this vertex's free parameters.
    position : optional
        Convenience vertex-position distribution, appended to
        `distributions` if given.
    physical : list, optional
        Physical-side distributions kept on the Vertex for the Weighter
        to read; not consumed by Vertex.compile() itself.
    kinematics : channels.Mixture, channels.Channel, _directed.Directed, optional
        Phase-space biasing description, compiled once per signature
        enumerated from `interactions`. None means no phase space is
        registered (the engine's default sampler is used).
    weighting : VertexWeightingMode, optional
        Defaults to Propagated(). Accepts the engine enum value however
        it is currently exposed (siren.injection.VertexWeightingMode
        today; siren.dataclasses.VertexWeightingMode / siren.Propagated /
        siren.Fixed once re-exported).
    expand : tuple, optional
        Expansion rules from expand.child() / expand.depth_below().
    continue_if : callable, optional
        (tree, parent, i) -> bool predicate consulted by
        expand.compile_expansion().
    label : str, optional
        Human-readable label, not consumed by the engine.

    Everything after `interactions` is keyword-only. Field storage uses
    __slots__: assigning an unknown attribute raises AttributeError
    instead of silently creating a new one, so a typo'd field name (a
    common failure mode of a duck-typed config surface) is caught
    immediately rather than silently ignored.
    """

    __slots__ = (
        "particle",
        "_resolved_particle",
        "interactions",
        "distributions",
        "physical",
        "kinematics",
        "weighting",
        "expand",
        "continue_if",
        "label",
    )

    def __init__(self, particle, interactions, *, distributions=None,
                 position=None, physical=None, kinematics=None,
                 weighting=None, expand=(), continue_if=None, label=None):
        self.particle = particle
        self._resolved_particle = _particles.resolve(particle)
        self.interactions = _as_list(interactions)
        if not self.interactions:
            raise ConfigurationError(
                "Vertex(particle={!r}): interactions must be a Decay/"
                "CrossSection or a non-empty list thereof".format(particle))

        dists = list(distributions) if distributions is not None else []
        if position is not None:
            dists.append(position)
        self.distributions = dists

        self.physical = list(physical) if physical is not None else []
        self.kinematics = kinematics
        self.weighting = weighting
        self.expand = tuple(expand)
        self.continue_if = continue_if
        self.label = label

    # ------------------------------------------------------------------ #
    #  Introspection                                                       #
    # ------------------------------------------------------------------ #

    def is_primary(self, *, is_primary):
        """Return whether this Vertex would compile as a primary process."""
        return bool(is_primary)

    def is_secondary(self, *, is_primary):
        """Return whether this Vertex would compile as a secondary process."""
        return not is_primary

    def _all_signatures(self):
        """(model, signature) pairs for every model's possible signatures."""
        pairs = []
        for model in self.interactions:
            for sig in _model_signatures(model):
                pairs.append((model, sig))
        return pairs

    def as_vertex_spec(self):
        """Build the expand.VertexSpec for this Vertex.

        `secondary_types` is the union, across every signature this
        vertex's models can produce, of the secondary particle types --
        consulted only by expand._validation's "secondaries with no
        expansion declaration" check, not by compile_expansion itself.
        """
        seen = []
        for _model, sig in self._all_signatures():
            for t in sig.secondary_types:
                if t not in seen:
                    seen.append(t)
        return _expand.VertexSpec(
            particle=self._resolved_particle,
            expand=self.expand,
            continue_if=self.continue_if,
            secondary_types=seen,
        )

    # ------------------------------------------------------------------ #
    #  Compilation                                                         #
    # ------------------------------------------------------------------ #

    def compile(self, *, is_primary, detector=None):
        """Compile this Vertex into an engine process.

        Returns a siren.injection.PrimaryInjectionProcess when
        `is_primary` is True, else a SecondaryInjectionProcess. Builds
        the InteractionCollection from `self.interactions`, assigns
        `self.distributions`, registers a phase space for every
        signature this vertex's models can produce (when `kinematics`
        is not None), and sets `weighting_mode`.
        """
        siren = _siren()

        interaction_collection = siren.interactions.InteractionCollection(
            self._resolved_particle, self.interactions)

        if is_primary:
            process = siren.injection.PrimaryInjectionProcess(
                self._resolved_particle, interaction_collection)
        else:
            process = siren.injection.SecondaryInjectionProcess(
                self._resolved_particle, interaction_collection)

        process.distributions = self.distributions

        if self.kinematics is not None:
            compilable = _as_compilable(self.kinematics)
            for model, sig in self._all_signatures():
                mcps = compilable.compile(
                    sig, detector=detector, models=[model])
                process.SetPhaseSpace(sig, mcps)

        process.weighting_mode = _resolve_weighting_mode(self.weighting)

        # `self.interactions`/`self.distributions`/`self.kinematics` are the
        # strong references keeping every Python model, distribution, and
        # mixture alive for as long as this Vertex is alive (the pybind
        # process objects above hold no __dict__ of their own to pin
        # keepalives on -- they are not built with dynamic_attr()). The
        # compiled process's Python-side trampolines (a DarkNews model, a
        # custom distribution) therefore survive gc.collect() exactly as
        # long as the Vertex that compiled them does; callers that want the
        # process to outlive the Vertex must keep the Vertex around too.
        return process
