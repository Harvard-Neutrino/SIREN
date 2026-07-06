"""Directed: sugar over channels.py for the common directed+fallback mix.

Directed(...) describes "point this daughter at a target volume, with a
physical/isotropic fallback so events that miss the target are not lost
from the estimator." It never introduces a new engine primitive: every
Directed compiles to exactly the channels.Mixture the channels algebra
would build by hand. Directed.by_signature({signature: Directed|Mixture})
covers the case where different signatures of the same vertex need
different directed treatments (e.g. a 2-body and a 3-body decay channel
of the same Decay model).
"""

from __future__ import annotations

from typing import Dict, Union

from . import channels as _channels
from .errors import ConfigurationError

__all__ = ["Directed"]


def _siren():
    import siren
    return siren


def _fallback(daughter):
    """physical() fallback, falling back further to isotropic() on failure.

    physical() late-binds a PhysicalDecayChannel/PhysicalCrossSectionChannel
    against the vertex's own interaction model at compile time; Directed
    prefers it as the fallback (it samples the real physical distribution
    rather than an isotropic proxy) but does not require a model be
    supplied -- Mixture.compile()/validate() will raise ConfigurationError
    from inside physical()'s factory if none is available at compile time.
    """
    return _channels.physical()


class Directed:
    """Sugar over channels.py for "directed daughter, with a fallback".

    Parameters
    ----------
    particle : str or ParticleType
        The daughter this Directed points at `toward`. Named, not
        indexed, so the same Directed can be reused across signatures
        where the daughter's position among secondaries differs.
    toward : Geometry
        The target volume passed to the underlying toward()/
        toward_3body()/scatter_toward() factory.
    fraction : float, optional
        Weight given to the directed term; `1 - fraction` goes to the
        fallback (physical() by default). Defaults to 0.98. Callers
        wanting the legacy 0.01/0.99 split pass fraction=0.99 explicitly
        -- this class never hardcodes that value.
    mode : optional
        DirectedMode passed through to the channel factory. Defaults to
        DirectedMode.Volume (channels.py's own default).
    spectator : str or ParticleType, optional
        Selects toward_3body()'s recursive strategy: the spectator
        daughter of a 3-body directed channel.
    pair_mass : channels.PairMass, optional
        Invariant-mass sampling law for a 3-body directed channel's pair.
    variable : optional
        ScatteringVariable selecting scatter_toward()'s sampled variable;
        its presence (along with q2_mode) is what selects the scattering
        factory over the 2-body/3-body ones.
    q2_mode : optional
        ScatteringQ2Mode passed through to scatter_toward().
    tiling : int, optional
        Grid size passed through to channels.tile() in place of a bare
        toward(); selects the tiled factory over the 2-body one.
    """

    def __init__(self, particle, toward, fraction: float = 0.98, *,
                 mode=None, spectator=None, pair_mass=None,
                 variable=None, q2_mode=None, tiling=None):
        if not (0.0 <= fraction <= 1.0):
            raise ConfigurationError(
                "Directed(fraction={!r}): fraction must be in [0, 1]".format(fraction))
        self.particle = particle
        self.toward = toward
        self.fraction = float(fraction)
        self.mode = mode
        self.spectator = spectator
        self.pair_mass = pair_mass
        self.variable = variable
        self.q2_mode = q2_mode
        self.tiling = tiling

    def _directed_channel(self):
        """Build the bare directed Channel/Tiling (no fallback mixed in)."""
        if self.tiling is not None:
            return _channels.tile(
                self.toward, self.tiling, daughter=self.particle,
                mode=self.mode)
        if self.variable is not None or self.q2_mode is not None:
            return _channels.scatter_toward(
                self.toward, variable=self.variable, q2_mode=self.q2_mode,
                mode=self.mode)
        if self.spectator is not None or self.pair_mass is not None:
            return _channels.toward_3body(
                self.particle, self.toward, spectator=self.spectator,
                pair_mass=self.pair_mass, mode=self.mode)
        return _channels.toward(self.particle, self.toward, mode=self.mode)

    def to_mixture(self, signature=None) -> "_channels.Mixture":
        """Build the channels.Mixture this Directed describes.

        `signature` is accepted (and ignored) so a plain Directed and a
        Directed.by_signature() mapping share the same call shape at the
        Vertex.compile() call site.
        """
        directed = self._directed_channel()
        if self.fraction >= 1.0:
            return _channels.Mixture([(1.0, directed)])
        fallback = _fallback(self.particle)
        return (self.fraction * directed) + ((1.0 - self.fraction) * fallback)

    def compile(self, signature, *, detector=None, models=None):
        """Compile this Directed's mixture for `signature`.

        Delegates to channels.Mixture.compile() after resolving
        to_mixture(); this is the same call shape Vertex.compile() uses
        for a bare channels.Mixture, so a Vertex's `kinematics=` accepts
        a Directed interchangeably with a Mixture.
        """
        return self.to_mixture(signature).compile(
            signature, detector=detector, models=models)

    def validate(self, signature=None):
        """Data-free configuration check, delegating to Mixture.validate()."""
        return self.to_mixture(signature).validate(signature)

    def describe(self) -> str:
        """The g(x) expression for this Directed's mixture."""
        return self.to_mixture().describe()

    @classmethod
    def by_signature(cls, mapping: Dict[object, Union["Directed", "_channels.Mixture"]]) -> "_DirectedBySignature":
        """A Directed-like object compiling a different mixture per signature.

        `mapping` keys are InteractionSignature (or anything else that
        compares equal to one); values are a Directed or a bare
        channels.Mixture. Looked up by exact signature match at compile
        time -- there is no partial/fallback matching.
        """
        return _DirectedBySignature(mapping)


class _DirectedBySignature:
    """Directed.by_signature(...) result: dispatches to_mixture() by signature."""

    __slots__ = ("_mapping",)

    def __init__(self, mapping):
        self._mapping = dict(mapping)

    def _entry_for(self, signature):
        if signature not in self._mapping:
            raise ConfigurationError(
                "Directed.by_signature(): no entry for signature {!r}; "
                "known signatures: {}".format(
                    signature, list(self._mapping.keys())))
        return self._mapping[signature]

    def to_mixture(self, signature) -> "_channels.Mixture":
        entry = self._entry_for(signature)
        if isinstance(entry, Directed):
            return entry.to_mixture(signature)
        return entry

    def compile(self, signature, *, detector=None, models=None):
        return self.to_mixture(signature).compile(
            signature, detector=detector, models=models)

    def validate(self, signature=None):
        if signature is not None:
            return self.to_mixture(signature).validate(signature)
        for sig, entry in self._mapping.items():
            mixture = entry.to_mixture(sig) if isinstance(entry, Directed) else entry
            mixture.validate(sig)
        return None

    def describe(self) -> str:
        parts = []
        for sig, entry in self._mapping.items():
            mixture = entry.to_mixture(sig) if isinstance(entry, Directed) else entry
            parts.append("{!r}: {}".format(sig, mixture.describe()))
        return "; ".join(parts)
