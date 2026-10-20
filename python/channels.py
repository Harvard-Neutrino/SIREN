"""Channel algebra compiling to siren.injection.MultiChannelPhaseSpace.

Public names: Channel, Mixture, isotropic, physical.

A Channel wraps one engine PhaseSpaceChannel factory that is resolved
against a concrete InteractionSignature at compile time (daughters and
spectators are named by particle, not by index).  Channels combine with
`*` (rescale) and `+` (mix) into a Mixture, whose weights are normalized
to sum to 1 at every construction.  Mixture.compile() resolves every
Channel's factory against a signature and builds the validated
MultiChannelPhaseSpace; Mixture.validate() runs the same compilation plus
a detector-free ConvertDensity probe against a template record, so a bad
an inconvertible measure combination fails
at configuration time rather than at first weight evaluation.
"""

from __future__ import annotations

import math
from typing import Callable, List, Optional, Sequence, Tuple, Union

from . import particles as _particles
from .errors import ConfigurationError, MeasureCompatibilityError

__all__ = ["Channel", "Mixture", "isotropic", "physical"]


def _siren():
    import siren
    return siren


def _resolve_index(signature, name_or_index, role):
    """Resolve a daughter/spectator spec to a secondary index.

    `name_or_index` may be an int (used directly) or a particle name/
    ParticleType resolved against `signature.secondary_types`.  Raises
    ConfigurationError naming the particle and the signature on an
    absent or ambiguous match; never silently skips.
    """
    if isinstance(name_or_index, int):
        return name_or_index
    try:
        ptype = _particles.resolve(name_or_index)
    except (ValueError, TypeError) as e:
        raise ConfigurationError(
            "{}: {}".format(role, e))
    secondaries = list(signature.secondary_types)
    matches = [i for i, t in enumerate(secondaries) if t == ptype]
    if not matches:
        raise ConfigurationError(
            "{}: particle {!r} not found among secondary_types {} of "
            "signature {}".format(
                role, name_or_index,
                [str(t) for t in secondaries], signature))
    if len(matches) > 1:
        raise ConfigurationError(
            "{}: particle {!r} is ambiguous among secondary_types {} of "
            "signature {} (matched indices {})".format(
                role, name_or_index,
                [str(t) for t in secondaries], signature, matches))
    return matches[0]


# ------------------------------------------------------------------ #
#  Shared mixture-node algebra                                        #
# ------------------------------------------------------------------ #

class _MixtureNode:
    """Common `*`/`+` algebra for Channel nodes.

    A node wraps a deferred factory(signature) -> engine PhaseSpaceChannel
    plus a weight; `_build` resolves it against a concrete signature only
    at compile/validate time.  Subclasses set `_factory`, `weight`,
    `label` and implement `_copy_with_weight` so `*` returns the same
    subclass.
    """

    _factory: Callable
    weight: float
    label: str

    def _build(self, signature, *, detector=None, models=None):
        return self._factory(signature, detector=detector, models=models)

    def _copy_with_weight(self, weight):
        raise NotImplementedError

    def __rmul__(self, scalar):
        if not isinstance(scalar, (int, float)):
            return NotImplemented
        return self._copy_with_weight(self.weight * float(scalar))

    def __mul__(self, scalar):
        return self.__rmul__(scalar)

    def __add__(self, other):
        if isinstance(other, _MixtureNode):
            return Mixture([(self.weight, self), (other.weight, other)])
        if isinstance(other, Mixture):
            return Mixture([(self.weight, self)] + other._scaled_entries())
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, Mixture):
            return Mixture(other._scaled_entries() + [(self.weight, self)])
        return NotImplemented

    def describe(self) -> str:
        return self.label


class Channel(_MixtureNode):
    """One weighted spec node wrapping a deferred engine-channel factory.

    `factory(signature)` builds the concrete engine PhaseSpaceChannel for
    that signature; it is called only at compile/validate time, never at
    construction, so a Channel can be built before any signature is known.
    """

    def __init__(self, factory: Callable, weight: float = 1.0,
                 label: str = "channel"):
        self._factory = factory
        self.weight = float(weight)
        self.label = label

    def _copy_with_weight(self, weight):
        return Channel(self._factory, weight, self.label)


# ------------------------------------------------------------------ #
#  Mixture                                                             #
# ------------------------------------------------------------------ #

class Mixture:
    """A normalized weighted list of Channels.

    Always normalized so the stored weights sum to 1: at construction and
    after every `+`.  Rejects an empty, negative, or non-finite weight set
    with ConfigurationError.
    """

    def __init__(self, entries: Sequence[Tuple[float, "Channel"]]):
        entries = list(entries)
        if not entries:
            raise ConfigurationError("Mixture: cannot construct from an empty list of channels")
        for w, ch in entries:
            if not math.isfinite(w):
                raise ConfigurationError(
                    "Mixture: non-finite weight {!r} for channel {!r}".format(w, ch.describe()))
            if w < 0.0:
                raise ConfigurationError(
                    "Mixture: negative weight {!r} for channel {!r}".format(w, ch.describe()))
        total = sum(w for w, _ in entries)
        if not (math.isfinite(total) and total > 0.0):
            raise ConfigurationError(
                "Mixture: weights sum to {!r}, must be finite and positive".format(total))
        self._entries: List[Tuple[float, "Channel"]] = [
            (w / total, ch) for w, ch in entries]
        self._scale = total

    def _scaled_entries(self):
        """The entries as (scale * relative_weight, channel) pairs."""
        return [(self._scale * w, ch) for w, ch in self._entries]

    def __add__(self, other):
        if isinstance(other, _MixtureNode):
            return Mixture(self._scaled_entries() + [(other.weight, other)])
        if isinstance(other, Mixture):
            return Mixture(self._scaled_entries() + other._scaled_entries())
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, _MixtureNode):
            return Mixture([(other.weight, other)] + self._scaled_entries())
        return NotImplemented

    def __rmul__(self, scalar):
        if not isinstance(scalar, (int, float)):
            return NotImplemented
        # Compose the scale multiplicatively so nested rescaling accumulates:
        # `p * (q * mixture)` carries p*q, matching _MixtureNode.__rmul__.
        scaled = Mixture(list(self._entries))
        scaled._scale = self._scale * float(scalar)
        return scaled

    __mul__ = __rmul__

    def describe(self) -> str:
        parts = ["{:.2f} * {}".format(w, ch.describe()) for w, ch in self._entries]
        return " + ".join(parts)

    # -------------------------------------------------------------- #
    #  Compilation                                                    #
    # -------------------------------------------------------------- #

    def _flatten(self, signature, *, detector=None, models=None):
        """Resolve every entry to (weight, engine PhaseSpaceChannel)."""
        engine_channels = []
        weights = []
        for w, ch in self._entries:
            built = ch._build(signature, detector=detector, models=models)
            engine_channels.append(built)
            weights.append(w)
        return engine_channels, weights

    def compile(self, signature, *, detector=None, models=None):
        """Build the siren.injection.MultiChannelPhaseSpace for `signature`.

        Uses the validating constructor: it normalizes the weights and
        throws MeasureCompatibilityError on a fatal topology/measure
        incompatibility.
        """
        siren = _siren()
        engine_channels, weights = self._flatten(
            signature, detector=detector, models=models)
        return siren.injection.MultiChannelPhaseSpace(engine_channels, weights)

    def validate(self, signature=None, *, models=None):
        """Data-free configuration check.

        Runs (a) compilation against a template signature via the
        validating MultiChannelPhaseSpace ctor, which normalizes weights
        and raises MeasureCompatibilityError on a fatal incompatibility,
        and (b) a ConvertDensity probe on a template record built from the
        signature's masses, surfacing incompatible measures at
        configuration time.  Does not run the detector-dependent
        ValidateChannelDensities probe (that runs in Injector._build).

        `models=` is forwarded to compile() exactly as compile() accepts
        it; required whenever the mixture holds a physical() channel (see
        physical()'s docstring).
        """
        siren = _siren()
        sig = signature if signature is not None else _template_signature()
        mcps = self.compile(sig, models=models)  # raises MeasureCompatibilityError if fatal

        record = _template_record(sig)
        topology = mcps.CommonTopology()
        common_measure = mcps.CommonMeasure()
        for ch in mcps.channels:
            measure = ch.Measure()
            if measure == common_measure:
                continue
            try:
                siren.injection.ConvertDensity(
                    1.0, measure, common_measure, topology, record)
            except MeasureCompatibilityError:
                # SolidAngleLab conversions outside Decay2Body topology are
                # unconvertible by design; any other channel combination that
                # cannot convert is a genuine configuration error.
                if topology != siren.injection.PhaseSpaceTopology.Decay2Body and (
                        measure.type == siren.injection.PhaseSpaceMeasureType.SolidAngleLab
                        or common_measure.type == siren.injection.PhaseSpaceMeasureType.SolidAngleLab):
                    continue
                raise
        return mcps


# ------------------------------------------------------------------ #
#  Template signature / record for detector-free validation           #
# ------------------------------------------------------------------ #

def _template_signature():
    siren = _siren()
    sig = siren.dataclasses.InteractionSignature()
    sig.primary_type = siren.particles.N4 if hasattr(siren.particles, "N4") else siren.particles.NuMu
    sig.target_type = siren.dataclasses.ParticleType.Decay
    sig.secondary_types = [siren.particles.NuLight, siren.particles.Gamma]
    return sig


def _template_record(signature, primary_mass=0.02, energy=0.05):
    """A minimal InteractionRecord for the given signature.

    Mirrors the idiom in tests/python/test_density_breakdown.py: plausible
    small masses, a forward-moving primary, and a vertex at the origin.
    Secondary masses are set to a small common value (no physical model is
    consulted -- this record only feeds the detector-free ConvertDensity
    probe, not any Sample()).
    """
    siren = _siren()
    n_sec = len(signature.secondary_types)
    m_sec = 0.001
    pz = math.sqrt(max(energy * energy - primary_mass * primary_mass, 0.0))
    rec = siren.dataclasses.InteractionRecord()
    rec.signature.primary_type = signature.primary_type
    rec.signature.target_type = signature.target_type
    rec.signature.secondary_types = list(signature.secondary_types)
    rec.primary_mass = primary_mass
    rec.primary_momentum = [energy, 0.0, 0.0, pz]
    rec.secondary_masses = [m_sec] * n_sec
    rec.secondary_momenta = [[0.0, 0.0, 0.0, 0.0]] * n_sec
    rec.secondary_helicities = [0] * n_sec
    rec.interaction_vertex = [0.0, 0.0, 0.0]
    rec.primary_initial_position = [0.0, 0.0, 0.0]
    return rec


# ------------------------------------------------------------------ #
#  PairMass                                                            #
# ------------------------------------------------------------------ #


# ------------------------------------------------------------------ #
#  Factories                                                           #
# ------------------------------------------------------------------ #

def isotropic(daughter: Union[int, str] = 0) -> Channel:
    """Isotropic2BodyChannel wrapped as a Channel.

    `daughter` is either a secondary index or a particle name/ParticleType
    resolved against the signature at compile time.
    """
    def factory(signature, *, detector=None, models=None):
        siren = _siren()
        idx = _resolve_index(signature, daughter, "isotropic(daughter=...)")
        return siren.injection.Isotropic2BodyChannel(idx)

    label = "isotropic" if daughter == 0 else "isotropic({!r})".format(daughter)
    return Channel(factory, weight=1.0, label=label)


def physical() -> Channel:
    """Late-binds PhysicalDecayChannel or PhysicalCrossSectionChannel.

    Resolution happens at compile time and needs the vertex's interaction
    model(s), passed as `models=[...]` to Mixture.compile()/validate().
    Exactly one model matching the signature (by GetPossibleSignatures, if
    available, else taken as-is) must be supplied; otherwise raises
    ConfigurationError with a one-line fix.
    """
    def factory(signature, *, detector=None, models=None):
        siren = _siren()
        if not models:
            raise ConfigurationError(
                "physical(): no interaction model available to late-bind; "
                "pass models=[<Decay or CrossSection>] to Mixture.compile()/"
                "validate()")
        candidates = list(models)
        if len(candidates) > 1:
            raise ConfigurationError(
                "physical(): {} candidate models supplied, expected exactly "
                "one; pass models=[<the single Decay or CrossSection for "
                "this vertex>]".format(len(candidates)))
        model = candidates[0]
        if isinstance(model, siren.interactions.Decay):
            return siren.injection.PhysicalDecayChannel(model, signature)
        if isinstance(model, siren.interactions.CrossSection):
            return siren.injection.PhysicalCrossSectionChannel(model, signature)
        raise ConfigurationError(
            "physical(): model {!r} is neither a Decay nor a CrossSection".format(model))

    return Channel(factory, weight=1.0, label="physical")
