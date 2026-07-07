"""Channel algebra compiling to siren.injection.MultiChannelPhaseSpace.

Public names: Channel, Mixture, PairMass, Tiling, isotropic, toward,
toward_3body, scatter_toward, physical, tile.

A Channel wraps one engine PhaseSpaceChannel factory that is resolved
against a concrete InteractionSignature at compile time (daughters and
spectators are named by particle, not by index).  Channels combine with
`*` (rescale) and `+` (mix) into a Mixture, whose weights are normalized
to sum to 1 at every construction.  Mixture.compile() resolves every
Channel's factory against a signature and builds the validated
MultiChannelPhaseSpace; Mixture.validate() runs the same compilation plus
a detector-free ConvertDensity probe against a template record, so a bad
3-body factorization index or an inconvertible measure combination fails
at configuration time rather than at first weight evaluation.
"""

from __future__ import annotations

import math
from typing import Callable, List, Optional, Sequence, Tuple, Union

from . import directed_tiling as _tiling
from . import particles as _particles
from .errors import ConfigurationError, MeasureCompatibilityError

__all__ = [
    "Channel",
    "Mixture",
    "PairMass",
    "Tiling",
    "isotropic",
    "toward",
    "toward_3body",
    "scatter_toward",
    "physical",
    "tile",
]


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
    """Common `*`/`+` algebra for Channel and Tiling.

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
    """A normalized weighted list of Channels (and Tilings).

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
        # A composition scale carried by `*`: `p * mixture` sets it so that
        # `p * mix_a + q * mix_b` gives each sub-channel its p (or q) share of
        # the combined mixture. It never affects a standalone mixture's
        # normalized `_entries`.
        self._scale = 1.0

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

    def validate(self, signature=None):
        """Data-free configuration check.

        Runs (a) compilation against a template signature via the
        validating MultiChannelPhaseSpace ctor, which normalizes weights
        and raises MeasureCompatibilityError on a fatal incompatibility,
        and (b) a ConvertDensity probe on a template record built from the
        signature's masses, surfacing bad 3-body factorization indices at
        configuration time.  Does not run the detector-dependent
        ValidateChannelDensities probe (that runs in Injector._build).
        """
        siren = _siren()
        sig = signature if signature is not None else _template_signature()
        mcps = self.compile(sig)  # raises MeasureCompatibilityError if fatal

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

class PairMass:
    """Invariant-mass sampling law for a 3-body directed channel's pair.

    Carries the fields toward_3body maps onto the
    DetectorDirected3BodyChannel ctor: mass_mode, resonance_mass,
    resonance_width, power_law_nu, power_law_offset, mass_cdf_nodes,
    mass_cdf_values.
    """

    def __init__(self, mass_mode, resonance_mass=0.0, resonance_width=0.0,
                 power_law_nu=0.8, power_law_offset=0.0,
                 mass_cdf_nodes=(), mass_cdf_values=()):
        self.mass_mode = mass_mode
        self.resonance_mass = resonance_mass
        self.resonance_width = resonance_width
        self.power_law_nu = power_law_nu
        self.power_law_offset = power_law_offset
        self.mass_cdf_nodes = list(mass_cdf_nodes)
        self.mass_cdf_values = list(mass_cdf_values)

    @classmethod
    def uniform(cls):
        siren = _siren()
        return cls(siren.injection.InvariantMassMode.Uniform)

    @classmethod
    def breit_wigner(cls, m, w):
        siren = _siren()
        return cls(siren.injection.InvariantMassMode.BreitWigner,
                    resonance_mass=m, resonance_width=w)

    @classmethod
    def power_law(cls, nu, offset=0.0):
        siren = _siren()
        return cls(siren.injection.InvariantMassMode.PowerLaw,
                    power_law_nu=nu, power_law_offset=offset)

    @classmethod
    def tabulated(cls, nodes, values):
        siren = _siren()
        return cls(siren.injection.InvariantMassMode.Tabulated,
                    mass_cdf_nodes=nodes, mass_cdf_values=values)


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


def toward(daughter: Union[int, str], target, *, mode=None, volume=-1.0,
           fraction: Optional[float] = None) -> Union[Channel, Mixture]:
    """DetectorDirected2BodyChannel toward `target`, daughter named by particle.

    If `fraction` is given and < 1, returns the mixture
    `fraction*toward(...) + (1-fraction)*isotropic(...)` (a physical
    fallback for events that miss the target); otherwise a bare directed
    Channel.
    """
    siren = _siren()
    if mode is None:
        mode = siren.injection.DirectedMode.Volume

    def factory(signature, *, detector=None, models=None):
        idx = _resolve_index(signature, daughter, "toward(daughter=...)")
        return siren.injection.DetectorDirected2BodyChannel(
            target, idx, mode, volume)

    label = "toward({!r})".format(daughter)
    directed = Channel(factory, weight=1.0, label=label)

    if fraction is None or fraction >= 1.0:
        return directed
    if not (0.0 <= fraction < 1.0):
        raise ConfigurationError(
            "toward(fraction={!r}): fraction must be in [0, 1)".format(fraction))
    return fraction * directed + (1.0 - fraction) * isotropic(daughter)


def toward_3body(directed: Union[int, str], target, *,
                  spectator: Optional[Union[int, str]] = None,
                  strategy: str = "recursive",
                  pair_mass: Optional[PairMass] = None,
                  mode=None) -> Channel:
    """DetectorDirected3BodyChannel toward `target`.

    strategy='direct' uses the direct ctor (directed_index only, other two
    indices inferred ascending). strategy='recursive' uses the recursive
    ctor (spectator/pair indices resolved per signature; `spectator`
    selects the spectator daughter, the other two secondaries become the
    pair, with `directed` as the directed member of the pair).

    The channel topology is inferred per signature: a real target type
    tags the channel Scatter2to3 (2->3 scattering), the Decay marker tags
    it Decay3Body.
    """
    siren = _siren()
    if mode is None:
        mode = siren.injection.DirectedMode.Volume
    if pair_mass is None:
        pair_mass = PairMass.uniform()
    if strategy not in ("recursive", "direct"):
        raise ConfigurationError(
            "toward_3body(strategy={!r}): must be 'recursive' or 'direct'".format(strategy))

    def factory(signature, *, detector=None, models=None):
        directed_idx = _resolve_index(signature, directed, "toward_3body(directed=...)")
        decay_marker = siren.dataclasses.Particle.ParticleType.Decay
        topology = (siren.injection.PhaseSpaceTopology.Decay3Body
                    if signature.target_type == decay_marker
                    else siren.injection.PhaseSpaceTopology.Scatter2to3)
        common_kwargs = dict(
            topology=topology,
            mass_mode=pair_mass.mass_mode,
            resonance_mass=pair_mass.resonance_mass,
            resonance_width=pair_mass.resonance_width,
            power_law_nu=pair_mass.power_law_nu,
            power_law_offset=pair_mass.power_law_offset,
            mode=mode,
            mass_cdf_nodes=pair_mass.mass_cdf_nodes,
            mass_cdf_values=pair_mass.mass_cdf_values,
        )
        if strategy == "direct":
            return siren.injection.DetectorDirected3BodyChannel(
                factorization=siren.injection.ThreeBodyMode.Direct,
                target=target, directed_index=directed_idx, **common_kwargs)

        n_sec = len(signature.secondary_types)
        if spectator is None:
            raise ConfigurationError(
                "toward_3body(strategy='recursive'): spectator= is required "
                "to resolve the pair indices")
        spectator_idx = _resolve_index(signature, spectator, "toward_3body(spectator=...)")
        others = [i for i in range(n_sec) if i != spectator_idx]
        if directed_idx not in others:
            raise ConfigurationError(
                "toward_3body: directed index {} coincides with spectator "
                "index {} for signature {}".format(directed_idx, spectator_idx, signature))
        pair_first_idx, pair_second_idx = others[0], others[1]
        return siren.injection.DetectorDirected3BodyChannel(
            factorization=siren.injection.ThreeBodyMode.Recursive,
            target=target, spectator_index=spectator_idx,
            pair_first_index=pair_first_idx, pair_second_index=pair_second_idx,
            directed_index=directed_idx, **common_kwargs)

    label = "toward_3body({!r}, strategy={!r})".format(directed, strategy)
    return Channel(factory, weight=1.0, label=label)


def scatter_toward(target, *, variable=None, q2_mode=None,
                    mediator_mass=0.0, mode=None) -> Channel:
    """DetectorDirectedScatteringChannel toward `target`."""
    siren = _siren()
    if mode is None:
        mode = siren.injection.DirectedMode.Volume
    if variable is None:
        variable = siren.injection.ScatteringVariable.Q2
    if q2_mode is None:
        q2_mode = siren.injection.ScatteringQ2Mode.Geometry

    def factory(signature, *, detector=None, models=None):
        return siren.injection.DetectorDirectedScatteringChannel(
            target, 0, variable, mode, q2_mode, mediator_mass)

    return Channel(factory, weight=1.0, label="scatter_toward")


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


# ------------------------------------------------------------------ #
#  Tiling                                                              #
# ------------------------------------------------------------------ #

class Tiling(_MixtureNode):
    """A grid/disjoint tiling of a geometry, treated as one Channel-like unit.

    Wraps directed_tiling.build_grid_tiling/disjointify; the optimizer
    tunes the tiling's inner weights as one nested level (see
    NestedMixtureChannel), never touching the outer mixture's weight for
    it directly.  Participates in Mixture algebra the same way a Channel
    does.
    """

    def __init__(self, factory: Callable, weight: float = 1.0,
                 label: str = "tile"):
        self._factory = factory
        self.weight = float(weight)
        self.label = label

    def _copy_with_weight(self, weight):
        return Tiling(self._factory, weight, self.label)


def tile(geometry, grid: int, method: str = "disjoint",
         daughter: Union[int, str] = 0, mode=None) -> Tiling:
    """Tile `geometry` into a grid of directed sub-channels.

    `grid` is the number of cells per axis (an axis-aligned Box region is
    required, matching directed_tiling.grid_cells).  method='disjoint'
    routes through directed_tiling.disjointify on the grid cells (already
    disjoint by construction, so disjointify is a no-op composition step
    kept for a uniform call path); any other method raises
    ConfigurationError.
    """
    if method != "disjoint":
        raise ConfigurationError(
            "tile(method={!r}): only 'disjoint' is supported".format(method))

    def factory(signature, *, detector=None, models=None):
        siren = _siren()
        if mode is None:
            resolved_mode = siren.injection.DirectedMode.Volume
        else:
            resolved_mode = mode
        idx = _resolve_index(signature, daughter, "tile(daughter=...)")
        cells = _tiling.grid_cells(geometry, grid)
        disjoint_cells = _tiling.disjointify(cells)
        chan_factory = _tiling.default_2body_factory(idx, resolved_mode)
        return _tiling._wrap(disjoint_cells, chan_factory, True, 100000)

    return Tiling(factory, weight=1.0, label="tile(grid={})".format(grid))
