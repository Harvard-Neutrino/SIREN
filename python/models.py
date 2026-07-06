"""Authoring bases for Python-implemented interaction models.

``DecayModel`` and ``CrossSectionModel`` collapse the pyDecay/pyCrossSection
pure-virtual sets to a handful of physics hooks. A subclass declares its
particles and measure once; the base derives the signature methods, the width
or cross-section overload pairs, ``FinalStateProbability`` as differential over
total, ``DensityVariables`` (normalized to a list of strings), and ``Topology``
from the final-state arity.

The default ``SampleFinalState`` samples the DECLARED measure through a
self-contained engine channel, so the Sample==Density closure holds by
construction whenever ``sample()`` is not overridden. It NEVER delegates to
``PhysicalDecayChannel``/``PhysicalCrossSectionChannel``: those channels call
the wrapped model's ``SampleFinalState`` (PhysicalChannelAdapters.cxx), which
dispatches back through the trampoline into this default and recurses forever.
Only measures with a self-contained channel are supported by the default; any
other declared measure requires an explicit ``sample()`` override, enforced at
model registration by ``_validation.audit_overrides``.

``decay_model_base(base=...)`` / ``cross_section_model_base(base=...)`` produce
these bases over a chosen pybind C++ base, so the DarkNews port can build on
``siren.interactions.DarkNewsDecay`` and keep its trampoline while still gaining
the authoring surface.
"""

import difflib

from . import interactions as _interactions
from . import dataclasses as _dataclasses
from . import injection as _injection
from . import particles as _particles
from .errors import ConfigurationError

Particle = _dataclasses.Particle
InteractionSignature = _dataclasses.InteractionSignature
InteractionRecord = _dataclasses.InteractionRecord

PhaseSpaceTopology = _injection.PhaseSpaceTopology
PhaseSpaceMeasure = _injection.PhaseSpaceMeasure
PhaseSpaceMeasureType = _injection.PhaseSpaceMeasureType

# Re-exported so authors spell measures as siren.Measure.SolidAngleRest() and
# topologies as siren.Topology.Decay2Body without reaching into siren.injection.
Measure = PhaseSpaceMeasure
Topology = PhaseSpaceTopology


# The pyDecay / pyCrossSection pure virtuals a subclass may need to supply. A
# name resolving only to the abstract pybind base (no Python override) is the
# silent pure-virtual-abort case audit_overrides turns into a loud error.
_DECAY_VIRTUALS = (
    "TotalDecayWidth",
    "TotalDecayWidthAllFinalStates",
    "DifferentialDecayWidth",
    "SampleFinalState",
    "GetPossibleSignatures",
    "GetPossibleSignaturesFromParent",
    "FinalStateProbability",
    "DensityVariables",
)

_CROSS_SECTION_VIRTUALS = (
    "TotalCrossSection",
    "DifferentialCrossSection",
    "InteractionThreshold",
    "SampleFinalState",
    "GetPossibleTargets",
    "GetPossibleTargetsFromPrimary",
    "GetPossiblePrimaries",
    "GetPossibleSignatures",
    "GetPossibleSignaturesFromParents",
    "FinalStateProbability",
    "DensityVariables",
)


def _resolve_type(name_or_type):
    """Resolve a particle name or ParticleType to a ParticleType."""
    if isinstance(name_or_type, Particle.ParticleType):
        return name_or_type
    return _particles.resolve(name_or_type)


def _normalize_density_variables(value):
    """Coerce a DensityVariables result to a list of strings.

    Accepts a list/tuple of names or a single comma-or-space-separated string.
    """
    if value is None:
        return []
    if isinstance(value, str):
        parts = value.replace(",", " ").split()
        return list(parts)
    return [str(v) for v in value]


def _topology_from_arity(n_finals):
    """Map a final-state count to a decay topology."""
    if n_finals == 2:
        return PhaseSpaceTopology.Decay2Body
    if n_finals == 3:
        return PhaseSpaceTopology.Decay3Body
    return PhaseSpaceTopology.DecayNBody


# A declared measure maps to a self-contained engine channel (one that writes an
# InteractionRecord directly and does not call back into SampleFinalState).
# SolidAngleRest 2-body -> Isotropic2BodyChannel, whose DetectorModel argument is
# unused (None is safe). Anything absent here has no closure-by-construction
# default; the author must override sample().
def _self_contained_channel(measure, n_finals, daughter_index):
    if n_finals == 2 and measure.type == PhaseSpaceMeasureType.SolidAngleRest:
        return _injection.Isotropic2BodyChannel(daughter_index)
    return None


class _ModelSubclassCheck:
    """Shared __init_subclass__ near-miss override rejection.

    difflib the names a subclass defines against the pybind virtual set of the
    chosen base; a close-but-wrong name (e.g. total_witdh, DifferentialXS) is a
    typo that would silently leave the real hook unimplemented, so it is
    rejected at class definition.
    """

    _virtual_names = ()
    _hook_names = ()

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        valid = set(cls._virtual_names) | set(cls._hook_names)
        own = {
            name for name, value in vars(cls).items()
            if callable(value) and not name.startswith("_")
        }
        for name in own:
            if name in valid:
                continue
            close = difflib.get_close_matches(name, valid, n=1, cutoff=0.8)
            if close:
                raise ConfigurationError(
                    "%s.%s looks like a typo of the model hook '%s'; "
                    "rename it or remove it" % (cls.__name__, name, close[0]))


def decay_model_base(base=None):
    """Build a DecayModel authoring base over a pybind Decay C++ base.

    ``base`` defaults to ``siren.interactions.Decay``. Pass a trampoline-derived
    base such as ``siren.interactions.DarkNewsDecay`` to keep that trampoline
    while gaining the authoring surface.
    """
    if base is None:
        base = _interactions.Decay

    class DecayModel(base, _ModelSubclassCheck):
        """Authoring base for a Python decay model.

        Declare ``parent`` and ``daughters`` (names or ParticleTypes) and a
        ``measure``; supply ``total_width()`` and ``differential_width(record)``.
        Override ``sample(record, random)`` only for a measure without a
        self-contained default; overriding it is the trigger to run
        ``siren.check_closure``.
        """

        parent = None
        daughters = ()
        measure = None
        daughter_index = 0

        _virtual_names = _DECAY_VIRTUALS
        _hook_names = (
            "total_width", "differential_width", "density_variables", "sample",
        )

        def __init__(self, *args, **kwargs):
            base.__init__(self)

        # ---- declared metadata -> signature methods ----

        def _signature(self):
            sig = InteractionSignature()
            sig.primary_type = _resolve_type(self.parent)
            sig.target_type = Particle.ParticleType.Decay
            sig.secondary_types = [_resolve_type(d) for d in self.daughters]
            return sig

        def GetPossibleSignatures(self):
            return [self._signature()]

        def GetPossibleSignaturesFromParent(self, primary_type):
            if _resolve_type(self.parent) == primary_type:
                return [self._signature()]
            return []

        def Topology(self):
            return _topology_from_arity(len(self.daughters))

        def Measure(self):
            if self.measure is None:
                return PhaseSpaceMeasure.Unspecified()
            return self.measure

        def DensityVariables(self):
            return _normalize_density_variables(self.density_variables())

        # ---- width overload pair from the single physics hook ----

        def _resolved_total_width(self):
            return float(self.total_width())

        def TotalDecayWidth(self, arg):
            if isinstance(arg, InteractionRecord):
                sig = self._signature()
                rsig = arg.signature
                if (rsig.primary_type != sig.primary_type
                        or rsig.target_type != sig.target_type
                        or list(rsig.secondary_types) != list(sig.secondary_types)):
                    return 0.0
                return self._resolved_total_width()
            if _resolve_type(self.parent) != arg:
                return 0.0
            return self._resolved_total_width()

        def TotalDecayWidthAllFinalStates(self, arg):
            if isinstance(arg, InteractionRecord):
                primary = arg.signature.primary_type
            else:
                primary = arg
            if _resolve_type(self.parent) != primary:
                return 0.0
            return self._resolved_total_width()

        def DifferentialDecayWidth(self, record):
            return float(self.differential_width(record))

        def FinalStateProbability(self, record):
            total = self._resolved_total_width()
            if total <= 0.0:
                return 0.0
            return float(self.differential_width(record)) / total

        # ---- closure-by-construction default sampler ----

        def SampleFinalState(self, record, random):
            # The engine calls this; it dispatches to sample() so an override
            # takes effect. The default sample() draws the declared measure.
            self.sample(record, random)

        def _sample_via_channel(self, channel, csdr, random):
            """Bridge a CrossSectionDistributionRecord through an engine channel.

            The channel writes secondary momenta into an InteractionRecord; copy
            them back onto the CSDR's secondary particle records. The DetectorModel
            argument is unused by self-contained channels, so None is passed.
            """
            ir = csdr.record
            channel.Sample(random, None, ir)
            secondaries = csdr.get_secondary_particle_records()
            for i, spr in enumerate(secondaries):
                spr.four_momentum = ir.secondary_momenta[i]
                spr.mass = ir.secondary_masses[i]

        # ---- physics hooks (subclass supplies these) ----

        def total_width(self):
            raise NotImplementedError(
                "%s must implement total_width()" % type(self).__name__)

        def differential_width(self, record):
            raise NotImplementedError(
                "%s must implement differential_width(record)"
                % type(self).__name__)

        def density_variables(self):
            """Names of the sampled kinematic variables (list or comma string)."""
            return []

        def sample(self, record, random):
            """Sample the declared measure via a self-contained engine channel.

            Override for a measure with no self-contained channel. This never
            delegates to PhysicalDecayChannel, whose Sample dispatches back into
            SampleFinalState and would recurse.
            """
            channel = _self_contained_channel(
                self.Measure(), len(self.daughters), self.daughter_index)
            if channel is None:
                raise ConfigurationError(
                    "%s declares measure %r with no self-contained sampling "
                    "channel; override sample(record, random)"
                    % (type(self).__name__, self.Measure()))
            self._sample_via_channel(channel, record, random)

    return DecayModel


def cross_section_model_base(base=None):
    """Build a CrossSectionModel authoring base over a pybind CrossSection base.

    ``base`` defaults to ``siren.interactions.CrossSection``. Pass a
    trampoline-derived base such as ``siren.interactions.DarkNewsCrossSection``
    to keep that trampoline.
    """
    if base is None:
        base = _interactions.CrossSection

    class CrossSectionModel(base, _ModelSubclassCheck):
        """Authoring base for a Python cross-section model.

        Declare ``primary``, ``target`` and ``finals`` (names or ParticleTypes)
        and a ``measure``; supply ``total_xs(record)`` and
        ``differential_xs(record)``. Override ``sample(record, random)`` only for
        a measure without a self-contained default.
        """

        primary = None
        target = None
        finals = ()
        measure = None
        threshold = 0.0

        _virtual_names = _CROSS_SECTION_VIRTUALS
        _hook_names = (
            "total_xs", "differential_xs", "density_variables", "sample",
        )

        def __init__(self, *args, **kwargs):
            base.__init__(self)

        # ---- declared metadata -> target/signature methods ----

        def _signature(self):
            sig = InteractionSignature()
            sig.primary_type = _resolve_type(self.primary)
            sig.target_type = _resolve_type(self.target)
            sig.secondary_types = [_resolve_type(f) for f in self.finals]
            return sig

        def GetPossiblePrimaries(self):
            return [_resolve_type(self.primary)]

        def GetPossibleTargets(self):
            return [_resolve_type(self.target)]

        def GetPossibleTargetsFromPrimary(self, primary_type):
            if _resolve_type(self.primary) == primary_type:
                return [_resolve_type(self.target)]
            return []

        def GetPossibleSignatures(self):
            return [self._signature()]

        def GetPossibleSignaturesFromParents(self, primary_type, target_type):
            if (_resolve_type(self.primary) == primary_type
                    and _resolve_type(self.target) == target_type):
                return [self._signature()]
            return []

        def Topology(self):
            return (PhaseSpaceTopology.Scatter2to2 if len(self.finals) == 2
                    else PhaseSpaceTopology.Scatter2to3)

        def Measure(self):
            if self.measure is None:
                return PhaseSpaceMeasure.Unspecified()
            return self.measure

        def DensityVariables(self):
            return _normalize_density_variables(self.density_variables())

        def InteractionThreshold(self, record):
            return float(self.threshold)

        # ---- cross-section hooks ----

        def TotalCrossSection(self, record):
            return float(self.total_xs(record))

        def DifferentialCrossSection(self, record):
            return float(self.differential_xs(record))

        def FinalStateProbability(self, record):
            total = float(self.total_xs(record))
            if total <= 0.0:
                return 0.0
            return float(self.differential_xs(record)) / total

        # ---- default sampler (same recursion-safe contract as decays) ----

        def SampleFinalState(self, record, random):
            self.sample(record, random)

        # ---- physics hooks ----

        def total_xs(self, record):
            raise NotImplementedError(
                "%s must implement total_xs(record)" % type(self).__name__)

        def differential_xs(self, record):
            raise NotImplementedError(
                "%s must implement differential_xs(record)" % type(self).__name__)

        def density_variables(self):
            return []

        def sample(self, record, random):
            """Sample the declared measure via a self-contained engine channel.

            Override for a measure with no self-contained channel; never
            delegates to PhysicalCrossSectionChannel, which would recurse.
            """
            channel = _self_contained_channel(
                self.Measure(), len(self.finals), 0)
            if channel is None:
                raise ConfigurationError(
                    "%s declares measure %r with no self-contained sampling "
                    "channel; override sample(record, random)"
                    % (type(self).__name__, self.Measure()))
            ir = record.record
            channel.Sample(random, None, ir)
            secondaries = record.get_secondary_particle_records()
            for i, spr in enumerate(secondaries):
                spr.four_momentum = ir.secondary_momenta[i]
                spr.mass = ir.secondary_masses[i]

    return CrossSectionModel


DecayModel = decay_model_base()
CrossSectionModel = cross_section_model_base()
