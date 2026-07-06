"""
Input validation helpers for SIREN's Python interface.

Uses the C++ distribution methods ``SetVariables()``,
``RequiredVariables()``, and ``DensityVariables()`` to validate
distribution lists.  No hardcoded type lists -- adding a new
distribution in C++ with the correct overrides makes it work here
automatically.
"""

from . import distributions as _d
from . import particles as _particles
from .errors import ConfigurationError

DV = _d.DistributionVariable

_REQUIRED_INJECTION_VARIABLES = {
    DV.PrimaryEnergy,
    DV.PrimaryDirection,
}

_REQUIRED_INJECTION_VERTEX = {
    DV.InitialPosition,
    DV.InteractionVertex,
}

_REQUIRED_PHYSICAL_VARIABLES = {
    DV.PrimaryEnergy,
    DV.PrimaryDirection,
}


def collect_set_variables(distributions):
    """Return the union of all SetVariables across a list of distributions."""
    result = set()
    for d in distributions:
        if hasattr(d, "SetVariables"):
            result |= d.SetVariables()
    return result


def collect_density_variables(distributions):
    """Return the union of all DensityVariables across a list of distributions."""
    result = set()
    for d in distributions:
        result.update(d.DensityVariables())
    return result


def validate_ordering(distributions):
    """Check that each distribution's RequiredVariables are satisfied
    by the SetVariables of preceding distributions.

    Parameters
    ----------
    distributions : list
        Ordered list of PrimaryInjectionDistribution objects.

    Raises
    ------
    ValueError
        If a distribution requires variables that haven't been set yet,
        with a message naming the distribution and the missing variables.
    """
    available = set()
    for i, dist in enumerate(distributions):
        if not hasattr(dist, "RequiredVariables"):
            continue
        required = dist.RequiredVariables()
        missing = required - available
        if missing:
            missing_names = sorted(v.name for v in missing)
            raise ValueError(
                f"Distribution {type(dist).__name__} (index {i}) requires "
                f"variables {missing_names} to be set first, but they are "
                f"not covered by preceding distributions. "
                f"Available so far: {sorted(v.name for v in available)}"
            )
        if hasattr(dist, "SetVariables"):
            available |= dist.SetVariables()


def validate_injection_distributions(distributions):
    """Check that injection distributions are complete and correctly ordered.

    Validates:
    1. All required variable roles are covered (energy, direction, vertex).
    2. Ordering constraints are satisfied (RequiredVariables).

    Parameters
    ----------
    distributions : list
        Ordered list of PrimaryInjectionDistribution objects.

    Raises
    ------
    ValueError
        If required distributions are missing or ordering is wrong.
    """
    covered = collect_set_variables(distributions)

    missing = []
    for var in _REQUIRED_INJECTION_VARIABLES:
        if var not in covered:
            missing.append(var.name)

    has_vertex = bool(covered & _REQUIRED_INJECTION_VERTEX)
    if not has_vertex:
        missing.append("InitialPosition/InteractionVertex "
                       "(any VertexPositionDistribution subclass)")

    if missing:
        raise ValueError(
            "Missing injection distributions for:\n"
            + "\n".join(f"  - {m}" for m in missing)
        )

    validate_ordering(distributions)


def validate_physical_distributions(distributions):
    """Check that physical distributions cover required roles.

    Physical distributions need energy and direction but NOT position
    or mass.

    Parameters
    ----------
    distributions : list
        List of WeightableDistribution objects.

    Raises
    ------
    ValueError
        If required distribution types are missing.
    """
    covered = collect_set_variables(distributions)

    missing = []
    for var in _REQUIRED_PHYSICAL_VARIABLES:
        if var not in covered:
            missing.append(var.name)

    if missing:
        raise ValueError(
            "Missing physical distributions for:\n"
            + "\n".join(f"  - {m}" for m in missing)
        )


def validate_reweighting_compatibility(injection_distributions, physical_distributions):
    """Check that physical distributions can reweight injection distributions.

    Physical DensityVariables must be a subset of injection DensityVariables.
    Variables that are delta functions in both injection and physical are OK.

    Parameters
    ----------
    injection_distributions : list
        Injection distributions.
    physical_distributions : list
        Physical distributions.

    Raises
    ------
    ValueError
        If physical distributions are differential in variables that the
        injection distributions don't cover.
    """
    inj_density = collect_density_variables(injection_distributions)
    phys_density = collect_density_variables(physical_distributions)

    extra = phys_density - inj_density
    if extra:
        raise ValueError(
            f"Physical distributions are differential in variables "
            f"{sorted(extra)} that are not covered by injection distributions. "
            f"Injection density variables: {sorted(inj_density)}, "
            f"Physical density variables: {sorted(phys_density)}"
        )


def _expand_rule_names(rule):
    """Return the particle name(s) a child()/depth_below() rule references.

    child() rules name a particle; depth_below() rules name none.
    """
    name = getattr(rule, "name", None)
    return [name] if name is not None else []


def validate_expansion_wiring(vertices):
    """Check that a set of Vertex expand declarations forms a closed graph.

    Parameters
    ----------
    vertices : iterable
        Objects with ``.particle``, ``.expand`` (tuple of rules from
        ``siren.expand.child``/``depth_below``), and ``.continue_if``.
        An object may additionally expose ``.secondary_types`` (an
        iterable of ParticleType/name the vertex's interactions can
        produce); when absent, case (c) below is skipped for that
        vertex since there is nothing to cross-check against.

    Raises
    ------
    ConfigurationError
        (a) An expand rule names a child particle with no registered
            Vertex for that particle type.
        (b) A registered secondary Vertex (any vertex other than a root,
            i.e. one that some other vertex's expand list could reach)
            is unreachable from every other vertex's expand list.
        (c) A vertex declares secondaries (via ``.secondary_types``) but
            has no expansion declaration at all (empty ``.expand``).
    """
    vertices = list(vertices)
    by_particle = {}
    for v in vertices:
        ptype = _particles.resolve(v.particle)
        by_particle[ptype] = v

    named_children = set()
    for v in vertices:
        for rule in v.expand:
            for name in _expand_rule_names(rule):
                ptype = _particles.resolve(name)
                named_children.add(ptype)
                # (a) expand rule names a child with no registered Vertex.
                if ptype not in by_particle:
                    raise ConfigurationError(
                        f"Vertex {v.particle!r} declares expand rule "
                        f"child({name!r}) but no Vertex is registered for "
                        f"particle {name!r}. Fix: register a Vertex for "
                        f"{name!r}, or remove this expand rule."
                    )
            # A depth_below rule (no named particle) expands every secondary,
            # so all of this vertex's secondary types become reachable.
            if not _expand_rule_names(rule):
                for t in (getattr(v, "secondary_types", None) or []):
                    named_children.add(_particles.resolve(t))

    # (b) a registered secondary Vertex unreachable from any expand list.
    root = _particles.resolve(vertices[0].particle) if vertices else None
    for v in vertices:
        ptype = _particles.resolve(v.particle)
        if ptype == root:
            continue
        if ptype not in named_children:
            raise ConfigurationError(
                f"Vertex {v.particle!r} is registered but unreachable: no "
                f"other Vertex's expand list names it. Fix: add "
                f"child({v.particle!r}) to the expand list of the vertex "
                f"that produces it, or remove this Vertex."
            )

    # (c) secondaries present with no expansion declaration at all.
    for v in vertices:
        secondary_types = getattr(v, "secondary_types", None)
        if not secondary_types:
            continue
        if not v.expand:
            raise ConfigurationError(
                f"Vertex {v.particle!r} has secondaries "
                f"{list(secondary_types)!r} but declares no expand rules. "
                f"Fix: add expand=[siren.expand.child(...)] naming which "
                f"secondaries should recurse."
            )


def check_expand_vs_legacy_stopping(has_legacy_stopping, has_expand_or_continue_if):
    """Reject configurations mixing legacy stopping_condition with expand.

    Parameters
    ----------
    has_legacy_stopping : bool
        True if a legacy ``stopping_condition`` callable was supplied.
    has_expand_or_continue_if : bool
        True if any Vertex supplies ``expand`` and/or ``continue_if``.

    Raises
    ------
    ConfigurationError
        If both are set: the legacy callable and the declarative
        expand/continue_if fields are mutually exclusive.
    """
    if has_legacy_stopping and has_expand_or_continue_if:
        raise ConfigurationError(
            "Both a legacy stopping_condition callable and Vertex "
            "expand/continue_if fields are set. Fix: use one control "
            "surface only -- drop stopping_condition and express the "
            "chain with Vertex(expand=[...], continue_if=...)."
        )


def validate_secondary_keys(*dicts_with_labels):
    """Check that no secondary-keyed dict carries an unknown key.

    Each argument is a (label, mapping) pair. The first pair is the reference
    key set (the secondary interactions); every other mapping's keys must be a
    subset of it. An EXTRA key -- a phase-space or weighting-mode entry for a
    particle with no registered interaction -- raises ConfigurationError, since
    it is almost always a typo. A MISSING key is allowed: not every secondary
    needs a custom phase space or weighting mode.
    """
    pairs = [(label, mapping) for label, mapping in dicts_with_labels
             if mapping is not None]
    if not pairs:
        return
    ref_label, ref_map = pairs[0]
    ref_keys = set(ref_map.keys())
    for label, mapping in pairs[1:]:
        keys = set(mapping.keys())
        extra = keys - ref_keys
        if extra:
            raise ConfigurationError(
                "{} has secondary keys {} not present in {} (keys {}). "
                "Fix: the key is likely a typo or a stray entry; it must "
                "match a registered secondary type.".format(
                    label, sorted(str(k) for k in extra),
                    ref_label, sorted(str(k) for k in ref_keys)))


def build_template_record(signature, *, primary_mass=0.02, energy=0.05):
    """A minimal InteractionRecord for a signature, for config-time probes.

    Fills plausible small masses and a forward-moving primary; consulted only
    by measure/density probes, never by a real Sample().
    """
    import math
    import siren
    n_sec = len(signature.secondary_types)
    pz = math.sqrt(max(energy * energy - primary_mass * primary_mass, 0.0))
    rec = siren.dataclasses.InteractionRecord()
    rec.signature.primary_type = signature.primary_type
    rec.signature.target_type = signature.target_type
    rec.signature.secondary_types = list(signature.secondary_types)
    rec.primary_mass = primary_mass
    rec.primary_momentum = [energy, 0.0, 0.0, pz]
    rec.secondary_masses = [0.001] * n_sec
    rec.secondary_momenta = [[0.0, 0.0, 0.0, 0.0]] * n_sec
    rec.secondary_helicities = [0] * n_sec
    rec.interaction_vertex = [0.0, 0.0, 0.0]
    rec.primary_initial_position = [0.0, 0.0, 0.0]
    return rec


def probe_channel_densities(mcps, signature, detector_model, random,
                            *, samples_per_channel=100):
    """Run the detector-dependent ValidateChannelDensities probe on a mixture.

    Draws samples through each channel and checks the sampled density matches
    the analytic density, surfacing a sampler/density mismatch at build time.
    Raises whatever the engine raises (MeasureCompatibilityError,
    WeightCalculationError) on failure; a passing probe returns its result.
    """
    template = build_template_record(signature)
    return mcps.ValidateChannelDensities(
        random, detector_model, template, samples_per_channel)


def _pybind_base_of(model):
    """Return the pybind C++ interaction base a Python model derives from.

    The authoring bases (and legacy direct subclasses) sit above one of
    Decay/CrossSection/DarkNewsDecay/DarkNewsCrossSection; return the most
    derived such base found on the MRO, or None if the model is not
    trampoline-derived.
    """
    from . import interactions as _interactions
    candidates = []
    for name in ("DarkNewsDecay", "DarkNewsCrossSection", "Decay", "CrossSection"):
        base = getattr(_interactions, name, None)
        if base is not None and isinstance(model, base):
            candidates.append((name, base))
    if not candidates:
        return None
    # DarkNews bases derive from Decay/CrossSection, so prefer them when present.
    for name in ("DarkNewsDecay", "DarkNewsCrossSection"):
        for cand_name, base in candidates:
            if cand_name == name:
                return base
    return candidates[0][1]


# The required virtual method names per pybind base. A model whose MRO resolves
# one of these only to the abstract pybind base (no Python override) would abort
# in C++ with an opaque pure-virtual message; audit_overrides turns that into a
# named ConfigurationError.
_REQUIRED_DECAY_METHODS = (
    "GetPossibleSignatures", "GetPossibleSignaturesFromParent",
    "TotalDecayWidth", "TotalDecayWidthAllFinalStates",
    "DifferentialDecayWidth", "FinalStateProbability",
    "DensityVariables", "SampleFinalState",
)
_REQUIRED_CROSS_SECTION_METHODS = (
    "GetPossiblePrimaries", "GetPossibleTargets",
    "GetPossibleTargetsFromPrimary", "GetPossibleSignatures",
    "GetPossibleSignaturesFromParents", "TotalCrossSection",
    "DifferentialCrossSection", "InteractionThreshold",
    "FinalStateProbability", "DensityVariables", "SampleFinalState",
)


def _is_pybind_registered(klass):
    """True for a class registered directly by a pybind module.

    A subclass authored in Python still carries the pybind11 metaclass once it
    inherits a bound C++ base, so the metaclass alone cannot separate the two.
    A directly-registered C++ type has no method __dict__ authored in Python:
    its methods are builtin_function_or_method / instancemethod wrappers with a
    __module__ under the compiled interactions module, and it defines no
    ordinary Python function in its own namespace.
    """
    import types
    for value in vars(klass).values():
        if isinstance(value, types.FunctionType):
            return False
    return True


def _is_trampoline_derived(model, pybind_base):
    """True when a Python-authored class sits above the pybind base in the MRO.

    C++-native models resolve their virtuals in C++, so no class between
    type(model) and the pybind base defines an ordinary Python function; only
    genuinely Python-authored models do.
    """
    for klass in type(model).__mro__:
        if klass is pybind_base:
            return False
        if not _is_pybind_registered(klass):
            return True
    return False


def _resolves_before_root(model, method_name, abstract_root):
    """True if method_name is implemented above the abstract Decay/CrossSection.

    Walks type(model).__mro__ and reports the method satisfied when any class
    before the abstract root defines it -- a Python override or a concrete C++
    intermediate (e.g. DarkNewsDecay, whose FinalStateProbability delegates to
    the Python differential hook). Only a method resolving no earlier than the
    abstract root, which is pure there, is unimplemented.
    """
    for klass in type(model).__mro__:
        if klass is abstract_root:
            return False
        if method_name in vars(klass):
            return True
    return False


def audit_overrides(interactions):
    """Fail loudly on trampoline models with an unimplemented required virtual.

    ``interactions`` is any iterable yielding interaction objects (or an
    InteractionCollection). For every trampoline-derived model, walk its MRO
    for the required virtuals of Decay/CrossSection; a method that resolves only
    to the abstract root (pure there, with no override or concrete intermediate)
    raises ConfigurationError naming the class and the missing method. Also
    enforces the default-sampler measure contract: a model relying on the base
    SampleFinalState for a measure with no self-contained channel must override
    sample().
    """
    from . import interactions as _interactions

    cross_section_root = getattr(_interactions, "CrossSection")
    decay_root = getattr(_interactions, "Decay")

    models = _collect_models(interactions)
    for model in models:
        pybind_base = _pybind_base_of(model)
        if pybind_base is None:
            continue
        # C++-native models implement their virtuals in C++; only audit models
        # authored in Python via the trampoline.
        if not _is_trampoline_derived(model, pybind_base):
            continue
        if isinstance(model, cross_section_root):
            required = _REQUIRED_CROSS_SECTION_METHODS
            abstract_root = cross_section_root
        else:
            required = _REQUIRED_DECAY_METHODS
            abstract_root = decay_root
        for method_name in required:
            if not _resolves_before_root(model, method_name, abstract_root):
                raise ConfigurationError(
                    "%s does not implement the required method '%s' of its "
                    "interaction base; define or correctly spell it"
                    % (type(model).__name__, method_name))
        _audit_default_sampler(model)


def _audit_default_sampler(model):
    """Enforce the recursion-safe default-sampler contract for authoring bases.

    A model that leaves sample() at the authoring-base default must declare a
    measure with a self-contained channel; otherwise the default would need
    PhysicalChannelAdapters, which recurses. Skip models that override sample()
    or are not authoring-base derived.
    """
    from . import models as _models

    base_class = _models_base_sample(model)
    if base_class is None:
        return
    base_sample = getattr(base_class, "sample", None)
    if base_sample is None or type(model).sample is not base_sample:
        return
    measure = model.Measure()
    finals = len(model.GetPossibleSignatures()[0].secondary_types)
    if _models._self_contained_channel(measure, finals, 0) is None:
        raise ConfigurationError(
            "%s uses the default sampler for measure %r, which has no "
            "self-contained channel; override sample(record, random)"
            % (type(model).__name__, measure))


def _models_base_sample(model):
    """Find the authoring-base class carrying the default SampleFinalState."""
    from . import models as _models
    for klass in type(model).__mro__:
        if klass.__name__ in ("DecayModel", "CrossSectionModel") \
                and klass.__module__ == _models.__name__:
            return klass
    return None


def _collect_models(interactions):
    """Yield the individual cross-section and decay models from a container."""
    if interactions is None:
        return []
    if hasattr(interactions, "GetCrossSections") and hasattr(interactions, "GetDecays"):
        return list(interactions.GetCrossSections()) + list(interactions.GetDecays())
    return list(interactions)
