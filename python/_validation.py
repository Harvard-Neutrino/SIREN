"""
Input validation helpers for SIREN's Python interface.

Uses the C++ distribution methods ``ProvidedVariables()``,
``RequiredVariables()``, and ``DensityVariables()`` to validate
distribution lists.  No hardcoded type lists -- adding a new
distribution in C++ with the correct overrides makes it work here
automatically.
"""

from collections import Counter

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

_PRIMARY_VERTEX_DENSITY = "InteractionVertexPosition"
_PRIMARY_AREA_DENSITY = "PointOfClosestApproach"
_PRIMARY_TRANSVERSE_DENSITIES = (
    "PrimaryPositionTransverse[0]",
    "PrimaryPositionTransverse[1]",
)
_PRIMARY_LONGITUDINAL_DENSITY = "PrimaryPositionLongitudinal"


def collect_set_variables(distributions):
    """Return the union of all ProvidedVariables across a list of distributions."""
    result = set()
    for d in distributions:
        if hasattr(d, "ProvidedVariables"):
            result |= d.ProvidedVariables()
    return result


def collect_density_variables(distributions):
    """Return the union of all DensityVariables across a list of distributions."""
    result = set()
    for d in distributions:
        result.update(d.DensityVariables())
    return result


def validate_ordering(distributions):
    """Check that each distribution's RequiredVariables are satisfied
    by the ProvidedVariables of preceding distributions.

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
        if hasattr(dist, "ProvidedVariables"):
            newly_set = dist.ProvidedVariables()
            overlap = newly_set & available
            if overlap:
                raise ValueError(
                    f"Distribution {type(dist).__name__} (index {i}) sets "
                    f"variables {sorted(v.name for v in overlap)} that have "
                    f"already been set by preceding distributions. "
                    f"Available so far: {sorted(v.name for v in available)}"
                )
            available |= newly_set


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


def validate_creation_position(process, *, is_primary):
    """Require an origin-preserving sampler for creation survival.

    Unknown/custom position samplers need a separately verified origin contract.
    PointSource is excluded because clipping can move its source to world entry.
    """
    if not process.GetWeightingMode().survival_from_creation:
        return
    if is_primary:
        base = _d.VertexPositionDistribution
        supported = (
            _d.PrimaryBoundedVertexDistribution,
            _d.PrimaryPhysicalVertexDistribution,
            _d.FixedTargetPositionDistribution,
        )
    else:
        base = _d.SecondaryVertexPositionDistribution
        supported = (
            _d.SecondaryBoundedVertexDistribution,
            _d.SecondaryPhysicalVertexDistribution,
        )
    positions = [d for d in process.distributions if isinstance(d, base)
                 and (not is_primary or DV.InteractionVertex in d.ProvidedVariables())]
    if len(positions) != 1 or type(positions[0]) not in supported:
        raise ConfigurationError(
            "PropagatedFromCreation requires a position sampler that retains "
            "the creation point: use a bounded/physical vertex sampler "
            "(or a primary FixedTargetPosition). Volume, column-depth, range, "
            "decay-range and clipped point-source samplers replace that point."
        )
    if is_primary and type(positions[0]) is _d.FixedTargetPositionDistribution:
        if any(DV.InitialPosition in d.ProvidedVariables()
               for d in process.distributions if d is not positions[0]):
            raise ConfigurationError(
                "FixedTargetPosition defines its own creation point; it cannot "
                "preserve an external origin for PropagatedFromCreation."
            )


def _effective_density_variables(distributions, physical=False):
    """Return differential components, preserving repeated factors.

    On the physical side a distribution's ``PhysicalDensityVariables`` are
    used when it exposes them (a segment-mode external table samples a
    longitudinal position on the injection side only; its physical density
    leaves that factor to the weighter's normalized position density).
    """
    result = Counter()
    for distribution in distributions:
        getter = distribution.DensityVariables
        if physical:
            getter = getattr(distribution, "PhysicalDensityVariables", getter)
        for variable in getter():
            if variable == _PRIMARY_VERTEX_DENSITY:
                result.update(_PRIMARY_TRANSVERSE_DENSITIES)
                result[_PRIMARY_LONGITUDINAL_DENSITY] += 1
            elif variable == _PRIMARY_AREA_DENSITY:
                result.update(_PRIMARY_TRANSVERSE_DENSITIES)
            else:
                result[variable] += 1
    return result


def _format_density_variables(variables):
    return sorted(
        name if count == 1 else f"{name} (x{count})"
        for name, count in variables.items()
        if count > 0
    )


def validate_reweighting_compatibility(
    injection_distributions,
    physical_distributions,
    *,
    compute_interaction_probability=True,
    compute_position_probability=True,
):
    """Check that physical distributions can reweight injection distributions.

    Effective physical density variables must be a subset of the injection
    density variables.  Variables that are delta functions in both injection
    and physical are OK.

    Primary position measures are compared in transverse/longitudinal
    coordinates.  A ``VertexPositionDistribution`` supplies all three.  A
    ``PrimaryAreaDistribution`` supplies the two transverse coordinates, and
    the downstream normalized-position factor supplies the longitudinal
    coordinate when ``compute_position_probability`` is true.  Consequently,
    ordinary propagated weighting may leave two injection coordinates for an
    external flux-per-area factor, while a fixed-target area distribution plus
    normalized position density fully matches a fixed-target vertex
    distribution.  The interaction-probability factor is dimensionless, so
    ``compute_interaction_probability`` does not add a density variable.

    Parameters
    ----------
    injection_distributions : list
        Injection distributions.
    physical_distributions : list
        Physical distributions.
    compute_interaction_probability : bool, optional
        Whether downstream weighting includes the dimensionless interaction
        probability.  Accepted as part of the complete weighting-mode context;
        it does not change differential-variable compatibility.
    compute_position_probability : bool, optional
        Whether downstream weighting includes normalized position density.

    Raises
    ------
    ValueError
        If physical distributions are differential in variables that the
        injection distributions don't cover.
    """
    inj_density = _effective_density_variables(injection_distributions)
    phys_density = _effective_density_variables(physical_distributions, physical=True)

    if compute_position_probability:
        phys_density[_PRIMARY_LONGITUDINAL_DENSITY] += 1

    extra = phys_density - inj_density
    if extra:
        raise ValueError(
            "Physical weighting is differential in variables "
            f"{_format_density_variables(extra)} that are not covered by "
            "injection distributions. "
            "Effective injection density variables: "
            f"{_format_density_variables(inj_density)}, effective physical "
            "density variables (including downstream weighting factors): "
            f"{_format_density_variables(phys_density)}"
        )


def validate_expansion_wiring(vertices):
    """Validate declarations and reachability from the primary vertex.

    Empty expansion terminates a vertex. Predicates can veto a reachable
    daughter at runtime; they cannot make an undeclared daughter reachable.
    """
    from .expand import _ChildRule, _DepthBelowRule
    vertices = list(vertices)
    if not vertices:
        return
    secondary = {}
    for index, vertex in enumerate(vertices[1:], 1):
        particle = _particles.resolve(vertex.particle)
        if particle in secondary:
            raise ConfigurationError("Duplicate secondary vertex for {!r}".format(particle))
        secondary[particle] = index

    for vertex in vertices:
        if vertex.continue_if is not None and not callable(vertex.continue_if):
            raise ConfigurationError("continue_if must be callable or None")
        daughters = getattr(vertex, "secondary_types", None)
        if daughters is not None:
            daughters = {_particles.resolve(t) for t in daughters}
        for rule in vertex.expand:
            if isinstance(rule, _ChildRule):
                particle = _particles.resolve(rule.name)
                if daughters is not None and particle not in daughters:
                    raise ConfigurationError(
                        "Vertex {!r} does not produce the declared daughter {!r}".format(
                            vertex.particle, rule.name))
                if particle not in secondary:
                    raise ConfigurationError(
                        "Vertex {!r} names child {!r}, but no secondary Vertex is registered".format(
                            vertex.particle, rule.name))
            elif not isinstance(rule, _DepthBelowRule):
                raise ConfigurationError("expand entries must be child() or depth_below() rules")
        if len(vertices) == 1 and (vertex.continue_if is not None or any(
                isinstance(rule, _DepthBelowRule) and rule.n > 0 for rule in vertex.expand)):
            raise ConfigurationError("Expansion requires a registered secondary Vertex")

    pending = [(0, 0)]
    reached = {0}
    for index, depth in pending:
        vertex = vertices[index]
        daughters = getattr(vertex, "secondary_types", None)
        candidates = (set(secondary) if daughters is None
                      else {_particles.resolve(t) for t in daughters})
        signatures = getattr(vertex, "secondary_signatures", None)
        for rule in vertex.expand:
            if isinstance(rule, _ChildRule):
                particle = _particles.resolve(rule.name)
                if signatures is None:
                    matches = (rule.index is None or
                               (0 <= rule.index < 2**63 and rule.index == int(rule.index)))
                else:
                    matches = any(rule._matches(t, i) for signature in signatures
                                  for i, t in enumerate(signature))
                permitted = {particle} if matches else set()
            else:
                permitted = candidates if rule.n > depth else set()
            for particle in permitted:
                child_index = secondary.get(particle)
                if child_index is not None and child_index not in reached:
                    reached.add(child_index)
                    pending.append((child_index, depth + 1))
    for index, vertex in enumerate(vertices[1:], 1):
        if index not in reached:
            raise ConfigurationError(
                "Vertex {!r} is unreachable from the primary expand rules".format(vertex.particle))


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


def model_secondary_masses(model, secondary_types):
    """Resolve an author-supplied mass list with a named contract error."""
    types = list(secondary_types)
    label = type(model).__name__ + '.SecondaryMasses'
    try:
        masses = [float(mass) for mass in model.SecondaryMasses(types)]
    except (TypeError, ValueError, OverflowError) as exc:
        raise ConfigurationError(
            '%s must return numeric masses for %d secondary types: %s'
            % (label, len(types), exc)) from exc
    if len(masses) != len(types):
        raise ConfigurationError(
            '%s returned %d masses for %d secondary types'
            % (label, len(masses), len(types)))
    return masses


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


_REQUIRED_DECAY_METHODS = (
    "equal",
    "GetPossibleSignatures", "GetPossibleSignaturesFromParent",
    "TotalDecayWidth", "TotalDecayWidthAllFinalStates",
    "DifferentialDecayWidth", "FinalStateProbability",
    "DensityVariables", "SampleFinalState",
)
_REQUIRED_CROSS_SECTION_METHODS = (
    "equal",
    "GetPossiblePrimaries", "GetPossibleTargets",
    "GetPossibleTargetsFromPrimary", "GetPossibleSignatures",
    "GetPossibleSignaturesFromParents", "TotalCrossSection",
    "DifferentialCrossSection", "InteractionThreshold",
    "FinalStateProbability", "DensityVariables", "SampleFinalState",
)


def is_trampoline(obj):
    """Whether an interaction or distribution has a Python-authored type.

    Compare type identity with the native extension exports. Python subclasses
    retain the pybind metaclass and may define no methods of their own.
    """
    from . import interactions as _interactions
    from . import injection as _injection
    return not any(type(obj) is value
                   for module in (_interactions, _d, _injection)
                   for value in vars(module).values())


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
    requires authoring models to supply a sampler matching their density.

    For an authoring-base-derived model, also requires its declared physics
    hook -- total_width/differential_width, or total_xs/differential_xs -- to
    be overridden rather than left at the authoring base's own
    NotImplementedError-raising default.
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
        if not is_trampoline(model):
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
        _audit_required_physics_hooks(model)


def _audit_default_sampler(model):
    """Require an explicit sampler on authoring models."""
    base_class = _models_base_sample(model)
    if base_class is None:
        return
    for name in ('SampleFinalState', 'sample'):
        implementation = getattr(model, name, None)
        if (callable(implementation)
                and getattr(implementation, '__func__', implementation)
                is not getattr(base_class, name)):
            return
    raise ConfigurationError(
        "%s must implement sample(record, random); a measure does not "
        "determine the sampling density" % type(model).__name__)


def _audit_required_physics_hooks(model):
    """Fail loudly when an authoring-base subclass never supplies a physics hook.

    ``total_width``/``differential_width`` (DecayModel) and ``total_xs``/
    ``differential_xs`` (CrossSectionModel) default to a NotImplementedError-
    raising stub on the authoring base itself. The required-virtual walk
    above never catches their omission: TotalDecayWidth/TotalCrossSection and
    the rest of the pybind virtual set are always implemented on the
    authoring base, so they resolve before the abstract root regardless of
    whether the subclass ever supplies a working physics hook. Walk the MRO
    the same way ``_audit_default_sampler`` does via ``_models_base_sample``:
    a hook resolving no earlier than the authoring base itself is still the
    NotImplementedError stub. Skip models that are not authoring-base
    derived.
    """
    base_class = _models_base_sample(model)
    if base_class is None:
        return
    if base_class.__name__ == "DecayModel":
        hook_names = ("total_width", "differential_width")
    else:
        hook_names = ("total_xs", "differential_xs")
    for hook_name in hook_names:
        if not _resolves_before_root(model, hook_name, base_class):
            raise ConfigurationError(
                "%s does not implement the required physics hook '%s'; "
                "define it" % (type(model).__name__, hook_name))


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
