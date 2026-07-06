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

    # (b) a registered secondary Vertex unreachable from any expand list.
    root = vertices[0].particle if vertices else None
    for v in vertices:
        if v.particle == root:
            continue
        ptype = _particles.resolve(v.particle)
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
    """Check that every secondary-keyed dict shares the interactions keys.

    Each argument is a (label, mapping) pair. The first pair is the reference
    key set (the secondary interactions); every other mapping must key on
    exactly the same ParticleTypes. A typo'd or missing key raises
    ConfigurationError naming the offending mapping and the key difference,
    so a phase-space or weighting-mode dict cannot silently miss a secondary.
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
