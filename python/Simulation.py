"""
High-level simulation interface that unifies event injection and weighting.

This class is model-agnostic.  All physics-model-specific logic
(DarkNews, HNL splines, etc.) lives in the process loaders under
``resources/processes/``.  The Simulation class only knows about the
generic interfaces: distributions, cross sections, and decays.

Usage::

    import siren

    sim = siren.Simulation(
        events=100_000,
        detector="IceCube",
        primary="NuMu",
        interactions="CSMSDISSplines",
        targets="Nucleon",
        process="CC",
        energy=siren.dist.PowerLaw(2, 1e3, 1e6),
        direction=siren.dist.IsotropicDirection(),
        position=siren.dist.ColumnDepth(600, 600.0),
    )
    results = sim.run()
    results.save("output/IceCube")
"""

from . import particles as _particles
from . import distributions as _distributions
from . import interactions as _interactions
from . import injection as _injection
from . import utilities as _utilities
from . import detector as _detector
from . import math as _math
from .Results import Results
from ._validation import (
    validate_physical_distributions,
    validate_reweighting_compatibility,
    collect_set_variables,
)
from ._validation import validate_physical_distributions
from .errors import ConfigurationError, MeasureCompatibilityError


def _is_vertex(obj):
    """Whether `obj` is a Vertex spec (duck-typed).

    A Vertex carries a resolved particle plus interaction and distribution
    lists; a plain particle name or ParticleType does not. Detection is by
    attribute presence so the check does not import the Vertex class.
    """
    return (
        not isinstance(obj, str)
        and hasattr(obj, "_resolved_particle")
        and hasattr(obj, "interactions")
        and hasattr(obj, "distributions")
    )


def _is_biasing_object(obj):
    """Whether `obj` is a channel-algebra biasing spec.

    Directed, channels.Mixture, and Directed.by_signature() results all
    expose ``to_mixture`` and ``compile``; the dict / Geometry forms handled
    by ``bias_targets`` do not.
    """
    return hasattr(obj, "to_mixture") or hasattr(obj, "compile")


class Simulation:
    """Single entry point for the inject-weight pipeline.

    Accepts all configuration as constructor keyword arguments.  Shared
    state (detector, particle type, interactions) is specified once and
    used for both injection and weighting internally.

    This class is model-agnostic: it does not contain any
    physics-model-specific logic.  Model-specific behavior (table
    management, secondary process setup, etc.) belongs in the
    process loaders invoked by ``load_processes``.

    Parameters
    ----------
    events : int
        Number of successful events to generate.  The deprecated alias
        ``n_events`` is still accepted but emits a ``DeprecationWarning``.
    detector : str or DetectorModel
        Detector name (e.g. ``"IceCube"``) or a pre-loaded DetectorModel.
    primary : str or ParticleType
        Primary particle type.  Accepts enum values or string names
        (e.g. ``"NuMu"``).
    interactions : str or list
        Interaction model name (e.g. ``"CSMSDISSplines"``) or a list of
        CrossSection/Decay objects.  When a string is provided, additional
        kwargs (``targets``, ``process``, ``isoscalar``, and any extra
        ``**process_kwargs``) are forwarded to ``load_processes``.
    energy : distribution, optional
        Energy distribution, used for both injection and weighting.
        Mutually exclusive with ``injection_energy``/``physical_energy``.
    direction : distribution, optional
        Direction distribution, used for both injection and weighting.
        Mutually exclusive with ``injection_direction``/``physical_direction``.
    position : distribution, optional
        Position distribution (injection only -- never used for weighting).
    injection_energy : distribution, optional
        Energy distribution for injection only.
    physical_energy : distribution, optional
        Energy distribution for weighting only.
    injection_direction : distribution, optional
        Direction distribution for injection only.
    physical_direction : distribution, optional
        Direction distribution for weighting only.
    injection_distributions : list, optional
        Explicit list of primary injection distributions, appended to the
        list built from the named slots.  Supplies an arbitrary set of
        primary distributions -- for example a custom
        VertexPositionDistribution -- directly, without the named slots or
        the Vertex API.  The named slots become optional as long as the
        combined list covers energy, direction, and a vertex role; ordering
        and completeness are validated the same way either route is used.
    physical_distributions : list, optional
        Explicit list of primary physical (weighting) distributions,
        appended to those built from the physical named slots.  Kept
        separate from ``injection_distributions``: a distribution placed in
        the injection list is not automatically reused for weighting.
    flux : distribution, optional
        Alias for ``physical_energy``.  When set alongside ``energy``,
        ``energy`` becomes the injection energy and ``flux`` the physical
        energy.
    seed : int, optional
        Random number seed.
    targets : list or str or ParticleType, optional
        Target particle type(s) for ``load_processes``.
    process : str, optional
        Process type filter for ``load_processes`` (e.g. ``"CC"``).
    isoscalar : bool, optional
        Whether to use isoscalar targets.  Forwarded to the interaction
        loader only when given, so models without the option (e.g.
        ``"MarleyCrossSection"``) load with their own defaults.
    secondary_interactions : dict, optional
        ``{ParticleType: [CrossSection/Decay]}`` for secondary processes.
    secondary_position : distribution or dict, optional
        Position distribution(s) for secondary vertices.  If a single
        distribution, it is used for all secondary types.  If a dict,
        keys are ParticleType and values are distributions or lists.
    stopping_condition : callable, optional
        ``f(datum, i) -> bool`` controlling secondary generation.
    bias_targets : Geometry or dict, optional
        Enable direction biasing for secondary decays/scattering.
        If a single Geometry (e.g. a fiducial volume), creates
        biased channels for all secondary signatures, directing
        the ``bias_daughter`` toward that geometry.  If a dict,
        maps ``{ParticleType: Geometry}`` for per-type targets, or
        ``{InteractionSignature: MultiChannelPhaseSpace}`` for full
        control.
    bias_daughter : str or ParticleType, optional
        Which daughter particle to bias toward the target.  Required
        when ``bias_targets`` is set.  Specified by ParticleType
        (e.g. ``"EMinus"``).  Must appear exactly once in each
        signature for the index to be unambiguous.
    bias_spectator : str or ParticleType, optional
        For 3-body decays, which daughter is the spectator (produced
        in the first 2-body step of the recursive decomposition).
        If not specified, defaults to the first daughter that is
        neither the ``bias_daughter`` nor its antiparticle.
    mass : float, optional
        Primary particle mass in GeV (default: 0 for neutrinos).
    **process_kwargs
        Additional keyword arguments forwarded to ``load_processes``
        when ``interactions`` is a string.  This allows model-specific
        parameters (e.g. DarkNews kwargs) to flow through without
        the Simulation class knowing about them.
    """

    def __init__(
        self,
        *,
        events=None,
        detector,
        primary,
        interactions=None,
        # Shared distributions (used for both injection and physical)
        energy=None,
        direction=None,
        position=None,
        # Split distributions (when injection != physical)
        injection_energy=None,
        physical_energy=None,
        injection_direction=None,
        physical_direction=None,
        # Explicit distribution lists (arbitrary primary distributions)
        injection_distributions=None,
        physical_distributions=None,
        # Flux convenience
        flux=None,
        # Process loading helpers
        seed=None,
        targets=None,
        process=None,
        isoscalar=None,
        # Secondary / BSM
        secondary_interactions=None,
        secondary_position=None,
        stopping_condition=None,
        # Secondary biasing
        bias_targets=None,
        bias_daughter=None,
        bias_spectator=None,
        # Particle mass
        mass=None,
        # Spec-form entry points (compile to the same internal structures)
        weighting=None,
        biasing=None,
        secondaries=(),
        # Forward to load_processes
        **process_kwargs,
    ):
        # ---- Resolve the events count (n_events is a deprecated alias) ----
        if "n_events" in process_kwargs:
            n_events = process_kwargs.pop("n_events")
            if events is not None:
                raise TypeError(
                    "Simulation() got both 'events' and 'n_events'. "
                    "'n_events' is a deprecated alias for 'events'; pass "
                    "only 'events'."
                )
            import warnings
            warnings.warn(
                "Simulation(n_events=...) is deprecated; use "
                "Simulation(events=...) instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            events = n_events
        if events is None:
            raise TypeError(
                "Simulation() missing required keyword-only argument: 'events'"
            )

        # ---- Store the per-vertex weighting mode (applied at build time) ----
        self._weighting_mode = weighting

        # ---- Resolve detector ----
        if isinstance(detector, str):
            self._detector_model = _utilities.load_detector(detector)
            self._detector_name = detector
        else:
            self._detector_model = detector
            self._detector_name = None

        # ---- Resolve the primary vertex ----
        # `primary` may be a full Vertex spec (carrying its own interactions
        # and distributions) or a plain particle name / ParticleType. A Vertex
        # is duck-typed by the fields Vertex exposes.
        primary_vertex = primary if _is_vertex(primary) else None
        if primary_vertex is not None:
            self._primary_type = primary_vertex._resolved_particle
        else:
            self._primary_type = _particles.resolve(primary)
        if targets is None:
            self._target_types = None
        else:
            try:
                self._target_types = [_particles.resolve(targets)]
            except TypeError:
                self._target_types = [_particles.resolve(t) if t is not None else None for t in targets]

        # ---- Resolve interactions ----
        self._secondary_processes = {}
        if primary_vertex is not None:
            if interactions is not None:
                raise TypeError(
                    "Simulation() got both a Vertex 'primary' and "
                    "'interactions'. A Vertex carries its own interactions; "
                    "pass one or the other.")
            self._primary_interactions = list(primary_vertex.interactions)
            self._process_metadata = ()
        else:
            if interactions is None:
                raise TypeError(
                    "Simulation() missing required keyword-only argument: "
                    "'interactions'")
            self._init_interactions(
                interactions, process, isoscalar, process_kwargs)

        # ---- Resolve secondary interactions ----
        if secondary_interactions is not None:
            resolved = {}
            for k, v in secondary_interactions.items():
                resolved[_particles.resolve(k) if isinstance(k, str) else k] = v
            self._secondary_processes = resolved

        # ---- Resolve distributions ----
        # A Vertex primary supplies its own injection distributions; the
        # explicit energy/direction/position kwargs otherwise drive the
        # injection and physical lists.
        if primary_vertex is not None and primary_vertex.distributions:
            self._resolve_distributions_from_vertex(
                primary_vertex,
                energy=energy,
                direction=direction,
                physical_energy=physical_energy,
                physical_direction=physical_direction,
                flux=flux,
                injection_distributions=injection_distributions,
                physical_distributions=physical_distributions,
            )
        else:
            self._resolve_distributions(
                energy=energy,
                direction=direction,
                position=position,
                injection_energy=injection_energy,
                physical_energy=physical_energy,
                injection_direction=injection_direction,
                physical_direction=physical_direction,
                flux=flux,
                mass=mass,
                injection_distributions=injection_distributions,
                physical_distributions=physical_distributions,
            )

        # ---- Resolve secondary interaction vertices (spec form) ----
        # A non-empty `secondaries` tuple of Vertex specs populates the
        # secondary processes, their injection distributions, and their
        # phase spaces; it is an alternative to the
        # secondary_interactions/secondary_position kwargs.
        self._secondary_injection_distributions = {}
        self._secondary_phase_spaces = {}  # {InteractionSignature: MultiChannelPhaseSpace}
        secondary_specs = list(secondaries) if secondaries else []
        if secondary_specs:
            if secondary_interactions is not None:
                raise ValueError(
                    "Simulation() got both 'secondaries' Vertex specs and "
                    "'secondary_interactions'. Use one secondary-configuration "
                    "surface only.")
            self._resolve_secondaries_from_specs(secondary_specs)

        # ---- Resolve secondary position ----
        if secondary_position is not None:
            self._resolve_secondary_position(secondary_position)
        elif self._secondary_processes and not self._secondary_injection_distributions:
            raise ValueError(
                "Secondary interaction processes were configured "
                f"(for {list(self._secondary_processes.keys())}), but no "
                "'secondary_position' distribution was provided. Each secondary "
                "process needs a secondary vertex distribution; pass "
                "secondary_position= to Simulation().")

        # ---- Resolve secondary biasing ----
        # `biasing` (a Directed / Mixture / Directed.by_signature object) and
        # the bias_targets/bias_daughter kwargs are mutually exclusive routes
        # to the same per-signature phase-space map.
        if biasing is not None and bias_targets is not None:
            raise ConfigurationError(
                "Simulation() got both 'biasing' and 'bias_targets'. These "
                "configure the same per-signature phase-space biasing; pass "
                "one or the other.")
        if biasing is not None:
            self._resolve_biasing_object(biasing)
        elif bias_targets is not None:
            if bias_daughter is None:
                raise ValueError(
                    "'bias_daughter' is required when 'bias_targets' is set. "
                    "Specify which daughter particle to direct toward the "
                    "target (e.g. bias_daughter='EMinus')."
                )
            self._resolve_bias_targets(
                bias_targets,
                _particles.resolve(bias_daughter),
                _particles.resolve(bias_spectator) if bias_spectator else None,
            )

        # ---- Store remaining config ----
        self._events = events
        self._seed = seed
        self._stopping_condition = stopping_condition

        # ---- Built lazily by run() ----
        self._injector = None
        self._weighter = None
        self._last_events = None
        self._last_gen_times = None

    # ------------------------------------------------------------------ #
    #  Interaction resolution                                              #
    # ------------------------------------------------------------------ #

    def _init_interactions(self, interactions, process, isoscalar, extra_kwargs):
        """Resolve interactions from a model name string or explicit list."""
        if isinstance(interactions, str):
            load_kwargs = dict(
                primary_types=[self._primary_type],
            )
            if isoscalar is not None:
                load_kwargs["isoscalar"] = isoscalar
            if self._target_types is not None:
                load_kwargs["target_types"] = self._target_types
            if process is not None:
                load_kwargs["process_types"] = [process]
            load_kwargs.update(extra_kwargs)

            bundle = _utilities.load_processes(interactions, **load_kwargs)

            if self._primary_type not in bundle.primary:
                available = list(bundle.primary.keys())
                raise ValueError(
                    f"Interaction model {interactions!r} did not produce "
                    f"processes for {self._primary_type}. "
                    f"Available: {available}"
                )
            self._primary_interactions = bundle.primary[self._primary_type]
            self._process_metadata = bundle.metadata
            if bundle.secondary:
                self._secondary_processes.update(bundle.secondary)
        elif isinstance(interactions, list):
            self._primary_interactions = interactions
            self._process_metadata = ()
        else:
            raise TypeError(
                f"'interactions' must be a model name string or a list of "
                f"CrossSection/Decay objects, got {type(interactions).__name__}"
            )

    # ------------------------------------------------------------------ #
    #  Distribution resolution                                             #
    # ------------------------------------------------------------------ #

    def _resolve_distributions(
        self,
        energy,
        direction,
        position,
        injection_energy,
        physical_energy,
        injection_direction,
        physical_direction,
        flux,
        mass,
        injection_distributions,
        physical_distributions,
    ):
        """Resolve the injection and physical distribution lists.

        The named keyword slots (``energy``/``direction``/``position`` and
        their ``injection_``/``physical_``/``flux`` variants) are a
        convenience that builds a distribution list.  Distributions passed
        explicitly through ``injection_distributions`` /
        ``physical_distributions`` are appended to that list, so an
        arbitrary set of primary distributions -- for example a custom
        VertexPositionDistribution -- can be supplied directly, without the
        named slots or the Vertex API.  Completeness (energy, direction,
        and a vertex role for injection) is checked against the combined
        list, so each role may be satisfied by either a named slot or an
        explicit distribution.

        Named-slot rules:
        - ``energy`` sets both injection and physical energy.
        - ``injection_energy`` / ``physical_energy`` override individually.
        - ``flux`` is an alias for ``physical_energy``.  When combined with
          ``energy``, ``energy`` becomes injection-only.
        - ``direction`` sets both.  ``injection_direction`` /
          ``physical_direction`` override individually.
        - ``position`` is injection-only (never in physical).
        - ``mass`` is injection-only (auto-inferred as 0 when omitted and no
          mass distribution is supplied explicitly).

        The explicit ``injection_distributions`` and
        ``physical_distributions`` are kept separate: an energy
        distribution placed in ``injection_distributions`` is used for
        injection only and does not populate the physical side (pass it in
        ``physical_distributions``, or use the named ``energy``/``flux``
        slots, for that).
        """
        DV = _distributions.DistributionVariable
        extra_injection = list(injection_distributions or [])
        extra_physical = list(physical_distributions or [])

        # ---- Energy ----
        inj_energy = injection_energy
        phys_energy = physical_energy

        if flux is not None and physical_energy is not None:
            raise ValueError(
                "Cannot specify both 'flux' and 'physical_energy'. "
                "'flux' is an alias for 'physical_energy'."
            )
        if flux is not None:
            phys_energy = flux

        if energy is not None:
            if inj_energy is not None:
                raise ValueError(
                    "Cannot specify both 'energy' and 'injection_energy'. "
                    "Use 'energy' for both, or the split form."
                )
            inj_energy = energy
            if phys_energy is None:
                # energy alone sets both; energy + flux keeps flux physical.
                phys_energy = energy
        if phys_energy is None and inj_energy is not None:
            phys_energy = inj_energy

        # ---- Direction ----
        inj_dir = injection_direction
        phys_dir = physical_direction

        if direction is not None:
            if inj_dir is not None:
                raise ValueError(
                    "Cannot specify both 'direction' and 'injection_direction'."
                )
            if phys_dir is not None:
                raise ValueError(
                    "Cannot specify both 'direction' and 'physical_direction'."
                )
            inj_dir = direction
            phys_dir = direction
        if phys_dir is None and inj_dir is not None:
            phys_dir = inj_dir

        # ---- Assemble injection list: named slots, then explicit extras ----
        injection = []
        if inj_energy is not None:
            injection.append(inj_energy)
        if inj_dir is not None:
            injection.append(inj_dir)
        if position is not None:
            injection.append(position)
        injection += extra_injection

        # ---- Completeness (a role may come from a slot OR an explicit dist) ----
        covered = collect_set_variables(injection)
        if DV.PrimaryEnergy not in covered:
            raise ValueError(
                "No energy distribution specified. Provide 'energy', "
                "'injection_energy' + 'physical_energy', or include a "
                "PrimaryEnergyDistribution in 'injection_distributions'."
            )
        if DV.PrimaryDirection not in covered:
            raise ValueError(
                "No direction distribution specified. Provide 'direction', "
                "'injection_direction' + 'physical_direction', or include a "
                "PrimaryDirectionDistribution in 'injection_distributions'."
            )
        if not (covered & {DV.InitialPosition, DV.InteractionVertex}):
            raise ValueError(
                "No position distribution specified. Provide 'position', or "
                "include a VertexPositionDistribution in "
                "'injection_distributions'."
            )

        # ---- Mass (injection only, auto-inferred unless supplied) ----
        if DV.PrimaryMass not in covered:
            mass_val = mass if mass is not None else 0.0
            injection.insert(0, _distributions.PrimaryMass(mass_val))
        elif mass is not None:
            raise ValueError(
                "Cannot specify both 'mass' and a mass distribution supplied "
                "through 'injection_distributions'."
            )

        # ---- Assemble physical list: named slots, then explicit extras ----
        physical = []
        if phys_energy is not None:
            physical.append(phys_energy)
        if phys_dir is not None:
            physical.append(phys_dir)
        physical += extra_physical

        self._injection_distributions = injection
        self._physical_distributions = physical

    def _resolve_distributions_from_vertex(
        self, vertex, energy, direction, physical_energy, physical_direction,
        flux, injection_distributions, physical_distributions,
    ):
        """Resolve distributions when the primary is a Vertex spec.

        The Vertex's own distribution list is the injection list (it already
        carries mass/energy/direction/position samplers); any explicit
        ``injection_distributions`` are appended to it, matching how the named
        keyword path composes them. The physical list is taken from the
        Vertex's ``physical`` field when present, else built from the explicit
        ``physical_energy``/``flux`` and ``physical_direction`` (falling back to
        the shared ``energy``/``direction`` kwargs, then to the energy- and
        direction-carrying members of the injection list). Explicit
        ``physical_distributions`` are appended to whichever physical list
        results.
        """
        self._injection_distributions = (
            list(vertex.distributions) + list(injection_distributions or []))

        if flux is not None and physical_energy is not None:
            raise ValueError(
                "Cannot specify both 'flux' and 'physical_energy'. "
                "'flux' is an alias for 'physical_energy'."
            )

        if getattr(vertex, "physical", None):
            self._physical_distributions = (
                list(vertex.physical) + list(physical_distributions or []))
            return

        phys_energy = physical_energy if physical_energy is not None else flux
        if phys_energy is None:
            phys_energy = energy
        phys_dir = physical_direction
        if phys_dir is None:
            phys_dir = direction

        phys = []
        if phys_energy is not None:
            phys.append(phys_energy)
        else:
            phys.extend(
                d for d in self._injection_distributions
                if isinstance(d, _distributions.PrimaryEnergyDistribution))
        if phys_dir is not None:
            phys.append(phys_dir)
        else:
            phys.extend(
                d for d in self._injection_distributions
                if isinstance(d, _distributions.PrimaryDirectionDistribution))
        phys.extend(physical_distributions or [])
        self._physical_distributions = phys

    def _resolve_secondaries_from_specs(self, secondaries):
        """Populate secondary structures from a list of Vertex specs.

        Each Vertex maps to one secondary type: its interaction list becomes
        the secondary process, its distribution list becomes the secondary
        injection distributions, and its ``kinematics`` (a Directed / Mixture)
        compiles to a per-signature phase space for every signature the
        vertex's models can produce.
        """
        for vertex in secondaries:
            if not _is_vertex(vertex):
                raise ConfigurationError(
                    "Simulation(secondaries=...): each entry must be a Vertex "
                    "spec, got {!r}".format(type(vertex).__name__))
            sec_type = vertex._resolved_particle
            self._secondary_processes[sec_type] = list(vertex.interactions)
            self._secondary_injection_distributions[sec_type] = list(
                vertex.distributions)
            if vertex.kinematics is not None:
                self._compile_biasing_for_models(
                    vertex.kinematics, vertex.interactions)

    def _resolve_biasing_object(self, biasing):
        """Build per-signature phase spaces from a channel-algebra spec.

        `biasing` is a Directed, a channels.Mixture, or a
        Directed.by_signature() result. It is compiled once per signature
        against the interaction model that owns that signature: the secondary
        models when secondaries are configured, else the primary models.
        """
        if not _is_biasing_object(biasing):
            raise ConfigurationError(
                "Simulation(biasing=...): expected a Directed, a "
                "channels.Mixture, or a Directed.by_signature() result, got "
                "{!r}".format(type(biasing).__name__))
        if self._secondary_processes:
            for models in self._secondary_processes.values():
                self._compile_biasing_for_models(biasing, models)
        else:
            self._compile_biasing_for_models(
                biasing, self._primary_interactions)

    def _compile_biasing_for_models(self, biasing, models):
        """Compile `biasing` for every signature of every model in `models`.

        Registers each compiled MultiChannelPhaseSpace under its signature in
        ``self._secondary_phase_spaces``. A Directed / Directed.by_signature()
        object exposes ``to_mixture(sig)``; a bare channels.Mixture is compiled
        directly. A per-signature compilation that raises a measure/
        configuration error is warned and skipped, matching the geometry-based
        bias path.
        """
        import warnings
        for model in models:
            for sig in model.GetPossibleSignatures():
                try:
                    if hasattr(biasing, "to_mixture"):
                        compilable = biasing.to_mixture(sig)
                    else:
                        compilable = biasing
                    mcps = compilable.compile(
                        sig, detector=self._detector_model, models=[model])
                except (MeasureCompatibilityError, ConfigurationError) as exc:
                    warnings.warn(
                        f"Skipping biasing for signature {sig}: {exc}",
                        stacklevel=3,
                    )
                    continue
                self._secondary_phase_spaces[sig] = mcps

    @staticmethod
    def _find_daughter_index(signature, daughter_type):
        """Find the index of a daughter type in a signature.

        Raises ValueError if the type is absent or ambiguous (appears
        more than once).
        """
        indices = [
            i for i, t in enumerate(signature.secondary_types)
            if t == daughter_type
        ]
        if len(indices) == 0:
            raise ValueError(
                f"Daughter type {daughter_type} not found in signature "
                f"{signature}. Available: {list(signature.secondary_types)}"
            )
        if len(indices) > 1:
            raise ValueError(
                f"Daughter type {daughter_type} appears {len(indices)} times "
                f"in signature {signature}. Use the per-signature API with "
                f"explicit indices instead."
            )
        return indices[0]

    @staticmethod
    def _build_phase_space_for_signature(
        target, sig, interaction, daughter_type, spectator_type=None,
        fraction=0.99,
    ):
        """Build a MultiChannelPhaseSpace for one signature.

        The physical fallback channel is appended FIRST at weight
        ``1 - fraction``; the directed channel is appended SECOND at weight
        ``fraction``. Channel selection walks cumulative weights in stored
        order, so this order is part of the mixture's identity: at
        ``fraction=0.99`` the compiled mixture is ``[physical@0.01,
        directed@0.99]``.

        Parameters
        ----------
        target : Geometry
            Target geometry to direct the daughter toward.
        sig : InteractionSignature
            The specific final-state signature.
        interaction : Decay or CrossSection
            The physics model for this interaction.
        daughter_type : ParticleType
            Which daughter to bias toward the target.
        spectator_type : ParticleType, optional
            For 3-body: which daughter is the spectator.
        fraction : float
            Weight assigned to the directed channel; ``1 - fraction`` goes
            to the physical fallback.

        Returns
        -------
        MultiChannelPhaseSpace
        """
        n_sec = len(sig.secondary_types)
        daughter_idx = Simulation._find_daughter_index(sig, daughter_type)

        channels = []
        weights = []

        # Physical channel (fallback for events that miss the target),
        # appended first so it occupies index 0 of the stored mixture.
        if isinstance(interaction, _interactions.Decay):
            channels.append(_injection.PhysicalDecayChannel(interaction, sig))
        else:
            channels.append(
                _injection.PhysicalCrossSectionChannel(interaction, sig))
        weights.append(1.0 - fraction)

        if n_sec == 2:
            channels.append(
                _injection.DetectorDirected2BodyChannel(
                    target, daughter_idx))
            weights.append(fraction)

        elif n_sec >= 3:
            # Determine spectator and pair indices
            if spectator_type is not None:
                spectator_idx = Simulation._find_daughter_index(
                    sig, spectator_type)
            else:
                # Default: first daughter that isn't the biased one
                spectator_idx = next(
                    i for i in range(n_sec) if i != daughter_idx
                )
            pair_indices = [
                i for i in range(n_sec)
                if i != spectator_idx
            ]
            directed_pair_idx = (
                daughter_idx if daughter_idx in pair_indices
                else pair_indices[0]
            )

            channels.append(
                _injection.DetectorDirected3BodyChannel(
                    factorization=_injection.ThreeBodyMode.Recursive,
                    target=target,
                    spectator_index=spectator_idx,
                    pair_first_index=pair_indices[0],
                    pair_second_index=pair_indices[1],
                    directed_index=directed_pair_idx,
                ))
            weights.append(fraction)

        # Normalize
        total = sum(weights)
        weights = [w / total for w in weights]

        mc = _injection.MultiChannelPhaseSpace()
        mc.channels = channels
        mc.weights = weights
        return mc

    def _resolve_bias_targets(self, bias_targets, daughter_type,
                              spectator_type, fraction=0.99):
        """Build per-signature phase spaces for direction biasing.

        Parameters
        ----------
        bias_targets : Geometry or dict
            If a Geometry, biases all secondary signatures toward it.
            If a dict, maps ``{ParticleType: Geometry}`` to direct that type's
            daughters toward its geometry, or
            ``{InteractionSignature: MultiChannelPhaseSpace}`` for full control.
        daughter_type : ParticleType
            Which daughter to direct toward the target.
        spectator_type : ParticleType or None
            For 3-body: which daughter is the spectator.
        fraction : float
            Weight assigned to the directed channel in each built mixture;
            ``1 - fraction`` goes to the physical fallback.
        """
        if isinstance(bias_targets, dict):
            for k, v in bias_targets.items():
                if isinstance(v, _injection.MultiChannelPhaseSpace):
                    # {InteractionSignature: MultiChannelPhaseSpace}: registered
                    # verbatim under the signature key.
                    self._secondary_phase_spaces[k] = v
                else:
                    # {ParticleType: Geometry}: direct that type's daughters
                    # toward the geometry, per signature that type's models can
                    # produce.
                    sec_type = _particles.resolve(k) if isinstance(k, str) else k
                    interactions_list = self._secondary_processes.get(sec_type)
                    if interactions_list is None:
                        raise ConfigurationError(
                            "bias_targets: key {!r} is not a configured "
                            "secondary particle type".format(sec_type))
                    self._bias_type_toward(
                        interactions_list, v, daughter_type, spectator_type,
                        fraction)
            return

        # Single geometry: build per-signature for all secondaries
        target = bias_targets
        for sec_type, interactions_list in self._secondary_processes.items():
            self._bias_type_toward(
                interactions_list, target, daughter_type, spectator_type,
                fraction)

    def _bias_type_toward(self, interactions_list, target, daughter_type,
                          spectator_type, fraction):
        """Register directed phase spaces for one secondary type's signatures.

        Builds one MultiChannelPhaseSpace per signature every model in
        ``interactions_list`` can produce, directing ``daughter_type`` toward
        ``target``. A signature in which the daughter is absent or ambiguous is
        skipped with a warning rather than aborting the whole configuration.
        """
        for interaction in interactions_list:
            for sig in interaction.GetPossibleSignatures():
                try:
                    mc = self._build_phase_space_for_signature(
                        target, sig, interaction,
                        daughter_type, spectator_type, fraction,
                    )
                    self._secondary_phase_spaces[sig] = mc
                except ValueError:
                    import warnings
                    warnings.warn(
                        f"Skipping biasing for signature {sig}: "
                        f"bias_daughter {daughter_type} not found "
                        f"or ambiguous in this channel.",
                        stacklevel=3,
                    )

    def _resolve_secondary_position(self, secondary_position):
        """Resolve secondary vertex position distributions."""
        if isinstance(secondary_position, dict):
            for k, v in secondary_position.items():
                key = _particles.resolve(k) if isinstance(k, str) else k
                self._secondary_injection_distributions[key] = (
                    v if isinstance(v, list) else [v]
                )
        else:
            for sec_type in self._secondary_processes.keys():
                self._secondary_injection_distributions[sec_type] = [
                    secondary_position
                ]

    # ------------------------------------------------------------------ #
    #  Build injector / weighter                                           #
    # ------------------------------------------------------------------ #

    def _build_injector(self):
        """Construct and return an Injector from stored config.

        The Layer-2 Injector wrapper validates the injection distributions in
        its own ``_build``, so the Simulation does not re-run
        ``validate_injection_distributions`` here.
        """
        # The per-vertex weighting mode is a constructor-only field on the
        # wrapper; when None (the default) it leaves the process at
        # Propagated(), so the engine configuration is unchanged.
        injector = _injection.Injector(
            primary_weighting_mode=self._weighting_mode)
        injector.number_of_events = self._events
        injector.detector_model = self._detector_model
        injector.primary_type = self._primary_type
        injector.primary_interactions = self._primary_interactions
        injector.primary_injection_distributions = self._injection_distributions

        if self._seed is not None:
            injector.seed = self._seed

        if self._secondary_processes:
            injector.secondary_interactions = self._secondary_processes
            injector.secondary_injection_distributions = (
                self._secondary_injection_distributions
            )
            if self._secondary_phase_spaces:
                # Group per-signature phase spaces by primary type
                # for the Injector wrapper's {ParticleType: {sig: mc}} format
                by_type = {}
                for sig, mc in self._secondary_phase_spaces.items():
                    pt = sig.primary_type
                    if pt not in by_type:
                        by_type[pt] = {}
                    by_type[pt][sig] = mc
                injector.secondary_phase_spaces = by_type

        if self._stopping_condition is not None:
            injector.stopping_condition = self._stopping_condition

        return injector

    def _build_weighter(self, injector):
        """Construct and return a Weighter from stored config."""
        validate_physical_distributions(self._physical_distributions)
        validate_reweighting_compatibility(
            self._injection_distributions,
            self._physical_distributions,
        )

        weighter = _injection.Weighter()
        weighter.injectors = [injector]
        weighter.detector_model = self._detector_model
        weighter.primary_type = self._primary_type
        weighter.primary_interactions = self._primary_interactions
        weighter.primary_physical_distributions = self._physical_distributions

        if self._secondary_processes:
            weighter.secondary_interactions = self._secondary_processes
            weighter.secondary_physical_distributions = {}

        return weighter

    # ------------------------------------------------------------------ #
    #  Run                                                                 #
    # ------------------------------------------------------------------ #

    def run(self, *, optimize=False, on_shortfall="warn"):
        """Generate events and compute weights.

        Parameters
        ----------
        optimize : bool or Plan, optional
            When truthy, tune the injection chain's multi-channel phase-space
            weights before generating the real events. ``True`` uses
            ``tune()``'s defaults; a :class:`siren.tune.Plan` supplies its
            keyword arguments. Tuning warms up the channel weights; the
            injector's counters are reset afterward so the tuning batch does
            not count against the requested event total.
        on_shortfall : str, optional
            Policy when fewer than ``events`` successful trees are generated:
            ``'warn'`` (default), ``'raise'``, or ``'ignore'``.

        Returns
        -------
        Results
            A :class:`Results` object containing events, weights, and
            per-event generation times. The generation times are the uniform
            run-average (total generation wall-clock time divided evenly
            across the events), not individually measured per-event times.
        """
        import time

        injector = self._build_injector()
        weighter = self._build_weighter(injector)

        if optimize:
            self._optimize(injector, weighter, optimize)
            # Clear the tuning batch's injected/attempt counters so the real
            # generation starts from a full budget.
            injector.reset()

        # Forward generation to the Layer-2 injector, which counts successes,
        # retries failed attempts, and applies the shortfall policy itself.
        t0 = time.time()
        events = injector.generate(self._events, on_shortfall=on_shortfall)
        elapsed = time.time() - t0
        # Generation is timed for the whole run, not per event, so record the
        # uniform run-average (total wall-clock time split evenly across the
        # events) rather than a measured per-event time.
        per_event = (elapsed / len(events)) if events else 0.0
        gen_times = [per_event] * len(events)

        weights = [weighter(event) for event in events]

        self._injector = injector
        self._weighter = weighter
        self._last_events = events
        self._last_gen_times = gen_times

        # Snapshot the requested count explicitly: generate() raises the
        # injector's engine attempt quota to events * 1000, so its live
        # number_of_events no longer reports the request.
        return Results(events, weights, gen_times, weighter, injector,
                       requested=self._events)

    def _optimize(self, injector, weighter, optimize):
        """Tune the injection chain's channel weights in place.

        ``optimize`` is either ``True`` (tune with defaults) or a
        :class:`siren.tune.Plan` carrying the ``tune()`` keyword arguments.
        """
        from .tune import tune as _tune, Plan as _Plan
        if isinstance(optimize, _Plan):
            plan = optimize
            _tune(
                injector, weighter,
                events=plan.events,
                rounds=plan.rounds,
                rule=plan.rule,
                metric=plan.metric,
                damping=plan.damping,
                min_weight=plan.min_weight,
                failure_mode=plan.failure_mode,
                group_directed=plan.group_directed,
            )
        else:
            _tune(injector, weighter)

    # ------------------------------------------------------------------ #
    #  Reweight                                                            #
    # ------------------------------------------------------------------ #

    def reweight(
        self,
        *,
        physical_energy=None,
        physical_direction=None,
        physical_distributions=None,
        interactions=None,
        secondary_interactions=None,
    ):
        """Reweight existing events with new physical parameters.

        Builds a new Weighter with the specified changes and applies it
        to the events from the last :meth:`run` call.  Does not
        re-generate events.  Weighter construction is cheap (pointer
        copies only), so this can be called many times for parameter
        scans.

        Any parameter not specified is kept from the original simulation.

        Parameters
        ----------
        physical_energy : distribution, optional
            New physical energy distribution.
        physical_direction : distribution, optional
            New physical direction distribution.
        physical_distributions : list, optional
            Explicit replacement for the entire physical distribution list.
            Mutually exclusive with ``physical_energy``/``physical_direction``
            (which patch individual roles instead).
        interactions : list, optional
            New primary interaction list (CrossSection/Decay objects).
        secondary_interactions : dict, optional
            New secondary interactions ``{ParticleType: [CrossSection/Decay]}``.

        Returns
        -------
        Results
            A new :class:`Results` with the same events but new weights.

        Raises
        ------
        RuntimeError
            If :meth:`run` has not been called yet.

        Examples
        --------
        Flux reweighting::

            results = sim.run()
            new_results = sim.reweight(
                physical_energy=siren.load_flux("BNB", tag="FHC_numu",
                                                physically_normalized=True),
            )

        Coupling scan::

            results = sim.run()
            for mu in [1e-6, 1e-7, 1e-8]:
                bundle = siren.load_processes("DarkNewsTables",
                                              mu_tr_mu4=mu, ...)
                new_results = sim.reweight(
                    interactions=bundle.primary[siren.particles.NuMu],
                    secondary_interactions=bundle.secondary,
                )
        """
        if self._injector is None:
            raise RuntimeError(
                "Must call run() before reweight(). "
                "No events have been generated yet."
            )

        # Build new physical distribution list.  An explicit
        # physical_distributions list replaces the whole set; otherwise the
        # named slots patch the individual energy/direction roles.
        if physical_distributions is not None:
            if physical_energy is not None or physical_direction is not None:
                raise ValueError(
                    "Cannot combine 'physical_distributions' with "
                    "'physical_energy'/'physical_direction'. The explicit "
                    "list replaces the entire physical distribution set."
                )
            phys_dists = list(physical_distributions)
        else:
            phys_dists = list(self._physical_distributions)
            if physical_energy is not None:
                phys_dists = [
                    d for d in phys_dists
                    if not isinstance(d, _distributions.PrimaryEnergyDistribution)
                ]
                phys_dists.append(physical_energy)
            if physical_direction is not None:
                phys_dists = [
                    d for d in phys_dists
                    if not isinstance(d, _distributions.PrimaryDirectionDistribution)
                ]
                phys_dists.append(physical_direction)

        validate_physical_distributions(phys_dists)
        validate_reweighting_compatibility(
            self._injection_distributions,
            phys_dists,
        )

        # Resolve interactions
        primary_interactions = (
            interactions if interactions is not None
            else self._primary_interactions
        )
        secondary = (
            secondary_interactions if secondary_interactions is not None
            else self._secondary_processes
        )

        # Build a new Weighter (cheap -- just pointer copies)
        weighter = _injection.Weighter()
        weighter.injectors = [self._injector]
        weighter.detector_model = self._detector_model
        weighter.primary_type = self._primary_type
        weighter.primary_interactions = primary_interactions
        weighter.primary_physical_distributions = phys_dists

        if secondary:
            weighter.secondary_interactions = secondary
            weighter.secondary_physical_distributions = {}

        weights = [weighter(event) for event in self._last_events]

        return Results(
            self._last_events, weights, self._last_gen_times,
            weighter, self._injector, requested=self._events,
        )

    # ------------------------------------------------------------------ #
    #  Properties for inspection                                           #
    # ------------------------------------------------------------------ #

    @property
    def detector_model(self):
        """The resolved DetectorModel."""
        return self._detector_model

    @property
    def primary_type(self):
        """The resolved primary ParticleType."""
        return self._primary_type

    @property
    def injection_distributions(self):
        """List of injection distributions (read-only)."""
        return list(self._injection_distributions)

    @property
    def physical_distributions(self):
        """List of physical distributions (read-only)."""
        return list(self._physical_distributions)

    @property
    def phase_spaces(self):
        """Per-signature phase space configurations (read-only).

        Returns a dict mapping ``InteractionSignature`` to
        ``MultiChannelPhaseSpace`` for all signatures that have
        biased phase space channels registered.  Empty if no biasing
        is configured.
        """
        return dict(self._secondary_phase_spaces)

    @property
    def process_metadata(self):
        """Extra values returned by load_processes beyond (primary, secondary).

        Model-specific loaders may return additional metadata (e.g. keys
        for DarkNews table saving).  This property provides access without
        the Simulation class needing to know what the metadata means.
        """
        return self._process_metadata

    @property
    def injector(self):
        """The Layer-2 :class:`siren.Injector` wrapper, built on first access.

        After :meth:`run` returns this is the injector that run built and
        used; run() always builds fresh and replaces any cached instance, so
        a pre-run access does not configure the run. Exposes ``report()``,
        ``export_tuning()``, ``apply_tuning()``, and ``engine``.
        """
        if self._injector is None:
            self._injector = self._build_injector()
        return self._injector

    @property
    def weighter(self):
        """The Layer-2 :class:`siren.Weighter` wrapper, built on first access.

        Built against :attr:`injector`; after :meth:`run` returns this is the
        weighter that run built and used (run() replaces any cached
        instance). Exposes ``explain(tree)`` and ``engine``.
        """
        if self._weighter is None:
            self._weighter = self._build_weighter(self.injector)
        return self._weighter

    def describe(self):
        """Return a human-readable multi-line summary of the pipeline.

        Reports the requested event count, primary type, detector, interaction
        and distribution kinds, any secondary processes, the per-signature
        biasing count with its directed fraction, and the seed.
        """
        lines = ["Simulation:"]
        lines.append("  events requested: {}".format(self._events))
        lines.append("  primary type:     {}".format(self._primary_type))
        detector = ("<loaded>" if self._detector_name is None
                    else self._detector_name)
        lines.append("  detector:         {}".format(detector))

        n_primary = len(self._primary_interactions)
        lines.append("  primary interactions: {} model(s)".format(n_primary))

        inj_kinds = ", ".join(
            sorted(type(d).__name__ for d in self._injection_distributions))
        phys_kinds = ", ".join(
            sorted(type(d).__name__ for d in self._physical_distributions))
        lines.append("  injection distributions: {}".format(inj_kinds or "none"))
        lines.append("  physical distributions:  {}".format(phys_kinds or "none"))

        if self._secondary_processes:
            lines.append("  secondary processes:")
            for sec_type, models in self._secondary_processes.items():
                lines.append("    {}: {} model(s)".format(
                    sec_type, len(models)))
        else:
            lines.append("  secondary processes: none")

        n_biased = len(self._secondary_phase_spaces)
        if n_biased:
            lines.append(
                "  biasing: {} signature(s) with directed channels".format(
                    n_biased))
        else:
            lines.append("  biasing: none")

        if self._weighting_mode is not None:
            lines.append("  weighting mode:   {}".format(self._weighting_mode))
        lines.append("  seed:             {}".format(
            self._seed if self._seed is not None else "<auto>"))
        return "\n".join(lines)

    def __repr__(self):
        parts = [
            f"events={self._events}",
            f"primary={self._primary_type}",
            f"detector={'<loaded>' if self._detector_name is None else self._detector_name!r}",
        ]
        return f"Simulation({', '.join(parts)})"
