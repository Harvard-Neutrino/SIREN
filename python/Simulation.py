"""
High-level simulation interface that unifies event injection and weighting.

This class is model-agnostic.  All physics-model-specific logic
(DarkNews, HNL splines, etc.) lives in the process loaders under
``resources/processes/``.  The Simulation class only knows about the
generic interfaces: distributions, cross sections, and decays.

Usage::

    import siren

    sim = siren.Simulation(
        n_events=100_000,
        detector="IceCube",
        primary="NuMu",
        interactions="CSMSDISSplines",
        target="Nucleon",
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
from ._util import GenerateEvents as _GenerateEvents
from .Results import Results
from ._validation import (
    validate_injection_distributions,
    validate_physical_distributions,
)


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
    n_events : int
        Number of events to generate.
    detector : str or DetectorModel
        Detector name (e.g. ``"IceCube"``) or a pre-loaded DetectorModel.
    primary : str or ParticleType
        Primary particle type.  Accepts enum values or string names
        (e.g. ``"NuMu"``).
    interactions : str or list
        Interaction model name (e.g. ``"CSMSDISSplines"``) or a list of
        CrossSection/Decay objects.  When a string is provided, additional
        kwargs (``target``, ``process``, ``isoscalar``, and any extra
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
    flux : distribution, optional
        Alias for ``physical_energy``.  When set alongside ``energy``,
        ``energy`` becomes the injection energy and ``flux`` the physical
        energy.
    seed : int, optional
        Random number seed.
    target : str or ParticleType, optional
        Target particle type for ``load_processes``.
    process : str, optional
        Process type filter for ``load_processes`` (e.g. ``"CC"``).
    isoscalar : bool
        Whether to use isoscalar targets (default ``True``).
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
        n_events,
        detector,
        primary,
        interactions,
        # Shared distributions (used for both injection and physical)
        energy=None,
        direction=None,
        position=None,
        # Split distributions (when injection != physical)
        injection_energy=None,
        physical_energy=None,
        injection_direction=None,
        physical_direction=None,
        # Flux convenience
        flux=None,
        # Process loading helpers
        seed=None,
        target=None,
        process=None,
        isoscalar=True,
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
        # Forward to load_processes
        **process_kwargs,
    ):
        # ---- Resolve particle types ----
        self._primary_type = _particles.resolve(primary)
        self._target_type = (
            _particles.resolve(target) if target is not None else None
        )

        # ---- Resolve detector ----
        if isinstance(detector, str):
            self._detector_model = _utilities.load_detector(detector)
            self._detector_name = detector
        else:
            self._detector_model = detector
            self._detector_name = None

        # ---- Resolve interactions ----
        self._secondary_processes = {}
        self._init_interactions(interactions, process, isoscalar, process_kwargs)

        # ---- Resolve secondary interactions ----
        if secondary_interactions is not None:
            resolved = {}
            for k, v in secondary_interactions.items():
                resolved[_particles.resolve(k) if isinstance(k, str) else k] = v
            self._secondary_processes = resolved

        # ---- Resolve distributions ----
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
        )

        # ---- Resolve secondary position ----
        self._secondary_injection_distributions = {}
        if secondary_position is not None:
            self._resolve_secondary_position(secondary_position)

        # ---- Resolve secondary biasing ----
        self._secondary_phase_spaces = {}  # {InteractionSignature: MultiChannelPhaseSpace}
        if bias_targets is not None:
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
        self._n_events = n_events
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
                isoscalar=isoscalar,
            )
            if self._target_type is not None:
                load_kwargs["target_types"] = [self._target_type]
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
    ):
        """Resolve the injection and physical distribution lists.

        Rules:
        - ``energy`` sets both injection and physical energy.
        - ``injection_energy`` / ``physical_energy`` override individually.
        - ``flux`` is an alias for ``physical_energy``.  When combined with
          ``energy``, ``energy`` becomes injection-only.
        - ``direction`` sets both.  ``injection_direction`` /
          ``physical_direction`` override individually.
        - ``position`` is injection-only (never in physical).
        - ``mass`` is injection-only (auto-inferred as 0 if omitted).
        """
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
            if phys_energy is not None:
                # energy + flux: energy is injection, flux is physical
                inj_energy = energy
            else:
                inj_energy = energy
                phys_energy = energy

        if inj_energy is None:
            raise ValueError(
                "No energy distribution specified. Provide 'energy', "
                "or 'injection_energy' + 'physical_energy'."
            )
        if phys_energy is None:
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

        if inj_dir is None:
            raise ValueError(
                "No direction distribution specified. Provide 'direction', "
                "or 'injection_direction' + 'physical_direction'."
            )
        if phys_dir is None:
            phys_dir = inj_dir

        # ---- Position (injection only) ----
        if position is None:
            raise ValueError(
                "No position distribution specified. Provide 'position' "
                "(any subclass of VertexPositionDistribution)."
            )

        # ---- Mass (injection only, auto-inferred) ----
        mass_val = mass if mass is not None else 0.0
        mass_dist = _distributions.PrimaryMass(mass_val)

        # ---- Assemble lists ----
        self._injection_distributions = [mass_dist, inj_energy, inj_dir, position]
        self._physical_distributions = [phys_energy, phys_dir]

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
        target, sig, interaction, daughter_type, spectator_type=None
    ):
        """Build a MultiChannelPhaseSpace for one signature.

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

        Returns
        -------
        MultiChannelPhaseSpace
        """
        n_sec = len(sig.secondary_types)
        daughter_idx = Simulation._find_daughter_index(sig, daughter_type)

        channels = []
        weights = []

        # Physical channel (fallback for events that miss the target)
        if isinstance(interaction, _interactions.Decay):
            channels.append(_injection.PhysicalDecayChannel(interaction))
        else:
            channels.append(
                _injection.PhysicalCrossSectionChannel(interaction))
        weights.append(0.01)

        if n_sec == 2:
            channels.append(
                _injection.DetectorDirected2BodyChannel(
                    target, daughter_idx))
            weights.append(0.99)

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
                    target,
                    spectator_idx,
                    pair_indices[0],
                    pair_indices[1],
                    directed_pair_idx,
                ))
            weights.append(0.99)

        elif isinstance(interaction, _interactions.CrossSection):
            channels.append(
                _injection.DetectorDirectedScatteringChannel(
                    target, daughter_idx))
            weights.append(0.99)

        # Normalize
        total = sum(weights)
        weights = [w / total for w in weights]

        mc = _injection.MultiChannelPhaseSpace()
        mc.channels = channels
        mc.weights = weights
        return mc

    def _resolve_bias_targets(self, bias_targets, daughter_type,
                              spectator_type):
        """Build per-signature phase spaces for direction biasing.

        Parameters
        ----------
        bias_targets : Geometry or dict
            If a Geometry, biases all secondary signatures toward it.
            If a dict, maps ``{InteractionSignature: MultiChannelPhaseSpace}``
            for full control.
        daughter_type : ParticleType
            Which daughter to direct toward the target.
        spectator_type : ParticleType or None
            For 3-body: which daughter is the spectator.
        """
        if isinstance(bias_targets, dict):
            for k, v in bias_targets.items():
                if isinstance(v, _injection.MultiChannelPhaseSpace):
                    self._secondary_phase_spaces[k] = v
                # else: v is a Geometry, k is a signature
                # (not yet supported, could add)
            return

        # Single geometry: build per-signature for all secondaries
        target = bias_targets
        for sec_type, interactions_list in self._secondary_processes.items():
            decays = [x for x in interactions_list
                      if isinstance(x, _interactions.Decay)]
            cross_sections = [x for x in interactions_list
                              if isinstance(x, _interactions.CrossSection)]

            for decay in decays:
                for sig in decay.GetPossibleSignatures():
                    try:
                        mc = self._build_phase_space_for_signature(
                            target, sig, decay,
                            daughter_type, spectator_type,
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

            for xs in cross_sections:
                for sig in xs.GetPossibleSignatures():
                    try:
                        mc = self._build_phase_space_for_signature(
                            target, sig, xs,
                            daughter_type, spectator_type,
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
        """Construct and return an Injector from stored config."""
        validate_injection_distributions(self._injection_distributions)

        injector = _injection.Injector()
        injector.number_of_events = self._n_events
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

    def run(self):
        """Generate events and compute weights.

        Returns
        -------
        Results
            A :class:`Results` object containing events, weights,
            and generation times.
        """
        injector = self._build_injector()
        events, gen_times = _GenerateEvents(injector)
        weighter = self._build_weighter(injector)
        weights = [weighter(event) for event in events]

        self._injector = injector
        self._weighter = weighter
        self._last_events = events
        self._last_gen_times = gen_times

        return Results(events, weights, gen_times, weighter, injector)

    # ------------------------------------------------------------------ #
    #  Reweight                                                            #
    # ------------------------------------------------------------------ #

    def reweight(
        self,
        *,
        physical_energy=None,
        physical_direction=None,
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

        # Build new physical distribution list
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
            weighter, self._injector,
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

    def __repr__(self):
        parts = [
            f"n_events={self._n_events}",
            f"primary={self._primary_type}",
            f"detector={'<loaded>' if self._detector_name is None else self._detector_name!r}",
        ]
        return f"Simulation({', '.join(parts)})"
