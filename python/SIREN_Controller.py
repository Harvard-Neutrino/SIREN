import os
import h5py
import numpy as np
import awkward as ak
import time

from . import utilities as _utilities
from . import detector as _detector
from . import injection as _injection
from . import distributions as _distributions
from . import dataclasses as _dataclasses
from . import interactions as _interactions
from . import geometry as _geometry
from . import math as _math

from . import _util

from . import darknews_version
if darknews_version() is not None:
    from .SIREN_DarkNews import PyDarkNewsInteractionCollection

# Helper functions

# attempts to merge multiple interaction collections
def MergeInteractionCollections(primary_type,int_col_list):
    cross_sections = []
    decays = []
    record = _dataclasses.InteractionRecord()
    record.signature.primary_type = primary_type
    for int_col in int_col_list:
        assert(int_col.MatchesPrimary(record))
        if int_col.HasCrossSections():
            cross_sections += list(int_col.GetCrossSections())
        if int_col.HasDecays():
            decays += list(int_col.GetDecays())
    return _interactions.InteractionCollection(primary_type, cross_sections, decays)



# Parent python class for handling event generation and weighting
class SIREN_Controller:

    def __init__(self, events_to_inject, experiment=None, detector_model_file=None, materials_model_file=None, seed=0, detector_model=None):
        """
        SIREN controller class constructor.
        :param int event_to_inject: number of events to generate
        :param str experiment: experiment name in string (default None)
        :param str detector_model_file: path to the detector model file (default None)
        :param str materials_model_file: path to the materials model file (default None)
        :param int seed: Optional random number generator seed (default 0)
        :param detector_model: a pre-built ``DetectorModel`` to use directly (e.g.
            a GDML composite from ``load_detector("SBN", detector=...)``); when
            provided it overrides the experiment/file loading (default None)
        """

        self.global_start = time.time()

        self.resources_dir = _util.resource_package_dir()

        # Initialize a random number generator
        self.random = _utilities.SIREN_random(seed)

        # Save number of events to inject
        self.events_to_inject = events_to_inject
        self.experiment = experiment

        # Empty list for our interaction trees
        self.events = []


        self.detector_model_file = detector_model_file
        self.materials_model_file = materials_model_file
        if detector_model is not None:
            # Use a pre-built DetectorModel directly (e.g. a GDML composite from
            # load_detector("SBN", detector=...)). It supersedes any
            # experiment/file configuration, so clear the file paths -- there
            # are no densities/materials files for file-derived helpers like
            # GetFiducialVolume() to consult.
            self.detector_model = detector_model
            self.detector_model_file = None
            self.materials_model_file = None
        else:
            if experiment is not None:
                # Find the density and materials files
                detector_dir = _util.get_detector_model_path(experiment)
                self.materials_model_file = os.path.join(detector_dir, "materials.dat")
                self.detector_model_file = os.path.join(detector_dir, "densities.dat")
            elif (self.detector_model_file is None or self.materials_model_file is None):
                raise ValueError("Must provide either an experiment name, both a detector model file and materials model file, or a pre-built detector_model")

            self.detector_model = _detector.DetectorModel()
            self.detector_model.LoadMaterialModel(self.materials_model_file)
            self.detector_model.LoadDetectorModel(self.detector_model_file)

        # Define the primary injection and physical process
        self.primary_injection_process = _injection.PrimaryInjectionProcess()
        self.primary_physical_process = _injection.PhysicalProcess()

        # Define lists for the secondary injection and physical processes
        self.secondary_injection_processes = []
        self.secondary_physical_processes = []

        # Set the fiducial volume
        self.fid_vol = self.GetFiducialVolume()

    def GetDetectorSectorGeometry(self, sector_name):
        for sector in self.detector_model.Sectors:
            if sector.name==sector_name:
                return sector.geo
        return None

    def SetInjectionProcesses(
        self,
        primary_type,
        primary_injection_distributions,
        secondary_types=[],
        secondary_injection_distributions=[],
        fid_vol_secondary=True,
    ):
        """
        SIREN injection process setter.
        :param ParticleType primary_type: The primary particle being generated
        :param dict<str,InjectionDistribution> primary_injection_distributions: The dict of injection distributions for the primary process
        :param list<ParticleType> secondary_types: The secondary particles being generated
        :param list<dict<str,InjectionDistribution> secondary_injection_distributions: List of dict of injection distributions for each secondary process
        :param bool fid_vol_secondary: whether to restrict secondary interactions to fiducial volume
        """

        # Define the primary injection process primary type
        self.primary_injection_process.primary_type = primary_type

        # Default injection distributions. The pybind Process exposes a
        # `distributions` list property (the old Add*Distribution methods were
        # removed upstream), so assemble the full list and assign it once.
        primary_idist_list = []
        if "mass" not in primary_injection_distributions.keys():
            primary_idist_list.append(_distributions.PrimaryMass(0))
        if "helicity" not in primary_injection_distributions.keys():
            primary_idist_list.append(_distributions.PrimaryNeutrinoHelicityDistribution())

        # Add all injection distributions
        for _, idist in primary_injection_distributions.items():
            primary_idist_list.append(idist)
        self.primary_injection_process.distributions = primary_idist_list

        # Loop through possible secondary interactions
        for i_sec, secondary_type in enumerate(secondary_types):
            secondary_injection_process = _injection.SecondaryInjectionProcess()
            secondary_injection_process.primary_type = secondary_type

            sec_idist_list = list(secondary_injection_distributions[i_sec])
            # Add the position distribution
            if fid_vol_secondary and self.fid_vol is not None:
                sec_idist_list.append(_distributions.SecondaryBoundedVertexDistribution(self.fid_vol))
            else:
                sec_idist_list.append(_distributions.SecondaryPhysicalVertexDistribution())
            secondary_injection_process.distributions = sec_idist_list

            self.secondary_injection_processes.append(secondary_injection_process)

    def SetPhysicalProcesses(
        self,
        primary_type,
        primary_physical_distributions,
        secondary_types=[],
        secondary_physical_distributions=[],
    ):
        """
        SIREN physical process setter.
        :param ParticleType primary_type: The primary particle being generated
        :param dict<str,PhysicalDistribution> primary_physical_distributions: The dict of physical distributions for the primary process
        :param list<ParticleType> secondary_types: The secondary particles being generated
        :param list<dict<str,PhysicalDistribution> secondary_physical_distributions: List of dict of physical distributions for each secondary process
        """

        # Define the primary physical process primary type
        self.primary_physical_process.primary_type = primary_type

        # Default physical distributions (assign the `distributions` list, as for
        # the injection processes above).
        primary_pdist_list = []
        if "mass" not in primary_physical_distributions.keys():
            primary_pdist_list.append(_distributions.PrimaryMass(0))
        if "helicity" not in primary_physical_distributions.keys():
            primary_pdist_list.append(_distributions.PrimaryNeutrinoHelicityDistribution())

        # Add all physical distributions
        for _, pdist in primary_physical_distributions.items():
            primary_pdist_list.append(pdist)
        self.primary_physical_process.distributions = primary_pdist_list

        # Loop through possible secondary interactions
        for i_sec, secondary_type in enumerate(secondary_types):
            secondary_physical_process = _injection.PhysicalProcess()
            secondary_physical_process.primary_type = secondary_type
            secondary_physical_process.distributions = list(secondary_physical_distributions[i_sec])
            self.secondary_physical_processes.append(secondary_physical_process)

    def SetProcesses(
        self,
        primary_type,
        primary_injection_distributions,
        primary_physical_distributions,
        secondary_types=[],
        secondary_injection_distributions=[],
        secondary_physical_distributions=[],
        fid_vol_secondary=True,
    ):
        """
        SIREN process setter.
        :param ParticleType primary_type: The primary particle being generated
        :param dict<str,InjectionDistribution> primary_injection_distributions: The dict of injection distributions for the primary process
        :param dict<str,PhysicalDistribution> primary_physical_distributions: The dict of physical distributions for the primary process
        :param list<ParticleType> secondary_types: The secondary particles being generated
        :param list<dict<str,InjectionDistribution> secondary_injection_distributions: List of dict of injection distributions for each secondary process
        :param list<dict<str,PhysicalDistribution> secondary_physical_distributions: List of dict of physical distributions for each secondary process
        :param bool fid_vol_secondary: whether to restrict secondary interactions to fiducial volume
        """
        self.SetInjectionProcesses(primary_type,primary_injection_distributions,secondary_types,secondary_injection_distributions,fid_vol_secondary)
        self.SetPhysicalProcesses(primary_type,primary_physical_distributions,secondary_types,secondary_physical_distributions)

    def InputDarkNewsModel(self, primary_type, table_dir, upscattering=True, decay=True, fill_tables_at_start=False, Emax=None, fid_vol_secondary=True, **kwargs):
        """
        Sets up the relevant processes and cross section/decay objects related to a provided DarkNews model dictionary.
        Will handle the primary cross section collection as well as the entire list of secondary processes

        :param _dataclasses.Particle.ParticleType primary_type: primary particle to be generated
        :param string table_dir: Directory for storing cross section and decay tables
        :param bool upscattering: Flag to set upscattering interactions
        :param bool decay: Flag to set decay interactions
        :param bool fill_tables_at_start: Flag to fill total/differential cross section tables upon initialization
        :param float Emax: maximum energy for cross section tables
        :param dict<str,val> kwargs: The dict of DarkNews model and cross section parameters
        """
        # Add nuclear targets to the model arguments
        kwargs["nuclear_targets"] = self.GetDetectorModelTargets()[1]
        # Initialize DarkNews cross sections and decays
        self.DN_processes = PyDarkNewsInteractionCollection(
            table_dir=table_dir, **kwargs
        )

        if fill_tables_at_start:
            if Emax is None:
                print("WARNING: Cannot fill cross section tables without specifying a maximum energy")
            else:
                self.DN_processes.FillCrossSectionTables(Emax=Emax)

        # Initialize primary InteractionCollection
        # Loop over available cross sections and save those which match primary type
        primary_cross_sections = []
        for cross_section in self.DN_processes.cross_sections:
            if primary_type == _dataclasses.Particle.ParticleType(
                cross_section.ups_case.nu_projectile.pdgid
            ):
                primary_cross_sections.append(cross_section)
        primary_interaction_collection = _interactions.InteractionCollection(
            primary_type, primary_cross_sections
        )

        # Initialize secondary processes and define secondary InteractionCollection objects
        secondary_decays = {}
        # Also keep track of the minimum decay width for defining the position distribution later
        self.DN_min_decay_width = np.inf
        # Loop over available decays, group by parent type
        for decay in self.DN_processes.decays:
            secondary_type = _dataclasses.Particle.ParticleType(
                decay.dec_case.nu_parent.pdgid
            )
            if secondary_type not in secondary_decays.keys():
                secondary_decays[secondary_type] = []
            secondary_decays[secondary_type].append(decay)
            total_decay_width = decay.TotalDecayWidth(secondary_type)
            if total_decay_width < self.DN_min_decay_width:
                self.DN_min_decay_width = total_decay_width
        # Now make the list of secondary cross section collections
        # Add new secondary injection and physical processes at the same time
        secondary_interaction_collections = []
        for secondary_type, decay_list in secondary_decays.items():

            # Define a secondary injection distribution if necessary
            inj_sec_defined = False
            phys_sec_defined = False
            existing_inj_process = None
            for secondary_injection_process in self.secondary_injection_processes:
                if secondary_injection_process.primary_type == secondary_type:
                    inj_sec_defined = True
                    existing_inj_process = secondary_injection_process
            for secondary_physical_process in self.secondary_physical_processes:
                if secondary_physical_process.primary_type == secondary_type:
                    phys_sec_defined = True

            if inj_sec_defined:
                secondary_injection_process = existing_inj_process
            else:
                secondary_injection_process = _injection.SecondaryInjectionProcess()
                secondary_injection_process.primary_type = secondary_type

            # Add the secondary position distribution (append to whatever the
            # process already carries; the pybind `distributions` is a list property).
            sec_dists = list(secondary_injection_process.distributions)
            if fid_vol_secondary and self.fid_vol is not None:
                sec_dists.append(_distributions.SecondaryBoundedVertexDistribution(self.fid_vol))
            else:
                sec_dists.append(_distributions.SecondaryPhysicalVertexDistribution())
            secondary_injection_process.distributions = sec_dists

            if not inj_sec_defined:
                self.secondary_injection_processes.append(secondary_injection_process)
            if not phys_sec_defined:
                secondary_physical_process = _injection.PhysicalProcess()
                secondary_physical_process.primary_type = secondary_type
                self.secondary_physical_processes.append(secondary_physical_process)

            secondary_interaction_collections.append(
                _interactions.InteractionCollection(secondary_type, decay_list)
            )

        # check whether to not include either process
        if not upscattering: primary_interaction_collection = None
        if not decay: secondary_interaction_collections = None
        self.SetInteractions(
            primary_interaction_collection=primary_interaction_collection,
            secondary_interaction_collections=secondary_interaction_collections
        )

    def InputDarkNewsDecay(self, primary_type, table_dir, **kwargs):
        """
        Sets up the relevant processes and decay objects related to a provided DarkNews model dictionary.
        Will handle the primary decay collection

        :param _dataclasses.Particle.ParticleType primary_type: primary particle to be generated
        :param string table_dir: Directory for storing cross section and decay tables
        :param dict<str,val> kwargs: The dict of DarkNews model and cross section parameters
        """
        # Add nuclear targets to the model arguments
        kwargs["nuclear_targets"] = self.GetDetectorModelTargets()[1]
        # Initialize DarkNews cross sections and decays
        self.DN_processes = PyDarkNewsInteractionCollection(
            table_dir=table_dir, **kwargs
        )

        # Initialize primary InteractionCollection
        # Loop over available cross sections and save those which match primary type
        primary_decays = []
        # Also keep track of the minimum decay width for defining the position distribution later
        self.DN_min_decay_width = np.inf
        for decay in self.DN_processes.decays:
            if primary_type == _dataclasses.Particle.ParticleType(
                decay.dec_case.nu_parent.pdgid
            ):
                primary_decays.append(decay)
                total_decay_width = decay.TotalDecayWidth(primary_type)
                if total_decay_width > 0 and total_decay_width < self.DN_min_decay_width:
                    self.DN_min_decay_width = total_decay_width
        primary_interaction_collection = _interactions.InteractionCollection(
            primary_type, primary_decays
        )

        self.SetInteractions(
            primary_interaction_collection=primary_interaction_collection
        )

    def GetFiducialVolume(self):
        """
        :return: identified fiducial volume for the experiment, None if not found
        """
        # A pre-built detector model (e.g. a GDML composite) has no
        # densities.dat to parse a fiducial line from.
        if self.detector_model_file is None:
            return None
        with open(self.detector_model_file) as file:
            fiducial_line = None
            detector_line = None
            for line in file:
                data = line.split()
                if len(data) <= 0:
                    continue
                elif data[0] == "fiducial":
                    fiducial_line = line
                elif data[0] == "detector":
                    detector_line = line
            if fiducial_line is None or detector_line is None:
                return None
            return _detector.DetectorModel.ParseFiducialVolume(fiducial_line, detector_line)
        return None

    def GetVolumePositionDistributionFromSector(self, sector_name):
        geo = self.GetDetectorSectorGeometry(sector_name)
        if geo is None:
            raise ValueError("Sector %s not found" % sector_name)
        # the position is in geometry coordinates
        # must update to detector coordintes
        det_position = self.detector_model.GeoPositionToDetPosition(_detector.GeometryPosition(geo.placement.Position))
        det_rotation = geo.placement.Quaternion
        det_placement = _geometry.Placement(det_position.get(), det_rotation)
        if type(geo)==_geometry.Cylinder:
            cylinder = _geometry.Cylinder(det_placement,geo.Radius,geo.InnerRadius,geo.Z)
            return _distributions.CylinderVolumePositionDistribution(cylinder)
        elif type(geo)==_geometry.Sphere:
            sphere = _geometry.Sphere(det_placement,geo.Radius,geo.InnerRadius)
            return _distributions.SphereVolumePositionDistribution(sphere)
        else:
            raise TypeError("Geometry type %s not supported for position distribution" % str(type(geo)))

    def GetDetectorModelTargets(self):
        """
        Determines the targets that exist inside the detector model
        :return: lists of targets and strings
        :rtype: (list<ParticleType>, list<str>)
        """
        count = 0
        targets = []
        target_strs = []
        while self.detector_model.Materials.HasMaterial(count):
            for _target in self.detector_model.Materials.GetMaterialTargets(count):
                if _target not in targets:
                    targets.append(_target)
                if str(_target).find("Nucleus") == -1:
                    continue
                else:
                    target_str = str(_target)[
                        str(_target).find("Type") + 5 : str(_target).find("Nucleus")
                    ]
                    if target_str == "H":
                        target_str = "H1"
                    if target_str not in target_strs:
                        target_strs.append(target_str)
            count += 1
        return targets, target_strs

    def SetInteractions(
        self, primary_interaction_collection=None, secondary_interaction_collections=None, injection=True, physical=True
    ):
        """
        Set cross sections for the primary and secondary processes
        If cross sections already exist for either, attempts to merge the interaction collections
        :param InteractionCollection primary_interaction_collection: The cross section collection for the primary process
        :param list<InteractionCollection> secondary_interaction_collections: The list of cross section collections for the primary process
        :param bool injection: whether to apply these interaction collections to the injection processes
        :param bool physical: whether to apply these interaction collections to the physical processes
        """

        if primary_interaction_collection is not None:
            # Set primary cross sections
            if injection:
                # set the primary type if doesn't exist
                if self.primary_injection_process.primary_type == _dataclasses.Particle.ParticleType.unknown:
                    self.primary_injection_process.primary_type = primary_interaction_collection.GetPrimaryType()
                else:
                    assert(self.primary_injection_process.primary_type == primary_interaction_collection.GetPrimaryType())
                if self.primary_injection_process.interactions is None:
                    self.primary_injection_process.interactions = primary_interaction_collection
                else:
                    self.primary_injection_process.interactions = MergeInteractionCollections(self.primary_injection_process.primary_type,
                                                                                            [self.primary_injection_process.interactions, primary_interaction_collection])
            if physical:
                if self.primary_physical_process.primary_type == _dataclasses.Particle.ParticleType.unknown:
                    self.primary_physical_process.primary_type = primary_interaction_collection.GetPrimaryType()
                else:
                    assert(self.primary_physical_process.primary_type == primary_interaction_collection.GetPrimaryType())
                if self.primary_physical_process.interactions is None:
                    self.primary_physical_process.interactions = primary_interaction_collection
                else:
                    self.primary_physical_process.interactions = MergeInteractionCollections(self.primary_physical_process.primary_type,
                                                                                            [self.primary_physical_process.interactions, primary_interaction_collection])

        if secondary_interaction_collections is not None:
            # Loop through secondary processes
            for sec_inj, sec_phys in zip(
                self.secondary_injection_processes, self.secondary_physical_processes
            ):
                assert sec_inj.primary_type == sec_phys.primary_type
                record = _dataclasses.InteractionRecord()
                record.signature.primary_type = sec_inj.primary_type
                found_collection = False
                # Loop through possible seconday cross sections
                for sec_ints in secondary_interaction_collections:
                    # Match cross section collection on the primary type
                    if sec_ints.MatchesPrimary(record):
                        # Set secondary cross sections
                        if injection:
                            if sec_inj.interactions is None:
                                sec_inj.interactions = sec_ints
                            else:
                                sec_inj.interactions = MergeInteractionCollections(sec_inj.primary_type,
                                                                                [sec_inj.interactions, sec_ints])
                        if physical:
                            if sec_phys.interactions is None:
                                sec_phys.interactions = sec_ints
                            else:
                                sec_phys.interactions = MergeInteractionCollections(sec_phys.primary_type,
                                                                                    [sec_phys.interactions, sec_ints])
                        found_collection = True
                if not found_collection and(sec_inj.interactions is None or sec_phys.interactions is None):
                    raise RuntimeError(
                        "Couldn't find cross section collection for secondary particle %s"
                        % record.signature.primary_type
                    )

    # set the stopping condition of the injector with a python function
    # must accept two arguments, assumes first is datum and the second is the index of the secondary particle
    def SetInjectorStoppingCondition(self, stopping_condition):
        self.injector.SetStoppingCondition(stopping_condition)

    # Initialize the injector, either from an existing .siren_injector file or from controller injection objects
    def InitializeInjector(self, filenames=None):
        if type(filenames) == str:
            if os.path.isfile(filenames):
                filenames = [filenames]
        self.injectors=[]
        filenames = None
        if filenames is None:
            assert(self.primary_injection_process.primary_type is not None)
            # Use controller injection objects
            self.injectors.append(
                _injection._Injector(
                    self.events_to_inject,
                    self.detector_model,
                    self.primary_injection_process,
                    self.secondary_injection_processes,
                    self.random,
                )
            )
        else:
            # Try initilalizing with the provided filenames
            assert(len(filenames)>0) # require at least one injector filename
            for filename in filenames:
                self.injectors.append(
                    _injection._Injector(
                        self.events_to_inject,
                        filename,
                        self.random,
                    )
                )
                self.injectors[-1].ResetInjectedEvents()
        self.injector = self.injectors[0] # presume that injection happens with only the first provided injector

    # Initialize the weighter, either from an existing .siren_weighter file or from controller injection objects
    def InitializeWeighter(self,filename=None):
        if filename is None:
            assert(self.primary_physical_process.primary_type is not None)
            # Use controller physical objects
            self.weighter = _injection._Weighter(
                self.injectors,
                self.detector_model,
                self.primary_physical_process,
                self.secondary_physical_processes,
            )
        else:
            # Try initilalizing with the provided filename
            self.weighter = _injection._Weighter(
                self.injectors,
                filename
            )

    # Initialize the injector and weighter objects
    # Use existing .siren_injector and/or .siren_weighter files if they exist
    def Initialize(self, injection_filenames=None, weighter_filename=None):

        # Define the injector object(s)
        self.InitializeInjector(filenames=injection_filenames)

        # Define the weighter object
        self.InitializeWeighter(filename=weighter_filename)

    # Generate events using the self.injector object
    def GenerateEvents(self, N=None, fill_tables_at_exit=True, verbose=True):
        if N is None:
            N = self.events_to_inject
        count = 0
        self.gen_times,self.global_times = [],[]
        prev_time = time.time()
        while (self.injector.InjectedEvents() < self.events_to_inject) and (count < N):
            if verbose: print("Injecting Event %d/%d  " % (count, N), end="\r")
            event = self.injector.GenerateEvent()
            self.events.append(event)
            t = time.time()
            self.gen_times.append(t-prev_time)
            self.global_times.append(t-self.global_start)
            prev_time = t
            count += 1
        if hasattr(self, "DN_processes"):
            self.DN_processes.SaveCrossSectionTables(fill_tables_at_exit=fill_tables_at_exit)
        return self.events

    # Load events from the custom SIREN event format
    def LoadEvents(self, filename):
        self.events = _dataclasses.LoadInteractionTrees(filename)
        self.gen_times = np.zeros_like(self.events)
        self.global_times = np.zeros_like(self.events)

    # Load events from a HepMC3 file written by SIREN
    def LoadEventsFromHepMC3(self, filename):
        from . import hepmc3 as _hepmc3
        self.events = _hepmc3.LoadInteractionTreesFromHepMC3(filename)
        self.gen_times = np.zeros_like(self.events)
        self.global_times = np.zeros_like(self.events)

    # Save events to hdf5, parquet, and/or custom SIREN filetypes
    # if the weighter exists, calculate the event weight too
    def SaveEvents(self, filename, fill_tables_at_exit=True,
                   hdf5=True, parquet=True, siren_events=True, # filetypes to save events
                   save_int_probs=False,save_int_params=False,save_survival_probs=False,
                   verbose=True, hepmc3=False, hepmc3_weights="auto", hepmc3_gzip=False):

        # Resolve the HepMC3 weight policy up front so per-event central values are
        # computed exactly once (BEFORE any file is written) and reused by the native
        # file, the HepMC3 file, and the tabular datasets below. Under "auto" the
        # controller's weighter (when present) populates the tree headers.
        hepmc3_cv = None
        hepmc3_state = "unweighted"
        if hepmc3 or siren_events:
            hepmc3_cv, hepmc3_state = _util.resolve_hepmc3_weight_policy(
                self.events, hepmc3_weights, getattr(self, "weighter", None))

        if siren_events:
            _dataclasses.SaveInteractionTrees(self.events, filename)
        if hepmc3:
            from . import hepmc3 as _hepmc3
            opts = _hepmc3.HepMC3WriterOptions()
            opts.weights_state = hepmc3_state
            opts.gzip = bool(hepmc3_gzip)
            # Store the generation counts (from the injector) as run metadata; used
            # to normalize the flux-averaged cross section. Absent after LoadEvents.
            if getattr(self, "injector", None) is not None:
                opts.attempted_events = int(self.injector.InjectionAttempts())
                opts.accepted_events = int(self.injector.InjectedEvents())
            out = filename + ".hepmc3"
            if hepmc3_gzip and not out.endswith(".gz"):
                out = out + ".gz"
            _hepmc3.SaveInteractionTreesAsHepMC3(self.events, out, opts)
        # A dictionary containing each dataset we'd like to save
        datasets = {
            "event_weight":[], # weight of entire event
            "event_gen_time":[], # generation time of each event
            "event_weight_time":[], # weight calculation time of each event
            "event_global_time":[], # global time of each event
            "num_interactions":[], # number of interactions per event
            "vertex":[], # vertex of each interaction of an event
            "vertex_time":[], # time of each interaction vertex of an event
            "primary_initial_position":[], # initial position of primary in each interaction of an event
            "primary_initial_time":[], # initial time of primary in each interaction of an event
            "in_fiducial":[], # whether or not each vertex is in the fiducial volume
            "primary_type":[], # primary type of each interaction
            "target_type":[], # target type of each interaction
            "num_secondaries":[], # number of secondary particles of each interaction
            "secondary_types":[], # secondary type of each interaction
            "primary_momentum":[], # primary momentum of each interaction
            "secondary_momenta":[], # secondary momentum of each interaction
            "secondary_times":[], # production time of each secondary of each interaction
            "parent_idx":[], # index of the parent interaction
            "num_daughters":[], # number of daughter interactions
        }
        if save_int_probs: datasets["int_probs"] = []
        if save_survival_probs: datasets["survival_probs"] = []
        for ie, event in enumerate(self.events):
            if verbose: print("Saving Event %d/%d  " % (ie, len(self.events)), end="\r")
            t0 = time.time()
            if hepmc3_cv is not None:
                datasets["event_weight"].append(hepmc3_cv[ie])  # reuse the CV computed above
            else:
                datasets["event_weight"].append(self.weighter.EventWeight(event) if hasattr(self,"weighter") else 0)
            if save_int_probs:
                datasets["int_probs"].append(self.weighter.GetInteractionProbabilities(event)if hasattr(self,"weighter") else [])
            if save_survival_probs:
                datasets["survival_probs"].append(self.weighter.GetSurvivalProbabilities(event)if hasattr(self,"weighter") else [])
            datasets["event_weight_time"].append(time.time()-t0)
            datasets["event_gen_time"].append(self.gen_times[ie])
            datasets["event_global_time"].append(self.global_times[ie])
            # add empty lists for each per interaction dataset
            for k in ["vertex",
                      "vertex_time",
                      "primary_initial_position",
                      "primary_initial_time",
                      "in_fiducial",
                      "primary_type",
                      "target_type",
                      "num_secondaries",
                      "secondary_types",
                      "primary_momentum",
                      "secondary_momenta",
                      "secondary_times",
                      "parent_idx",
                      "num_daughters"]:
                datasets[k].append([])
            if save_int_params:
                datasets.setdefault("int_params", [])
                datasets["int_params"].append({})
            # parent index of each interaction, taken from the tree's parent edges
            parent_indices = _util.get_parent_indices(event.tree)
            # loop over interactions
            for id, datum in enumerate(event.tree):
                if save_int_params:
                    for param_name, param_value in datum.record.interaction_parameters.items():
                        if param_name not in datasets["int_params"][-1]:
                            datasets["int_params"][-1][param_name] = []
                        datasets["int_params"][-1][param_name].append(param_value)
                datasets["vertex"][-1].append(np.array(datum.record.interaction_vertex,dtype=float))
                datasets["vertex_time"][-1].append(datum.record.interaction_time)
                datasets["primary_initial_position"][-1].append(np.array(datum.record.primary_initial_position,dtype=float))
                datasets["primary_initial_time"][-1].append(datum.record.primary_initial_time)

                 # primary particle stuff
                datasets["primary_type"][-1].append(int(datum.record.signature.primary_type))
                datasets["primary_momentum"][-1].append(np.array(datum.record.primary_momentum, dtype=float))
                datasets["num_daughters"][-1].append(len(datum.daughter_indices))

                # parent interaction index (from the tree's parent/daughter edges)
                datasets["parent_idx"][-1].append(parent_indices[id])

                if self.fid_vol is not None:
                    pos = _math.Vector3D(datasets["vertex"][-1][-1])
                    dir = _math.Vector3D(datasets["primary_momentum"][-1][-1][1:])
                    dir.normalize()
                    datasets["in_fiducial"][-1].append(self.fid_vol.IsInside(pos,dir))
                else:
                    datasets["in_fiducial"][-1].append(False)

                # target particle stuff
                datasets["target_type"][-1].append(int(datum.record.signature.target_type))

                # secondary particle stuff
                datasets["secondary_types"][-1].append([])
                datasets["secondary_momenta"][-1].append([])
                datasets["secondary_times"][-1].append([])
                # events loaded from old files predate secondary_times;
                # fall back to the vertex time for them
                sec_times = datum.record.secondary_times
                if len(sec_times) != len(datum.record.secondary_momenta):
                    sec_times = [datum.record.interaction_time] * len(datum.record.secondary_momenta)
                for isec, (sec_type, sec_momenta, sec_time) in enumerate(zip(datum.record.signature.secondary_types,
                                                                             datum.record.secondary_momenta,
                                                                             sec_times)):
                    datasets["secondary_types"][-1][-1].append(int(sec_type))
                    datasets["secondary_momenta"][-1][-1].append(np.array(sec_momenta,dtype=float))
                    datasets["secondary_times"][-1][-1].append(sec_time)
                datasets["num_secondaries"][-1].append(isec+1)
            datasets["num_interactions"].append(id+1)

        # save injector and weighter (writes <filename>.siren_injector and
        # <filename>.siren_weighter alongside the event file). These are the
        # pybind _Injector/_Weighter objects, so use their C++ serialization
        # methods (SaveInjector writes the literal path; SaveWeighter appends
        # the .siren_weighter suffix). The weighter is optional.
        self.injector.SaveInjector(filename + ".siren_injector")
        if hasattr(self, "weighter"):
            self.weighter.SaveWeighter(filename)

        # save events
        ak_array = ak.Array(datasets)
        if hdf5:
            fout = h5py.File(filename+".hdf5", "w")
            group = fout.create_group("Events")
            form, length, container = ak.to_buffers(ak.to_packed(ak_array), container=group)
            group.attrs["form"] = form.to_json()
            group.attrs["length"] = length
            fout.close()
        if parquet:
            ak.to_parquet(ak_array,filename+".parquet")

        # save darknews cross section tables
        if hasattr(self, "DN_processes"):
            self.DN_processes.SaveCrossSectionTables(fill_tables_at_exit=fill_tables_at_exit)
