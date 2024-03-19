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

from .LIDarkNews import PyDarkNewsInteractionCollection



# Parent python class for handling event generation
class LIController:
    def __init__(self, events_to_inject, experiment, seed=0):
        """
        LI class constructor.
        :param int event_to_inject: number of events to generate
        :param str experiment: experiment name in string
        :param int seed: Optional random number generator seed
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

        # Find the density and materials files
        materials_file = _util.get_material_model_path(experiment)
        detector_model_file = _util.get_detector_model_path(experiment)

        self.detector_model = _detector.DetectorModel()
        self.detector_model.LoadMaterialModel(materials_file)
        self.detector_model.LoadDetectorModel(detector_model_file)

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

    def SetProcesses(
        self,
        primary_type,
        primary_injection_distributions,
        primary_physical_distributions,
        secondary_types=[],
        secondary_injection_distributions=[],
        secondary_physical_distributions=[],
    ):
        """
        LI process setter.
        :param ParticleType primary_type: The primary particle being generated
        :param dict<str,InjectionDistribution> primary_injection_distributions: The dict of injection distributions for the primary process
        :param dict<str,PhysicalDistribution> primary_physical_distributions: The dict of physical distributions for the primary process
        :param list<ParticleType> secondary_types: The secondary particles being generated
        :param list<dict<str,InjectionDistribution> secondary_injection_distributions: List of dict of injection distributions for each secondary process
        :param list<dict<str,PhysicalDistribution> secondary_physical_distributions: List of dict of physical distributions for each secondary process
        """

        # Define the primary injection and physical process
        self.primary_injection_process.primary_type = primary_type
        self.primary_physical_process.primary_type = primary_type

        # Default injection distributions
        if "mass" not in primary_injection_distributions.keys():
            self.primary_injection_process.AddPrimaryInjectionDistribution(
                _distributions.PrimaryMass(0)
            )

        # Default injection distributions
        if "mass" not in primary_physical_distributions.keys():
            self.primary_physical_process.AddPhysicalDistribution(
                _distributions.PrimaryMass(0)
            )

        # Default injection distributions
        if "helicity" not in primary_injection_distributions.keys():
            self.primary_injection_process.AddPrimaryInjectionDistribution(
                _distributions.PrimaryNeutrinoHelicityDistribution()
            )

        # Default physical distributions
        if "helicity" not in primary_physical_distributions.keys():
            self.primary_physical_process.AddPhysicalDistribution(
                _distributions.PrimaryNeutrinoHelicityDistribution()
            )

        # Add all injection distributions
        for _, idist in primary_injection_distributions.items():
            self.primary_injection_process.AddPrimaryInjectionDistribution(idist)
        # Add all physical distributions
        for _, pdist in primary_physical_distributions.items():
            self.primary_physical_process.AddPhysicalDistribution(pdist)

        # Loop through possible secondary interactions
        for i_sec, secondary_type in enumerate(secondary_types):
            secondary_injection_process = _injection.SecondaryInjectionProcess()
            secondary_physical_process = _injection.PhysicalProcess()
            secondary_injection_process.primary_type = secondary_type
            secondary_physical_process.primary_type = secondary_type

            # Add all injection distributions
            for idist in secondary_injection_distributions[i_sec]:
                secondary_injection_process.AddSecondaryInjectionDistribution(idist)
            # Add all physical distributions
            for pdist in secondary_physical_distributions[i_sec]:
                secondary_physical_process.AddPhysicalDistribution(pdist)

            # Add the position distribution
            if self.fid_vol is not None:
                secondary_injection_process.AddSecondaryInjectionDistribution(
                    _distributions.SecondaryBoundedVertexDistribution(self.fid_vol)
                )
            else:
                secondary_injection_process.AddSecondaryInjectionDistribution(
                    _distributions.SecondaryPhysicalVertexDistribution()
                )

            self.secondary_injection_processes.append(secondary_injection_process)
            self.secondary_physical_processes.append(secondary_physical_process)

    

    
    def InputDarkNewsModel(self, primary_type, table_dir, fill_tables_at_start=False, Emax=None, **kwargs):
        """
        Sets up the relevant processes and cross section/decay objects related to a provided DarkNews model dictionary.
        Will handle the primary cross section collection as well as the entire list of secondary processes

        :param _dataclasses.Particle.ParticleType primary_type: primary particle to be generated
        :param string table_dir: Directory for storing cross section and decay tables
        :param string fill_tables_at_start: Flag to fill total/differential cross section tables upon initialization
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
            # Define a sedcondary injection distribution
            secondary_injection_process = _injection.SecondaryInjectionProcess()
            secondary_physical_process = _injection.PhysicalProcess()
            secondary_injection_process.primary_type = secondary_type
            secondary_physical_process.primary_type = secondary_type

            # Add the secondary position distribution
            if self.fid_vol is not None:
                secondary_injection_process.AddSecondaryInjectionDistribution(
                    _distributions.SecondaryBoundedVertexDistribution(self.fid_vol)
                )
            else:
                secondary_injection_process.AddSecondaryInjectionDistribution(
                    _distributions.SecondaryPhysicalVertexDistribution()
                )

            self.secondary_injection_processes.append(secondary_injection_process)
            self.secondary_physical_processes.append(secondary_physical_process)

            secondary_interaction_collections.append(
                _interactions.InteractionCollection(secondary_type, decay_list)
            )

        self.SetInteractions(
            primary_interaction_collection, secondary_interaction_collections
        )

    def GetFiducialVolume(self):
        """
        :return: identified fiducial volume for the experiment, None if not found
        """
        detector_model_file = _util.get_detector_model_path(self.experiment)
        with open(detector_model_file) as file:
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
        self, primary_interaction_collection, secondary_interaction_collections=None
    ):
        """
        Set cross sections for the primary and secondary processes
        :param InteractionCollection primary_interaction_collection: The cross section collection for the primary process
        :param list<InteractionCollection> secondary_interaction_collections: The list of cross section collections for the primary process
        """
        if secondary_interaction_collections is None:
            secondary_interaction_collections = []

        # Set primary cross sections
        self.primary_injection_process.interactions = primary_interaction_collection
        self.primary_physical_process.interactions = primary_interaction_collection

        # Loop through secondary processes
        for sec_inj, sec_phys in zip(
            self.secondary_injection_processes, self.secondary_physical_processes
        ):
            assert sec_inj.primary_type == sec_phys.primary_type
            record = _dataclasses.InteractionRecord()
            record.signature.primary_type = sec_inj.primary_type
            found_collection = False
            # Loop through possible seconday cross sections
            for sec_xs in secondary_interaction_collections:
                # Match cross section collection on  the primary type
                if sec_xs.MatchesPrimary(record):
                    sec_inj.interactions = sec_xs
                    sec_phys.interactions = sec_xs
                    found_collection = True
            if not found_collection:
                print(
                    "Couldn't find cross section collection for secondary particle %s; Exiting"
                    % record.primary_type
                )
                exit(0)

    def Initialize(self):
        # Define stopping condition
        # TODO: make this more general
        def StoppingCondition(datum, i):
            return True

        # Define the injector object
        self.injector = _injection.Injector(
            self.events_to_inject,
            self.detector_model,
            self.primary_injection_process,
            self.secondary_injection_processes,
            self.random,
        )

        self.injector.SetStoppingCondition(StoppingCondition)

        # Define the weighter object
        self.weighter = _injection.LeptonTreeWeighter(
            [self.injector],
            self.detector_model,
            self.primary_physical_process,
            self.secondary_physical_processes,
        )
    
    def GetCylinderVolumePositionDistributionFromSector(self, sector_name):
        geo = self.GetDetectorSectorGeometry(sector_name)
        if geo is None:
            print("Sector %s not found. Exiting"%sector_name)
            exit(0)
        # the position of this cylinder is in geometry coordinates
        # must update to detector coordintes
        det_position = self.detector_model.GeoPositionToDetPosition(_detector.GeometryPosition(geo.placement.Position))
        det_rotation = geo.placement.Quaternion
        det_placement = _geometry.Placement(det_position.get(), det_rotation)
        cylinder = _geometry.Cylinder(det_placement,geo.Radius,geo.InnerRadius,geo.Z)
        return _distributions.CylinderVolumePositionDistribution(cylinder)

    def GenerateEvents(self, N=None, fill_tables_at_exit=True):
        if N is None:
            N = self.events_to_inject
        count = 0
        self.gen_times,self.global_times = [],[]
        prev_time = time.time()
        while (self.injector.InjectedEvents() < self.events_to_inject) and (count < N):
            print("Injecting Event %d/%d  " % (count, N), end="\r")
            tree = self.injector.GenerateEvent()
            self.events.append(tree)
            t = time.time()
            self.gen_times.append(t-prev_time)
            self.global_times.append(t-self.global_start)
            prev_time = t
            count += 1
        if hasattr(self, "DN_processes"):
            self.DN_processes.SaveCrossSectionTables(fill_tables_at_exit=fill_tables_at_exit)
        return self.events

    def SaveEvents(self, filename, fill_tables_at_exit=True, hdf5=True, parquet=True):
        
        # A dictionary containing each dataset we'd like to save
        datasets = {
            "event_weight":[], # weight of entire event
            "event_gen_time":[], # generation time of each event
            "event_global_time":[], # global time of each event
            "num_interactions":[], # number of interactions per event
            "vertex":[], # vertex of each interaction in an event
            "in_fiducial":[], # whether or not each vertex is in the fiducial volume
            "primary_type":[], # primary type of each interaction
            "target_type":[], # target type of each interaction
            "num_secondaries":[], # number of secondary particles of each interaction
            "secondary_types":[], # secondary type of each interaction
            "primary_momentum":[], # primary momentum of each interaction
            "secondary_momenta":[], # secondary momentum of each interaction
        }
        for ie, event in enumerate(self.events):
            print("Saving Event %d/%d  " % (ie, len(self.events)), end="\r")
            datasets["event_weight"].append(self.weighter.EventWeight(event))
            datasets["event_gen_time"].append(self.gen_times[ie])
            datasets["event_global_time"].append(self.global_times[ie])
            # add empty lists for each per interaction dataset
            for k in ["vertex",
                      "in_fiducial",
                      "primary_type",
                      "target_type",
                      "num_secondaries",
                      "secondary_types",
                      "primary_momentum",
                      "secondary_momenta"]:
                datasets[k].append([])
            # loop over interactions
            for id, datum in enumerate(event.tree):
                
                datasets["vertex"][-1].append(np.array(datum.record.interaction_vertex,dtype=float))

                 # primary particle stuff
                datasets["primary_type"][-1].append(str(datum.record.signature.primary_type))
                datasets["primary_momentum"][-1].append(np.array(datum.record.primary_momentum, dtype=float))
                
                if self.fid_vol is not None:
                    pos = _math.Vector3D(datasets["vertex"][-1][-1])
                    dir = _math.Vector3D(datasets["primary_momentum"][-1][-1][1:])
                    dir.normalize()
                    datasets["in_fiducial"][-1].append(self.fid_vol.IsInside(pos,dir))
                else:
                    datasets["in_fiducial"][-1].append(False)

                # target particle stuff
                datasets["target_type"][-1].append(str(datum.record.signature.target_type))
                
                # secondary particle stuff
                datasets["secondary_types"][-1].append([])
                datasets["secondary_momenta"][-1].append([])
                for isec, (sec_type, sec_momenta) in enumerate(zip(datum.record.signature.secondary_types,
                                                                   datum.record.secondary_momenta)):
                    datasets["secondary_types"][-1][-1].append(str(sec_type))
                    datasets["secondary_momenta"][-1][-1].append(np.array(sec_momenta,dtype=float))
                datasets["num_secondaries"][-1].append(isec)
            datasets["num_interactions"].append(id)

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
        if hasattr(self, "DN_processes"):
            self.DN_processes.SaveCrossSectionTables(fill_tables_at_exit=fill_tables_at_exit)
