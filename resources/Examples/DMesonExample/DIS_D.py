import os

import siren
from siren.SIREN_Controller import SIREN_Controller

# Number of events to inject
events_to_inject = 5

# Expeirment to run
experiment = "IceCube"

# Define the controller
controller = SIREN_Controller(events_to_inject, experiment, seed = 1)

# Particle to inject
primary_type = siren.dataclasses.Particle.ParticleType.NuMu

cross_section_model = "CSMSDISSplines"

xsfiledir = siren.utilities.get_cross_section_model_path(cross_section_model)

# Cross Section Model
target_type = siren.dataclasses.Particle.ParticleType.Nucleon

DIS_xs = siren.interactions.CharmDISFromSpline(
    os.path.join(xsfiledir, "dsdxdy_nu_CC_iso.fits"),
    os.path.join(xsfiledir, "sigma_nu_CC_iso.fits"),
    [primary_type],
    [target_type], "m"
)

primary_xs = siren.interactions.InteractionCollection(primary_type, [DIS_xs])
controller.SetInteractions(primary_xs)

# Primary distributions
primary_injection_distributions = {}
primary_physical_distributions = {}

mass_dist = siren.distributions.PrimaryMass(0)
primary_injection_distributions["mass"] = mass_dist
primary_physical_distributions["mass"] = mass_dist

# energy distribution
edist = siren.distributions.PowerLaw(2, 1e5, 1e10)
primary_injection_distributions["energy"] = edist
primary_physical_distributions["energy"] = edist

# direction distribution
direction_distribution = siren.distributions.IsotropicDirection()
primary_injection_distributions["direction"] = direction_distribution
primary_physical_distributions["direction"] = direction_distribution

# position distribution
muon_range_func = siren.distributions.LeptonDepthFunction()
position_distribution = siren.distributions.ColumnDepthPositionDistribution(
    600, 600.0, muon_range_func, set(controller.GetDetectorModelTargets()[0])
)
primary_injection_distributions["position"] = position_distribution

# SetProcesses
controller.SetProcesses(
    primary_type, primary_injection_distributions, primary_physical_distributions
)

def add_secondary_to_controller(controller, secondary_type, secondary_xsecs, secondary_decays = None):
    if secondary_decays is not None:
        secondary_collection = siren.interactions.InteractionCollection(secondary_type, [secondary_xsecs], [secondary_decays])
    else:
        secondary_collection = siren.interactions.InteractionCollection(secondary_type, [secondary_xsecs])
    secondary_injection_process = siren.injection.SecondaryInjectionProcess()
    secondary_physical_process = siren.injection.PhysicalProcess()
    secondary_injection_process.primary_type = secondary_type
    secondary_physical_process.primary_type = secondary_type
    secondary_injection_process.AddSecondaryInjectionDistribution(siren.distributions.SecondaryPhysicalVertexDistribution())
    controller.secondary_injection_processes.append(secondary_injection_process)
    controller.secondary_physical_processes.append(secondary_physical_process)

    return secondary_collection

# secondary interactions
charms = siren.dataclasses.Particle.ParticleType.Charm
DPlus = siren.dataclasses.Particle.ParticleType.DPlus
D0 = siren.dataclasses.Particle.ParticleType.D0
charm_hadronization = siren.interactions.CharmHadronization()
DPlus_decay = siren.interactions.CharmMesonDecay(primary_type = DPlus)
D0_decay = siren.interactions.CharmMesonDecay(primary_type = D0)
D_energy_loss = siren.interactions.DMesonELoss()

secondary_charm_collection = add_secondary_to_controller(controller, charms, charm_hadronization)
secondary_DPlus_collection = add_secondary_to_controller(controller, DPlus, D_energy_loss, DPlus_decay)
secondary_D0_collection = add_secondary_to_controller(controller, D0, D_energy_loss, D0_decay)

controller.SetInteractions(primary_xs, [secondary_charm_collection, secondary_D0_collection, secondary_DPlus_collection])

controller.Initialize()

# def stop(datum, i):
#     secondary_type = datum.record.signature.secondary_types[i]
#     return ((secondary_type != siren.dataclasses.Particle.ParticleType.Charm) and (secondary_type != siren.dataclasses.Particle.ParticleType.DPlus))

def stop(datum, i):
    return False

controller.SetInjectorStoppingCondition(stop)

events = controller.GenerateEvents()

os.makedirs("output", exist_ok=True)

controller.SaveEvents("output/FullSim")
