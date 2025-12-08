import os
import numpy as np
import siren
from siren.SIREN_Controller import SIREN_Controller
import nuflux

# Number of events to inject
events_to_inject = 10000

# Expeirment to run
experiment = "IceCube"

# physical flux model to use
physical_flux = "atmos"

# Define the controller
controller = SIREN_Controller(events_to_inject, experiment, seed = 1)

# Particle to inject
primary_type = siren.dataclasses.Particle.ParticleType.NuMu

xs_option = "" # current choices are the empty string and cutoff-""
xsfiledir = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/miaochenjin/CharmXS/xsec_splines/M_Muon-{}105MeV".format(xs_option)

if xs_option == "":
    spline_option = "M_Muon-105MeV" # this is to account for the fact that I put the output names of this particular spline incorrectly
else:
    spline_option  = ""

# Cross Section Model
target_type = siren.dataclasses.Particle.ParticleType.Nucleon

DIS_xs = siren.interactions.CharmDISFromSpline(
    os.path.join(xsfiledir, "{}dsdxdynu-N-cc-HERAPDF15LO_EIG_central.fits".format(spline_option)),
    os.path.join(xsfiledir, "{}sigmanu-N-cc-HERAPDF15LO_EIG_central.fits".format(spline_option)),
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
edist = siren.distributions.PowerLaw(2, 1e4, 1e10)
primary_injection_distributions["energy"] = edist

# make an atmospheric flux
flux = nuflux.makeFlux('honda2006')
nu_type=nuflux.NuMu
erange = np.logspace(2,6,100)
erange_atmo = np.logspace(2,6,100)
cosrange = np.linspace(0,1,100)
atmo_flux_tables = {}
particle = nuflux.NuMu
siren_key = siren.dataclasses.Particle.ParticleType(int(particle))
atmo_flux_tables[siren_key] = np.zeros(len(erange))
for i,e in enumerate(erange):
    f = flux.getFlux(particle,e,cosrange)
    atmo_flux_tables[siren_key][i] += 0.01*np.sum(f) * 1e4 * 2 * np.pi

# this is for weighting the events to the astrophysical flux
if physical_flux == "astro":
    edist_astro = siren.distributions.PowerLaw(2, 1e4, 1e10)
    norm = 1e-18 * 1e4 * 4 * np.pi # GeV^-1 m^-2 s^-1
    edist_astro.SetNormalizationAtEnergy(norm,1e5)
    primary_physical_distributions["energy"] = edist_astro
elif physical_flux == "atmos":
    edist_atmo = siren.distributions.TabulatedFluxDistribution(erange_atmo,atmo_flux_tables[primary_type],True)
    primary_physical_distributions["energy"] = edist_atmo
else:
    primary_injection_distributions["energy"] = edist


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

# # secondary interactions
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

# secondary_DPlus_collection = add_secondary_to_controller(controller, DPlus, DPlus_decay)
# secondary_D0_collection = add_secondary_to_controller(controller, D0, D0_decay)

controller.SetInteractions(primary_xs, [secondary_charm_collection, secondary_D0_collection, secondary_DPlus_collection])
# controller.SetInteractions(primary_xs, [secondary_charm_collection])
# controller.SetInteractions(primary_xs, [])

controller.Initialize()

def stop(datum, i):
    return False

controller.SetInjectorStoppingCondition(stop)

events = controller.GenerateEvents()

print("finished generating events")

outdir = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/miaochenjin/DBSearch/SIREN_outputs"
expname = "0819_LE_debug_CharmHadron_atmos"
savedir = os.path.join(outdir, expname)

os.makedirs(savedir, exist_ok=True)

controller.SaveEvents("{}/{}_".format(savedir, expname))
