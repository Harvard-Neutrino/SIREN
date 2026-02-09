import os
import numpy as np
import siren
from siren import utilities
from siren._util import GenerateEvents, SaveEvents, get_processes_model_path

SaveDarkNewsProcesses = siren.resources.processes.DarkNewsTables.SaveDarkNewsProcesses

# Define a DarkNews model
model_kwargs = {
    "m4": 0.47,  # 0.140,
    "mu_tr_mu4": 2.50e-6,  # 1e-6, # GeV^-1
    "UD4": 0,
    "Umu4": 0,
    "epsilon": 0.0,
    "gD": 0.0,
    "decay_product": "photon",
    "noHC": True,
    "HNLtype": "dirac",
}

# Number of events to inject
events_to_inject = 100000

# Expeirment to run
experiment = "MINERvA"
detector_model = utilities.load_detector(experiment)

# Particle to inject
primary_type = siren.dataclasses.Particle.ParticleType.NuMu

table_name = f"DarkNewsTables-v{siren.utilities.darknews_version()}/"
table_name += "Dipole_M%2.2e_mu%2.2e"%(model_kwargs["m4"],model_kwargs["mu_tr_mu4"])
table_dir = os.path.join(get_processes_model_path("DarkNewsTables"),table_name)
os.makedirs(table_dir,exist_ok=True)

# Load DarkNews processes
primary_processes, secondary_processes, primary_ups_keys, secondary_dec_keys = utilities.load_processes(
    "DarkNewsTables",
    primary_type=primary_type,
    detector_model = detector_model,
    table_name = table_name,
    **model_kwargs,
)

print(primary_processes)

# Find minimum decay width for secondary decays, to set the decay range for the position distribution
DN_min_decay_width = np.inf
# Loop over available decays, group by parent type
for secondary_type, decays in secondary_processes.items():
    for decay in decays:
        is_decay = False
        for signature in decay.GetPossibleSignatures():
            is_decay |= signature.target_type == siren.dataclasses.Particle.ParticleType.Decay
        if not is_decay:
            continue
        decay_width = decay.TotalDecayWidth(secondary_type)
        DN_min_decay_width = min(DN_min_decay_width, decay_width)

print("Minimum decay width for secondary decays: ", DN_min_decay_width)

# Primary distributions
primary_injection_distributions = {}
primary_physical_distributions = {}

mass_ddist = siren.distributions.PrimaryMass(0)

# energy distribution
edist = siren.utilities.load_flux("NUMI", tag="FHC_ME_numu", physically_normalized=True)
edist_gen = siren.utilities.load_flux("NUMI", tag="FHC_ME_numu", min_energy=model_kwargs["m4"], max_energy=20, physically_normalized=False)

# direction distribution
direction_distribution = siren.distributions.FixedDirection(siren.math.Vector3D(0, 0, 1.0))

# position distribution
decay_range_func = siren.distributions.DecayRangeFunction(
    model_kwargs["m4"], DN_min_decay_width, 3, 240
)
position_distribution = siren.distributions.DecayRangePositionDistribution(
    1.24,
    5.0,
    decay_range_func,
)

primary_injection_distributions = [
    mass_ddist,
    edist_gen,
    direction_distribution,
    position_distribution,
]

primary_physical_distributions = [
    edist,
    direction_distribution,
]

fiducial_volume = siren.utilities.get_fiducial_volume(experiment)
secondary_injection_distributions = {}
for secondary_type in secondary_processes.keys():
    secondary_injection_distributions[secondary_type] = [
        siren.distributions.SecondaryBoundedVertexDistribution(fiducial_volume)
    ]

def stop(datum, i):
    secondary_type = datum.record.signature.secondary_types[i]
    return secondary_type != siren.dataclasses.Particle.ParticleType.N4

injector = siren.injection.Injector()
injector.number_of_events = events_to_inject
injector.detector_model = detector_model
injector.primary_type = primary_type
injector.primary_interactions = primary_processes[primary_type]
injector.primary_injection_distributions = primary_injection_distributions
injector.secondary_interactions = secondary_processes
injector.secondary_injection_distributions = secondary_injection_distributions
injector.stopping_condition = stop

events,gen_times = GenerateEvents(injector)

# Output the events
os.makedirs("output", exist_ok=True)

weighter = siren.injection.Weighter()
weighter.injectors = [injector]
weighter.detector_model = detector_model
weighter.primary_type = primary_type
weighter.primary_interactions = primary_processes[primary_type]
weighter.secondary_interactions = secondary_processes
weighter.primary_physical_distributions = primary_physical_distributions
weighter.secondary_physical_distributions = {}

SaveEvents(events,weighter,gen_times,output_filename="output/CCM_Dipole")

# save cross section tables
SaveDarkNewsProcesses(table_dir,
                      primary_processes,
                      primary_ups_keys,
                      secondary_processes,
                      secondary_dec_keys)



weights = [weighter(event) for event in events]


