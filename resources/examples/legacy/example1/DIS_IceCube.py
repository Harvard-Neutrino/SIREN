import os
import siren
from siren._util import GenerateEvents,SaveEvents

# Number of events to inject
events_to_inject = int(1e5)

# Experiment to run
experiment = "IceCube"
detector_model = siren.utilities.load_detector(experiment)

# Particle to inject
primary_type = siren.dataclasses.Particle.ParticleType.NuMu

# Cross-section model to use
cross_section_model = "CSMSDISSplines"

# Load the cross-section model
primary_processes, _ = siren.utilities.load_processes(
    cross_section_model,
    primary_types=[primary_type],
    target_types=[siren.dataclasses.Particle.ParticleType.Nucleon],
    isoscalar=True,
    process_types=["CC"]
)

# Extract the primary cross-sections for the primary type
primary_cross_sections = primary_processes[primary_type]

# Set up the Injector
injector = siren.injection.Injector()
injector.number_of_events = events_to_inject
injector.detector_model = detector_model
injector.primary_type = primary_type
injector.primary_interactions = primary_cross_sections

# Directly set the distributions
injector.primary_injection_distributions = [
    siren.distributions.PrimaryMass(0),  # Mass distribution
    siren.distributions.PowerLaw(2, 1e3, 1e6),  # Energy distribution
    siren.distributions.IsotropicDirection(),  # Direction distribution
    siren.distributions.ColumnDepthPositionDistribution(600, 600.0, siren.distributions.LeptonDepthFunction())  # Position distribution
]

# Generate events
events,gen_times = GenerateEvents(injector)

# Set up the Weighter for event weighting (without position distribution)
weighter = siren.injection.Weighter()
weighter.injectors = [injector]
weighter.detector_model = detector_model
weighter.primary_type = primary_type
weighter.primary_interactions = primary_cross_sections
weighter.primary_physical_distributions = [
    siren.distributions.PowerLaw(2, 1e3, 1e6),  # Energy distribution
    siren.distributions.IsotropicDirection()  # Direction distribution
]


# Output events and weights
os.makedirs("output", exist_ok=True)
SaveEvents(events,weighter,gen_times,output_filename="output/IceCube")
