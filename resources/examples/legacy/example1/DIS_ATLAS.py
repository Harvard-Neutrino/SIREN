import os
import siren
from siren._util import GenerateEvents,SaveEvents

seed = 99

# Number of events to inject
events_to_inject = int(1e5)

# Expeirment to run
experiment = "ATLAS"

# Load the detector model
detector_model = siren.utilities.load_detector(experiment)

# Particle to inject
primary_type = siren.dataclasses.Particle.ParticleType.NuMu

cross_section_model = "CSMSDISSplines"

# Cross Section Model
target_type = siren.dataclasses.Particle.ParticleType.Nucleon

primary_processes, secondary_processes = siren.utilities.load_processes(
    cross_section_model,
    primary_types = [primary_type],
    target_types = [target_type],
    isoscalar = True,
    process_types = ["CC"],
)

# Choose the direction we will inject from
# let's just inject upwards
injection_dir = siren.math.Vector3D(0, 0, 1)
injection_dir.normalize()

# Build the position distribution using properties of the geometry
fiducial_volume_name = "tilecal"
geo = None
for sector in detector_model.Sectors:
    if sector.name == fiducial_volume_name:
        geo = sector.geo
        break

det_position = detector_model.GeoPositionToDetPosition(siren.detector.GeometryPosition(geo.placement.Position))
det_rotation = geo.placement.Quaternion
det_placement = siren.geometry.Placement(det_position.get(), det_rotation)
cylinder = siren.geometry.Cylinder(det_placement, geo.Radius, geo.InnerRadius, geo.Z)

position_distribution = siren.distributions.CylinderVolumePositionDistribution(cylinder)

primary_injection_distributions = [
    siren.distributions.PrimaryMass(0),
    # energy distribution
    # HE SN flux from ATLAS paper
    siren.utilities.load_flux(
        "HE_SN",
        tag="numu",
        min_energy=100,
        max_energy=1e6,
        physically_normalized=True),
    siren.distributions.FixedDirection(injection_dir),
    position_distribution,
]

primary_physical_distributions = [
    # energy distribution
    # HE SN flux from ATLAS paper
    siren.utilities.load_flux(
        "HE_SN",
        tag="numu",
        min_energy=100,
        max_energy=1e6,
        physically_normalized=False),
    siren.distributions.FixedDirection(injection_dir),
]

injector = siren.injection.Injector()
injector.seed = seed
injector.number_of_events = events_to_inject
injector.detector_model = detector_model
injector.primary_type = primary_type
injector.primary_interactions = primary_processes[primary_type]
injector.primary_injection_distributions = primary_injection_distributions

events,gen_times = GenerateEvents(injector)

weighter = siren.injection.Weighter()
weighter.injectors = [injector]
weighter.detector_model = detector_model
weighter.primary_type = primary_type
weighter.primary_interactions = primary_processes[primary_type]
weighter.primary_physical_distributions = primary_physical_distributions

# Output events and weights
os.makedirs("output", exist_ok=True)
print("Saving events")
SaveEvents(events,weighter,gen_times,output_filename="output/ATLAS")

