import os
import numpy as np
import siren
from siren import utilities
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--target", type=int, help = '0 for upper tungsten target, 1 for lower tungsten target', default = 0)
args = parser.parse_args()
target = args.target


# Number of events to inject
events_to_inject = 10000

# Experiment to run
experiment = "CCM"
detector_model = utilities.load_detector(experiment)

# Fiducial volume
fiducial_volume = siren.utilities.get_fiducial_volume(experiment)

# Particle to inject (NuE for Marley)
primary_type = siren.dataclasses.Particle.ParticleType.NuE

# Load MARLEY processes
primary_processes, secondary_processes = siren.resources.load_processes(
        "MarleyCrossSection",
        primary_types=[primary_type],
        process_types=["CC"],
)

print("Primary processes loaded:")
for key, val in primary_processes.items():
    print(f"Type: {key}, Process: {val}")

# Mass distribution
mass_ddist = siren.distributions.PrimaryMass(0)

# Primary distributions
edist = siren.distributions.PiDARNuEDistribution()

# Cone direction distribution from tungsten target
opening_angle = np.arcsin(1.21 / 23.0)

#Darcy's target
if (target == 0):
    ### upper target!
    target_origin = siren.math.Vector3D(0, 0, 0.1375)
else:
    ### lower target!
    target_origin = siren.math.Vector3D(0, 0, -0.241)

detector_origin = siren.math.Vector3D(23, 0, -0.65)

target_dir = detector_origin - target_origin
target_dir.normalize()
target_inj_ddist = siren.distributions.Cone(target_dir, opening_angle)

phys_ddist = siren.distributions.IsotropicDirection()

# Position distribution
max_dist = 25
# Position distribution: use fixed target distribution that considers physical volume of the target
if (target == 0):
    ### upper target!
    target_cylinder = siren.geometry.Cylinder(siren.geometry.Placement(target_origin - detector_origin), 0.05, 0.0, 0.091)
else:
    ### lower target!
    target_cylinder = siren.geometry.Cylinder(siren.geometry.Placement(target_origin - detector_origin), 0.05, 0.0, 0.298)

target_pos_dist = siren.distributions.FixedTargetPositionDistribution(target_cylinder, fiducial_volume, max_dist)

primary_injection_distributions = [
    mass_ddist,  # Mass distribution
    edist,  # Energy distribution
    target_inj_ddist,  # Direction distribution
    target_pos_dist  # Position distribution
]

primary_physical_distributions = [
    edist,  # Energy distribution
    phys_ddist,  # Direction distribution
]

injector = siren.injection.Injector()
injector.number_of_events = events_to_inject
injector.detector_model = detector_model
injector.primary_type = primary_type
injector.primary_interactions = primary_processes[primary_type]
injector.primary_injection_distributions = primary_injection_distributions

# Generate events
print("Injecting events")
events = [injector.generate_event() for _ in range(events_to_inject)]

weighter = siren.injection.Weighter()
weighter.injectors = [injector]
weighter.detector_model = detector_model
weighter.primary_type = primary_type
weighter.primary_interactions = primary_processes[primary_type]
weighter.secondary_interactions = secondary_processes
weighter.primary_physical_distributions = primary_physical_distributions
weighter.secondary_physical_distributions = {}

print("Weighting events")
weights = [weighter(event) for event in events]

successful_events = [event for event in events if len(event.tree) > 0]
failed_events = [event for event in events if len(event.tree) == 0]
print(f"Successful events: {len(successful_events)}")
print(f"Failed events: {len(failed_events)}")


#Get the attributes of the InteractionRecord and print them
for i, event in enumerate(successful_events):
    print(f"Event {i+1}:")
    for j, datum in enumerate(event.tree):
        print(f"  Vertex {j+1}:")
        if hasattr(datum, 'record'):
            record = datum.record
            # List all attributes of the record to see what's available
            print(f"    Attributes of InteractionRecord:")
            for attr in dir(record):
                if not attr.startswith("_"):  # Ignore private or internal methods
                    print(f"      {attr}: {getattr(record, attr)}")


# Output the events
os.makedirs("output", exist_ok=True)

# Extract and save interaction vertex and energy of the primary particle from each event
event_data = []

for i, (event, weight) in enumerate(zip(events, weights)):
    if len(event.tree) == 0:
        continue
    event_info = {
        'event_id': i,
        'vertices': [],
        'energies': [],
        "weight": weight,
    }
    for datum in event.tree:
        if hasattr(datum, 'record'):
            record = datum.record
            event_info['vertices'].append(record.interaction_vertex)
            print(f"Interaction Vertex: {record.interaction_vertex}")
            event_info['energies'].append(record.primary_momentum[0])
            print(f"Energy: {record.primary_momentum[0]}")
            print(f"Weight: {weight}")
    event_data.append(event_info)

with open('output/successful_events_summary.pkl', 'wb') as f:
    pickle.dump(event_data, f)

