import pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load the data from the pickle file
with open('output/successful_events_summary.pkl', 'rb') as f:
    event_data = pickle.load(f)

# Extract interaction vertices and primary energies
vertices = []
energies = []
for event in event_data:
    if 'energies' in event and event['energies']:
        energies.append(event['energies'][0])
    if 'vertices' in event:
        for vertex in event['vertices']:
            vertices.append(vertex)

# Convert vertices to a numpy array for plotting
vertices = np.array(vertices)

# Plot the interaction vertices in 3D
if vertices.ndim == 2 and vertices.shape[1] == 3:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], c='b', marker='o')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')
    ax.set_title('Interaction Vertices within the Detector')
    plt.savefig('output/interaction_vertices_3d.png')
else:
    print("Error: vertices do not have the correct shape for plotting.")

# Plot the primary energies as a histogram
plt.figure()
plt.hist(energies, bins=100, color='b', alpha=0.7)
plt.xlabel('Primary Energy [GeV]')
plt.ylabel('Events')
plt.title(f'Distribution of Primary Energies (Total Events: {len(energies)})')
plt.savefig('output/primary_energy_distribution.png')
