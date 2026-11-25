import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load the events saved by CC_MARLEY_CCM.py
events = ak.from_parquet("output/CCM_MARLEY.parquet")

# Keep events with at least one interaction (failed injection attempts are
# saved with an empty interaction list).
events = events[ak.num(events["vertex"]) > 0]

# Extract interaction vertices and primary energies. Each event stores one
# entry per interaction; the first interaction is the primary CC vertex.
vertices = np.asarray(ak.flatten(events["vertex"]))
energies = np.asarray(events["primary_momentum"][:, 0, 0])

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
