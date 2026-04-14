# Benjamin Smithers
# benjamin.smithers@mavs.uta.edu

# this example script ...
#   + imports the SIREN libraries 
#   + creates two injectors, and their operator 
#   + tells the operator to execute the process 

import SIREN as LI
from math import pi
import os 

# use this if you are on the Cobalt testbed 
#xs_folder = "/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms_differential_v1.0"

# for now, just use the example cross sections that come pre-installed. 
# this looks in the parent folder to this file's containing folder 
xs_folder = os.path.join( os.path.dirname(__file__), '..' )

# Now, we'll make a new injector for muon tracks 
n_events    = 55000
diff_xs     = xs_folder + "/test_xs.fits"
total_xs    = xs_folder + "/test_xs_total.fits"
is_ranged   = True
final_1     = siren.Particle.MuMinus
final_2     = siren.Particle.Hadrons
the_injector = siren.Injector( n_events , final_1, final_2, diff_xs, total_xs, is_ranged)



deg = pi/180.

# define some defaults 
minE        = 1000.     # [GeV]
maxE        = 100000.   # [GeV]
gamma       = 2. 
minZenith   = 80.*deg
maxZenith   = 180.*deg
minAzimuth  = 0.*deg
maxAzimuth  = 180.*deg

# construct the controller 
controller  = siren.Controller( the_injector, minE, maxE, gamma, minAzimuth, maxAzimuth, minZenith, maxZenith)  

# specify the output
controller.NameOutfile("./data_output.h5")
controller.NameLicFile("./config.lic")

# run the simulation
controller.Execute()
