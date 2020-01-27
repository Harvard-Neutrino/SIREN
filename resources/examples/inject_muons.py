# Benjamin Smithers
# benjamin.smithers@mavs.uta.edu

# this example script ...
#   + imports the Nuzooka libraries 
#   + creates a muon Nuzooka, and its operator 
#   + tells the operator to fire 

import pylepton_injector as LI
from math import pi

# define the injector 
n_events    = 50000
diff_xs     = "/home/benito/software/cross_sections/dsdxdy_nubar_CC_iso.fits"
total_xs    = "/home/benito/software/cross_sections/sigma_nubar_CC_iso.fits"
is_ranged   = False
final_1     = LI.Particle.NuEBar
final_2     = LI.Particle.Hadrons

# the above describes a NC interaction, so let's use Volume Mode
# create the injector list using the above parameters
the_injector = LI.injector( n_events, final_1, final_2, diff_xs, total_xs, is_ranged)

# Now, we'll make a new injector for muon tracks 
n_events    = 55000
diff_xs     = "/home/benito/software/cross_sections/dsdxdy_nubar_NC_iso.fits"
total_xs    = "/home/benito/software/cross_sections/sigma_nubar_NC_iso.fits"
is_ranged   = True
final_1     = LI.Particle.MuMinus
final_2     = LI.Particle.Hadrons
the_next_injector = LI.injector( n_events , final_1, final_2, diff_xs, total_xs, is_ranged)



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
controller  = LI.Controller( the_injector, minE, maxE, gamma, minAzimuth, maxAzimuth, minZenith, maxZenith)  
controller.AddInjector( the_next_injector )
# injection radius and endcap length are left as defaults

# specify the output
controller.Output("./data_output.h5")
controller.LICFile("./config.lic")

# run the simulation
controller.Execute()
