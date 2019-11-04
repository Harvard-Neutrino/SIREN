# Benjamin Smithers
# benjamin.smithers@mavs.uta.edu

# this example script ...
#   + imports the Nuzooka libraries 
#   + creates a muon Nuzooka, and its operator 
#   + tells the operator to fire 

import LeptonInjector as LI

# define the injector 
final_type1 = LI.ParticleType.MuPlus
final_type2 = LI.ParticleType.Hadrons
n_events    = 100
diff_xs     = "/path/to/example/diff_xs.fits"
total_xs    = "/path/to/example/total_xs.fits"
is_ranged   = True # use the ranged mode injection
# because the final state had a charged muon and hadrons, the event is the result of an antimuon neutrino undergoing a charged current interaction

# create the injector list using the above parameters
LI_injector = [LI.injector( n_events, final_type1, final_type2, diff_xs, total_xs, is_ranged)]

minE        = 1000*LI.Constants.GeV
maxE        = 100000*LI.Constants.GeV
gamma       = 2. #unitless
minZenith   = 80.*LI.Constants.Degrees
maxZenith   = 180.*LI.Constants.Degrees
minAzimuth  = 0.*LI.Constants.Degrees
maxAzimuth  = 180.*LI.Constants.Degrees

# construct the controller 
controller  = LI.Controller( LI_injector, minE, maxE, gamma, minAzimuth, maxAzimuth, minZenith, maxZenith)  
# injection radius and endcap length are left as defaults

# specify the output
Controller.Output("./data_output.h5")

# run the simulation
Controller.Execute()
