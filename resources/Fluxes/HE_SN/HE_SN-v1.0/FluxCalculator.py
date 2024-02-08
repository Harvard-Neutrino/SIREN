import os

def MakeFluxFile(tag, abs_flux_dir):
    '''
    tag is not used here
    '''
    return os.path.join(abs_flux_dir,"dN_dE_SNe_2n_D1_0_s20_t100d_NuMu_d10kpc.txt")