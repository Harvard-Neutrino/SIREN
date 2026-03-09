from scipy.interpolate import interp1d
import numpy as np
import os

def bar_scaling(abs_flux_dir, m_n_min=0.002, m_n_max=0.4, m_n_step=25):

    
    # Energy in GeV
    m_n, ratio = np.loadtxt(os.path.join(abs_flux_dir,'ratio.dat'), usecols=(0, 1), unpack=True)
    
    extended_m_n = np.arange(m_n_min, m_n_max + m_n_step, m_n_step)
    
    # Perform interpolation only where necessary
    mask = (extended_m_n < m_n[0]) | (extended_m_n > m_n[-1])
    extended_min = np.interp(extended_m_n, m_n, ratio)
    
    interp_func = interp1d(m_n, ratio, kind='linear', bounds_error=False, fill_value=(ratio[0], ratio[-1]))
    extended_min[~mask] = interp_func(extended_m_n[~mask])
    
    # Create final interpolation function
    bar_flux_scaling = interp1d(extended_m_n, extended_min, kind='linear', bounds_error=False, fill_value='extrapolate')
    
    return bar_flux_scaling

def MakeFluxFile(tag, abs_flux_dir, mass_cut=None):

    '''
    Accepts the following tags:
        {nue,nuebar,numu,numubar}
    '''
    
    particle, enhance = tag.split("_")
        
    if enhance == 'MINUS':
        bar = True
        bar_scale = bar_scaling(abs_flux_dir)
        
    elif enhance == 'PLUS':
        bar = False
    else:
        print(f'Enhance tag "{tag}" is invalid')
        exit(0)

    if particle not in ["numu", "numubar"]:
        print("%s particle specified in tag %s is not valid"%(particle))
        exit(0)
        
    input_flux_file = os.path.join(abs_flux_dir,
                                    "kaon-flux-data.dat")

    output_flux_file = os.path.join(abs_flux_dir, 
                                    f"kaon-flux-{particle}_bar.dat") if bar else os.path.join(abs_flux_dir, f"kaon-flux-{particle}.dat")

        
    with open(input_flux_file,"r") as fin:
        all_lines = fin.readlines()
        data = [line.strip().split() for line in all_lines[:]]
                    
        with open(output_flux_file,"w") as fout:
            for row in data:
                E = float(row[0])
                flux = 10**float(row[1]) * bar_scale(float(row[0])) if bar else 10**float(row[1])
                flux*=2e-16 # put flux in units of nu/m^2/GeV/POT
                print(E, flux, file=fout)
                
    return output_flux_file 
