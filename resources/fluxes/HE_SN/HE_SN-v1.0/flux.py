import os
import siren

def load_flux(tag=None, min_energy=None, max_energy=None, physically_normalized=True):
    '''
    only supported tag is "numu"
    '''
    if tag!="numu":
        raise ValueError("Tag %s not supported for HE SN" % tag)

    has_energy_range = min_energy is not None

    abs_flux_dir = os.path.dirname(__file__)
    input_flux_file = os.path.join(abs_flux_dir,
                                   "dN_dE_SNe_2n_D1_0_s20_t100d_NuMu_d10kpc.txt")

    all_lines = open(input_flux_file, "r").readlines()
    headers = all_lines[0].strip().split()
    data = [line.strip().split() for line in all_lines[1:]]
    e_idx = 0
    flux_idx = 1

    energies = [float(row[e_idx]) for row in data]
    flux = [float(row[flux_idx]) * 1e4 for row in data] # put flux in units of nu/m^2/GeV/100d

    table = None
    if has_energy_range:
        table = siren.distributions.TabulatedFluxDistribution(min_energy, max_energy, energies, flux, physically_normalized)
    else:
        table = siren.distributions.TabulatedFluxDistribution(energies, flux, physically_normalized)

    return table

