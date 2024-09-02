import os
import siren


def load_flux(tag=None, min_energy=None, max_energy=None, physically_normalized=True):
    '''
    Accepts the following tags:
        {FHC,RHC}_{LE,ME}_{nue,nuebar,numu,numubar}
    '''

    if tag is None:
        raise TypeError("\"tag\" is a required argument")
    try:
        tag = str(tag)
    except:
        raise RuntimeError("\"tag\" must convert to a str")
    if min_energy is None != max_energy is None:
        raise RuntimeError("Neither or both \"min_energy\" and \"max_energy\" must be provided")
    has_energy_range = min_energy is not None

    mode, energy, particle = tag.split("_")
    if mode not in ["FHC","RHC"]:
        raise ValueError("%s beam mode specified in tag %s is not valid"%(mode,tag))
    if energy not in ["LE","ME"]:
        raise ValueError("%s energy mode specified in tag %s is not valid"%(energy,tag))
    if particle not in ["nue","numu","nuebar","numubar"]:
        raise ValueError("%s particle specified in tag %s is not valid"%(particle,tag))

    abs_flux_dir = os.path.dirname(__file__)
    input_flux_file = os.path.join(abs_flux_dir,
                                   "NUMI_%s_%s.dat" % (mode, energy))

    all_lines = open(input_flux_file, "r").readlines()
    headers = all_lines[0].strip().split()
    data = [line.strip().split() for line in all_lines[1:]]
    pid = headers.index(particle)
    e_idx = 0

    energies = [float(row[e_idx]) for row in data]
    flux = [float(row[pid]) * 1e4 for row in data] # put flux in units of nu/m^2/GeV/POT

    table = None
    if has_energy_range:
        table = siren.distributions.TabulatedFluxDistribution(min_energy, max_energy, energies, flux, physically_normalized)
    else:
        table = siren.distributions.TabulatedFluxDistribution(energies, flux, physically_normalized)

    return table
