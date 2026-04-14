import os
import siren


def load_flux(tag=None, min_energy=None, max_energy=None, physically_normalized=True):
    """
    Accepts the following tags:
        {FHC,RHC}_{nue,nuebar,numu,numubar}
    """

    if tag is None:
        raise TypeError("\"tag\" is a required argument")
    try:
        tag = str(tag)
    except:
        raise RuntimeError("\"tag\" must convert to a str")
    if min_energy is None != max_energy is None:
        raise RuntimeError("Neither or both \"min_energy\" and \"max_energy\" must be provided")
    has_energy_range = min_energy is not None

    mode, particle = tag.split("_")
    if mode not in ["FHC", "RHC"]:
        raise ValueError("%s beam mode specified in tag %s is not valid" % (mode, tag))
    if particle not in ["nue", "numu", "nuebar", "numubar"]:
        raise ValueError("%s particle specified in tag %s is not valid" % (particle, tag))

    abs_flux_dir = os.path.dirname(__file__)
    input_flux_file = os.path.join(abs_flux_dir, "BNB_%s.dat" % mode)

    all_lines = open(input_flux_file, "r").readlines()
    headers = all_lines[0].strip().split()
    data = [line.strip().split() for line in all_lines[1:]]
    pid = headers.index(particle)
    e_low_idx = 0
    e_high_idx = 1

    energies = [(float(row[e_low_idx]) + float(row[e_high_idx]))/2.0 for row in data]
    flux = [float(row[pid]) / 50 * 1000 * 1e4 for row in data] # put flux in units of nu/m^2/GeV/POT

    table = None
    if has_energy_range:
        table = siren.distributions.TabulatedFluxDistribution(min_energy, max_energy, energies, flux, physically_normalized)
    else:
        table = siren.distributions.TabulatedFluxDistribution(energies, flux, physically_normalized)

    return table

