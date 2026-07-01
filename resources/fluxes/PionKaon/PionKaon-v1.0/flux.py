import os
import siren

def load_flux(tag=None, min_energy=None, max_energy=None, physically_normalized=True):
    if tag is None:
        raise TypeError('"tag" is a required argument')
    tag = str(tag)

    has_energy_range = (min_energy is not None)
    if (min_energy is None) != (max_energy is None):
        raise RuntimeError('Neither or both "min_energy" and "max_energy" must be provided')

    parts = tag.split("_")
    if len(parts) != 2:
        raise ValueError(f'Tag "{tag}" must be of the form <mode>_<particle>')

    mode, particle = parts
    valid_modes     = ["FHC", "pion", "kaon"]
    valid_particles = ["numu", "numubar", "nue", "nuebar"]

    if mode not in valid_modes:
        raise ValueError(f'Beam mode "{mode}" not valid. Choose from {valid_modes}')
    if particle not in valid_particles:
        raise ValueError(f'Particle "{particle}" not valid. Choose from {valid_particles}')

    dat_name = {
        "FHC":  "PionKaon_FHC_all.dat",
        "pion": "PionKaon_FHC_pion.dat",
        "kaon": "PionKaon_FHC_kaon.dat",
    }[mode]

    abs_flux_dir    = os.path.dirname(__file__)
    input_flux_file = os.path.join(abs_flux_dir, dat_name)

    all_lines = open(input_flux_file, "r").readlines()
    headers   = all_lines[0].strip().split()
    data      = [line.strip().split() for line in all_lines[1:] if line.strip()]

    pid = headers.index(particle)
    energies = [(float(row[0]) + float(row[1])) / 2.0 for row in data]
    flux     = [float(row[pid]) / (50 * 1000 * 1e4) for row in data]

    if has_energy_range:
        table = siren.distributions.TabulatedFluxDistribution(
            min_energy, max_energy, energies, flux, physically_normalized)
    else:
        table = siren.distributions.TabulatedFluxDistribution(
            energies, flux, physically_normalized)
    return table
