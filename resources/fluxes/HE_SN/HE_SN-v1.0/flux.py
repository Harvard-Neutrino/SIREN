import os
import siren
from siren.download import ensure_files, writable_data_dir, resolve_data_path

_INSTALL_DIR = os.path.dirname(os.path.abspath(__file__))
_DATA_BASE = "https://raw.githubusercontent.com/SIREN-Generator/SIREN-data/main/fluxes/HE_SN/HE_SN-v1.0"

_ABS_DIR = None

def _get_abs_dir():
    global _ABS_DIR
    if _ABS_DIR is None:
        _ABS_DIR = writable_data_dir(_INSTALL_DIR)
    return _ABS_DIR


def fetch_data():
    abs_dir = _get_abs_dir()
    ensure_files([
        {"path": os.path.join(abs_dir, "dN_dE_SNe_2n_D1_0_s20_t100d_NuMu_d10kpc.txt"),
         "url": f"{_DATA_BASE}/dN_dE_SNe_2n_D1_0_s20_t100d_NuMu_d10kpc.txt",
         "sha256": "b6dc394dbcccbfb19a09273656b9305d0611f625e1f523a9f04e3b51cbbfb00c"},
    ])


def load_flux(tag=None, min_energy=None, max_energy=None, physically_normalized=True):
    '''
    only supported tag is "numu"
    '''
    fetch_data()

    if tag!="numu":
        raise ValueError("Tag %s not supported for HE SN" % tag)

    has_energy_range = min_energy is not None

    input_flux_file = resolve_data_path(_INSTALL_DIR, _get_abs_dir(),
                                       "dN_dE_SNe_2n_D1_0_s20_t100d_NuMu_d10kpc.txt")

    with open(input_flux_file, "r") as f:
        all_lines = f.readlines()
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
