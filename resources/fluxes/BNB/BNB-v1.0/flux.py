import os
import siren
from siren.download import ensure_files, writable_data_dir

_ABS_DIR = writable_data_dir(os.path.dirname(os.path.abspath(__file__)))
_DATA_BASE = "https://raw.githubusercontent.com/SIREN-Generator/SIREN-data/main/fluxes/BNB/BNB-v1.0"

_DATA_FILES = [
    {"path": os.path.join(_ABS_DIR, "BNB_FHC.dat"), "url": f"{_DATA_BASE}/BNB_FHC.dat",
     "sha256": "093a46a07c66c170ad08278ebded5ae6c314c2a73c7696884623181dd03dd887"},
    {"path": os.path.join(_ABS_DIR, "BNB_RHC.dat"), "url": f"{_DATA_BASE}/BNB_RHC.dat",
     "sha256": "3b16ada8932faae4ef6412dcefad2a29ef3d450e8af60f925eaf6f2ca1de1369"},
]


def fetch_data():
    ensure_files(_DATA_FILES)


def load_flux(tag=None, min_energy=None, max_energy=None, physically_normalized=True):
    """
    Accepts the following tags:
        {FHC,RHC}_{nue,nuebar,numu,numubar}
    """

    fetch_data()

    if tag is None:
        raise TypeError("\"tag\" is a required argument")
    tag = str(tag)
    if (min_energy is None) != (max_energy is None):
        raise RuntimeError("Neither or both \"min_energy\" and \"max_energy\" must be provided")
    has_energy_range = min_energy is not None

    mode, particle = tag.split("_")
    if mode not in ["FHC", "RHC"]:
        raise ValueError("%s beam mode specified in tag %s is not valid" % (mode, tag))
    if particle not in ["nue", "numu", "nuebar", "numubar"]:
        raise ValueError("%s particle specified in tag %s is not valid" % (particle, tag))

    input_flux_file = os.path.join(_ABS_DIR, "BNB_%s.dat" % mode)

    with open(input_flux_file, "r") as f:
        all_lines = f.readlines()
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
