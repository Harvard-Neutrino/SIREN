import os
import siren
from siren.download import ensure_files, writable_data_dir, resolve_data_path

_INSTALL_DIR = os.path.dirname(os.path.abspath(__file__))
_DATA_BASE = "https://raw.githubusercontent.com/SIREN-Generator/SIREN-data/main/fluxes/NUMI/NUMI-v1.0"

_ABS_DIR = None

def _get_abs_dir():
    global _ABS_DIR
    if _ABS_DIR is None:
        _ABS_DIR = writable_data_dir(_INSTALL_DIR)
    return _ABS_DIR

_DATA_ENTRIES = [
    ("FHC", "LE", "38fc38719d88dd96fdd0f849c9921fcee224a9faaa6a4f0d0a569d408bc3b372"),
    ("FHC", "ME", "2401354dc561a51a84de87e5f47308fdb736b2739316b275fd0959e1428c9860"),
    ("FHC", "ME_unofficial", "3c5e1b3afe0e992602442bac13b07f2359a290dfe31c45a6b0026bf1542befd1"),
    ("RHC", "LE", "e994bb199e441a0aa695ec4a6b1698dbbf4d876e17e0c650d252c9943cd7da2c"),
    ("RHC", "ME", "52c538efd182b03bcf45577bdf81172575fe4247eec5cc980a41764e08fbabf4"),
    ("RHC", "ME_unofficial", "63bfbff1b27303ee04c7aade3841481a22416b2f5b39cd807398a797c61c90d4"),
]


def fetch_data():
    abs_dir = _get_abs_dir()
    ensure_files([
        {"path": os.path.join(abs_dir, f"NUMI_{mode}_{energy}.dat"),
         "url": f"{_DATA_BASE}/NUMI_{mode}_{energy}.dat",
         "sha256": sha}
        for mode, energy, sha in _DATA_ENTRIES
    ])


def load_flux(tag=None, min_energy=None, max_energy=None, physically_normalized=True):
    '''
    Accepts the following tags:
        {FHC,RHC}_{LE,ME}_{nue,nuebar,numu,numubar}
    '''

    fetch_data()

    if tag is None:
        raise TypeError("\"tag\" is a required argument")
    tag = str(tag)
    if (min_energy is None) != (max_energy is None):
        raise RuntimeError("Neither or both \"min_energy\" and \"max_energy\" must be provided")
    has_energy_range = min_energy is not None

    parts = tag.split("_")
    if len(parts) < 3:
        raise ValueError(
            "Tag %r is not valid. Expected {FHC,RHC}_{LE,ME,ME_unofficial}_{particle}" % tag)
    mode = parts[0]
    particle = parts[-1]
    energy = "_".join(parts[1:-1])
    if mode not in ["FHC","RHC"]:
        raise ValueError("%s beam mode specified in tag %s is not valid"%(mode,tag))
    if energy not in ["LE","ME","ME_unofficial"]:
        raise ValueError("%s energy mode specified in tag %s is not valid"%(energy,tag))
    if particle not in ["nue","numu","nuebar","numubar"]:
        raise ValueError("%s particle specified in tag %s is not valid"%(particle,tag))

    input_flux_file = resolve_data_path(_INSTALL_DIR, _get_abs_dir(),
                                       "NUMI_%s_%s.dat" % (mode, energy))

    with open(input_flux_file, "r") as f:
        all_lines = f.readlines()
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
