import os
import glob
from typing import List, Optional
import siren
import collections
from siren.download import ensure_tar_archive, writable_data_dir

# MarleyCrossSection-v2.0: MARLEY v2 hybrid model resources.
#
# Differences with respect to v1.0:
#  - Self-contained after download: the complete, SHA-256-verified data bundle
#    comes from SIREN-data (react files, CRPA response tables, nuclear charge
#    radii, v2-format structure files with per-level half-lives, masses, and
#    parities). No dependency on $PREFIX or an external MARLEY installation.
#  - The CC process is the MARLEY v2 recommended hybrid: measured discrete
#    levels (Bhattacharya 2009) + HF-CRPA continuum, loaded together as a
#    single MarleyCrossSection object (multi-react constructor).
basepath = os.path.dirname(os.path.abspath(__file__))

_DATA_ARCHIVE = "MarleyCrossSection-v2.0.tar.xz"
_DATA_URL = (
    "https://raw.githubusercontent.com/SIREN-Generator/SIREN-data/main/"
    "processes/MarleyCrossSection/MarleyCrossSection-v2.0/"
    + _DATA_ARCHIVE
)
_DATA_SHA256 = "eb64ee2b330001205c96118d5dbdd2021c41d4f08b506363cef8a700433a6aa8"
_DATA_DIR = None


def _search_path(data_dir):
    return ':'.join([
        data_dir,
        os.path.join(data_dir, 'react'),
        os.path.join(data_dir, 'structure')
    ])


def _get_data_dir():
    global _DATA_DIR
    if _DATA_DIR is None:
        _DATA_DIR = writable_data_dir(basepath)
    return _DATA_DIR


def fetch_data():
    """Download and extract the MARLEY v2 input bundle from SIREN-data."""
    data_dir = _get_data_dir()
    ensure_tar_archive(
        _DATA_URL, _DATA_ARCHIVE, data_dir, sha256=_DATA_SHA256)
    return data_dir


default_marley_search_path = _search_path(_get_data_dir())

#lists of neutrinos and antineutrinos
neutrinos = [
        siren.dataclasses.Particle.ParticleType.NuE,
        siren.dataclasses.Particle.ParticleType.NuMu,
        siren.dataclasses.Particle.ParticleType.NuTau,
]
antineutrinos = [
        siren.dataclasses.Particle.ParticleType.NuEBar,
        siren.dataclasses.Particle.ParticleType.NuMuBar,
        siren.dataclasses.Particle.ParticleType.NuTauBar,
]

#list of processes and the react files that build each one.
#CC is the v2 hybrid: BOTH files form a single cross-section object.
processes = ["CC", "CEvNS", "ES"]
reactions_by_process = {
    "CC": ["ve40ArCC_HF-CRPA.react", "ve40ArCC_Bhattacharya2009-Discrete.react"],
    "CEvNS": ["CEvNS40Ar.react"],
    "ES": ["ES.react"],
}

# Auxiliary data files needed by the v2 engine, with the relative name each
# one must keep on the MARLEY search path (the HF-CRPA manifest references the
# response tables as "crpa/<table>.dat"). Keep this list explicit because the
# data files are not present until fetch_data() extracts the archive.
_crpa_tables = [f"responses_ar40_crpa_J{j}_G3.dat" for j in range(6)]
_logger_config = [("config/logger.js", "data/config/logger.js")]
aux_by_process = {
    "CC": [(os.path.join("react", "crpa", f), os.path.join("crpa", f))
           for f in _crpa_tables]
          + [("nuclear_charge_radii.js", "nuclear_charge_radii.js")]
          + _logger_config,
    "CEvNS": [("nuclear_charge_radii.js", "nuclear_charge_radii.js")]
             + _logger_config,
    "ES": _logger_config,
}

#mapping of processes to primary particles
primaries_by_process = {
    "CC": [siren.dataclasses.Particle.ParticleType.NuE, siren.dataclasses.Particle.ParticleType.NuEBar],
    "ES": neutrinos + antineutrinos,
    "CEvNS": neutrinos + antineutrinos,
}

def _find_file(search_path, fname):
    for path in search_path.split(':'):
        full_path = os.path.join(path, fname)
        if os.path.exists(full_path):
            return full_path
    raise FileNotFoundError(f"Could not find file {fname} in search path {search_path}")

#Get primary particles, setting default to NuE if none are provided, only supports NuE and NuEBar
def _get_primary_types(primary_types):
    if primary_types is None:
        primary_types = [
                siren.dataclasses.Particle.ParticleType.NuE,
        ]

    supported_primaries = [siren.dataclasses.Particle.ParticleType.NuE, siren.dataclasses.Particle.ParticleType.NuEBar]
    for i, p in enumerate(primary_types):
        if p not in supported_primaries:
            raise ValueError(f"primary_types[{i}] \"{p}\" not supported. Allowed primary_types are: {supported_primaries}")

    if len(primary_types) == 0:
        print("Warning: len(primary_types) == 0")

    return primary_types

#Get process types, setting default to CC, CEvNS, and ES if none are provided
def _get_process_types(process_types):
    if process_types is None:
        process_types = ["CC", "CEvNS", "ES"]

    for i, p in enumerate(process_types):
        if p not in processes:
            raise ValueError(f"process_types[{i}] \"{p}\" not supported. Allowed processes are: {processes}")

    if len(process_types) == 0:
        print("Warning: len(process_types) == 0")

    return process_types

#Load processes based on primary particles and process types
def load_processes(
    primary_types: Optional[List[siren.dataclasses.Particle.ParticleType]] = None,
    process_types: Optional[List[str]] = None,
    marley_search_path: Optional[str] = None,
    ):

    if marley_search_path is None:
        marley_search_path = _search_path(fetch_data())

    primary_types = _get_primary_types(primary_types)
    process_types = _get_process_types(process_types)

    primary_processes_dict = collections.defaultdict(list)

    for process in process_types:
        #Resolve every react file of this process (CC = hybrid pair)
        react_fnames = [_find_file(marley_search_path, r) for r in reactions_by_process[process]]

        nuclide_index_fname = _find_file(marley_search_path, "nuclide_index.txt")
        nuclide_path = os.path.dirname(nuclide_index_fname)
        nuclide_fnames = sorted(glob.glob(os.path.join(nuclide_path, "*.dat")))
        masses_fname = _find_file(marley_search_path, "mass_table.js")
        gs_parity_fname = _find_file(marley_search_path, "gs_spin_parity_table.txt")

        #Auxiliary files: (path inside the bundle, relative name for MARLEY)
        aux_files = [_find_file(marley_search_path, rel_src)
                     for rel_src, _ in aux_by_process[process]]
        aux_names = [name for _, name in aux_by_process[process]]

        xs = siren.interactions.MarleyCrossSection(
            react_fnames, nuclide_index_fname, nuclide_fnames,
            masses_fname, gs_parity_fname, aux_files, aux_names)

        reaction_primary_types = set(primaries_by_process[process]) & set(primary_types)

        #Add cross section to primary_processes_dict
        for primary_type in reaction_primary_types:
            primary_processes_dict[primary_type].append(xs)

    return dict(primary_processes_dict), {}
