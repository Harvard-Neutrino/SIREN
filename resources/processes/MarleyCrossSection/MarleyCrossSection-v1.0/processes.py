import os
import glob
from typing import List, Optional
import siren
import collections
from siren.download import ensure_tar_archive, writable_data_dir

# MarleyCrossSection-v1.0: MARLEY v1-era nuclear data (measured discrete levels
# + Cheoun QRPA) evaluated with the MARLEY 2.0.0 engine. Kept primarily as a
# validation control for the single-react code path.
#
# Data comes from the same SHA-256-verified SIREN-data bundle used by
# MarleyCrossSection-v2.0; the v1-format react files live under react/v1/.
basepath = os.path.dirname(os.path.abspath(__file__))

# Keep these constants in sync with MarleyCrossSection-v2.0/processes.py
# (both bundles download the same archive).
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
    """Download and extract the MARLEY input bundle from SIREN-data."""
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

#list of processes and mapping of reaction files to processes
processes = ["CC", "CEvNS", "ES"]
process_by_reaction = {
    "v1/ve40ArCC_Bhattacharya1998.react": "CC",
    "v1/ve40ArCC_Liu1998.react": "CC",
    "v1/ve40ArCC_Bhattacharya2009.react": "CC",
    "CEvNS40Ar.react": "CEvNS",
    "ES.react": "ES",
}

# Runtime data files required by the MARLEY 2.0.0 engine, with the relative
# name each one must keep on the MARLEY search path.
_runtime_aux = [
    ("config/logger.js", "data/config/logger.js"),
    ("nuclear_charge_radii.js", "nuclear_charge_radii.js"),
]

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
    reaction_name: Optional[str] = None,
    ):

    if marley_search_path is None:
        marley_search_path = _search_path(fetch_data())

    primary_types = _get_primary_types(primary_types)
    process_types = _get_process_types(process_types)
    reaction_names = []

    if reaction_name is not None:
        reaction_names = [reaction_name]
    else:
        if "CC" in process_types and (siren.dataclasses.Particle.ParticleType.NuE in primary_types or siren.dataclasses.Particle.ParticleType.NuEBar in primary_types):
            reaction_names.append("v1/ve40ArCC_Bhattacharya2009.react")
        if "CEvNS" in process_types:
            reaction_names.append("CEvNS40Ar.react")
        if "ES" in process_types:
            reaction_names.append("ES.react")

    primary_processes_dict = collections.defaultdict(list)

    for reaction_name in reaction_names:
        react_fname = _find_file(marley_search_path, reaction_name)
        nuclide_index_fname = _find_file(marley_search_path, "nuclide_index.txt")
        nuclide_path = os.path.dirname(nuclide_index_fname)
        nuclide_fnames = sorted(glob.glob(os.path.join(nuclide_path, "*.dat")))
        masses_fname = _find_file(marley_search_path, "mass_table.js")
        gs_parity_fname = _find_file(marley_search_path, "gs_spin_parity_table.txt")

        aux_files = [_find_file(marley_search_path, rel_src)
                     for rel_src, _ in _runtime_aux]
        aux_names = [name for _, name in _runtime_aux]

        xs = siren.interactions.MarleyCrossSection(
            [react_fname], nuclide_index_fname, nuclide_fnames,
            masses_fname, gs_parity_fname, aux_files, aux_names)

        reaction_primary_types = set(primaries_by_process[process_by_reaction[reaction_name]])
        reaction_primary_types = reaction_primary_types & set(primary_types)

        #Add cross section to primary_processes_dict
        for primary_type in reaction_primary_types:
            primary_processes_dict[primary_type].append(xs)

    return dict(primary_processes_dict), {}
