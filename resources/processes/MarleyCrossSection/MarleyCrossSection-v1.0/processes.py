import os
import glob
from typing import Tuple, List, Any, Optional
import siren
import collections

#basepath and search path for marley reaction files
basepath = os.environ['PREFIX']
basepath = os.path.join(basepath, 'share/marley')
default_marley_search_path = ':'.join([
    basepath,
    os.path.join(basepath, 'react'),
    os.path.join(basepath, 'structure')
    ])

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
    "ve40ArCC_Bhattacharya1998.react": "CC",
    "ve40ArCC_Liu1998.react": "CC",
    "ve40ArCC_Bhattacharya2009.react": "CC",
    "CEvNS40Ar.react": "CEvNS",
    "ES.react": "ES",
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
    reaction_name: Optional[str] = None,
    ):

    if marley_search_path is None:
        marley_search_path = default_marley_search_path

    primary_types = _get_primary_types(primary_types)
    process_types = _get_process_types(process_types)
    reaction_names = []

    if reaction_name is not None:
        reaction_names = [reaction_name]
    else:
        if "CC" in process_types and (siren.dataclasses.Particle.ParticleType.NuE in primary_types or siren.dataclasses.Particle.ParticleType.NuEBar in primary_types):
            reaction_names.append("ve40ArCC_Bhattacharya2009.react")
        if "CEvNS" in process_types:
            reaction_names.append("CEvNS40Ar.react")
        if "ES" in process_types:
            reaction_names.append("ES.react")

    primary_processes_dict = collections.defaultdict(list)

    for reaction_name in reaction_names:
        #Load cross section from marley reaction file
        react_fname = _find_file(marley_search_path, reaction_name)
        nuclide_index_fname = _find_file(marley_search_path, "nuclide_index.txt")
        nuclide_path = os.path.dirname(nuclide_index_fname)
        nuclide_fnames = glob.glob(os.path.join(nuclide_path, "*.dat"))
        masses_fname = _find_file(marley_search_path, "mass_table.js")
        gs_parity_fname = _find_file(marley_search_path, "gs_spin_parity_table.txt")

        xs = siren.interactions.MarleyCrossSection(react_fname, nuclide_index_fname, nuclide_fnames, masses_fname, gs_parity_fname)
        reaction_primary_types = set(primaries_by_process[process_by_reaction[reaction_name]])
        reaction_primary_types = reaction_primary_types & set(primary_types)

        #Add cross section to primary_processes_dict
        for primary_type in reaction_primary_types:
            primary_processes_dict[primary_type].append(xs)

    return dict(primary_processes_dict), {}

