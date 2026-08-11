import os
import glob
from typing import List, Optional
import siren
import collections

# MarleyCrossSection-v2.0: MARLEY v2 hybrid model resources.
#
# Differences with respect to v1.0:
#  - Self-contained: every data file ships INSIDE this bundle (react files,
#    CRPA response tables, nuclear charge radii, v2-format structure files
#    with per-level half-lives, masses, parities). No dependency on $PREFIX
#    or on an external MARLEY installation.
#  - The CC process is the MARLEY v2 recommended hybrid: measured discrete
#    levels (Bhattacharya 2009) + HF-CRPA continuum, loaded together as a
#    single MarleyCrossSection object (multi-react constructor).
basepath = os.path.dirname(os.path.abspath(__file__))
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

#list of processes and the react files that build each one.
#CC is the v2 hybrid: BOTH files form a single cross-section object.
processes = ["CC", "CEvNS", "ES"]
reactions_by_process = {
    "CC": ["ve40ArCC_HF-CRPA.react", "ve40ArCC_Bhattacharya2009-Discrete.react"],
    "CEvNS": ["CEvNS40Ar.react"],
    "ES": ["ES.react"],
}

#auxiliary data files needed by the v2 engine, with the relative name each
#one must keep on the MARLEY search path (the HF-CRPA manifest references the
#response tables as "crpa/<table>.dat")
aux_by_process = {
    "CC": [(os.path.join("react", "crpa", os.path.basename(f)), os.path.join("crpa", os.path.basename(f)))
           for f in sorted(glob.glob(os.path.join(basepath, "react", "crpa", "*.dat")))]
          + [("nuclear_charge_radii.js", "nuclear_charge_radii.js")],
    "CEvNS": [("nuclear_charge_radii.js", "nuclear_charge_radii.js")],
    "ES": [],
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
        marley_search_path = default_marley_search_path

    primary_types = _get_primary_types(primary_types)
    process_types = _get_process_types(process_types)

    primary_processes_dict = collections.defaultdict(list)

    for process in process_types:
        #Resolve every react file of this process (CC = hybrid pair)
        react_fnames = [_find_file(marley_search_path, r) for r in reactions_by_process[process]]

        nuclide_index_fname = _find_file(marley_search_path, "nuclide_index.txt")
        nuclide_path = os.path.dirname(nuclide_index_fname)
        nuclide_fnames = glob.glob(os.path.join(nuclide_path, "*.dat"))
        masses_fname = _find_file(marley_search_path, "mass_table.js")
        gs_parity_fname = _find_file(marley_search_path, "gs_spin_parity_table.txt")

        #Auxiliary files: (path inside the bundle, relative name for MARLEY)
        aux_files = [os.path.join(basepath, rel_src) for rel_src, _ in aux_by_process[process]]
        aux_names = [name for _, name in aux_by_process[process]]

        xs = siren.interactions.MarleyCrossSection(
            react_fnames, nuclide_index_fname, nuclide_fnames,
            masses_fname, gs_parity_fname, aux_files, aux_names)

        reaction_primary_types = set(primaries_by_process[process]) & set(primary_types)

        #Add cross section to primary_processes_dict
        for primary_type in reaction_primary_types:
            primary_processes_dict[primary_type].append(xs)

    return dict(primary_processes_dict), {}
