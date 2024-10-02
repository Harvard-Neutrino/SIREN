import os
from typing import Tuple, List, Any, Optional
import siren
import collections


basepath = os.environ['PREFIX']
basepath = os.path.join(basepath, 'share/marley')
marley_search_path = ':'.join([
    basepath,
    os.path.join(basepath, 'react'),
    os.path.join(basepath, 'structure')
    ])

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

processes = ["CC", "CEvNS", "ES"]
process_by_reaction = {
    "ve40ArCC_Bhattacharya1998.react": "CC",
    "ve40ArCC_Liu1998.react": "CC",
    "ve40ArCC_Bhattacharya2009.react": "CC",
    "CEvNS40Ar.react": "CEvNS",
    "ES.react": "ES",
}

primaries_by_process = {
    "CC": [siren.dataclasses.Particle.ParticleType.NuE, siren.dataclasses.Particle.ParticleType.NuEBar],
    "ES": neutrinos + antineutrinos,
    "CEvNS": neutrinos + antineutrinos,
}

def _get_primary_types(primary_types):
    if primary_types is None:
        primary_types = [
                siren.dataclasses.Particle.ParticleType.NuE,
        ]

    supported_primaries = [siren.dataclasses.Particle.ParticleType.NuE, siren.dataclasses.Particle.ParticleType.NuEBar]
    for i, p in enumerate(primary_types):
        if p not in supported_primaries:
            raise ValueError(f"primary_types[{i}] \"{p}\" not supported")

    if len(primary_types) == 0:
        print("Warning: len(primary_types) == 0")

    return primary_types

def _get_process_types(process_types):
    if process_types is None:
        process_types = ["CC", "CEvNS", "ES"]

    for i, p in enumerate(process_types):
        if p not in processes:
            raise ValueError(f"process_types[{i}] \"{p}\" not supported. Allowed processes are: {processes}")

    if len(process_types) == 0:
        print("Warning: len(process_types) == 0")

    return process_types

def load_processes(
    primary_types: Optional[List[siren.dataclasses.Particle.ParticleType]] = None,
    process_types: Optional[List[str]] = None,
    marley_search_path: Optional[str] = None,
    reaction_name: Optional[str] = None,
    ):

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

        xs = siren.interactions.MarleyCrossSection(f"{{reactions: [\"{reaction_name}\"]}}", marley_search_path)
        reaction_primary_types = set(primaries_by_process[process_by_reaction[reaction_name]])
        reaction_primary_types = reaction_primary_types & set(primary_types)

        for primary_type in reaction_primary_types:
            primary_processes_dict[primary_type].append(xs)

    return dict(primary_processes_dict), {}

