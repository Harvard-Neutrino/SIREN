import os
from typing import Tuple, List, Any, Optional
import siren
import collections

base_path = os.path.dirname(os.path.abspath(__file__))

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
processes = ["CC", "NC"]


def _get_primary_types(primary_types):
    if primary_types is None:
        primary_types = [
                siren.dataclasses.Particle.ParticleType.NuE,
                siren.dataclasses.Particle.ParticleType.NuMu,
                siren.dataclasses.Particle.ParticleType.NuTau,
                siren.dataclasses.Particle.ParticleType.NuEBar,
                siren.dataclasses.Particle.ParticleType.NuMuBar,
                siren.dataclasses.Particle.ParticleType.NuTauBar,
        ]

    supported_primaries = neutrinos + antineutrinos
    for i, p in enumerate(primary_types):
        if p not in supported_primaries:
            raise ValueError(f"primary_types[{i}] \"{p}\" not supported")

    if len(primary_types) == 0:
        print("Warning: len(primary_types) == 0")

    return primary_types

def _get_isoscalar(isoscalar):
    if isoscalar is None:
        isoscalar = True

    if not isoscalar:
        raise ValueError("Non-isoscalar splines are not supported for CSMSDISSplines-v1.0")

    return isoscalar


def _get_target_types(isoscalar, target_types):
    if target_types is None:
        if isoscalar:
            target_types = [siren.dataclasses.Particle.ParticleType.Nucleon]
        else:
            target_types = [siren.dataclasses.Particle.ParticleType.PPlus, siren.dataclasses.Particle.ParticleType.Neutron]

    if len(target_types) == 0:
        print("Warning: len(target_types) == 0")

    return target_types

def _get_process_types(process_types):
    if process_types is None:
        process_types = ["CC", "NC"]

    for i, p in enumerate(process_types):
        if p not in processes:
            raise ValueError(f"process_types[{i}] \"{p}\" not supported")

    if len(process_types) == 0:
        print("Warning: len(process_types) == 0")

    return process_types


def load_processes(
    primary_types: Optional[List[siren.dataclasses.Particle.ParticleType]] = None,
    target_types: Optional[List[siren.dataclasses.Particle.ParticleType]] = None,
    isoscalar: Optional[bool] = None,
    process_types: Optional[List[str]] = None,
    ):

    primary_types = _get_primary_types(primary_types)
    isoscalar = _get_isoscalar(isoscalar)
    target_types = _get_target_types(isoscalar, target_types)
    process_types = _get_process_types(process_types)

    neutrino_types = [t for t in primary_types if t in neutrinos]
    antineutrino_types = [t for t in primary_types if t in antineutrinos]

    primary_processes = []
    primary_processes_dict = collections.defaultdict(list)

    for process_type in process_types:
        for primaries, nunubar in [[neutrinos, "nu"], [antineutrinos, "nubar"]]:
            if isoscalar:
                dxs_file = os.path.join(base_path, f"dsdxdy_{nunubar}_{process_type}_iso.fits")
                xs_file = os.path.join(base_path, f"sigma_{nunubar}_{process_type}_iso.fits")
                xs = siren.interactions.DISFromSpline(dxs_file, xs_file, primaries, target_types, "m")
                primary_processes.append(xs)
                for primary_type in primaries:
                    primary_processes_dict[primary_type].append(xs)
            else:
                pass

    return dict(primary_processes_dict), {}

