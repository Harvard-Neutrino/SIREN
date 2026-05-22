from scipy.interpolate import interp1d
import numpy as np
import os
from siren.download import ensure_files, writable_data_dir

_ABS_DIR = writable_data_dir(os.path.dirname(os.path.abspath(__file__)))
_DATA_BASE = "https://raw.githubusercontent.com/SIREN-Generator/SIREN-data/main/fluxes/T2K_Kaons/T2K_Kaons-v1.0"

_DATA_FILES = [
    {"path": os.path.join(_ABS_DIR, "kaon-flux-data.dat"), "url": f"{_DATA_BASE}/kaon-flux-data.dat",
     "sha256": "32e06bb7f454547fa30ac4f6f34caf50700c53683ec40963a904c1fb91f87d56"},
    {"path": os.path.join(_ABS_DIR, "ratio.dat"), "url": f"{_DATA_BASE}/ratio.dat",
     "sha256": "ae85fd610073a6158065b4b24ab351d42a353ad7f704d1dbb58d3c6cdf911e98"},
    {"path": os.path.join(_ABS_DIR, "TOT_PLUS_NUMU.dat"), "url": f"{_DATA_BASE}/TOT_PLUS_NUMU.dat",
     "sha256": "fdce382924d17d4841c1b2e8da6dcc7932d7a3dc9c5ee68b0a438035083a8b20"},
    {"path": os.path.join(_ABS_DIR, "TOT_MINUS_NUMUBAR.dat"), "url": f"{_DATA_BASE}/TOT_MINUS_NUMUBAR.dat",
     "sha256": "347a94e69261a25b3fb5174c193e0abbae167e817f5d5e50c9da213270e26531"},
]


def fetch_data():
    ensure_files(_DATA_FILES)


def bar_scaling():
    energy, ratio = np.loadtxt(os.path.join(_ABS_DIR, 'ratio.dat'), usecols=(0, 1), unpack=True)
    return interp1d(energy, ratio, kind='linear', bounds_error=False, fill_value=(ratio[0], ratio[-1]))

def MakeFluxFile(tag, output_dir=None):
    if output_dir is None:
        output_dir = _ABS_DIR

    parts = tag.split("_")
    if len(parts) != 2:
        raise ValueError(f"Tag must be '{{particle}}_{{PLUS|MINUS}}', got '{tag}'")
    particle, enhance = parts

    if enhance == 'MINUS':
        bar = True
        bar_scale = bar_scaling()
    elif enhance == 'PLUS':
        bar = False
    else:
        raise ValueError(f'Enhance tag "{enhance}" is invalid, expected PLUS or MINUS')

    if particle not in ["numu", "numubar"]:
        raise ValueError(f"\"{particle}\" particle specified in tag \"{tag}\" is not valid")

    input_flux_file = os.path.join(_ABS_DIR, "kaon-flux-data.dat")

    if bar:
        output_flux_file = os.path.join(output_dir, f"kaon-flux-{particle}_bar.dat")
    else:
        output_flux_file = os.path.join(output_dir, f"kaon-flux-{particle}.dat")

    with open(input_flux_file,"r") as fin:
        all_lines = fin.readlines()
        data = [line.strip().split() for line in all_lines[:]]

        with open(output_flux_file,"w") as fout:
            for row in data:
                E = float(row[0])
                flux = 10**float(row[1]) * bar_scale(float(row[0])) if bar else 10**float(row[1])
                flux*=2e-16 # put flux in units of nu/m^2/GeV/POT
                print(E, flux, file=fout)

    return output_flux_file

def load_flux(tag, abs_flux_dir=None):
    fetch_data()
    return MakeFluxFile(tag, output_dir=abs_flux_dir)
