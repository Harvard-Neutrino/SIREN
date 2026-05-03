from scipy.interpolate import interp1d
import numpy as np
import os
import pathlib

def bar_scaling(abs_flux_dir):
    energy, ratio = np.loadtxt(os.path.join(abs_flux_dir, 'ratio.dat'), usecols=(0, 1), unpack=True)
    return interp1d(energy, ratio, kind='linear', bounds_error=False, fill_value=(ratio[0], ratio[-1]))

def MakeFluxFile(tag, abs_flux_dir):
    parts = tag.split("_")
    if len(parts) != 2:
        raise ValueError(f"Tag must be '{{particle}}_{{PLUS|MINUS}}', got '{tag}'")
    particle, enhance = parts

    if enhance == 'MINUS':
        bar = True
        bar_scale = bar_scaling(abs_flux_dir)
    elif enhance == 'PLUS':
        bar = False
    else:
        raise ValueError(f'Enhance tag "{enhance}" is invalid, expected PLUS or MINUS')

    if particle not in ["numu", "numubar"]:
        raise ValueError(f"\"{particle}\" particle specified in tag \"{tag}\" is not valid")

    input_flux_file = os.path.join(abs_flux_dir,
                                    "kaon-flux-data.dat")

    output_flux_file = os.path.join(abs_flux_dir,
                                    f"kaon-flux-{particle}_bar.dat") if bar else os.path.join(abs_flux_dir, f"kaon-flux-{particle}.dat")


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
    if abs_flux_dir is None:
        abs_flux_dir = str(pathlib.Path(__file__).resolve().parent)
    return MakeFluxFile(tag, abs_flux_dir)
