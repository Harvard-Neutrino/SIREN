import os
from siren.download import ensure_files, writable_data_dir

_ABS_DIR = writable_data_dir(os.path.dirname(os.path.abspath(__file__)))
_DATA_BASE = "https://raw.githubusercontent.com/SIREN-Generator/SIREN-data/main/fluxes/T2K_NEAR/T2K_NEAR-v1.0"

_DATA_FILES = [
    {"path": os.path.join(_ABS_DIR, "T2K_PLUS_250kA.dat"), "url": f"{_DATA_BASE}/T2K_PLUS_250kA.dat",
     "sha256": "bc7958cda04788976d5c0e75a5efca075b62d63d82cce0a67f7fea539bb33eab"},
    {"path": os.path.join(_ABS_DIR, "T2K_MINUS_250kA.dat"), "url": f"{_DATA_BASE}/T2K_MINUS_250kA.dat",
     "sha256": "e2d6bedd56f4c30ad765ec0320970565ad2d63faf5cb4cdd9004d4f6cd965ffe"},
]


def fetch_data():
    ensure_files(_DATA_FILES)


def load_flux(tag):

    '''
    Accepts the following tags:
        {PLUS, MINUS}_{nue,nuebar,numu,numubar}
    '''

    fetch_data()

    enhance, particle = tag.split("_")

    if enhance not in ["MINUS", "PLUS"]:
        raise ValueError("%s 250kA enhancement specified in tag %s is not valid" % (enhance, tag))
    if particle not in ["numu", "numubar", "nue", "nuebar"]:
        raise ValueError("%s particle specified in tag %s is not valid" % (particle, tag))

    input_flux_file = os.path.join(_ABS_DIR,
                                   "T2K_%s_250kA.dat"%(enhance))

    output_flux_file = os.path.join(_ABS_DIR,
                                    "T2KOUT_%s.dat"%(tag))

    with open(input_flux_file,"r") as fin:
        all_lines = fin.readlines()
        headers = all_lines[1].strip().split()
        data = [line.strip().split() for line in all_lines[3:]]
        pid = headers.index(particle)
        with open(output_flux_file,"w") as fout:
            for row in data:
                E, flux = (float(row[1])+float(row[3]))/2, float(row[pid+2])
                flux*=2e-16 # put flux in units of nu/m^2/GeV/POT
                print(E, flux, file=fout)
    return output_flux_file
