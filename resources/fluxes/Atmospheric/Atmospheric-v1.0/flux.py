import os
import numpy as np
from siren.download import ensure_files, writable_data_dir

_ABS_DIR = writable_data_dir(os.path.dirname(os.path.abspath(__file__)))

_DATA_BASE = (
    "https://raw.githubusercontent.com/SIREN-Generator/SIREN-data/"
    "main/fluxes/Atmospheric/Atmospheric-v1.0"
)

MODELS = ["bartol", "daemonflux", "H3a_SIBYLL21", "H3a_SIBYLL23C", "honda2006"]
PARTICLES = ["nue", "numu", "nuebar", "numubar", "nutau", "nutaubar"]
OSC_STATES = ["osc", "unosc"]

_NPZ_SHA256 = {
    "bartol": "af0806ea7a6dae9fd0e8196aa4d3896ca33f19638bc9a5e688921f9aa4c12ee4",
    "daemonflux": "40fcec18c4f66ff556b7c725b2ff84559f7e9a7cd6e150efe7fd120285af57a9",
    "H3a_SIBYLL21": "9f3280c76c1f3558a4522647a395131cb9f84020dfd18f128fc2b7bc612bfc96",
    "H3a_SIBYLL23C": "e291e782a2d149a7e979c7332b63e479130c801db6b355d7b140b06968fd84ae",
    "honda2006": "9afe3fe86e67e6653e95ae1e9629964a73d1063c5abaf6e39eabf5444c964559",
}

_NPZ_FILES = [
    {
        "path": os.path.join(_ABS_DIR, f"{model}.npz"),
        "url": f"{_DATA_BASE}/{model}.npz",
        "sha256": _NPZ_SHA256[model],
    }
    for model in MODELS
]


def fetch_data():
    """Download all atmospheric flux npz files."""
    ensure_files(_NPZ_FILES)


def load_flux(tag):
    '''
    Accepts the following tags:
        {model}_{osc_status}_{particle1}_{particle2}_...{particleN}
        where model is one of {bartol,daemonflux,H3a_SIBYLL21,H3a_SIBYLL23C,honda2006}
        osc_status is one of {osc,unosc}
        particlei is one of {nue,numu,nuebar,numubar,nutau,nutaubar}
    '''
    fields = tag.split("_")
    if len(fields) < 3:
        raise ValueError(
            f"Tag '{tag}' is not valid. Expected "
            "{{model}}_{{osc_status}}_{{particle1}}_{{particle2}}_...")
    model, osc_status = fields[0], fields[1]
    particles = fields[2:]
    if model not in MODELS:
        raise ValueError(f"Unknown model '{model}' in tag '{tag}'")
    if osc_status not in OSC_STATES:
        raise ValueError(f"Unknown osc_status '{osc_status}' in tag '{tag}'")
    for p in particles:
        if p not in PARTICLES:
            raise ValueError(f"Unknown particle '{p}' in tag '{tag}'")

    npz_spec = {"path": os.path.join(_ABS_DIR, f"{model}.npz"),
                "url": f"{_DATA_BASE}/{model}.npz",
                "sha256": _NPZ_SHA256.get(model, "")}
    ensure_files([npz_spec])

    output_flux_file = os.path.join(_ABS_DIR, f"Atmospheric_{tag}_flux.txt")
    if os.path.isfile(output_flux_file):
        return output_flux_file

    npz = np.load(npz_spec["path"])
    energy = npz["energy"]
    cos_theta = npz["cos_theta"]

    tot_flux = None
    for particle in particles:
        key = f"{particle}_{osc_status}"
        if key not in npz:
            raise ValueError(f"Variant '{key}' not found in {model}.npz")
        grid = npz[key].astype(np.float64)
        if tot_flux is None:
            tot_flux = grid.copy()
        else:
            tot_flux += grid

    with open(output_flux_file, "w") as fout:
        for i, e in enumerate(energy):
            for j, c in enumerate(cos_theta):
                print(e, c, tot_flux[i, j], file=fout)

    return output_flux_file
