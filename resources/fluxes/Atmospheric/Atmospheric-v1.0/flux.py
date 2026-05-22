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
    # Match model by longest known prefix (handles underscores in names
    # like H3a_SIBYLL21)
    model = None
    for m in sorted(MODELS, key=len, reverse=True):
        if tag.startswith(m + "_"):
            model = m
            break
    if model is None:
        raise ValueError(
            f"Tag '{tag}' does not start with a known model. "
            f"Expected one of: {', '.join(MODELS)}")
    rest = tag[len(model) + 1:]
    parts = rest.split("_")
    if len(parts) < 2:
        raise ValueError(
            f"Tag '{tag}' is not valid. Expected "
            "{{model}}_{{osc_status}}_{{particle1}}_{{particle2}}_...")
    osc_status = parts[0]
    particles = parts[1:]
    if osc_status not in OSC_STATES:
        raise ValueError(f"Unknown osc_status '{osc_status}' in tag '{tag}'")
    for p in particles:
        if p not in PARTICLES:
            raise ValueError(f"Unknown particle '{p}' in tag '{tag}'")

    fetch_data()

    output_flux_file = os.path.join(_ABS_DIR, f"Atmospheric_{tag}_flux.txt")
    if os.path.isfile(output_flux_file):
        return output_flux_file

    npz_path = os.path.join(_ABS_DIR, f"{model}.npz")
    with np.load(npz_path) as npz:
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

    ee, cc = np.meshgrid(energy, cos_theta, indexing='ij')
    np.savetxt(output_flux_file,
               np.column_stack([ee.ravel(), cc.ravel(), tot_flux.ravel()]))

    return output_flux_file
