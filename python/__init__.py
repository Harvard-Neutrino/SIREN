from . import utilities
from . import math
from . import dataclasses
from . import geometry
from . import detector
from . import interactions
from . import distributions
from . import injection

from . import _util

# set up some public-facing utilities functions
utilities.get_resource_package_dir = _util.resource_package_dir
utilities.get_detector_model_path = _util.get_detector_model_path
utilities.get_material_model_path = _util.get_material_model_path
utilities.get_cross_section_model_path = _util.get_cross_section_model_path
utilities.get_flux_model_path = _util.get_flux_model_path
utilities.load_flux = _util.load_flux

def darknews_version():
    try:
        import DarkNews
        return _util.normalize_version(DarkNews.__version__)
    except:
        print("WARNING: DarkNews is not installed in the local environment")
        return None
utilities.darknews_version = darknews_version

