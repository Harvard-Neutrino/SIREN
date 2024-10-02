from . import utilities
from . import math
from . import dataclasses
from . import geometry
from . import detector
from . import interactions
from . import distributions
from . import injection

from . import _util
from . import Injector
from . import Weighter
from . import resources

# Intropspect package version
import sys
if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    import importlib_metadata as metadata
__version__ = metadata.version(__package__)
del sys
del metadata

# set up some public-facing utilities functions
utilities.get_resource_package_dir = _util.resource_package_dir
utilities.get_detector_model_path = _util.get_detector_model_path
utilities.get_processes_model_path = _util.get_processes_model_path
utilities.get_flux_model_path = _util.get_flux_model_path
utilities.load_flux = _util.load_flux
utilities.load_detector = _util.load_detector
utilities.load_processes = _util.load_processes
utilities.get_fiducial_volume = _util.get_fiducial_volume

# Override the Injector with the python wrapper
injection._Injector = injection.Injector
injection.Injector = Injector.Injector
del Injector

# Override the Weighter with the python wrapper
injection._Weighter = injection.Weighter
injection.Weighter = Weighter.Weighter
del Weighter

dataclasses.Particle.ParticleType = dataclasses.ParticleType

def darknews_version():
    try:
        import DarkNews
        return _util.normalize_version(DarkNews.__version__)
    except:
        print("WARNING: DarkNews is not installed in the local environment")
        return None
utilities.darknews_version = darknews_version

