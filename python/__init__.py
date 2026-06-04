from . import utilities
from . import math
from . import dataclasses
from . import geometry
from . import detector
from . import interactions
from . import distributions
from . import injection

from . import optimize

from . import _util
from . import resources
from . import visualization

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
utilities.get_cross_section_model_path = _util.get_processes_model_path
utilities.get_flux_model_path = _util.get_flux_model_path
utilities.load_flux = _util.load_flux
utilities.load_detector = _util.load_detector
utilities.load_processes = _util.load_processes
utilities.get_fiducial_volume = _util.get_fiducial_volume

# Override the Injector with the python wrapper
injection._Injector = injection.Injector
del injection.Injector
from . import Injector
injection.Injector = Injector.Injector
del Injector

# Override the Weighter with the python wrapper
injection._Weighter = injection.Weighter
del injection.Weighter
from . import Weighter
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

# ====================================================================== #
#  Public high-level API (Phases 1-3 of interface redesign)               #
# ====================================================================== #

# ---- Sub-namespaces (tab-completable) ----
from . import particles
from . import dist

# ---- Simulation and Results ----
from .Simulation import Simulation
from .Results import Results

# ---- Top-level convenience functions ----
load_detector = utilities.load_detector
load_flux = utilities.load_flux
load_processes = utilities.load_processes
get_fiducial_volume = utilities.get_fiducial_volume

ProcessBundle = _util.ProcessBundle
GenerateEvents = _util.GenerateEvents
SaveEvents = _util.SaveEvents
LoadEvents = _util.LoadEvents
get_volume_position_distribution_from_sector = _util.get_volume_position_distribution_from_sector
