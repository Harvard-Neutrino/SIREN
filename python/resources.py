__all__ = ["load_flux", "load_detector", "load_processes", "fluxes", "detectors", "processes"]

from . import _util
import os
import sys
import importlib
import logging
from importlib.util import spec_from_loader

# Set up logging configuration
logger = logging.getLogger(__name__)
if not logger.handlers:
    #handler = logging.FileHandler("./resources.log", mode="w")
    handler = logging.StreamHandler()
    #formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
    #handler.setFormatter(formatter)
    handler.setFormatter(_util.CustomFormatter())
    logger.addHandler(handler)
    logger.setLevel(logging.WARN)

load_flux = _util.load_flux
load_detector = _util.load_detector
load_processes = _util.load_processes

class ResourceList:
    def __init__(self, resource_type, list_method, load_method):
        self.__resource_type = resource_type
        self.__list_method = list_method
        self.__load_method = load_method
        self.__path__ = []
        self.__name__ = __name__ + "." + resource_type
        self.__package__ = __name__ + "." + resource_type

    def __getattr__(self, name):
        logger.debug("ResourceList(%r).__getattr__(%r)", self.__resource_type, name)
        resources = self.__list_method()
        logger.debug("Resources: %r", resources)
        if name in resources:
            logger.debug("%r is listed in resources, attempting to load", name)
            load_res = self.__load_method(name)
            logger.debug(f"Loaded resource: {load_res}")
            return load_res
        else:
            # Default behaviour
            return object.__getattribute__(self, name)

    def __hasattr__(self, name):
        logger.debug("ResourceList(%r).__hasattr__(%r)", self.__resource_type, name)
        resources = self.__list_method()
        logger.debug("Resources: %r", resources)
        if name in resources:
            return True
        else:
            # Default behaviour
            try:
                object.__getattribute__(self, name)
                return True
            except AttributeError:
                return False

    def __dir__(self):
        dirs = dir(self.__class__)
        dirs += list(self.__dict__.keys())
        dirs += self.__list_method()
        return sorted(dirs)

fluxes = ResourceList('fluxes', _util.list_fluxes, _util._get_flux_loader)
detectors = ResourceList('detectors', _util.list_detectors, _util._get_detector_loader)
processes = ResourceList('processes', _util.list_processes, _util._get_process_loader)

class _SIRENResourcesMetaPathImporter(object):
    """
    A meta path importer to import siren.resources and its submodules.
    """

    def __init__(self, siren_module_name):
        logger.debug("Importer init")
        self.name = siren_module_name
        self.known_modules = {}

    def _add_module(self, mod, *fullnames):
        logger.debug("Adding module %r with fullnames %r", mod, fullnames)
        for fullname in fullnames:
            self.known_modules[self.name + "." + fullname] = mod

    def _get_module(self, fullname):
        logger.debug("Getting module %r", fullname)
        return self.__get_module(fullname)

    def find_module(self, fullname, path=None):
        logger.debug("Finding module %r with path %r", fullname, path)
        mod = self.__get_module(fullname)
        return self
        try:
            mod = self.__get_module(fullname)
            return self
        except ImportError:
            return None
        return None

    def find_spec(self, fullname, path, target=None):
        logger.debug("Finding spec for %r with path %r and target %r", fullname, path, target)
        mod = self.__get_module(fullname)
        return spec_from_loader(fullname, self)
        try:
            mod = self.__get_module(fullname)
            return spec_from_loader(fullname, self)
        except ImportError:
            return None

    def __get_module(self, fullname):
        logger.debug("Getting module %r", fullname)
        if fullname in self.known_modules:
            mod = self.known_modules[fullname]
            logger.debug("Found module %r in known_modules", fullname)
            return mod

        if not fullname.startswith(self.name + "."):
            logger.debug("%r does not start with %r", fullname, self.name + ".")
            raise ImportError("This loader does not know module " + fullname)

        subname = fullname[len(self.name) + 1:]
        split_name = subname.split(".")
        mod = None
        for i in reversed(range(len(split_name))):
            name = self.name + "." + ".".join(split_name[:i + 1])
            logger.debug("Testing name: %r", name)
            if name in self.known_modules:
                logger.debug("%r found in known_modules", name)
                parent_mod = self.known_modules[name]
                success = True
                for j in range(i + 1, len(split_name)):
                    logger.debug("parent_mod: %r", parent_mod)
                    logger.debug("parent_mod name: %r", self.name + '.' + '.'.join(split_name[:j]))
                    submod = split_name[j]
                    logger.debug("Looking for %r in parent_mod", submod)
                    logger.debug("parent_mod dir: %r", dir(parent_mod))
                    logger.debug("parent_mod has %r: %r", submod, hasattr(parent_mod, submod))
                    try:
                        parent_mod = getattr(parent_mod, submod)
                    except AttributeError:
                        logger.debug("parent_mod does not have %r", submod)
                        success = False
                        break
                    logger.debug("parent_mod does have %r", submod)
                    parent_mod.__path__ = []
                    parent_mod.__name__ = self.name + "." + ".".join(split_name[:j+1])
                    if os.path.isdir(os.path.join(os.path.dirname(sys.modules[self.name].__file__), *split_name[:j+1])):
                        parent_mod.__package__ = self.name + "." + ".".join(split_name[:j+1])
                    else:
                        parent_mod.__package__ = self.name + "." + ".".join(split_name[:j])
                    self.known_modules[parent_mod.__name__] = parent_mod
                if success:
                    logger.debug("Found the module: %r", parent_mod)
                    mod = parent_mod
                    self.known_modules[name] = mod
                    break
            else:
                logger.debug("%r not in known_modules", name)
        if mod is None:
            raise ImportError("This loader does not know module " + fullname)

        return mod

    def load_module(self, fullname):
        logger.debug("Loading module %r", fullname)
        try:
            # In case of a reload
            return sys.modules[fullname]
        except KeyError:
            pass
        mod = self.__get_module(fullname)
        mod.__loader__ = self
        sys.modules[fullname] = mod
        return mod

    def is_package(self, fullname):
        """
        Return True if the named module is a package.
        (Required for proper spec objects with Python 3.4+, see PEP451)
        """
        logger.debug("Is package: %r", fullname)
        return hasattr(self.__get_module(fullname), "__path__")

    def get_code(self, fullname):
        """Return None (required if is_package is implemented)"""
        logger.debug("Get code for %r", fullname)
        self.__get_module(fullname)  # May raise ImportError
        return None
    get_source = get_code  # Alias for get_code

    def create_module(self, spec):
        logger.debug("Create module with spec: %r", spec)
        return self.load_module(spec.name)

    def exec_module(self, module):
        logger.debug("Exec module: %r", module)
        pass

_importer = _SIRENResourcesMetaPathImporter(__name__)

_importer._add_module(fluxes, "fluxes")
_importer._add_module(detectors, "detectors")
_importer._add_module(processes, "processes")

logger.debug("Module name: %r", __name__)

# Complete the moves implementation.
# This code is at the end of this module to speed up module loading.
# Turn this module into a package.
__path__ = []  # required for PEP 302 and PEP 451
__package__ = __name__  # see PEP 366 @ReservedAssignment
if globals().get("__spec__") is not None:
    __spec__.submodule_search_locations = []  # PEP 451 @UndefinedVariable

# Remove other meta path importers of the same type to avoid conflicts.
if sys.meta_path:
    for i, importer in enumerate(sys.meta_path):
        # Since there could be multiple instances of this importer (from different six modules),
        # we check by name rather than isinstance().
        if (type(importer).__name__ == "_SIRENResourcesMetaPathImporter" and
                importer.name == __name__):
            del sys.meta_path[i]
            break
    del i, importer

# Finally, add the importer to the meta path import hook.
sys.meta_path.append(_importer)

del _util
del ResourceList

