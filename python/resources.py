__all__ = ["load_flux", "load_detector", "load_processes", "fluxes", "detectors", "processes"]

from . import _util
import os
import sys
import importlib
from importlib.util import spec_from_loader

output = sys.stdout

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
        print(f"ResourceList(\"{self.__resource_type}\").__getattr__(\"{name}\")")
        resources = self.__list_method()
        print(f"Resources: {resources}")
        if name in resources:
            return self.__load_method(name)
        else:
            # Default behaviour
            return object.__getattribute__(self, name)

    def __hasattr__(self, name):
        print(f"ResourceList(\"{self.__resource_type}\").__getattr__(\"{name}\")")
        resources = self.__list_method()
        print(f"Resources: {resources}")
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
        print("Importer init", file=output)
        self.name = siren_module_name
        self.known_modules = {}

    def _add_module(self, mod, *fullnames):
        print("Adding module", mod, fullnames, file=output)
        for fullname in fullnames:
            self.known_modules[self.name + "." + fullname] = mod

    def _get_module(self, fullname):
        print("Getting module", fullname, file=output)
        return self.__get_module(fullname)

    def find_module(self, fullname, path=None):
        print("Finding module", fullname, path, file=output)
        mod = self.__get_module(fullname)
        return self
        try:
            mod = self.__get_module(fullname)
            return self
        except ImportError:
            return None
        return None

    def find_spec(self, fullname, path, target=None):
        print("Finding spec", fullname, path, target, file=output)
        mod = self.__get_module(fullname)
        return spec_from_loader(fullname, self)
        try:
            mod = self.__get_module(fullname)
            return spec_from_loader(fullname, self)
        except ImportError:
            return None

    def __get_module(self, fullname):
        print("Getting module", fullname, file=output)
        try:
            return self.known_modules[fullname]
        except KeyError:
            if not fullname.startswith(self.name + "."):
                raise ImportError("This loader does not know module " + fullname)

            subname = fullname[len(self.name) + 1:]
            split_name = subname.split(".")
            mod = None
            for i in reversed(range(len(split_name))):
                name = self.name + "." + ".".join(split_name[:i + 1])
                print(f"Testing name: {name}")
                if name in self.known_modules:
                    print(f"{name} in known_modules")
                    parent_mod = self.known_modules[name]
                    success = True
                    for j in range(i + 1, len(split_name)):
                        print(f"parent_mod: {parent_mod}")
                        print(f"parent_mod name: {self.name + '.' + '.'.join(split_name[:j])}")
                        submod = split_name[j]
                        print(f"Looking for {submod} in parent_mod")
                        print(dir(parent_mod))
                        print(f"parent_mod has {submod}: {hasattr(parent_mod, submod)}")
                        try:
                            parent_mod = getattr(parent_mod, submod)
                        except AttributeError:
                            print(f"parent_mod does not have {submod}")
                            success = False
                            break
                        print(f"parent_mod does have {submod}")
                        parent_mod.__path__ = []
                        parent_mod.__name__ = self.name + "." + ".".join(split_name[:j+1])
                        if os.path.isdir(os.path.join(os.path.dirname(sys.modules[self.name].__file__), *split_name[:j+1])):
                            parent_mod.__package__ = self.name + "." + ".".join(split_name[:j+1])
                        else:
                            parent_mod.__package__ = self.name + "." + ".".join(split_name[:j])
                        self.known_modules[parent_mod.__name__] = parent_mod
                    if success:
                        print(f"Found the module! {parent_mod}")
                        mod = parent_mod
                        self.known_modules[name] = mod
                        break
                else:
                    print(f"{name} not in known_modules")
            if mod is None:
                raise ImportError("This loader does not know module " + fullname)

            return mod

    def load_module(self, fullname):
        print("Loading module", fullname, file=output)
        try:
            # in case of a reload
            return sys.modules[fullname]
        except KeyError:
            pass
        mod = self.__get_module(fullname)
        mod.__loader__ = self
        sys.modules[fullname] = mod
        return mod

    def is_package(self, fullname):
        """
        Return true, if the named module is a package.

        We need this method to get correct spec objects with
        Python 3.4 (see PEP451)
        """
        print("Is package", fullname, file=output)
        return hasattr(self.__get_module(fullname), "__path__")

    def get_code(self, fullname):
        """Return None

        Required, if is_package is implemented"""
        print("Get code", fullname, file=output)
        self.__get_module(fullname)  # eventually raises ImportError
        return None
    get_source = get_code  # same as get_code

    def create_module(self, spec):
        print("Create module", spec, file=output)
        return self.load_module(spec.name)

    def exec_module(self, module):
        print("Exec module", module, file=output)
        pass

_importer = _SIRENResourcesMetaPathImporter(__name__)

_importer._add_module(fluxes, "fluxes")
_importer._add_module(detectors, "detectors")
_importer._add_module(processes, "processes")

print(__name__)

# Complete the moves implementation.
# This code is at the end of this module to speed up module loading.
# Turn this module into a package.
__path__ = []  # required for PEP 302 and PEP 451
__package__ = __name__  # see PEP 366 @ReservedAssignment
if globals().get("__spec__") is not None:
    __spec__.submodule_search_locations = []  # PEP 451 @UndefinedVariable
# Remove other six meta path importers, since they cause problems. This can
# happen if six is removed from sys.modules and then reloaded. (Setuptools does
# this for some reason.)
if sys.meta_path:
    for i, importer in enumerate(sys.meta_path):
        # Here's some real nastiness: Another "instance" of the six module might
        # be floating around. Therefore, we can't use isinstance() to check for
        # the six meta path importer, since the other six instance will have
        # inserted an importer with different class.
        if (type(importer).__name__ == "_SIRENResourcesMetaPathImporter" and
                importer.name == __name__):
            del sys.meta_path[i]
            break
    del i, importer
# Finally, add the importer to the meta path import hook.
sys.meta_path.append(_importer)

del _util
del ResourceList
