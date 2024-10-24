__all__ = ["load_flux", "load_detector", "load_processes", "fluxes", "detectors", "processes"]

from . import _util

load_flux = _util.load_flux
load_detector = _util.load_detector
load_processes = _util.load_processes

class ResourceList:
    def __init__(self, resource_type, list_method, load_method):
        self.__resource_type = resource_type
        self.__list_method = list_method
        self.__load_method = load_method

    def __getattr__(self, name):
        resources = self.__list_method()
        if name in resources:
            return self.__load_method(name)
        else:
            # Default behaviour
            return object.__getattribute__(self, name)

    def __dir__(self):
        dirs = dir(self.__class__)
        dirs += list(self.__dict__.keys())
        dirs += self.__list_method()
        return sorted(dirs)

fluxes = ResourceList('fluxes', _util.list_fluxes, _util._get_flux_loader)
detectors = ResourceList('detectors', _util.list_detectors, _util._get_detector_loader)
processes = ResourceList('processes', _util.list_processes, _util._get_process_loader)

del _util
del ResourceList
