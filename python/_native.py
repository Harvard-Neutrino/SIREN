"""Native runtime identity checks, without importing extension modules."""

import ctypes
from pathlib import Path
import re
import sys


def linux_mapped_paths(maps):
    paths = []
    for line in maps.splitlines():
        fields = line.split(None, 5)
        if len(fields) == 6 and fields[5].startswith("/"):
            paths.append(Path(fields[5].replace(r"\012", "\n")))
    return paths


def loaded_libraries():
    if sys.platform == "darwin":
        loader = ctypes.CDLL(None)
        loader._dyld_image_count.restype = ctypes.c_uint32
        loader._dyld_get_image_name.argtypes = [ctypes.c_uint32]
        loader._dyld_get_image_name.restype = ctypes.c_char_p
        return [Path(loader._dyld_get_image_name(i).decode())
                for i in range(loader._dyld_image_count())]
    if sys.platform.startswith("linux"):
        return linux_mapped_paths(Path("/proc/self/maps").read_text())
    if sys.platform == "win32":
        loader = ctypes.WinDLL("kernel32", use_last_error=True)
        loader.GetCurrentProcess.restype = ctypes.c_void_p
        loader.K32EnumProcessModules.argtypes = [ctypes.c_void_p, ctypes.c_void_p,
                                               ctypes.c_uint32, ctypes.c_void_p]
        loader.K32EnumProcessModules.restype = ctypes.c_int
        loader.GetModuleFileNameW.argtypes = [ctypes.c_void_p, ctypes.c_wchar_p,
                                            ctypes.c_uint32]
        loader.GetModuleFileNameW.restype = ctypes.c_uint32
        count = 256
        while True:
            modules = (ctypes.c_void_p * count)()
            needed = ctypes.c_uint32()
            if not loader.K32EnumProcessModules(loader.GetCurrentProcess(), modules,
                                               ctypes.sizeof(modules), ctypes.byref(needed)):
                raise ctypes.WinError(ctypes.get_last_error())
            if needed.value <= ctypes.sizeof(modules):
                break
            count = needed.value // ctypes.sizeof(ctypes.c_void_p)
        paths = []
        for handle in modules[:needed.value // ctypes.sizeof(ctypes.c_void_p)]:
            path = ctypes.create_unicode_buffer(32768)
            if not loader.GetModuleFileNameW(handle, path, len(path)):
                raise ctypes.WinError(ctypes.get_last_error())
            paths.append(Path(path.value))
        return paths
    raise RuntimeError("Native runtime checks are unavailable on " + sys.platform)


def reject_standalone_runtime():
    """A native core and a wheel core must not mint independent particle IDs."""
    standalone = re.compile(r"^(lib)?siren(?:\.[0-9]+)*\.(dylib|so(?:\.[0-9]+)*|dll)$", re.I)
    conflicts = [str(path) for path in loaded_libraries()
                 if standalone.fullmatch(path.name[:-10] if path.name.endswith(" (deleted)") else path.name)]
    if conflicts:
        raise ImportError(
            "The standalone SIREN library is already loaded: " + ", ".join(conflicts)
            + ". Native SIREN and the Python wheel have separate process state; "
            "use separate processes instead of loading both cores."
        )
