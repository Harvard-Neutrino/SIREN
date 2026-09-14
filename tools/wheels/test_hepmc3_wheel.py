"""Check native wheel tags, the loaded core library, and HepMC3 round trips.

Run against an installed wheel from outside the source tree. Release validation
also makes the original source/build directories unavailable and clears loader
search-path overrides. Missing native support is a failure, never a skip.
"""

import ctypes
from email.parser import Parser
import hashlib
from importlib import metadata
from pathlib import Path
import sys
import tempfile

from packaging.tags import parse_tag, sys_tags

wheel = metadata.distribution("siren")
headers = Parser().parsestr(wheel.read_text("WHEEL"))
assert headers["Root-Is-Purelib"] == "false", headers
wheel_tags = set().union(*(parse_tag(value) for value in headers.get_all("Tag")))
assert wheel_tags.intersection(sys_tags()), wheel_tags
assert all(tag.abi != "none" and tag.platform != "any" for tag in wheel_tags), (
    wheel_tags
)

# The interpreter supplies Python; the wheel must not bundle a second runtime.
assert not any(
    path.name == "Python"
    or path.name.lower().startswith("libpython")
    or (path.name.lower().startswith("python") and path.suffix == ".dll")
    for path in wheel.files
), "Wheel bundles an embedding runtime"

import siren
from siren import hepmc3
from siren import dataclasses as d

if sys.platform == "darwin":
    loader = ctypes.CDLL(None)
    loader._dyld_image_count.restype = ctypes.c_uint32
    loader._dyld_get_image_name.argtypes = [ctypes.c_uint32]
    loader._dyld_get_image_name.restype = ctypes.c_char_p
    libraries = [
        Path(loader._dyld_get_image_name(i).decode())
        for i in range(loader._dyld_image_count())
    ]
elif sys.platform.startswith("linux"):
    libraries = [
        Path(line.split()[-1])
        for line in Path("/proc/self/maps").read_text().splitlines()
        if "/" in line
    ]
elif sys.platform == "win32":
    loader = ctypes.WinDLL("kernel32", use_last_error=True)
    loader.GetModuleHandleW.argtypes = [ctypes.c_wchar_p]
    loader.GetModuleHandleW.restype = ctypes.c_void_p
    loader.GetModuleFileNameW.argtypes = [
        ctypes.c_void_p,
        ctypes.c_wchar_p,
        ctypes.c_uint32,
    ]
    loader.GetModuleFileNameW.restype = ctypes.c_uint32
    libraries = []
    for name in {
        path.name
        for path in wheel.files
        if path.name.lower().startswith(("siren", "libsiren"))
        and path.suffix.lower() == ".dll"
    }:
        handle = loader.GetModuleHandleW(name)
        if handle:
            path = ctypes.create_unicode_buffer(32768)
            assert loader.GetModuleFileNameW(handle, path, len(path)), (
                ctypes.get_last_error()
            )
            libraries.append(Path(path.value))
else:
    raise AssertionError(
        f"Loaded-library verification is not implemented on {sys.platform}"
    )

core = {
    path.resolve()
    for path in libraries
    if path.name.lower().startswith("libsiren")
    or (path.name.lower().startswith("siren") and path.suffix.lower() == ".dll")
}
assert len(core) == 1, core
installed_files = {wheel.locate_file(path).resolve() for path in wheel.files}
assert core <= installed_files, f"Core library is not from this wheel: {core}"
for path in core:
    print(f"Loaded {path}: sha256={hashlib.sha256(path.read_bytes()).hexdigest()}")

rng = siren.utilities.SIREN_random(7)
assert 0 <= rng.Uniform(0, 1) <= 1
rec = d.InteractionRecord()
sig = rec.signature
sig.primary_type = d.ParticleType.NuMu
sig.secondary_types = [d.ParticleType.NuMu]
rec.signature = sig
rec.primary_momentum = [1.0, 0.0, 0.0, 1.0]
rec.secondary_ids = [d.ParticleID(1, 1)]
rec.secondary_momenta = [[1.0, 0.0, 0.0, 1.0]]
rec.secondary_masses = [0.0]
rec.secondary_helicities = [0.0]
tree = d.InteractionTree()
tree.add_entry(rec, None)

with tempfile.TemporaryDirectory() as directory:
    for suffix in (".hepmc3", ".hepmc3.gz"):
        path = str(Path(directory) / ("wheel_smoke" + suffix))
        options = hepmc3.HepMC3WriterOptions()
        options.gzip = suffix.endswith(".gz")
        hepmc3.SaveInteractionTreesAsHepMC3([tree], path, options)
        if options.gzip:
            assert Path(path).read_bytes()[:2] == b"\x1f\x8b"
        loaded = hepmc3.LoadInteractionTreesFromHepMC3(path)
        assert len(loaded) == 1, loaded
        assert len(loaded[0].tree) == 1, loaded[0].tree
        assert list(loaded[0].tree[0].record.primary_momentum) == list(
            rec.primary_momentum
        )

print("Native wheel tags, packaged core library, and HepMC3 round trips OK")
