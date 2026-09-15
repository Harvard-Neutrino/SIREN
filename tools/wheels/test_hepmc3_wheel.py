"""Check native wheel tags, the loaded core library, and HepMC3 round trips.

Run against an installed wheel from outside the source tree. Release validation
also makes the original source/build directories unavailable and clears loader
search-path overrides. Missing native support is a failure, never a skip.
Pass --standalone-library to additionally check the wheel with a native SIREN
installation on the loader search path, in a fresh interpreter.
"""

import argparse
from email.parser import Parser
import hashlib
from importlib import metadata
import os
import re
from pathlib import Path
import subprocess
import sys
import tempfile

from packaging.tags import parse_tag, sys_tags


def is_core_library(path):
    name = path.name.lower()
    if name.endswith(" (deleted)"):
        name = name[:-10]
    return (
        name.startswith("libsiren") and (name.endswith(".dylib") or ".so" in name)
    ) or (name.startswith("siren") and name.endswith(".dll"))


def library_identity(path):
    """Match version aliases without conflating distinct repair-tool hashes."""
    name = path.name.lower()
    if name.endswith(" (deleted)"):
        name = name[:-10]
    if name.endswith(".dylib"):
        return re.sub(r"(\.[0-9]+)(?:\.[0-9]+)*(?=\.dylib$)", r"\1", name)
    return re.sub(r"(?<=\.so)(\.[0-9]+)(?:\.[0-9]+)*$", r"\1", name)


def verify_bundled_libraries(wheel, libraries):
    packaged = {library_identity(path) for path in wheel.files
                if path.name.endswith((".dylib", ".dll")) or (
                    ".so" in path.name and ".cpython-" not in path.name
                    and ".abi3." not in path.name)}
    installed = {wheel.locate_file(path).resolve() for path in wheel.files}
    libraries = {path.resolve() for path in libraries}
    foreign = {path.resolve() for path in libraries
               if library_identity(path) in packaged and path.resolve() not in installed}
    if foreign:
        # Other wheels (e.g. SciPy on macOS) can load independent copies with
        # the same SONAME. Accept those only if SIREN's copy is also loaded.
        # An actual replacement, an unowned native-prefix copy, or a stale
        # SIREN installation still fails this check.
        own_loaded = {library_identity(path) for path in libraries & installed}
        other_files = set()
        for distribution in metadata.distributions():
            if distribution.metadata["Name"].lower() == "siren":
                continue
            other_files.update(distribution.locate_file(path).resolve()
                               for path in (distribution.files or ())
                               if library_identity(path) in own_loaded)
        foreign = {path for path in foreign if path not in other_files}
    assert not foreign, f"Bundled libraries loaded from outside this wheel: {foreign}"


def check_loader_override(variable, expected):
    assert os.environ.get(variable) == expected, (
        f"{variable} did not reach the child interpreter; "
        "this interpreter cannot validate the requested loader override"
    )


def main():
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

    packaged_core = [path for path in wheel.files if is_core_library(path)]
    assert len(packaged_core) == 1, (
        f"Expected one core library in the wheel: {packaged_core}"
    )

    import siren
    from siren import hepmc3
    from siren import dataclasses as d

    from siren._native import loaded_libraries
    libraries = loaded_libraries()

    core = {path.resolve() for path in libraries if is_core_library(path)}
    assert len(core) == 1, core
    installed_files = {wheel.locate_file(path).resolve() for path in wheel.files}
    assert core <= installed_files, f"Core library is not from this wheel: {core}"
    verify_bundled_libraries(wheel, libraries)
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

    print("Native wheel tags, bundled library origins, and HepMC3 round trips OK")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--standalone-library", type=Path)
    parser.add_argument("--clean-environment", action="store_true",
                        help="Run in a child without build/repair loader search paths")
    parser.add_argument("--expect-loader-override", nargs=2, metavar=("VARIABLE", "VALUE"),
                        help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.clean_environment:
        if args.expect_loader_override is not None:
            parser.error("Cannot clear an override that this child is meant to check")
        env = {key: value for key, value in os.environ.items() if key not in (
            "PYTHONPATH", "DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH", "LD_LIBRARY_PATH")}
        command = [sys.executable, str(Path(__file__).resolve())]
        if args.standalone_library is not None:
            command.extend(["--standalone-library", str(args.standalone_library)])
        sys.exit(subprocess.run(command, env=env).returncode)
    if args.expect_loader_override is not None:
        check_loader_override(*args.expect_loader_override)
    main()
    if args.standalone_library is not None:
        library = args.standalone_library.resolve()
        if not library.is_file():
            parser.error(f"Standalone library does not exist: {library}")
        if sys.platform not in ("darwin", "linux"):
            parser.error("The standalone loader override check requires macOS or Linux")
        variable = (
            "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
        )
        env = dict(os.environ)
        env[variable] = os.pathsep.join(
            [str(library.parent)] + ([env[variable]] if env.get(variable) else [])
        )
        print(f"Checking with {variable}={env[variable]}", flush=True)
        subprocess.run(
            [sys.executable, str(Path(__file__).resolve()),
             "--expect-loader-override", variable, env[variable]], env=env, check=True
        )
