"""Install a CMake wheel into a verified Python destination, honoring DESTDIR."""

import argparse
import csv
from email.parser import Parser
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile


def same_directory(left, right):
    try:
        return os.path.samefile(left, right)
    except (OSError, TypeError):
        return False


def installed_distributions(paths):
    from importlib import metadata

    distributions = {}
    for dist in metadata.distributions(path=list(dict.fromkeys(paths))):
        if (dist.metadata.get("Name") or "").lower() != "siren":
            continue
        record = dist.read_text("RECORD")
        if record is None:
            raise RuntimeError("Cannot replace SIREN without its installed file record")
        # purelib and platlib can be aliases of the same directory.
        # Resolve directory aliases but preserve a symlink at the recorded
        # filename: replacement must remove the link, not its external target.
        # metadata.files can omit missing paths, including dangling symlinks.
        # Read RECORD directly to retain every path pip was meant to remove.
        paths = (Path(dist.locate_file(row[0])) for row in csv.reader(io.StringIO(record)) if row)
        files = sorted(str(path.parent.resolve() / path.name) for path in paths)
        distributions[tuple(files)] = {"version": dist.version, "files": files}
    return list(distributions.values())


def inspect_layout(prefix, root):
    # Query the installing pip, including distro-specific schemes, instead of
    # reconstructing its paths from sys.prefix or a hardcoded site-packages.
    # An unavailable API fails this read-only probe before any installation.
    from pip._internal.locations import get_scheme
    from importlib import metadata

    active = get_scheme("siren", isolated=True)
    target = get_scheme("siren", prefix=str(prefix), root=root, isolated=True)
    active_paths = [active.purelib, active.platlib]
    target_paths = [target.purelib, target.platlib]
    try:
        selected = str(metadata.distribution("siren").locate_file(""))
    except metadata.PackageNotFoundError:
        selected = None
    return {"prefix": sys.prefix, "version": list(sys.version_info[:2]),
            "implementation": sys.implementation.name,
            "active": active_paths, "target": target_paths, "selected": selected,
            "destination": installed_distributions(target_paths)}


def isolated_environment():
    # -I excludes PYTHONPATH/user-site; --isolated ignores pip environment and
    # user settings. Disable config files too so they cannot redirect the write.
    return dict(os.environ, PIP_CONFIG_FILE=os.devnull)


def probe(python, prefix, root):
    command = [str(python), "-I", str(Path(__file__).resolve()), "--inspect",
               "--prefix", str(prefix)]
    if root:
        command += ["--root", root]
    return json.loads(subprocess.check_output(command, env=isolated_environment(), text=True))


def replaces_in_place(layout):
    return all(same_directory(target, active)
               for target, active in zip(layout["target"], layout["active"]))


def select_interpreter(prefix, root):
    python = sys.executable
    layout = probe(python, prefix, root)
    if root:
        # Even a spelling that maps back onto the live environment is staging,
        # never permission to uninstall or overwrite that environment.
        if any(same_directory(target, active)
               for target in layout["target"] for active in layout["active"]):
            raise RuntimeError("DESTDIR staging overlaps the live Python installation")
        return python, layout, False
    if replaces_in_place(layout):
        return python, layout, True
    if same_directory(prefix, layout["prefix"]):
        raise RuntimeError("pip --prefix uses different package directories from this interpreter; "
                           "configure Python_EXECUTABLE and CMAKE_INSTALL_PREFIX for a matching venv")
    for name in ("python", "python3"):
        candidate = prefix / "bin" / name
        if not candidate.is_file():
            continue
        target = probe(candidate, prefix, None)
        if (target["version"] != layout["version"]
                or target["implementation"] != layout["implementation"]):
            raise RuntimeError("Destination Python version differs from the build interpreter; "
                               "rebuild with its Python_EXECUTABLE")
        if not replaces_in_place(target):
            raise RuntimeError("Destination interpreter does not install into the requested package directories")
        return str(candidate), target, True
    if (prefix / "pyvenv.cfg").exists() or (prefix / "conda-meta").exists():
        raise RuntimeError("Destination Python environment has no usable interpreter")
    return python, layout, False


def install(prefix, wheel, root):
    python, layout, replace = select_interpreter(prefix, root)
    existing = layout["destination"]
    if not replace and existing:
        raise RuntimeError("Destination already contains SIREN; use its matching Python interpreter "
                           "or a fresh staging prefix")
    if replace:
        if len(existing) > 1:
            raise RuntimeError("Destination contains multiple SIREN installations; resolve them before replacing")
        if layout["selected"] is not None and not any(
                same_directory(layout["selected"], path) for path in layout["active"]):
            raise RuntimeError("Another SIREN installation shadows the destination; refusing an out-of-environment uninstall")
    with zipfile.ZipFile(wheel) as archive:
        metadata = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
        if len(metadata) != 1:
            raise RuntimeError("Expected one wheel METADATA file")
        headers = Parser().parsestr(archive.read(metadata[0]).decode())
        if (headers.get("Name") or "").lower() != "siren" or not headers.get("Version"):
            raise RuntimeError("Expected a SIREN wheel with version metadata")
        expected = headers["Version"]
    command = [str(python), "-I", "-m", "pip", "--isolated", "--disable-pip-version-check",
               "install", "--no-index", "--no-deps",
               "--force-reinstall" if replace else "--ignore-installed",
               "--prefix", str(prefix)]
    if root:
        command.extend(["--root", root])
    subprocess.run(command + [str(wheel)], env=isolated_environment(), check=True)
    after = probe(python, prefix, root)["destination"]
    if len(after) != 1 or after[0]["version"] != expected:
        raise RuntimeError("pip did not leave exactly the requested SIREN installation")
    old_files = {path for dist in existing for path in dist["files"]}
    obsolete = old_files - set(after[0]["files"])
    if any(Path(path).exists() or Path(path).is_symlink() for path in obsolete):
        raise RuntimeError("pip left obsolete SIREN files after replacement")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--inspect", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--root", help=argparse.SUPPRESS)
    parser.add_argument("wheel", type=Path, nargs="?")
    args = parser.parse_args()
    if args.inspect:
        print(json.dumps(inspect_layout(args.prefix, args.root)))
    else:
        if args.wheel is None:
            parser.error("wheel is required")
        install(args.prefix.absolute(), args.wheel, os.environ.get("DESTDIR"))


if __name__ == "__main__":
    main()
