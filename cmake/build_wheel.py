"""Build the CMake wheel from a fresh PythonWheel staging tree when its contents change."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


# Packaging metadata copied next to the staged package for the wheel backend.
PACKAGING_FILES = {
    "pyproject.toml": "pyproject.toml",
    "README.md": "README.md",
    "LICENSE": "LICENSE",
    "package/CMakeLists.txt": "CMakeLists.txt",
}


def file_digest(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def staged_digests(staging):
    """Content of every staged entry, keyed by its path inside the wheel.

    The backend packages what it reads, so links are followed: a file link
    contributes its target's bytes and a directory link its target's entries.
    A dangling link or a directory cycle is recorded as such instead.
    """
    digests = {}

    def walk(directory, ancestors):
        for path in sorted(directory.iterdir()):
            relative = path.relative_to(staging).as_posix()
            if path.is_dir():
                target = path.resolve()
                if target in ancestors:
                    digests[relative] = "cycle:" + os.readlink(str(path))
                else:
                    walk(path, ancestors | {target})
            elif path.is_file():
                digests[relative] = file_digest(path)
            else:
                digests[relative] = "dangling:" + os.readlink(str(path))

    walk(staging, {staging.resolve()})
    return digests


def changed_inputs(previous, current):
    for group in sorted(current.keys() | previous.keys()):
        before, after = previous.get(group, {}), current.get(group, {})
        if isinstance(before, dict) and isinstance(after, dict):
            for name in sorted(before.keys() | after.keys()):
                if name not in before or name not in after or before[name] != after[name]:
                    yield f"{group}: {name} ({before.get(name)} -> {after.get(name)})"
        elif before != after:
            yield f"{group}: {before} -> {after}"


def read_state(path):
    try:
        state = json.loads(path.read_text())
        if not isinstance(state, dict):
            raise ValueError("expected an object")
        return state
    except FileNotFoundError:
        return {}
    except (OSError, ValueError) as error:
        print(f"Wheel rebuild: ignoring invalid cache {path}: {error}", flush=True)
        return {}


def stage(args, staging):
    """Install the PythonWheel component into an empty staging tree.

    CMake's install rules are the only definition of the wheel contents: there
    is no separate source walk to keep consistent with them, and a file removed
    from the source tree is absent from a fresh tree.
    """
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir()
    command = [args.cmake, "--install", str(args.build), "--config", args.config,
               "--prefix", str(staging), "--component", "PythonWheel"]
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        print(result.stdout + result.stderr, flush=True)
        raise subprocess.CalledProcessError(result.returncode, command)
    for name, target in PACKAGING_FILES.items():
        shutil.copy2(args.source / name, staging / target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--library-dir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--cmake", required=True)
    parser.add_argument("--no-build-isolation", action="store_true")
    args = parser.parse_args()
    args.source = args.source.resolve()
    args.build = args.build.resolve()
    staging = args.build / "python_staging"
    wheels = args.build / "dist_wheels"
    state = args.build / ".wheel_inputs.json"

    stage(args, staging)
    # Compare the staged contents, not timestamps: CMake rewrites install
    # scripts on every generate and a restore can preserve size and mtime.
    current = {
        "arguments": {key: str(value) for key, value in vars(args).items()},
        "staged": staged_digests(staging),
        "interpreter": sys.executable,
        "environment": {key: os.environ.get(key)
                        for key in ("ARCHFLAGS", "MACOSX_DEPLOYMENT_TARGET")},
    }
    previous = read_state(state)
    wheel_files = list(wheels.glob("*.whl"))
    if previous == current and len(wheel_files) == 1:
        print("Staged wheel contents unchanged")
        return
    for change in changed_inputs(previous, current):
        print("Wheel rebuild: " + change, flush=True)
    if len(wheel_files) != 1:
        print(f"Wheel rebuild: expected one output wheel, found {len(wheel_files)}", flush=True)

    isolation = ["--no-build-isolation"] if args.no_build_isolation else []
    # Keep the last successful wheel and state if the backend fails. Publish
    # complete files with same-filesystem renames, including the cache itself.
    with tempfile.TemporaryDirectory(prefix=".wheel-build-", dir=args.build) as folder:
        pending = Path(folder)
        subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-deps", *isolation,
                        "--config-settings=cmake.define.SIREN_WHEEL_LIBRARY_DIR=" + args.library_dir,
                        "--wheel-dir", str(pending), str(staging)], check=True)
        outputs = list(pending.glob("*.whl"))
        if len(outputs) != 1:
            raise RuntimeError(f"Expected one built wheel, found {len(outputs)}")
        wheels.mkdir(exist_ok=True)
        outputs[0].replace(wheels / outputs[0].name)
        for obsolete in wheel_files:
            if obsolete.name != outputs[0].name:
                obsolete.unlink()
        try:
            new_state = pending / state.name
            new_state.write_text(json.dumps(current, indent=2, sort_keys=True) + "\n")
            new_state.replace(state)
        except OSError as error:
            # The wheel is complete. An unwritable cache (including a directory
            # at this filename) must not turn a successful package into failure.
            print(f"Wheel built, but input cache could not be saved: {error}; "
                  "the next invocation will rebuild", flush=True)


if __name__ == "__main__":
    main()
