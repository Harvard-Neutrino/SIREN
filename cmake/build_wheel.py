"""Build a staged wheel only when payload or install-rule contents change."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def file_digest(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
    except (ValueError, UnicodeError) as error:
        print(f"Wheel rebuild: ignoring invalid cache {path}: {error}", flush=True)
        return {}


def payload_digest(path, source):
    try:
        return file_digest(path)
    except FileNotFoundError:
        # A stale glob manifest can still list a removed Python/resource file.
        # Directory installation will omit it; a missing binary or build input
        # must remain a failure, never yield an incomplete but cached wheel.
        if any(source / name in path.parents for name in ("python", "resources")):
            return None
        raise FileNotFoundError(f"Required wheel input missing: {path}; rebuild its CMake target")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--library-dir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--cmake", required=True)
    parser.add_argument("--no-build-isolation", action="store_true")
    args = parser.parse_args()
    staging = args.build / "python_staging"
    wheels = args.build / "dist_wheels"
    state = args.build / ".wheel_inputs.json"

    # CMake rewrites install scripts on every generate, even without a change.
    # Compare contents, including subdirectory install rules and payload files.
    # A content-preserving relink should not rebuild a wheel, while a restore
    # with unchanged size/mtime must not leave stale Python/resources in it.
    current = {
        "arguments": {key: str(value) for key, value in vars(args).items()},
        "payload": {name: payload_digest(Path(name), args.source)
                    for name in args.inputs.read_text().splitlines()},
        "install rules": {str(path): file_digest(path)
                          for path in sorted(args.build.rglob("cmake_install.cmake"))
                          if staging not in path.parents},
        "interpreter": sys.executable,
        "environment": {key: os.environ.get(key)
                        for key in ("ARCHFLAGS", "MACOSX_DEPLOYMENT_TARGET")},
    }
    previous = read_state(state)
    wheel_files = list(wheels.glob("*.whl"))
    if previous == current and len(wheel_files) == 1:
        print("Wheel payload and install rules unchanged")
        return
    for change in changed_inputs(previous, current):
        print("Wheel rebuild: " + change, flush=True)
    if len(wheel_files) != 1:
        print(f"Wheel rebuild: expected one output wheel, found {len(wheel_files)}", flush=True)

    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir()
    subprocess.run([args.cmake, "--install", str(args.build), "--config", args.config,
                    "--prefix", str(staging), "--component", "PythonWheel"], check=True)
    for name in ("pyproject.toml", "README.md", "LICENSE"):
        shutil.copy2(args.source / name, staging / name)
    shutil.copy2(args.source / "package/CMakeLists.txt", staging / "CMakeLists.txt")
    isolation = ["--no-build-isolation"] if args.no_build_isolation else []
    # Keep the last successful wheel and state if staging or pip fails. Publish
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
        new_state = pending / state.name
        new_state.write_text(json.dumps(current, indent=2, sort_keys=True) + "\n")
        new_state.replace(state)


if __name__ == "__main__":
    main()
