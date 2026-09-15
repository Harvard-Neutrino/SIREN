"""Build a staged wheel only when payload or install-rule contents change."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


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
                if before.get(name) != after.get(name):
                    yield f"{group}: {name} ({before.get(name)} -> {after.get(name)})"
        elif before != after:
            yield f"{group}: {before} -> {after}"


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
        "payload": {name: file_digest(Path(name))
                    for name in args.inputs.read_text().splitlines()},
        "install rules": {str(path): file_digest(path)
                          for path in sorted(args.build.rglob("cmake_install.cmake"))
                          if staging not in path.parents},
        "interpreter": sys.executable,
        "environment": {key: os.environ.get(key)
                        for key in ("ARCHFLAGS", "MACOSX_DEPLOYMENT_TARGET")},
    }
    previous = json.loads(state.read_text()) if state.is_file() else {}
    wheel_files = list(wheels.glob("*.whl"))
    if previous == current and len(wheel_files) == 1:
        print("Wheel payload and install rules unchanged")
        return
    for change in changed_inputs(previous, current):
        print("Wheel rebuild: " + change, flush=True)
    if len(wheel_files) != 1:
        print(f"Wheel rebuild: expected one output wheel, found {len(wheel_files)}", flush=True)

    for path in (staging, wheels):
        if path.exists():
            shutil.rmtree(path)
        path.mkdir()
    subprocess.run([args.cmake, "--install", str(args.build), "--config", args.config,
                    "--prefix", str(staging), "--component", "PythonWheel"], check=True)
    for name in ("pyproject.toml", "README.md", "LICENSE"):
        shutil.copy2(args.source / name, staging / name)
    shutil.copy2(args.source / "package/CMakeLists.txt", staging / "CMakeLists.txt")
    isolation = ["--no-build-isolation"] if args.no_build_isolation else []
    subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-deps", *isolation,
                    "--config-settings=cmake.define.SIREN_WHEEL_LIBRARY_DIR=" + args.library_dir,
                    "--wheel-dir", str(wheels), str(staging)], check=True)
    state.write_text(json.dumps(current, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
