"""Build a staged wheel only when payload or install-rule contents change."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


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
    state = args.build / ".wheel_inputs.sha256"
    stamp = args.build / ".build_wheel"

    # CMake rewrites install scripts on every generate, even without a change.
    # Compare their contents, including subdirectory install rules. Payload
    # mtimes remain build-system inputs; additions/removals are in the manifest.
    inputs = []
    for name in args.inputs.read_text().splitlines():
        path = Path(name)
        stat = path.stat()
        inputs.append((name, stat.st_size, stat.st_mtime_ns))
    rules = [(str(path), hashlib.sha256(path.read_bytes()).hexdigest())
             for path in sorted(args.build.rglob("cmake_install.cmake"))
             if staging not in path.parents]
    signature = hashlib.sha256(json.dumps([
        vars(args), inputs, rules, sys.executable,
        {key: os.environ.get(key) for key in ("ARCHFLAGS", "MACOSX_DEPLOYMENT_TARGET")},
    ], default=str, sort_keys=True).encode()).hexdigest()
    if (state.is_file() and state.read_text() == signature
            and len(list(wheels.glob("*.whl"))) == 1):
        print("Wheel payload and install rules unchanged")
        stamp.touch()
        return

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
    state.write_text(signature)
    stamp.touch()


if __name__ == "__main__":
    main()
