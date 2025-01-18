#!/usr/bin/env python3

"""
Parse top-level pyproject.toml to extract [project] metadata,
then generate a minimal pyproject.toml using setuptools for the staging area.
"""

import sys
import os

# If you're on Python < 3.11, you can 'pip install tomli' and do `import tomli`.
# On Python >= 3.11, tomllib is in the standard library.
try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # fallback for older Python

try:
    import tomli_w as toml  # Use tomli_w to *write* TOML back out
except ModuleNotFoundError as e:
    print(f"[Error] Please install the 'tomli-w' package: `pip install toml-w`", file=sys.stderr)
    raise

def main():
    if len(sys.argv) < 3:
        print("Usage: parse_pyproject.py <input_pyproject> <output_pyproject>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    with open(input_file, "rb") as f:
        data = tomllib.load(f)

    # The main [project] metadata we want to reuse
    project_data = data.get("project", {})

    # Build a new structure for the staging pyproject
    # We'll always want a simple setuptools build-backend for "prebuilt .so" packaging.
    new_toml = {}
    new_toml["build-system"] = {
        "requires": ["setuptools>=42", "wheel"],
        "build-backend": "setuptools.build_meta"
    }

    # Copy project metadata verbatim (name, version, dependencies, etc.)
    # Possibly override or add fields if you want to differ from the top-level.
    new_toml["project"] = project_data

    # Optionally add [tool.setuptools] to handle package_data.
    # e.g., "include-package-data = true"
    new_toml["tool"] = {
        "setuptools": {
            "include-package-data": True,
            "package-data": {
                "siren": ["*.so", "*.dylib"],
            },
        }
    }

    # Write the merged TOML to output_file
    with open(output_file, "wb") as out:
        toml.dump(new_toml, out)

if __name__ == "__main__":
    main()
