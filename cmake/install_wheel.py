"""Install a CMake wheel into its requested prefix, honoring DESTDIR."""

import argparse
import os
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("wheel", type=Path)
    args = parser.parse_args()
    prefix = args.prefix.absolute()
    root = os.environ.get("DESTDIR")
    destination = Path(root) / prefix.relative_to(prefix.anchor) if root else prefix
    # A normal reinstall must remove obsolete modules and old dist-info. Pip
    # discovers packages in this interpreter, however, so it must not uninstall
    # them when writing to a separate prefix or a DESTDIR staging tree.
    replace = destination.resolve() == Path(sys.prefix).resolve()
    command = [sys.executable, "-m", "pip", "install", "--no-deps",
               "--force-reinstall" if replace else "--ignore-installed",
               "--prefix", str(prefix)]
    if root:
        command.extend(["--root", root])
    subprocess.run(command + [str(args.wheel)], check=True)


if __name__ == "__main__":
    main()
