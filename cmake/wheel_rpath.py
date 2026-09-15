"""Put the wheel sibling directory first in an installed Mach-O's RPATHs."""

from pathlib import Path
import re
import subprocess
import sys


def prefer_siblings(path):
    output = subprocess.check_output(["/usr/bin/otool", "-l", str(path)], text=True)
    # Universal binaries repeat the load commands for each architecture.
    per_arch = []
    for section in re.split(r"^\S.*:\n", output, flags=re.MULTILINE):
        paths = re.findall(r"cmd LC_RPATH\n\s+cmdsize \d+\n\s+path (.*?) \(offset \d+\)", section)
        if paths:
            per_arch.append(paths)
    if not per_arch or any(paths != per_arch[0] for paths in per_arch):
        raise RuntimeError(f"Expected consistent RPATHs across architectures: {path}: {per_arch}")
    paths = per_arch[0]
    desired = ["@loader_path"] + [entry for entry in paths if entry != "@loader_path"]
    if paths == desired:
        return
    # CMake can retain link-derived RPATHs before INSTALL_RPATH during its
    # install-time rewrite. Reorder only the installed PythonWheel copy.
    command = ["/usr/bin/install_name_tool"]
    for entry in paths:
        command.extend(["-delete_rpath", entry])
    subprocess.run(command + [str(path)], check=True)
    command = ["/usr/bin/install_name_tool"]
    for entry in desired:
        command.extend(["-add_rpath", entry])
    subprocess.run(command + [str(path)], check=True)
    subprocess.run(["/usr/bin/codesign", "--force", "--sign", "-", str(path)], check=True)


if __name__ == "__main__":
    prefer_siblings(Path(sys.argv[1]))
