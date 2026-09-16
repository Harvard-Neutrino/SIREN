"""Put the wheel sibling directory first in an installed Mach-O's RPATHs."""

from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


def prefer_siblings(path):
    architectures = subprocess.check_output(
        ["/usr/bin/lipo", "-archs", str(path)], text=True).split()
    if not architectures:
        raise RuntimeError(f"No Mach-O architectures found: {path}")
    # install_name_tool applies an edit to every slice, but a universal binary
    # can have different RPATHs (including none) in each architecture. Edit thin
    # copies and only replace the staged binary after all edits/signing succeed.
    with tempfile.TemporaryDirectory(prefix=".wheel-rpath-", dir=path.parent) as folder:
        folder = Path(folder)
        slices = []
        changed = False
        for architecture in architectures:
            thin = folder / architecture
            if len(architectures) > 1:
                subprocess.run(["/usr/bin/lipo", str(path), "-thin", architecture,
                                "-output", str(thin)], check=True)
            else:
                shutil.copy2(path, thin)
            slices.append(thin)
            output = subprocess.check_output(["/usr/bin/otool", "-l", str(thin)], text=True)
            paths = re.findall(
                r"cmd LC_RPATH\n\s+cmdsize \d+\n\s+path (.*?) \(offset \d+\)", output)
            desired = list(dict.fromkeys(["@loader_path", *paths]))
            if paths == desired:
                continue
            changed = True
            # Deleting then adding in separate invocations is required when a
            # path appears in both lists. An empty list has nothing to delete.
            if paths:
                command = ["/usr/bin/install_name_tool"]
                for entry in dict.fromkeys(paths):
                    command.extend(["-delete_rpath", entry])
                subprocess.run(command + [str(thin)], check=True)
            command = ["/usr/bin/install_name_tool"]
            for entry in desired:
                command.extend(["-add_rpath", entry])
            subprocess.run(command + [str(thin)], check=True)
        if not changed:
            return
        result = slices[0]
        if len(slices) > 1:
            result = folder / "universal"
            subprocess.run(["/usr/bin/lipo", "-create", *map(str, slices),
                            "-output", str(result)], check=True)
        subprocess.run(["/usr/bin/codesign", "--force", "--sign", "-", str(result)], check=True)
        shutil.copystat(path, result)
        result.replace(path)


if __name__ == "__main__":
    prefer_siblings(Path(sys.argv[1]))
