import sysconfig
import platform

SUPPORTED_PLATFORMS = [
    'linux-aarch64',
    'linux-x86_64',
    'musllinux-x86_64',
    'linux-i686',
    'linux-ppc64le',
    'linux-s390x',
    'win-amd64',
    'win-32',
    'macosx-x86_64',
    'macosx-arm64',
]

def get_plat():
    plat = sysconfig.get_platform()
    plat_split = plat.split("-")
    arch = plat_split[-1]
    if arch == "win32":
        plat = "win-32"
    elif arch in ["universal2", "intel"]:
        plat = "macosx-" + str(platform.uname().machine)
    elif len(plat_split) > 2:
        plat = str(plat_split[0]) + "-" + str(arch)
    assert plat in SUPPORTED_PLATFORMS, ('invalid platform ' + str(plat))
    return plat
