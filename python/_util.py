import os
import sys

THIS_DIR = os.path.abspath(os.path.dirname(__file__))

# From pyzolib/paths.py (https://bitbucket.org/pyzo/pyzolib/src/tip/paths.py)
def appdata_dir(appname=None, roaming=False):
    """appdata_dir(appname=None, roaming=False)

    Get the path to the application directory, where applications are allowed
    to write user specific files (e.g. configurations). For non-user specific
    data, consider using common_appdata_dir().
    If appname is given, a subdir is appended (and created if necessary).
    If roaming is True, will prefer a roaming directory (Windows Vista/7).
    """

    # Define default user directory
    userDir = os.getenv("LEPTONINJECTOR_USERDIR", None)
    if userDir is None:
        userDir = os.path.expanduser("~")
        if not os.path.isdir(userDir):  # pragma: no cover
            userDir = "/var/tmp"  # issue #54

    # Get system app data dir
    path = None
    if sys.platform.startswith("win"):
        path1, path2 = os.getenv("LOCALAPPDATA"), os.getenv("APPDATA")
        path = (path2 or path1) if roaming else (path1 or path2)
    elif sys.platform.startswith("darwin"):
        path = os.path.join(userDir, "Library", "Application Support")
    # On Linux and as fallback
    if not (path and os.path.isdir(path)):
        path = userDir

    # Maybe we should store things local to the executable (in case of a
    # portable distro or a frozen application that wants to be portable)
    prefix = sys.prefix
    if getattr(sys, "frozen", None):
        prefix = os.path.abspath(os.path.dirname(sys.executable))
    for reldir in ("settings", "../settings"):
        localpath = os.path.abspath(os.path.join(prefix, reldir))
        if os.path.isdir(localpath):  # pragma: no cover
            try:
                open(os.path.join(localpath, "test.write"), "wb").close()
                os.remove(os.path.join(localpath, "test.write"))
            except IOError:
                pass  # We cannot write in this directory
            else:
                path = localpath
                break

    # Get path specific for this app
    if appname:
        if path == userDir:
            appname = "." + appname.lstrip(".")  # Make it a hidden directory
        path = os.path.join(path, appname)
        if not os.path.isdir(path):  # pragma: no cover
            os.makedirs(path, exist_ok=True)

    # Done
    return path

# from imageio
# https://github.com/imageio/imageio/blob/65d79140018bb7c64c0692ea72cb4093e8d632a0/imageio/core/util.py
def resource_dirs():
    """resource_dirs()

    Get a list of directories where leptoninjector resources may be located.
    The first directory in this list is the "resources" directory in
    the package itself. The second directory is the appdata directory
    (~/.leptoninjector on Linux). The list further contains the application
    directory (for frozen apps), and may include additional directories
    in the future.
    """
    dirs = [resource_package_dir()]
    # Resource dir baked in the package.
    # Appdata directory
    try:
        dirs.append(appdata_dir("leptoninjector"))
    except Exception:  # pragma: no cover
        pass  # The home dir may not be writable
    # Directory where the app is located (mainly for frozen apps)
    if getattr(sys, "frozen", None):
        dirs.append(os.path.abspath(os.path.dirname(sys.executable)))
    elif sys.path and sys.path[0]:
        dirs.append(os.path.abspath(sys.path[0]))
    return dirs


# from imageio
# https://github.com/imageio/imageio/blob/65d79140018bb7c64c0692ea72cb4093e8d632a0/imageio/core/util.py
def resource_package_dir():
    """package_dir

    Get the resources directory in the leptoninjector package installation
    directory.

    Notes
    -----
    This is a convenience method that is used by `resource_dirs` and
    leptoninjector entry point scripts.
    """
    # Make pkg_resources optional if setuptools is not available
    try:
        # Avoid importing pkg_resources in the top level due to how slow it is
        # https://github.com/pypa/setuptools/issues/510
        import pkg_resources
    except ImportError:
        pkg_resources = None

    if pkg_resources:
        # The directory returned by `pkg_resources.resource_filename`
        # also works with eggs.
        pdir = pkg_resources.resource_filename("leptoninjector", "resources")
    else:
        # If setuptools is not available, use fallback
        pdir = os.path.abspath(os.path.join(THIS_DIR, "resources"))
    return pdir


# from imageio
# https://github.com/imageio/imageio/blob/65d79140018bb7c64c0692ea72cb4093e8d632a0/imageio/core/util.py
def get_platform():
    """get_platform()

    Get a string that specifies the platform more specific than
    sys.platform does. The result can be: linux32, linux64, win32,
    win64, osx32, osx64. Other platforms may be added in the future.
    """
    # Get platform
    if sys.platform.startswith("linux"):
        plat = "linux%i"
    elif sys.platform.startswith("win"):
        plat = "win%i"
    elif sys.platform.startswith("darwin"):
        plat = "osx%i"
    elif sys.platform.startswith("freebsd"):
        plat = "freebsd%i"
    else:  # pragma: no cover
        return None

    return plat % (struct.calcsize("P") * 8)  # 32 or 64 bits


# from imageio
# https://github.com/imageio/imageio/blob/65d79140018bb7c64c0692ea72cb4093e8d632a0/imageio/core/util.py
def has_module(module_name):
    """Check to see if a python module is available."""
    if sys.version_info > (3, 4):
        import importlib

        name_parts = module_name.split(".")
        for i in range(len(name_parts)):
            if importlib.util.find_spec(".".join(name_parts[: i + 1])) is None:
                return False
        return True
    else:  # pragma: no cover
        import imp

        try:
            imp.find_module(module_name)
        except ImportError:
            return False
        return True
