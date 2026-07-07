import os
import re
import sys
import uuid
import pydoc
import pathlib
import importlib
import functools

import time
import warnings
from siren import dataclasses as _dataclasses
from siren import math as _math
from siren.interactions import DarkNewsCrossSection,DarkNewsDecay
import numpy as np
import awkward as ak
import h5py
import pickle
import logging

# Set up logging configuration
logger = logging.getLogger(__name__)

class CustomFormatter(logging.Formatter):

    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s (%(pathname)s:%(lineno)d)"

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

if not logger.handlers:
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    #console_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))
    console_handler.setFormatter(CustomFormatter())
    logger.addHandler(console_handler)
    logger.setLevel(logging.WARN)

def log_newline(n=1):
    blank_fmt = logging.Formatter('')
    formatters = []
    for handler in logger.handlers:
        formatters.append(handler.formatter)
        handler.setFormatter(blank_fmt)

    for i in range(n):
        logger.warning(" \n")

    for handler, formatter in zip(logger.handlers, formatters):
        handler.setFormatter(formatter)

try:
    from DarkNews.nuclear_tools import NuclearTarget
except:
    pass

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

    Get a list of directories where siren resources may be located.
    The first directory in this list is the "resources" directory in
    the package itself. The second directory is the appdata directory
    (~/.siren on Linux). The list further contains the application
    directory (for frozen apps), and may include additional directories
    in the future.
    """
    dirs = [resource_package_dir()]
    # Resource dir baked in the package.
    # Appdata directory
    try:
        dirs.append(appdata_dir("siren"))
    except Exception:  # pragma: no cover
        pass  # The home dir may not be writable
    # Directory where the app is located (mainly for frozen apps)
    if getattr(sys, "frozen", None):
        dirs.append(os.path.abspath(os.path.dirname(sys.executable)))
    elif sys.path and sys.path[0]:
        dirs.append(os.path.abspath(sys.path[0]))
    return dirs


# from imageio
def resource_package_dir():
    """Get the resources directory inside the installed siren package.

    This is a convenience method used by ``resource_dirs`` and siren
    entry point scripts.
    """
    from importlib.resources import files
    return str(files("siren") / "resources")


# from imageio
# https://github.com/imageio/imageio/blob/65d79140018bb7c64c0692ea72cb4093e8d632a0/imageio/core/util.py
def get_platform():
    """get_platform()

    Get a string that specifies the platform more specific than
    sys.platform does. The result can be: linux32, linux64, win32,
    win64, osx32, osx64. Other platforms may be added in the future.
    """
    import struct
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


def load_module(name, path, persist=True):
    """Load a module with a specific name and path"""
    url = pathlib.Path(os.path.abspath(path)).as_uri()
    #module_name = f"{name}-{str(uuid.uuid5(uuid.NAMESPACE_URL, url))}"
    module_name = name
    if module_name in sys.modules:
        return sys.modules[module_name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception as e:
        del sys.modules[module_name]
        raise
    module = sys.modules[module_name]
    if not persist:
        del sys.modules[module_name]
    return module


_VERSION_PATTERN = r"""
    v?
    (?:
        (?:(?P<epoch>[0-9]+)!)?                           # epoch
        (?P<release>[0-9]+(?:\.[0-9]+)*)                  # release segment
        (?P<pre>                                          # pre-release
            [-_\.]?
            (?P<pre_l>(a|b|c|rc|alpha|beta|pre|preview))
            [-_\.]?
            (?P<pre_n>[0-9]+)?
        )?
        (?P<post>                                         # post release
            (?:-(?P<post_n1>[0-9]+))
            |
            (?:
                [-_\.]?
                (?P<post_l>post|rev|r)
                [-_\.]?
                (?P<post_n2>[0-9]+)?
            )
        )?
        (?P<dev>                                          # dev release
            [-_\.]?
            (?P<dev_l>dev)
            [-_\.]?
            (?P<dev_n>[0-9]+)?
        )?
    )
    (?:\+(?P<local>[a-z0-9]+(?:[-_\.][a-z0-9]+)*))?       # local version
"""

_version_regex = re.compile(
    r"^\s*" + _VERSION_PATTERN + r"\s*$",
    re.VERBOSE | re.IGNORECASE,
)

_UNVERSIONED_MODEL_PATTERN = (
    r"""
    (?P<model_name>
        (?:
            [a-zA-Z0-9]+
        )
        |
        (?:
            (?:[a-zA-Z0-9]+(?:[-_\.][a-zA-Z0-9]+)*(?:[-_\.][a-zA-Z]+[a-zA-Z0-9]*))?
        )
    )
    (?:
        -
        (?P<version>"""
    + _VERSION_PATTERN
    + r"))?"
)

_MODEL_PATTERN = (
    r"""
    (?P<model_name>
        [a-zA-Z0-9]+(?:(?:-(?!v?[0-9])|[_\.])[a-zA-Z0-9]+)*
    )
    (?:
        -
        (?P<version>"""
    + _VERSION_PATTERN
    + r"))?"
)

_model_regex = re.compile(
    r"^\s*" + _MODEL_PATTERN + r"\s*$",
    re.VERBOSE | re.IGNORECASE,
)

def decompose_version(version):
    # Break the version string into its components
    matches = _version_regex.match(version)
    if matches is None:
        return dict()
    else:
        return matches.groupdict()


def normalize_version(version):
    # Normalize the version string
    d = decompose_version(version)
    n_version = ""

    # Add epoch if present
    if d["epoch"] is not None:
        n_version += str(int(d["epoch"])) + "!"

    # Add release segment
    if d["release"] is not None:
        # Remove leading zeros from each segment
        n_version += ".".join(
            [str(int(s)) for s in d["release"].strip(". ").lower().split(".")]
        )

    # Add pre-release segment
    if d["pre"] is not None:
        if d["pre_l"] is not None:
            # Add pre-release label
            if d["pre_l"] in ["a", "alpha"]:
                n_version += "a"
            elif d["pre_l"] in ["b", "beta"]:
                n_version += "b"
            elif d["pre_l"] in ["c", "rc", "pre", "preview"]:
                n_version += "rc"
            # Add pre-release number
            if d["pre_n"] is not None:
                n_version += str(int(d["pre_n"]))

    # Add post-release segment
    if d["post"] is not None:
        n_version += ".post"
        # Add post-release number
        if d["post_n1"] is not None:
            n_version += str(int(d["post_n1"]))
        elif d["post_n2"] is not None:
            n_version += str(int(d["post_n2"]))

    # Add dev-release segment
    if d["dev"] is not None:
        n_version += ".dev"
        # Add dev-release number
        if d["dev_n"] is not None:
            n_version += str(int(d["dev_n"]))

    # Add local segment
    if d["local"] is not None:
        n_version += "+"
        segments = []
        # Local segment can contain letters and numbers
        for s in d["local"].lower().split("."):
            # Remove leading zeros from each segment if it is a number
            try:
                segments.append(str(int(s)))
            except:
                segments.append(s)
        # Join the segments with a dot
        n_version += ".".join(segments)
    return n_version


def tokenize_version(version):
    # Tokenize the version string
    d = decompose_version(normalize_version(version))
    tokens = []
    tokens = dict()

    if d["epoch"] is not None:
        # Add epoch if present
        tokens["epoch"] = int(d["epoch"])
    else:
        # Default epoch is 0
        tokens["epoch"] = 0

    # Add release segment
    if d["release"] is not None:
        # Remove leading zeros from each segment
        tokens["release"] = tuple(
            [int(s) for s in d["release"].strip(". ").lower().split(".")]
        )

    # Add pre-release segment
    if d["pre"] is not None:
        # Pre-release segment starts with 0 to indicate that it is present
        # Because a pre-release version comes before a release version
        pre_token = [0]

        if d["pre_l"] is not None:
            if d["pre_l"] in ["a", "alpha"]:
                pre_token.append(1)
            elif d["pre_l"] in ["b", "beta"]:
                pre_token.append(2)
            elif d["pre_l"] in ["c", "rc", "pre", "preview"]:
                pre_token.append(3)
            if d["pre_n"] is not None:
                pre_token.append(int(d["pre_n"]))
            else:
                pre_token.append(0)
            tokens["pre"] = tuple(pre_token)
        else:
            tokens["pre"] = tuple(pre_token)
    else:
        # Pre-release segment starts with 1 if not present to indicate that it is not present
        # Because a release version comes after a pre-release version
        tokens["pre"] = (1,)

    # Add post-release segment
    if d["post"] is not None:
        # If post-release segment is present, it starts with 1
        # Because a post-release version comes after its corresponding (pre-)release version
        post_token = [1]
        if d["post_n1"] is not None:
            post_token.append(int(d["post_n1"]))
        elif d["post_n2"] is not None:
            post_token.append(int(d["post_n2"]))
        else:
            post_token.append(0)
        tokens["post"] = tuple(post_token)
    else:
        # If post-release segment is not present, it starts with 0
        # Because a post-release version comes after its corresponding (pre-)release version
        tokens["post"] = (0,)

    # Add dev-release segment
    if d["dev"] is not None:
        # If dev-release segment is present, it starts with 0
        # Because a dev-release version comes before its corresponding (pre-,post-)release version
        dev_token = [0]
        if d["dev_n"] is not None:
            dev_token.append(int(d["dev_n"]))
        else:
            dev_token.append(0)
        tokens["dev"] = tuple(dev_token)
    else:
        # If dev-release segment is not present, it starts with 1
        # Because a dev-release version comes before its corresponding (pre-,post-)release version
        tokens["dev"] = (1,)

    # Add local segment
    if d["local"] is not None:
        # If local segment is present, it starts with 1
        # Because a local version comes after its corresponding (pre-,post-,dev-)release version
        local_token = [1]
        for s in d["local"].lower().split("."):
            try:
                local_token.append((1, int(s)))
            except:
                local_token.append((0, s))
        tokens["local"] = tuple(local_token)
    else:
        # If local segment is not present, it starts with 0
        # Because a local version comes after its corresponding (pre-,post-,dev-)release version
        tokens["local"] = (0,)

    # Return the tokenized version
    token_list = [
        tokens["epoch"],
        tokens["release"],
        tokens["pre"],
        tokens["post"],
        tokens["dev"],
        tokens["local"],
    ]

    return tuple(token_list)


def _get_base_directory(resources_dir, prefix):
    base_dir = resources_dir
    if prefix is not None:
        base_dir = os.path.join(base_dir, prefix)
    return base_dir

def _find_model_folder_and_file(base_dir, model_name, must_exist, specific_file=None):
    model_names = [
        f for f in os.listdir(base_dir) if not os.path.isfile(os.path.join(base_dir, f))
    ]

    exact_model_names = [f for f in model_names if f.lower() == model_name.lower()]

    if len(exact_model_names) == 0:
        model_names = [f for f in model_names if f.lower().startswith(model_name.lower())]
    else:
        model_names = exact_model_names

    if len(model_names) == 0 and must_exist:
        raise ValueError(f"No model folders found for {model_name}\nSearched in {base_dir}")
    elif len(model_names) == 0 and not must_exist:
        return model_name, False, None
    elif len(model_names) == 1:
        name = model_names[0]
        if specific_file is not None:
            specific_file_path = os.path.join(base_dir, name, specific_file)
            if os.path.isfile(specific_file_path):
                return name, True, specific_file_path
            else:
                return name, True, None
        else:
            return name, True, None
    else:
        raise ValueError(f"Multiple directories found for {model_name}\nSearched in {base_dir}")

def _get_model_files(base_dir, model_name, is_file, folder_exists, version=None):
    if folder_exists:
        if version:
            version_dir = os.path.join(base_dir, model_name, f"v{version}")
            if os.path.isdir(version_dir):
                return [
                    f for f in os.listdir(version_dir)
                    if is_file == os.path.isfile(os.path.join(version_dir, f))
                ]
        return [
            f for f in os.listdir(os.path.join(base_dir, model_name))
            if is_file == os.path.isfile(os.path.join(base_dir, model_name, f))
        ]
    return []

def _extract_model_versions(model_files, model_regex, model_name):
    model_versions = []
    for f in model_files:
        d = model_regex.match(f)
        if d is not None:
            if d.groupdict()["version"] is not None:
                model_versions.append(normalize_version(d.groupdict()["version"]))
            else:
                logging.warning(f"Input model file has no version: {f}")
        elif f.lower().startswith(model_name.lower()):
            logging.warning(f"Unable to parse version from {f}")
    return model_versions

def _get_model_file_name(version, model_versions, model_files, model_name, suffix, must_exist):
    if version is None and must_exist:
        version_idx, version = max(enumerate(model_versions), key=lambda x: tokenize_version(x[1]))
        return model_files[version_idx]
    elif version is None and not must_exist:
        version = "v1"
        return f"{model_name}-v{version}{suffix}"
    else:
        version = normalize_version(version)
        if must_exist:
            if version not in model_versions:
                raise ValueError(f"No model found for {model_name}-{version}")
            version_idx = model_versions.index(version)
            return model_files[version_idx]
        else:
            if version in model_versions:
                version_idx = model_versions.index(version)
                return model_files[version_idx]
            else:
                return f"{model_name}-v{version}{suffix}"


def _get_model_folder(base_dir, model_name, must_exist):
    model_names = [
        f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))
    ]

    exact_model_names = [f for f in model_names if f.lower() == model_name.lower()]

    if len(exact_model_names) == 0:
        model_names = [f for f in model_names if f.lower().startswith(model_name.lower())]
    else:
        model_names = exact_model_names

    if len(model_names) == 0 and must_exist:
        raise ValueError(f"No model folders found for {model_name}\nSearched in {base_dir}")
    elif len(model_names) == 0 and not must_exist:
        return model_name, False
    elif len(model_names) == 1:
        return model_names[0], True
    else:
        raise ValueError(f"Multiple directories found for {model_name}\nSearched in {base_dir}")

def _get_model_subfolders(base_dir, model_regex):
    model_subfolders = [
        f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))
    ]
    model_subfolders = [
        f for f in model_subfolders if model_regex.match(f) is not None
    ]
    return model_subfolders


def _get_model_path(model_name, prefix=None, suffix=None, is_file=True, must_exist=True, specific_file=None):
    model_regex = re.compile(
        r"^\s*" + _MODEL_PATTERN + ("" if suffix is None else r"(?:" + suffix + r")?") + r"\s*$",
        re.VERBOSE | re.IGNORECASE,
    )
    suffix = "" if suffix is None else suffix

    resources_dir = resource_package_dir()
    base_dir = _get_base_directory(resources_dir, prefix)

    d = model_regex.match(model_name)
    if d is None:
        raise ValueError(f"Invalid model name: {model_name}")
    d = d.groupdict()
    model_search_name, version = d["model_name"], d["version"]

    if version is not None:
        version = normalize_version(version)

    found_model_name, folder_exists = _get_model_folder(base_dir, model_search_name, must_exist)

    model_dir = os.path.join(base_dir, found_model_name)

    if not must_exist and not folder_exists:
        if version is None:
            model_dir = os.path.join(model_dir, f"{found_model_name}-v1")
        else:
            model_dir = os.path.join(model_dir, f"{found_model_name}-v{version}")

        return model_dir

    top_level_has_specific_file = specific_file is not None and os.path.isfile(os.path.join(model_dir, specific_file))

    if version is None and top_level_has_specific_file:
        return model_dir

    model_subfolders = _get_model_subfolders(model_dir, model_regex)

    if len(model_subfolders) == 0:
        if top_level_has_specific_file:
            return model_dir
        if must_exist:
            raise ValueError(f"No model folders found for {model_search_name}\nSearched in {model_dir}")
        else:
            if version is None:
                version = "v1"

            model_dir = os.path.join(model_dir, f"{found_model_name}-v{version}")
            return model_dir

    models_and_versions = []
    for f in model_subfolders:
        d = model_regex.match(f).groupdict()
        if d["version"] is not None:
            models_and_versions.append((f, normalize_version(d["version"])))

    matching_models = [(m, v) for m, v in models_and_versions if v == version]

    if len(matching_models) == 1:
        model_dir = os.path.join(model_dir, matching_models[0][0])
        return model_dir
    elif len(matching_models) > 1:
        raise ValueError(f"Multiple directories found for {model_search_name} with version {version}\nSearched in {model_dir}")

    if top_level_has_specific_file:
        return model_dir

    if len(matching_models) == 0:
        if must_exist and version is not None:
            raise ValueError(f"No model folders found for {model_search_name} with version {version}\nSearched in {model_dir}")

    found_model_subfolder, subfolder_version = max(models_and_versions, key=lambda x: tokenize_version(x[1]))

    return os.path.join(model_dir, found_model_subfolder)


_resource_folder_by_name = {
    "flux": "fluxes",
    "detector": "detectors",
    "processes": "processes",
}


def get_flux_model_path(model_name, must_exist=True):
    return _get_model_path(model_name, prefix=_resource_folder_by_name["flux"], is_file=False, must_exist=must_exist, specific_file=f"flux.py")


def get_detector_model_path(model_name, must_exist=True):
    return _get_model_path(model_name, prefix=_resource_folder_by_name["detector"], is_file=False, must_exist=must_exist, specific_file=f"detector.py")


def get_processes_model_path(model_name, must_exist=True):
    return _get_model_path(model_name, prefix=_resource_folder_by_name["processes"], is_file=False, must_exist=must_exist, specific_file="processes.py")

def import_resource(resource_type, resource_name):
    folder = _resource_folder_by_name[resource_type]
    specific_file = f"{resource_type}.py"

    abs_dir = _get_model_path(resource_name, prefix=folder, is_file=False, must_exist=True, specific_file=specific_file)

    fname = os.path.join(abs_dir, f"{resource_type}.py")
    if not os.path.isfile(fname):
        # No .py loader script; caller falls back to file-based loading.
        return None
    try:
        mod = load_module(f"siren-{resource_type}-{resource_name}", fname, persist=False)
    except Exception as e:
        log_newline()
        logger.warning(f"Encountered exception while loading '{resource_type}' '{resource_name}' from '{fname}':\n\t{e}")
        raise e
    return mod


def get_resource_loader(resource_type, resource_name):
    resource_module = import_resource(resource_type, resource_name)
    if resource_module is None:
        return None
    loader_name = f"load_{resource_type}"
    if not hasattr(resource_module, loader_name):
        logger.warning(f"from '{resource_module.__file__}' module '{resource_module.__name__}' has no attribute '{loader_name}'")
        raise AttributeError(f"from '{resource_module.__file__}' module '{resource_module.__name__}' has no attribute '{loader_name}'")
    loader = getattr(resource_module, loader_name)
    class Functor:
        def __init__(self, func):
            self.func = func
        def __call__(self, *args, **kwargs):
            return self.func(*args, **kwargs)
        def _repr_pretty_(self, p, cycle):
            p.pretty(self.func)
    functor = Functor(loader)
    for key in dir(resource_module):
        if key.startswith("__"):
            continue
        if key == f"load_{resource_type}":
            continue
        setattr(functor, key, getattr(resource_module, key))
    return functools.update_wrapper(functor, loader)


def get_tabulated_flux_model_path(model_name, must_exist=True):
    return _get_model_path(model_name, prefix=_resource_folder_by_name["flux"], is_file=False, must_exist=must_exist)


def get_tabulated_flux_file(model_name, tag, must_exist=True):
    return load_resource("flux", model_name, tag)


def load_resource(resource_type, resource_name, *args, **kwargs):
    loader = get_resource_loader(resource_type, resource_name)
    if loader is None:
        return None
    resource = loader(*args, **kwargs)
    return resource


def load_flux(model_name, *args, **kwargs):
    return load_resource("flux", model_name, *args, **kwargs)


def _detector_file_loader(model_name):
    resource_name = model_name
    folder = _resource_folder_by_name["detector"]

    abs_dir = _get_model_path(resource_name, prefix=folder, is_file=False, must_exist=True, specific_file=None)

    densities_fname = os.path.join(abs_dir, "densities.dat")
    materials_fname = os.path.join(abs_dir, "materials.dat")

    if os.path.isfile(densities_fname) and os.path.isfile(materials_fname):
        from . import detector as _detector
        detector_model = _detector.DetectorModel()
        detector_model.LoadMaterialModel(materials_fname)
        detector_model.LoadDetectorModel(densities_fname)
        return detector_model

    raise ValueError("Could not find detector loading script \"{script_fname}\" or densities and materials files \"{densities_fname}\", \"materials_fname\"")


def load_detector(model_name, *args, **kwargs):
    resource = load_resource("detector", model_name, *args, **kwargs)
    if resource is not None:
        return resource
    return _detector_file_loader(model_name)


class ProcessBundle:
    """Normalized return type from ``load_processes``.

    Behaves like a 2-tuple ``(primary, secondary)`` for unpacking::

        primary, secondary = siren.load_processes("CSMSDISSplines", ...)

    Any extra return values from the process loader are available as
    the ``metadata`` attribute (a tuple of additional return values).
    """

    def __init__(self, primary, secondary, *extra):
        self.primary = primary
        self.secondary = secondary
        self.metadata = extra

    def __iter__(self):
        yield self.primary
        yield self.secondary

    def __len__(self):
        return 2

    def __getitem__(self, idx):
        return (self.primary, self.secondary)[idx]

    def __repr__(self):
        n_primary = sum(len(v) for v in self.primary.values())
        n_secondary = sum(len(v) for v in self.secondary.values())
        extra = f", +{len(self.metadata)} metadata" if self.metadata else ""
        return f"ProcessBundle({n_primary} primary, {n_secondary} secondary{extra})"


def load_processes(model_name, *args, **kwargs):
    result = load_resource("processes", model_name, *args, **kwargs)
    if result is None:
        raise ValueError(
            f"No process loader found for model '{model_name}'. "
            "Check the model name and that its resource directory is installed.")
    if isinstance(result, tuple):
        if len(result) >= 2:
            return ProcessBundle(result[0], result[1], *result[2:])
        return ProcessBundle(result[0], {})
    return ProcessBundle(result, {})


def get_detector_model_targets(detector_model):
    """Return the set of target ParticleTypes (nuclei) present in *detector_model*.

    These are the material targets used by depth/range-based vertex
    distributions (e.g. RangePositionDistribution), not the primary
    projectile types.
    """
    targets = set()
    count = 0
    while detector_model.Materials.HasMaterial(count):
        targets.update(detector_model.Materials.GetMaterialTargets(count))
        count += 1
    return targets

def get_fiducial_volume(experiment):
    """
    :return: identified fiducial volume for the experiment, None if not found
    """
    detector_model_file = get_detector_model_path(experiment) + "/densities.dat"
    with open(detector_model_file) as file:
        fiducial_line = None
        detector_line = None
        for line in file:
            data = line.split()
            if len(data) <= 0:
                continue
            elif data[0] == "fiducial":
                fiducial_line = line
            elif data[0] == "detector":
                detector_line = line
        if fiducial_line is None or detector_line is None:
            return None
        from . import detector as _detector
        return _detector.DetectorModel.ParseFiducialVolume(fiducial_line, detector_line)
    return None

def get_volume_position_distribution_from_sector(detector_model, sector_name):
    """Create a position distribution from a named detector sector.

    Extracts the geometry from the sector, converts coordinates from
    geometry frame to detector frame, and returns the appropriate
    volume position distribution (Cylinder or Sphere).

    Parameters
    ----------
    detector_model : DetectorModel
        The loaded detector model.
    sector_name : str
        Name of the sector to use (e.g. "tilecal", "fiducial").

    Returns
    -------
    CylinderVolumePositionDistribution or SphereVolumePositionDistribution
    """
    from . import detector as _detector
    from . import geometry as _geometry
    from . import distributions as _distributions

    geo = None
    for sector in detector_model.Sectors:
        if sector.name == sector_name:
            geo = sector.geo
            break
    if geo is None:
        available = [s.name for s in detector_model.Sectors]
        raise ValueError(
            f"Sector {sector_name!r} not found. Available: {available}"
        )

    det_position = detector_model.GeoPositionToDetPosition(
        _detector.GeometryPosition(geo.placement.Position)
    )
    det_rotation = geo.placement.Quaternion
    det_placement = _geometry.Placement(det_position.get(), det_rotation)

    if isinstance(geo, _geometry.Cylinder):
        cylinder = _geometry.Cylinder(
            det_placement, geo.Radius, geo.InnerRadius, geo.Z
        )
        return _distributions.CylinderVolumePositionDistribution(cylinder)
    elif isinstance(geo, _geometry.Sphere):
        sphere = _geometry.Sphere(
            det_placement, geo.Radius, geo.InnerRadius
        )
        return _distributions.SphereVolumePositionDistribution(sphere)
    else:
        raise TypeError(
            f"Sector geometry type {type(geo).__name__} not supported "
            f"for position distribution"
        )


def list_fluxes():
    return sorted(_get_model_subfolders(_get_base_directory(resource_package_dir(), "fluxes"), _model_regex))

def list_detectors():
    dirs = sorted(_get_model_subfolders(_get_base_directory(resource_package_dir(), "detectors"), _model_regex))
    dirs = [d for d in dirs if d != "visuals"]
    return dirs

def list_processes():
    return sorted(_get_model_subfolders(_get_base_directory(resource_package_dir(), "processes"), _model_regex))

def flux_docs(flux_name):
    loader = get_resource_loader("flux", flux_name)
    if loader is None:
        raise ValueError(f"Could not find documentation for flux {flux_name}")
    return loader.__doc__

def detector_docs(detector_name):
    loader = get_resource_loader("detector", detector_name)
    if loader is not None:
        return loader.__doc__

    resource_name = detector_name
    folder = _resource_folder_by_name["detector"]

    abs_dir = _get_model_path(resource_name, prefix=folder, is_file=False, must_exist=True, specific_file=None)

    densities_fname = os.path.join(abs_dir, "densities.dat")
    materials_fname = os.path.join(abs_dir, "materials.dat")

    lines = []
    if os.path.isfile(densities_fname):
        with open(densities_fname) as file:
            new_lines = []
            for l in file.readlines():
                l = l.strip()
                if l.startswith("#"):
                    new_lines.append(l)
                else:
                    break
            if len(new_lines) > 0:
                lines.append(f"Detector definition: {densities_fname}")
            lines.extend(new_lines)

    if os.path.isfile(materials_fname):
        with open(materials_fname) as file:
            new_lines = []
            for l in file.readlines():
                l = l.strip()
                if l.startswith("#"):
                    new_lines.append(l)
                else:
                    break
            if len(lines) > 0 and len(new_lines) > 0:
                lines.append("")
                lines.append(f"Material definitions: {materials_fname}")
            lines.extend(new_lines)

    doc = "\n".join(lines)

    if len(lines) == 0:
        raise ValueError(f"Could not find documentation for detector {detector_name}")

    return doc


def process_docs(process_name):
    loader = get_resource_loader("processes", process_name)
    if loader is None:
        raise ValueError(f"Could not find documentation for process {process_name}")
    return loader.__doc__

def flux_help(flux_name):
    doc = flux_docs(flux_name)
    pydoc.pager(doc)

def detector_help(detector_name):
    doc = detector_docs(detector_name)
    pydoc.pager(doc)

def process_help(process_name):
    doc = process_docs(process_name)
    pydoc.pager(doc)

def _get_process_loader(process_name):
    return get_resource_loader("processes", process_name)

def _get_flux_loader(flux_name):
    return get_resource_loader("flux", flux_name)

def _get_detector_loader(detector_name):
    loader = get_resource_loader("detector", detector_name)
    if loader is not None:
        return loader

    def load_detector():
        return _detector_file_loader(detector_name)

    load_detector.__doc__ = detector_docs(detector_name)

    return load_detector



###### Injector helper functions #######

# Generate events using an injector object, returning (events, gen_times).
def GenerateEvents(injector, N=None):
    """Generate up to N events from an injector; return (events, gen_times).

    Deprecated: use ``siren.generate()`` or ``Injector.generate()`` instead.

    For the python-wrapper ``Injector``, generation is delegated to
    ``injector.generate(N, on_shortfall='ignore')``, which retries failed
    attempts internally and returns only successful (non-empty) trees --
    unlike the historical raw loop, a failed ``GenerateEvent`` attempt is
    never appended to ``events``. ``gen_times`` cannot be measured per event
    through that call, so it is synthesized: the whole call is timed once and
    the total elapsed time is divided evenly across the returned events (an
    aggregate estimate, not a per-event wall-clock measurement).

    For a raw (non-wrapper) injector -- one with no ``generate`` method, e.g.
    the engine's own ``_Injector`` -- the historical per-attempt loop over
    ``generate_event()``/``injected_events`` is used unchanged, including
    appending empty trees on failure and recording a true per-event
    wall-clock time for each attempt.
    """
    warnings.warn(
        "siren._util.GenerateEvents is deprecated; use siren.generate() or "
        "Injector.generate()",
        DeprecationWarning, stacklevel=2)

    from .Injector import Injector as _PyInjector

    if N is None:
        N = injector.number_of_events

    if isinstance(injector, _PyInjector):
        t0 = time.time()
        events = injector.generate(N, on_shortfall="ignore")
        elapsed = time.time() - t0
        n = len(events)
        per_event = (elapsed / n) if n > 0 else 0.0
        gen_times = [per_event] * n
        return events, gen_times

    count = 0
    gen_times = []
    prev_time = time.time()
    events = []
    while (injector.injected_events < injector.number_of_events) and (count < N):
        print("Injecting Event %d/%d  " % (count, N), end="\r")
        event = injector.generate_event()
        events.append(event)
        t = time.time()
        gen_times.append(t-prev_time)
        prev_time = t
        count += 1
    return events,gen_times

def get_parent_indices(tree):
    """Return the parent interaction index for each datum in an InteractionTree.

    ``tree`` is the ordered list of ``InteractionTreeDatum`` objects, i.e.
    ``event.tree``. The returned list holds, for each datum, the index (within
    ``tree``) of the interaction whose secondary particle became this datum's
    primary, or ``-1`` for a root/primary interaction that has no parent.

    Parentage is read directly from each datum's ``.parent_index``
    """
    parent_indices = []
    for datum in tree:
        if datum.is_root():
            parent_indices.append(-1)
        else:
            parent_indices.append(int(datum.parent_index))
    return parent_indices


def _event_weight(weighter, tree):
    """Central-value weight of one tree via a Weighter instance or a callable.
    EventWeight is preferred so a Weighter that also happens to be callable is
    still invoked through its documented method."""
    if hasattr(weighter, "EventWeight"):
        return float(weighter.EventWeight(tree))
    return float(weighter(tree))


def _set_header_cv(tree, w):
    """Overwrite only the CV slot (index 0) of a tree's header weights, keeping
    any additional named weights. Mirrors convert_siren_events_to_hepmc3."""
    weights = list(tree.header.weights) or [0.0]
    weights[0] = float(w)
    tree.header.weights = weights


def _trees_carry_trusted_header_weights(trees):
    """True if any tree has non-empty header weights AND its provenance does not
    mark them as the "unweighted" placeholder CV.

    A tree loaded back from an "unweighted" HepMC3 export carries a placeholder
    CV (commonly 1.0) in ``header.weights`` purely because HepMC3 always stores
    a weight vector; its ``header.provenance["siren.weights_state"]`` (restored
    by the reader from the run-level ``siren.weights_state`` attribute) says
    "unweighted". Such placeholders must never be promoted to a trusted
    "header" state -- only trees with no recorded provenance (e.g. built fresh
    in-process) or with provenance explicitly "header"/"computed" count.
    """
    for tree in trees:
        if len(tree.header.weights) == 0:
            continue
        state = tree.header.provenance.get("siren.weights_state")
        if state is not None and str(state).lower() == "unweighted":
            continue
        return True
    return False


def _header_cv_weights(events):
    """Read back the CV slot (index 0) of each tree's header weights.

    Used when the resolved policy trusts existing headers as-is (state
    "header") so the tabular event_weight column agrees with the CV already
    written to the native/HepMC3 headers, instead of being recomputed on a
    different path. A tree with no header weight defaults to 1.0, matching the
    placeholder the writer substitutes for an empty header, so the tabular
    column and the HepMC3 file agree on that tree too.
    """
    return [float(tree.header.weights[0]) if len(tree.header.weights) > 0 else 1.0
            for tree in events]


def resolve_hepmc3_weight_policy(events, hepmc3_weights, weighter):
    """Resolve the HepMC3 weight policy for a set of trees.

    ``hepmc3_weights`` selects how per-event central-value weights reach the
    HepMC3 (and native) output:

    - ``"auto"``: if a weighter is available, compute every event weight once and
      populate the tree headers (state ``"computed"``); else if any tree already
      carries non-empty header weights NOT themselves marked
      ``siren.weights_state == "unweighted"`` in provenance, trust them (state
      ``"header"``); else leave headers untouched (state ``"unweighted"``). This
      keeps a load-then-reexport of an unweighted HepMC3 file from laundering
      its placeholder CV into a trusted state.
    - ``"none"``: leave headers untouched (state ``"unweighted"``).
    - ``"header"``: trust the existing header weights (state ``"header"``).
    - a Weighter instance (has ``EventWeight``) or a callable: compute weights the
      same as auto-with-weighter (state ``"computed"``).
    - a sequence of floats: one CV weight per event, written to the headers
      (state ``"header"``); length must equal ``len(events)``.

    Headers are populated in place so both the native ``.siren_events`` file and
    the HepMC3 file carry the same CV. The tabular-convention sentinel
    (weighter-is-None -> weight 0) is NEVER written into a header: unweighted
    means untouched headers, never zeros.

    Returns ``(cv_weights, weights_state)`` where ``cv_weights`` is a list of the
    per-event CV weights -- computed here, or read back out of each tree's
    trusted header CV slot (index 0) under a ``"header"`` state -- so a caller
    (the tabular ``event_weight`` column) always agrees with what was written
    to the native/HepMC3 headers. ``cv_weights`` is ``None`` only under
    ``"unweighted"``, where there is no meaningful per-event weight at all.
    """
    # Explicit per-event sequence of floats.
    if hepmc3_weights is not None and not isinstance(hepmc3_weights, str) \
            and not callable(hepmc3_weights) and not hasattr(hepmc3_weights, "EventWeight"):
        try:
            seq = list(hepmc3_weights)
        except TypeError:
            seq = None
        if seq is not None:
            if len(seq) != len(events):
                raise ValueError(
                    "hepmc3_weights sequence length %d != number of events %d"
                    % (len(seq), len(events)))
            for tree, w in zip(events, seq):
                _set_header_cv(tree, w)
            return [float(w) for w in seq], "header"

    # A weighter passed directly (callable or Weighter instance) -> compute.
    direct_weighter = None
    if callable(hepmc3_weights) or hasattr(hepmc3_weights, "EventWeight"):
        direct_weighter = hepmc3_weights

    if hepmc3_weights is None or hepmc3_weights == "none":
        return None, "unweighted"

    if hepmc3_weights == "header":
        return _header_cv_weights(events), "header"

    # "auto" (or a directly-passed weighter): prefer a weighter, then existing
    # headers, then unweighted.
    active = direct_weighter if direct_weighter is not None else weighter
    if hepmc3_weights == "auto" or direct_weighter is not None:
        if active is not None:
            cv = [_event_weight(active, tree) for tree in events]
            for tree, w in zip(events, cv):
                _set_header_cv(tree, w)
            return cv, "computed"
        # No weighter: trust existing headers only if they are not themselves
        # the "unweighted" placeholder CV round-tripped from a prior export.
        if _trees_carry_trusted_header_weights(events):
            return _header_cv_weights(events), "header"
        return None, "unweighted"

    raise ValueError("unrecognized hepmc3_weights policy: %r" % (hepmc3_weights,))


def SaveEvents(events,
               weighter=None,
               gen_times=None,
               save_hdf5=True,
               save_parquet=True,
               save_siren_events=True,
               save_hepmc3=False,
               hepmc3_weights="auto",
               hepmc3_gzip=False,
               fid_vol=None,
               injector=None,
               output_filename=None,
               pot=None):

    # pot is recorded only as an HDF5 attribute; accepting it while HDF5
    # output is disabled would drop it silently.
    if pot is not None and not save_hdf5:
        raise ValueError(
            "pot is recorded in the HDF5 output; pass save_hdf5=True or omit pot")

    # An omitted gen_times must not crash the per-event loop below (which
    # unconditionally indexes it) after files have already been written;
    # default to one zero per event, matching the shape SIREN_Controller
    # supplies (a plain list of floats, one per event).
    if gen_times is None:
        gen_times = [0.0] * len(events)

    # Resolve the HepMC3 weight policy up front so per-event central values are
    # computed exactly once (BEFORE any file is written) and shared by the native
    # file, the HepMC3 file, and the tabular datasets below. Under "auto" the
    # available weighter populates the headers; the computed array is reused for
    # the event_weight column so the CV is never recomputed.
    hepmc3_cv = None
    hepmc3_state = "unweighted"
    if save_hepmc3 or save_siren_events:
        hepmc3_cv, hepmc3_state = resolve_hepmc3_weight_policy(
            events, hepmc3_weights, weighter)

    # Optionally save things. Headers are already populated (when the policy chose
    # to) so the native .siren_events file carries the same CV as the HepMC3 file.
    if save_siren_events: _dataclasses.SaveInteractionTrees(events, output_filename)
    if save_hepmc3:
        from . import hepmc3 as _hepmc3
        opts = _hepmc3.HepMC3WriterOptions()
        opts.weights_state = hepmc3_state
        opts.gzip = bool(hepmc3_gzip)
        # Generation counts (from the injector, when supplied) as run metadata:
        # attempted/accepted normalize the FATX; events_to_inject is the pooled-
        # weighting seed N_i. Absent when no injector is passed.
        if injector is not None:
            opts.attempted_events = int(injector.InjectionAttempts())
            opts.accepted_events = int(injector.InjectedEvents())
            opts.events_to_inject = int(injector.EventsToInject())
        out = output_filename + ".hepmc3"
        if hepmc3_gzip and not out.endswith(".gz"):
            out = out + ".gz"
        _hepmc3.SaveInteractionTreesAsHepMC3(events, out, opts)
    # A dictionary containing each dataset we'd like to save
    datasets = {
        "event_weight":[], # weight of entire event
        "event_gen_time":[], # generation time of each event
        "event_weight_time":[], # generation time of each event
        "num_interactions":[], # number of interactions per event
        "vertex":[], # vertex of each interaction in an event
        "in_fiducial":[], # whether or not each vertex is in the fiducial volume
        "primary_type":[], # primary type of each interaction
        "target_type":[], # target type of each interaction
        "num_secondaries":[], # number of secondary particles of each interaction
        "secondary_types":[], # secondary type of each interaction
        "primary_momentum":[], # primary momentum of each interaction
        "secondary_momenta":[], # secondary momentum of each interaction
        "parent_idx":[], # index of the parent interaction
    }
    for ie, event in enumerate(events):
        print("Saving Event %d/%d  " % (ie, len(events)), end="\r")
        t0 = time.time()
        if hepmc3_cv is not None:
            datasets["event_weight"].append(hepmc3_cv[ie])  # reuse the CV computed above
        elif weighter is None:
            datasets["event_weight"].append(0)
        elif callable(weighter):
            datasets["event_weight"].append(weighter(event))
        else:
            # A siren.injection.Weighter object is not callable; it exposes
            # EventWeight (same duck-typing resolve_hepmc3_weight_policy uses).
            datasets["event_weight"].append(weighter.EventWeight(event))
        datasets["event_weight_time"].append(time.time()-t0)
        datasets["event_gen_time"].append(gen_times[ie])
        # add empty lists for each per interaction dataset
        for k in ["vertex",
                  "in_fiducial",
                  "primary_type",
                  "target_type",
                  "num_secondaries",
                  "secondary_types",
                  "primary_momentum",
                  "secondary_momenta",
                  "parent_idx"]:
            datasets[k].append([])
        # parent index of each interaction, taken from the tree's parent edges
        parent_indices = get_parent_indices(event.tree)
        # loop over interactions
        for id, datum in enumerate(event.tree):
            datasets["vertex"][-1].append(np.array(datum.record.interaction_vertex,dtype=float))

             # primary particle stuff
            datasets["primary_type"][-1].append(int(datum.record.signature.primary_type))
            datasets["primary_momentum"][-1].append(np.array(datum.record.primary_momentum, dtype=float))

            # parent interaction index (from the tree's parent/daughter edges)
            datasets["parent_idx"][-1].append(parent_indices[id])

            if fid_vol is not None:
                pos = _math.Vector3D(datasets["vertex"][-1][-1])
                dir = _math.Vector3D(datasets["primary_momentum"][-1][-1][1:])
                dir.normalize()
                datasets["in_fiducial"][-1].append(fid_vol.IsInside(pos,dir))
            else:
                datasets["in_fiducial"][-1].append(False)

            # target particle stuff
            datasets["target_type"][-1].append(int(datum.record.signature.target_type))

            # secondary particle stuff
            datasets["secondary_types"][-1].append([])
            datasets["secondary_momenta"][-1].append([])
            for isec, (sec_type, sec_momenta) in enumerate(zip(datum.record.signature.secondary_types,
                                                               datum.record.secondary_momenta)):
                datasets["secondary_types"][-1][-1].append(int(sec_type))
                datasets["secondary_momenta"][-1][-1].append(np.array(sec_momenta,dtype=float))
            datasets["num_secondaries"][-1].append(isec+1)
        datasets["num_interactions"].append(id+1)

    # save events
    ak_array = ak.Array(datasets)
    if save_hdf5:
        fout = h5py.File(output_filename+".hdf5", "w")
        group = fout.create_group("Events")
        form, length, container = ak.to_buffers(ak.to_packed(ak_array), container=group)
        group.attrs["form"] = form.to_json()
        group.attrs["length"] = length
        # Protons-on-target normalization as run metadata, so a saved dataset
        # records the exposure it represents.
        if pot is not None:
            group.attrs["pot"] = float(pot)
        fout.close()
    if save_parquet:
        ak.to_parquet(ak_array, output_filename+".parquet")

# Load events from the custom SIREN event format
def LoadEvents(filename):
    return _dataclasses.LoadInteractionTrees(filename)


# Load events from a HepMC3 file written by SIREN
def LoadEventsFromHepMC3(filename, strict=True):
    from . import hepmc3 as _hepmc3
    return _hepmc3.LoadInteractionTreesFromHepMC3(filename, strict)


def convert_siren_events_to_hepmc3(in_path, out_path=None, weighter=None, options=None):
    """Convert a SIREN ``.siren_events`` file to a HepMC3 ascii file.

    Parameters
    ----------
    in_path : str
        Path to a ``.siren_events`` file (as written by ``SaveEvents``), or the
        base path without the suffix; both are accepted.
    out_path : str, optional
        Output path. Defaults to ``in_path`` with a trailing ``.siren_events``
        replaced by ``.hepmc3`` (or ``.hepmc3`` appended otherwise).
    weighter : callable or Weighter, optional
        When given, each tree's central-value weight (``tree.header.weights[0]``)
        is overwritten with the computed weight so the HepMC3 CV weight is set
        (weights_state ``"computed"``); a callable is invoked as ``weighter(tree)``,
        otherwise ``weighter.EventWeight(tree)`` is used. When ``None`` the existing
        header weights are kept (``"header"`` if any tree carries trusted weights,
        else ``"unweighted"``); a tree whose header weights are themselves a
        placeholder round-tripped from a prior "unweighted" HepMC3 export is not
        trusted (see ``_trees_carry_trusted_header_weights``).
    options : siren.hepmc3.HepMC3WriterOptions, optional
        Passed through to the writer (e.g. to request gzip or set FATX metadata).
        Its ``weights_state`` is set to reflect the resolved policy above.

    Returns
    -------
    str
        The output path written.

    Notes
    -----
    This function IS the weight-patch path. To update the central-value weights
    of a written HepMC3 file, re-run this function on the original
    ``.siren_events`` file with a new (or corrected) ``weighter``; the fresh CVs
    overwrite the tree headers and a new HepMC3 file is written. A HepMC3 file is
    NEVER a valid input for patching: the HepMC3 reader is run-level-lossy by
    design (the run-level injection/normalization metadata a Weighter needs to
    recompute weights is not round-tripped), so weights must always be recomputed
    from the native ``.siren_events`` + ``.siren_weighter`` pair, not from a
    HepMC3 round trip.
    """
    from . import hepmc3 as _hepmc3
    # LoadInteractionTrees takes a base name and appends ".siren_events" itself
    # (mirroring SaveInteractionTrees); strip a trailing suffix from the
    # documented ".siren_events" input so it is not doubled.
    base = _strip_known_suffix(in_path, ".siren_events")
    trees = _dataclasses.LoadInteractionTrees(base)
    if out_path is None:
        out_path = base + ".hepmc3"
    if options is None:
        options = _hepmc3.HepMC3WriterOptions()
        touched_fields = {}
    else:
        # HepMC3WriterOptions exposes no copy constructor to pybind, so it
        # cannot be copy.copy()'d; save the field this function overwrites
        # (weights_state) and restore it afterward instead, so a caller
        # reusing one options object across multiple calls does not see state
        # leak between them.
        touched_fields = {"weights_state": options.weights_state}
    # If gzip is requested, ensure a .gz suffix so the file self-describes (the reader
    # also sniffs the gzip magic bytes, but a .gz name is the clearer convention).
    if getattr(options, "gzip", False) and not out_path.endswith(".gz"):
        out_path = out_path + ".gz"
    if weighter is not None:
        for tree in trees:
            w = weighter(tree) if callable(weighter) else weighter.EventWeight(tree)
            _set_header_cv(tree, w)
        options.weights_state = "computed"
    elif _trees_carry_trusted_header_weights(trees):
        options.weights_state = "header"
    else:
        options.weights_state = "unweighted"
    try:
        _hepmc3.SaveInteractionTreesAsHepMC3(trees, out_path, options)
    finally:
        for name, value in touched_fields.items():
            setattr(options, name, value)
    return out_path


def _strip_known_suffix(path, suffix):
    """Return ``path`` with a trailing ``suffix`` removed if present."""
    if path is not None and path.endswith(suffix):
        return path[:-len(suffix)]
    return path


def _resolve_set_paths(entry):
    """Normalize one ``combine_and_export_hepmc3`` set entry.

    ``entry`` is either a mapping (dict) or a 2-tuple/list. A mapping recognizes
    the keys ``events`` (or ``events_path``), ``weighter`` (or ``weighter_path``),
    and the optional per-set count fields ``attempted``, ``accepted`` and
    ``events_to_inject``. A 2-tuple is ``(events_path, weighter_path)``. A bare
    string is treated as a shared base path (the naming convention
    ``SIREN_Controller.SaveEvents`` uses, where the events file is
    ``<base>.siren_events`` and the weighter file is ``<base>.siren_weighter``).

    A mapping may also carry ``weights_state``. Pooling recomputes every weight
    from the native pair, so it only makes sense for weighted sets; a set marked
    ``"unweighted"`` is rejected (see ``combine_and_export_hepmc3``), so its state
    is surfaced here.

    Returns ``(events_base, weighter_base, counts, weights_state)`` where the two
    bases are the suffix-free stems the C++ ``LoadInteractionTrees`` /
    ``_Weighter(..., stem)`` loaders expect, ``counts`` is a dict with any of the
    three count fields that were supplied (missing ones absent), and
    ``weights_state`` is the declared state or ``None`` when unspecified.
    """
    counts = {}
    weights_state = None
    if isinstance(entry, str):
        events_path = entry
        weighter_path = entry
    elif isinstance(entry, dict):
        events_path = entry.get("events", entry.get("events_path"))
        weighter_path = entry.get("weighter", entry.get("weighter_path"))
        if events_path is None:
            raise ValueError("combine set entry is missing an 'events' path")
        if weighter_path is None:
            # Fall back to the SaveEvents naming convention: the weighter shares
            # the events file's base path with the .siren_weighter suffix.
            weighter_path = _strip_known_suffix(events_path, ".siren_events")
        for key in ("attempted", "accepted", "events_to_inject"):
            if key in entry and entry[key] is not None:
                counts[key] = int(entry[key])
        weights_state = entry.get("weights_state")
    else:
        try:
            events_path, weighter_path = entry
        except (TypeError, ValueError):
            raise ValueError(
                "combine set entry must be a dict, an (events, weighter) tuple, "
                "or a base-path string; got %r" % (entry,))
    events_base = _strip_known_suffix(events_path, ".siren_events")
    weighter_base = _strip_known_suffix(weighter_path, ".siren_weighter")
    return events_base, weighter_base, counts, weights_state


def _process_interaction_signature(process):
    """Return a sorted, hashable marker of one process's interaction collection.

    Combines the collection's primary type with the (primary, target,
    secondaries) tuples of every cross section and decay it holds. Returns
    ``None`` if the signatures cannot be read (so a guard built on it degrades
    to a no-op rather than raising on an exotic process).
    """
    try:
        collection = process.interactions
        interactions = list(collection.GetCrossSections()) + list(collection.GetDecays())
        markers = [("primary_type", int(collection.GetPrimaryType()))]
        for inter in interactions:
            for sig in inter.GetPossibleSignatures():
                markers.append((
                    int(sig.primary_type),
                    int(sig.target_type),
                    tuple(int(s) for s in sig.secondary_types),
                ))
    except Exception:
        return None
    return tuple(sorted(markers, key=repr))


def _primary_process_signatures(injector):
    """Return a marker of an injector's INJECTION-side primary interaction
    signatures, for a cheap (best-effort, non-fatal) cross-set sanity check.

    This is the injection side, which is deliberately allowed to differ between
    pooled sets; it is only a fallback proxy used when the physical side cannot
    be read (see ``_physical_config_signature`` and the caller).
    """
    # Injectors deserialized from a weighter file are raw C++ Injector objects;
    # if a Python Injector wrapper slipped in, unwrap it.
    raw = getattr(injector, "_Injector__injector", injector)
    if raw is None:
        raw = injector
    try:
        process = raw.GetPrimaryProcess()
    except Exception:
        return None
    return _process_interaction_signature(process)


def _physical_config_signature(weighter):
    """Return a marker of a Weighter's PHYSICAL side: the interaction signatures
    of its primary physical process plus every secondary physical process.

    This is the config that MUST match across pooled sets (the injection side
    may differ). Returns ``None`` if any piece cannot be read -- e.g. an exotic
    Python-trampoline process, or a stale build whose ``Weighter`` pybind lacks
    the ``GetPrimaryPhysicalProcess`` accessor -- so the caller can fall back to
    the best-effort injection-side proxy instead of raising spuriously.
    """
    getter = getattr(weighter, "GetPrimaryPhysicalProcess", None)
    if getter is None:
        return None
    try:
        primary = getter()
    except Exception:
        return None
    if primary is None:
        return None
    parts = [_process_interaction_signature(primary)]
    try:
        secondaries = list(weighter.GetSecondaryPhysicalProcesses())
    except Exception:
        secondaries = []
    for sec in secondaries:
        parts.append(_process_interaction_signature(sec))
    if any(part is None for part in parts):
        return None
    return tuple(parts)


def combine_and_export_hepmc3(sets, out_path, options=None):
    """Pool multiple SIREN simulation sets and export ONE HepMC3 file.

    Combining simulation sets is not a matter of averaging stored per-set
    scalars: the C++ Weighter is natively pooled, so
    ``EventWeight(tree) = 1 / sum_i(N_i * G_i / P_i)`` over its injector list, and
    the per-injector physical probability ``P_i`` depends on injector ``i``'s
    injection bounds. The only correct combination is to build ONE Weighter whose
    injector list is the UNION of every set's DISTINCT injectors and re-run
    ``EventWeight`` on every event of every set. There is no shortcut through the
    per-set weights stored in the individual files. Sets whose weighter files
    resolve to the same path (a run split across multiple saves, or a set
    listed twice) contribute their injectors only once to the union, since
    including the same injector twice would double-count its term in the
    pooled denominator.

    Parameters
    ----------
    sets : list
        One entry per simulation set. Each entry provides a ``.siren_events``
        path and a ``.siren_weighter`` path, as either a dict (keys ``events`` /
        ``weighter``, plus optional per-set count fields ``attempted``,
        ``accepted``, ``events_to_inject``), an ``(events, weighter)`` tuple, or a
        single base-path string (the ``SIREN_Controller.SaveEvents`` naming
        convention where both files share a base path). This tooling operates on
        the native ``.siren_events`` + ``.siren_weighter`` pair ONLY; never feed
        it HepMC3 files, whose reader is run-level-lossy by design. A dict entry
        may also carry ``weights_state``; a set marked ``"unweighted"`` is
        rejected -- unweighted and computed sets must not be silently pooled.
    out_path : str
        Output HepMC3 path (``.hepmc3`` appended if absent; ``.gz`` appended when
        the passed ``options`` request gzip).
    options : siren.hepmc3.HepMC3WriterOptions, optional
        Passed through to the writer (a shallow copy is mutated; the caller's
        object is left untouched). Its ``weights_state`` is forced to
        ``"computed"`` (the weights written are the freshly pooled CVs). A
        pooled attempted / accepted / events_to_inject count is written only
        when EVERY set supplies that field -- a partial sum would understate
        the FATX normalization, which the writer always computes over ALL
        pooled events; a field missing from any set is omitted (with a
        warning) rather than partially summed.

    Returns
    -------
    (str, list of float)
        The output path written, and the flat list of pooled central-value
        weights in the order the events were exported (all of set 0, then set 1,
        ...), so callers can reuse them without re-reading the file.

    Notes
    -----
    Every set must share the same physical configuration (detector + physical
    processes); only the injection (generation) side may differ between sets.
    The pooled Weighter is built from the FIRST set's ``.siren_weighter`` file
    (which supplies the shared detector + physical processes) with the union of
    all sets' injectors passed live to override the file's injectors. Each set's
    physical-side interaction signatures (its primary and secondary physical
    processes, read via ``Weighter.GetPrimaryPhysicalProcess`` /
    ``GetSecondaryPhysicalProcesses``) are compared against the first set's, and
    a genuine mismatch RAISES ``ValueError`` -- pooling different physics would
    silently produce meaningless weights. The injection side is deliberately not
    checked (it is allowed to vary); it is only used as a best-effort warning
    fallback when the physical signature cannot be read.
    """
    from . import injection as _injection
    from . import hepmc3 as _hepmc3

    if not sets:
        raise ValueError("combine_and_export_hepmc3 requires at least one set")

    # a. Load each set's trees and Weighter.
    loaded_trees = []
    loaded_weighters = []
    per_set_counts = []
    first_weighter_base = None
    reference_phys_sig = None
    reference_sig = None
    all_injectors = []
    # Dedup injectors by the resolved absolute path of the .siren_weighter file
    # they were loaded from: two set entries pointing at the same weighter file
    # (a run split across multiple SaveEvents calls, or a set listed twice)
    # carry the SAME underlying injector(s). Pooling EventWeight = 1 /
    # sum_i(N_i G_i / P_i) would double-count that injector's term if its
    # weighter file contributed to all_injectors more than once, silently
    # halving every exported CV. There is no injector identity/equality check
    # exposed to Python, so path identity is the only reliable signal available
    # here.
    seen_weighter_paths = set()
    seen_events_paths = set()
    for i_set, entry in enumerate(sets):
        events_base, weighter_base, counts, weights_state = _resolve_set_paths(entry)
        # Refuse to silently pool an unweighted set: pooling recomputes weights
        # from the native pair, so a set that declares itself unweighted (no
        # meaningful weights) must not be mixed with weighted/computed sets.
        if weights_state is not None and str(weights_state).lower() == "unweighted":
            raise ValueError(
                "combine_and_export_hepmc3: set %d declares weights_state "
                "'unweighted'; unweighted and computed sets must not be pooled. "
                "Only weighted sets (with a usable .siren_weighter) can be "
                "combined." % i_set)
        # A usable weighter file is mandatory: without it there is no way to
        # recompute a pooled weight, and combining would produce garbage.
        weighter_path = weighter_base + ".siren_weighter"
        if not os.path.isfile(weighter_path):
            raise ValueError(
                "combine_and_export_hepmc3: set %d has no weighter file at "
                "'%s.siren_weighter'; pooling requires each set's weighter."
                % (i_set, weighter_base))
        # The events file is loaded further below via LoadInteractionTrees,
        # which fails deep inside cereal on a mistyped path; check it here so
        # the caller gets a clear message naming the path instead.
        events_path = events_base + ".siren_events"
        if not os.path.isfile(events_path):
            raise ValueError(
                "combine_and_export_hepmc3: set %d has no events file at "
                "'%s.siren_events'." % (i_set, events_base))
        # Skip a set whose events file was already pooled: the same events file
        # listed twice would otherwise export every event (and its FATX weight)
        # twice. Distinct event files that share one injector are the intended
        # split-run case and are NOT skipped here -- only their shared injector is
        # deduplicated below.
        resolved_events_path = os.path.abspath(events_path)
        if resolved_events_path in seen_events_paths:
            logger.warning(
                "combine_and_export_hepmc3: set %d's events file '%s' was already "
                "pooled by an earlier set; skipping the duplicate." % (i_set, events_path))
            continue
        seen_events_paths.add(resolved_events_path)
        if first_weighter_base is None:
            first_weighter_base = weighter_base

        trees = _dataclasses.LoadInteractionTrees(events_base)
        loaded_trees.append(trees)
        per_set_counts.append(counts)

        # Load this set's serialized weighter (no live injectors -> use the
        # file's own injectors, detector and physical processes).
        set_weighter = _injection._Weighter([], weighter_base)
        loaded_weighters.append(set_weighter)

        set_injectors = list(set_weighter.GetInjectors())
        # b-guard: the PHYSICAL config (detector + physical processes) is what
        # must match across sets; the injection side may legitimately differ
        # (that is the whole point of pooling distinct injectors). The Weighter
        # exposes its primary and secondary physical processes, so we build a
        # real physical-side signature and HARD-ERROR on a genuine mismatch --
        # pooling two sets with different physics would silently produce
        # meaningless weights (and would trip the Weighter's own
        # process/injector head-match assertion).
        phys_sig = _physical_config_signature(set_weighter)
        if phys_sig is not None:
            if reference_phys_sig is None:
                reference_phys_sig = phys_sig
            elif phys_sig != reference_phys_sig:
                raise ValueError(
                    "combine_and_export_hepmc3: set %d's physical configuration "
                    "(primary/secondary physical process interaction "
                    "signatures) differs from the first set's; pooled sets must "
                    "share the same physical process. Only the injection "
                    "(generation) side may differ between sets." % i_set)
        else:
            # Physical side unreadable (exotic Python-trampoline process, or a
            # stale build without the physical-process accessor). Fall back to
            # the best-effort injection-side proxy: it cannot validate the
            # physical config, so a difference is only a warning.
            for inj in set_injectors:
                sig = _primary_process_signatures(inj)
                if sig is None:
                    continue
                if reference_sig is None:
                    reference_sig = sig
                elif sig != reference_sig:
                    logger.warning(
                        "combine_and_export_hepmc3: set %d's injector has an "
                        "injection-side interaction signature that differs from "
                        "the first set's, and the physical config could not be "
                        "read to verify compatibility; ensure all pooled sets "
                        "share the same physical process." % i_set)

        # Dedup: only add this set's injectors to the pooled union the first
        # time its weighter file (by resolved absolute path) is seen. A repeat
        # (same set listed twice, or a run split across multiple SaveEvents
        # calls that all point at the same .siren_weighter) would otherwise
        # double-count that injector's term in the pooled denominator.
        resolved_weighter_path = os.path.abspath(weighter_path)
        if resolved_weighter_path in seen_weighter_paths:
            logger.warning(
                "combine_and_export_hepmc3: set %d's weighter file '%s' was "
                "already pooled by an earlier set; skipping its injectors to "
                "avoid double-counting them in the pooled weight."
                % (i_set, weighter_path))
        else:
            seen_weighter_paths.add(resolved_weighter_path)
            all_injectors.extend(set_injectors)

    # b. Pool: build ONE Weighter over the union of injectors. Passing the live
    # union overrides the file's injectors; the file supplies detector +
    # physical processes.
    pooled = _injection._Weighter(all_injectors, first_weighter_base)

    # c. Recompute every event's weight against the pooled Weighter and write it
    # into the CV slot (header.weights[0]) via the list-copy pattern.
    pooled_weights = []
    flat_trees = []
    for trees in loaded_trees:
        for tree in trees:
            w = pooled.EventWeight(tree)
            _set_header_cv(tree, w)
            pooled_weights.append(float(w))
            flat_trees.append(tree)

    # d. Export ONE HepMC3 file. Weights are the freshly pooled CVs.
    if options is None:
        options = _hepmc3.HepMC3WriterOptions()
        touched_fields = {}
    else:
        # HepMC3WriterOptions exposes no copy constructor to pybind, so it
        # cannot be copy.copy()'d; save the fields this function overwrites
        # and restore them afterward instead, so a caller reusing one options
        # object across multiple calls does not see state leak between them.
        touched_fields = {"weights_state": options.weights_state}
    options.weights_state = "computed"
    # Sum a per-set count field only when EVERY set supplied it: the FATX
    # weight_sum the writer pools always covers ALL sets, so a sum over a
    # strict subset would silently understate the normalization denominator.
    # A negative/omitted count means "not provided" to the writer, so omitting
    # it here (rather than writing a partial sum) is the safe fallback -- the
    # writer's own normalization degrades sensibly (or the caller can supply
    # the true combined count directly via `options`).
    for field, opt_name in (("attempted", "attempted_events"),
                            ("accepted", "accepted_events"),
                            ("events_to_inject", "events_to_inject")):
        provided = [c[field] for c in per_set_counts if field in c]
        if len(provided) == len(per_set_counts):
            if opt_name not in touched_fields:
                touched_fields[opt_name] = getattr(options, opt_name)
            setattr(options, opt_name, int(sum(provided)))
        elif provided:
            logger.warning(
                "combine_and_export_hepmc3: only %d/%d sets supplied a '%s' "
                "count; omitting the pooled '%s' (a partial sum would "
                "understate the FATX normalization, which covers all sets)."
                % (len(provided), len(per_set_counts), field, opt_name))

    if not out_path.endswith(".hepmc3") and not out_path.endswith(".hepmc3.gz"):
        out_path = out_path + ".hepmc3"
    if getattr(options, "gzip", False) and not out_path.endswith(".gz"):
        out_path = out_path + ".gz"
    try:
        _hepmc3.SaveInteractionTreesAsHepMC3(flat_trees, out_path, options)
    finally:
        for name, value in touched_fields.items():
            setattr(options, name, value)

    return out_path, pooled_weights

