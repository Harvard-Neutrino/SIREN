import os
import re
import sys
import uuid
import pathlib
import importlib

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
# https://github.com/imageio/imageio/blob/65d79140018bb7c64c0692ea72cb4093e8d632a0/imageio/core/util.py
def resource_package_dir():
    """package_dir

    Get the resources directory in the siren package installation
    directory.

    Notes
    -----
    This is a convenience method that is used by `resource_dirs` and
    siren entry point scripts.
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
        pdir = pkg_resources.resource_filename("siren", "resources")
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


def load_module(name, path, persist=True):
    """Load a module with a specific name and path"""
    url = pathlib.Path(os.path.abspath(path)).as_uri()
    module_name = f"{name}-{str(uuid.uuid5(uuid.NAMESPACE_URL, url))}"
    if module_name in sys.modules:
        return module
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module_name)
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

_MODEL_PATTERN = (
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


def _get_model_path(model_name, prefix=None, suffix=None, is_file=True, must_exist=True):
    # Get the path to the model file
    _model_regex = re.compile(
        r"^\s*" + _MODEL_PATTERN + ("" if suffix is None else r"(?:" + suffix + r")?") + r"\s*$",
        re.VERBOSE | re.IGNORECASE,
    )
    if suffix is None:
        suffix = ""
    # Get the path to the resources directory
    resources_dir = resource_package_dir()
    base_dir = resources_dir

    # Add prefix if present
    if prefix is not None:
        base_dir = os.path.join(base_dir, prefix)

    # Get the model name and version
    d = _model_regex.match(model_name)
    if d is None:
        raise ValueError("Invalid model name: {}".format(model_name))
    d = d.groupdict()
    model_name = d["model_name"]
    version = d["version"]

    # Search for the model folder in the resources directory
    model_names = [
        f for f in os.listdir(base_dir) if not os.path.isfile(os.path.join(base_dir, f))
    ]
    model_names = [f for f in model_names if f.lower().startswith(model_name.lower())]

    folder_exists = False

    if len(model_names) == 0 and must_exist:
        # Whoops, we didn't find the model folder!
        raise ValueError(
            "No model folders found for {}\nSearched in ".format(model_name, base_dir)
        )
    elif len(model_names) == 0 and not must_exist:
        # Let's use the provided model name as the folder name
        model_name = model_name
    elif len(model_names) == 1:
        # We found the model folder!
        folder_exists = True
        model_name = model_names[0]
    else:
        # Multiple model folders found, we cannot decide which one to use
        raise ValueError(
            "Multiple directories found for {}\nSearched in ".format(
                model_name, base_dir
            )
        )

    if folder_exists:
        # Search for the model file in the model folder
        model_files = [
            f
            for f in os.listdir(os.path.join(base_dir, model_name))
            if is_file == os.path.isfile(os.path.join(base_dir, model_name, f))
        ]
    else:
        model_files = []

    # From the found model files, extract the model versions
    model_versions = []
    for f in model_files:
        d = _model_regex.match(f)
        if d is not None:
            if d.groupdict()["version"] is not None:
                model_versions.append(normalize_version(d.groupdict()["version"]))
            else:
                print(ValueError(
                    "Input model file has no version: {}\nSearched in ".format(
                        f, os.path.join(base_dir, model_name)
                    )
                ))
        elif f.lower().startswith(model_name.lower()):
            print(ValueError(
                "Unable to parse version from {}\nFound in ".format(
                    f, os.path.join(base_dir, model_name)
                )
            ))

    # Raise an error if no model file is found and we require it to exist
    if len(model_versions) == 0 and must_exist:
        raise ValueError(
            "No model found for {}\nSearched in ".format(
                model_name, os.path.join(base_dir, model_name)
            )
        )

    if version is None and must_exist:
        # If no version is provided, use the latest version
        version_idx, version = max(
            enumerate(model_versions), key=lambda x: tokenize_version(x[1])
        )
        model_file_name = model_files[version_idx]
    elif version is None and not must_exist:
        # If no version is provided and we don't require it to exist, default to v1
        version = "v1"
        model_file_name = "{}-v{}{}".format(model_name, version, suffix)
    else:
        # A version is provided
        version = normalize_version(version)
        if must_exist:
            # If the version must exist, raise an error if it doesn't
            if version not in model_versions:
                raise ValueError(
                    "No model found for {}-{}\nSearched in ".format(
                        model_name, version, os.path.join(base_dir, model_name)
                    )
                )
            version_idx = model_versions.index(version)
            model_file_name = model_files[version_idx]
        else:
            # The version doesn't have to exist
            if version in model_versions:
                # If the version exists, use it
                version_idx = model_versions.index(version)
                model_file_name = model_files[version_idx]
            else:
                # Otherwise use the provided version
                model_file_name = "{}-v{}{}".format(model_name, version, suffix)

    return os.path.join(base_dir, model_name, model_file_name)


def get_detector_model_file_path(model_name, must_exist=True):
    return _get_model_path(model_name, prefix="Detectors/densities", suffix=".dat", is_file=True, must_exist=must_exist)


def get_material_model_file_path(model_name, must_exist=True):
    return _get_model_path(model_name, prefix="Detectors/materials", suffix=".dat", is_file=True, must_exist=must_exist)


_resource_folder_by_name = {
    "flux": "Fluxes",
    "detector": "Detectors",
    "processes": "Processes",
}


def get_flux_model_path(model_name, must_exist=True):
    return _get_model_path(model_name, prefix=_resource_folder_by_name["flux"], is_file=False, must_exist=must_exist)


def get_detector_model_path(model_name, must_exist=True):
    return _get_model_path(model_name, prefix=_resource_folder_by_name["detector"], is_file=False, must_exist=must_exist)


def get_processes_model_path(model_name, must_exist=True):
    return _get_model_path(model_name, prefix=_resource_folder_by_name["processes"], is_file=False, must_exist=must_exist)


def load_resource(resource_name, resource_type, *args, **kwargs):
    folder = _resource_folder_by_name[resource_type]

    abs_dir = _get_model_path(model_name, prefix=folder, is_file=False, must_exist=True)

    fname = os.path.join(abs_flux_dir, f"{resource_name}.py")
    assert(os.path.isfile(fname))
    resource_module = load_module(f"siren-{resource_type}-{model_name}", fname, persist=False)
    loader = getattr(resource_module, f"load_{resource_name}")
    resource = loader(*args, **kwargs)
    return resource


def load_flux(model_name, *args, **kwargs):
    return load_resource("flux", model_name, *args, **kwargs)


def load_detector(model_name, *args, **kwargs):
    return load_resource("detector", model_name, *args, **kwargs)


def load_processes(model_name, *args, **kwargs):
    return load_resource("processes", model_name, *args, **kwargs)
