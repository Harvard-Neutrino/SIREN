"""Smoke tests for the top-level siren package surface."""
import importlib
import subprocess
import sys
import textwrap

import pytest


def test_top_level_import():
    siren = importlib.import_module("siren")
    assert hasattr(siren, "__version__")
    assert isinstance(siren.__version__, str)


@pytest.mark.parametrize(
    "submodule",
    [
        "utilities",
        "math",
        "dataclasses",
        "geometry",
        "detector",
        "interactions",
        "distributions",
        "injection",
        "resources",
        "_util",
    ],
)
def test_submodule_imports(submodule):
    importlib.import_module(f"siren.{submodule}")


@pytest.mark.parametrize(
    "name",
    ["load_flux", "load_detector", "load_processes",
     "get_flux_model_path", "get_detector_model_path", "get_processes_model_path",
     "get_resource_package_dir", "get_fiducial_volume"],
)
def test_utilities_public_helpers(name):
    from siren import utilities
    assert hasattr(utilities, name), f"siren.utilities is missing public helper {name!r}"
    assert callable(getattr(utilities, name))


def test_resources_public_helpers():
    from siren import resources
    for name in ("load_flux", "load_detector", "load_processes"):
        assert callable(getattr(resources, name)), f"siren.resources.{name} not callable"
    for name in ("fluxes", "detectors", "processes"):
        assert hasattr(resources, name), f"siren.resources missing {name}"


# ======================================================================
# Import must not make the process multi-threaded
# ======================================================================

# Modules that start native worker threads on import. Pulling any of them in
# from siren's import path makes every SIREN process multi-threaded before user
# code runs, and fork() is only safe from a single-threaded process: fork
# duplicates the calling thread but inherits the whole memory image, so a lock
# another thread happened to hold is inherited locked and can never be
# released. That turns a multiprocessing Pool using the "fork" start method
# into an intermittent, silent deadlock. numexpr (16 threads, via
# DarkNews -> pandas) and awkward (2 threads) can both cause this issue.
_THREAD_STARTING_MODULES = ("numexpr", "pandas", "awkward")

def _run_probe(body):
    """Run *body* in a clean interpreter that has imported siren."""
    out = subprocess.run([sys.executable, "-c", textwrap.dedent(body)],
                         capture_output=True, text=True, timeout=300)
    assert out.returncode == 0, f"probe failed:\n{out.stdout}\n{out.stderr}"
    return out.stdout.strip()


def test_import_does_not_pull_in_thread_starting_modules():
    """siren must stay importable without starting background threads.

    These are lazily imported by the code that needs them (awkward by
    _util.SaveEvents, DarkNews by SIREN_DarkNews / DNModelContainer). Moving
    any of them back to module scope reintroduces the fork hazard.
    """
    loaded = _run_probe(f"""
        import sys
        import siren
        print(",".join(m for m in {_THREAD_STARTING_MODULES!r} if m in sys.modules))
    """)
    assert not loaded, (
        f"'import siren' pulled in {loaded}, which start native worker threads "
        "and make fork() unsafe; import them lazily where they are used")


@pytest.mark.skipif(not hasattr(__import__("os"), "fork"),
                    reason="requires POSIX fork()")
def test_import_leaves_process_forkable():
    """Forking after 'import siren' must not trip CPython's own warning.

    CPython >= 3.12 warns when fork() is called from a multi-threaded process.
    It counts OS-level threads, so unlike threading.active_count() this also
    catches native pools started by extension modules -- which is exactly the
    case that bit us.
    """
    result = _run_probe("""
        import os, warnings
        import siren
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            pid = os.fork()
            if pid == 0:
                os._exit(0)
            os.waitpid(pid, 0)
        print(";".join(str(w.message) for w in caught
                       if "multi-threaded" in str(w.message)))
    """)
    assert not result, (
        f"fork() after 'import siren' is unsafe: {result}. A multi-threaded "
        "process can inherit a held lock into the child, deadlocking "
        "multiprocessing Pools that use the 'fork' start method")


def test_three_body_mode_bound():
    """siren.injection.ThreeBodyMode exposes Direct and Recursive."""
    import siren
    assert hasattr(siren.injection.ThreeBodyMode, "Direct")
    assert hasattr(siren.injection.ThreeBodyMode, "Recursive")
