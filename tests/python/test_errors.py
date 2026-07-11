"""Typed-exception surface + cross-module translation.

These tests pin the contract that engine configuration failures raise typed
exceptions, that all typed exceptions are importable from siren.utilities and
derive from RuntimeError (existing pytest.raises(RuntimeError) call sites must
keep matching), and that the exception translators -- registered in the
utilities pybind module -- work even when a raise originates in a sibling
module.
"""
import subprocess
import sys
import textwrap

import pytest

siren = pytest.importorskip("siren")


# The full set of typed exceptions the contract requires on siren.utilities.
EXPECTED_EXCEPTIONS = [
    "ConfigurationError",
    "MeasureCompatibilityError",
    "WeightCalculationError",
    "AddProcessFailure",
    "SecondaryProcessFailure",
    "InjectionFailure",
    "PythonImplementationError",
]


# ------------------------------------------------------------------ #
#  Surface: each exception exists and is a RuntimeError subclass       #
# ------------------------------------------------------------------ #

@pytest.mark.parametrize("name", EXPECTED_EXCEPTIONS)
def test_exception_exists_on_utilities(name):
    assert hasattr(siren.utilities, name), (
        "siren.utilities is missing typed exception {!r}".format(name))


@pytest.mark.parametrize("name", EXPECTED_EXCEPTIONS)
def test_exception_is_runtimeerror_subclass(name):
    exc = getattr(siren.utilities, name)
    assert isinstance(exc, type)
    assert issubclass(exc, RuntimeError), (
        "{} must derive from RuntimeError so existing "
        "pytest.raises(RuntimeError) pins keep matching".format(name))


# ------------------------------------------------------------------ #
#  Cross-module registration: a fresh interpreter that imports the raw #
#  compiled sibling modules (geometry, utilities, injection) directly  #
#  -- NOT through the siren package __init__ -- must still translate    #
#  typed exceptions.  Registration lives in the utilities module, so    #
#  a raise from the injection sibling still crosses as the typed type.   #
# ------------------------------------------------------------------ #

_SUBPROCESS_SRC = textwrap.dedent(
    """
    import importlib.machinery
    import importlib.util
    import os
    import sys

    pkgdir = sys.argv[1]

    import glob

    def load(name):
        # Load the compiled module by file path, bypassing the siren package
        # __init__ (which would import the whole stack).  This models a consumer
        # that only pulls in a subset of sibling modules.  The module MUST be
        # loaded under its bare name so pybind's PyInit_<name> is found, and
        # registered in sys.modules so later modules see it.  Glob every
        # platform extension suffix (.so, .pyd, .abi3.so, ...) so this is
        # portable across macOS/Linux/Windows rather than assuming ".so".
        matches = []
        for suffix in importlib.machinery.EXTENSION_SUFFIXES:
            matches += glob.glob(os.path.join(pkgdir, name + "*" + suffix))
        assert matches, "no compiled module for " + name
        spec = importlib.util.spec_from_file_location(name, matches[0])
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        return mod

    # Import utilities FIRST (this is where the translators register), then a
    # couple of sibling modules.  The order below deliberately loads geometry
    # before injection to prove the registration is module-independent.
    utilities = load("utilities")
    geometry = load("geometry")
    math_mod = load("math")
    dataclasses = load("dataclasses")
    detector = load("detector")
    interactions = load("interactions")
    distributions = load("distributions")
    injection = load("injection")

    # Every typed exception must be present on the raw utilities module and be
    # a RuntimeError subclass.
    for nm in {names!r}:
        exc = getattr(utilities, nm, None)
        assert exc is not None, "utilities missing " + nm
        assert issubclass(exc, RuntimeError), nm + " not a RuntimeError"

    # A C++ raise from the injection sibling must cross as the typed exception,
    # even though we never imported the siren package facade.
    iso = injection.Isotropic2BodyChannel(0)
    raised = None
    try:
        injection.MultiChannelPhaseSpace([iso, iso], [1.0])
    except utilities.ConfigurationError as e:
        raised = "typed"
    except RuntimeError as e:
        raised = "runtimeerror-only"
    assert raised == "typed", (
        "cross-module raise did not translate to the typed exception; got: "
        + str(raised))
    print("OK")
    """
)


def test_cross_module_registration_subprocess():
    """In a fresh interpreter importing only the raw sibling .so modules (not
    the siren package __init__), a C++ raise from the injection module still
    translates to the typed ConfigurationError -- proving the translators,
    registered in the utilities module, are visible cross-module.
    """
    import os

    pkgdir = os.path.dirname(siren.__file__)
    src = _SUBPROCESS_SRC.format(names=EXPECTED_EXCEPTIONS)
    proc = subprocess.run(
        [sys.executable, "-c", src, pkgdir],
        capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, (
        "subprocess failed:\nSTDOUT:\n{}\nSTDERR:\n{}".format(
            proc.stdout, proc.stderr))
    assert "OK" in proc.stdout
