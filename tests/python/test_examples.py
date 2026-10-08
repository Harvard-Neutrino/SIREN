"""The maintained example scripts compile and avoid deprecated entry points.

The scripts under ``resources/examples/example{1,2,3}`` are the documented
spec-form usage of the package. This module keeps them importable-by-syntax
and free of the deprecated ``SIREN_Controller`` and ``GenerateEvents``
surfaces; the ``legacy/`` copies are excluded on purpose. It does not run
the examples, which need detector data, process tables, and (for MARLEY and
charm) optional native dependencies.
"""
from __future__ import annotations

import py_compile
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[2] / "resources" / "examples"
MAINTAINED = ("example1", "example2", "example3")
DEPRECATED = ("SIREN_Controller", "GenerateEvents")


def _maintained_scripts():
    return sorted(
        script
        for directory in MAINTAINED
        for script in (EXAMPLES / directory).glob("*.py")
    )


def _script_id(script):
    return str(script.relative_to(EXAMPLES))


@pytest.mark.parametrize("script", _maintained_scripts(), ids=_script_id)
def test_example_compiles(script, tmp_path):
    py_compile.compile(str(script), cfile=str(tmp_path / "compiled.pyc"),
                       doraise=True)


@pytest.mark.parametrize("script", _maintained_scripts(), ids=_script_id)
def test_example_avoids_deprecated_entry_points(script):
    text = script.read_text()
    offenders = [name for name in DEPRECATED if name in text]
    assert not offenders, f"{_script_id(script)} uses deprecated {offenders}"
