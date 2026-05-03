"""Shared pytest fixtures."""
from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def project_root() -> Path:
    here = Path(__file__).resolve()
    for parent in (here, *here.parents):
        if (parent / "pyproject.toml").exists() and (parent / "resources").exists():
            return parent
    raise RuntimeError(f"could not locate project root from {here}")


@pytest.fixture(scope="session")
def resources_dir(project_root: Path) -> Path:
    return project_root / "resources"


@pytest.fixture(scope="session")
def detectors_dir(resources_dir: Path) -> Path:
    return resources_dir / "detectors"


@pytest.fixture(scope="session")
def fluxes_dir(resources_dir: Path) -> Path:
    return resources_dir / "fluxes"


@pytest.fixture(scope="session")
def processes_dir(resources_dir: Path) -> Path:
    return resources_dir / "processes"
