"""Shared pytest fixtures."""
from __future__ import annotations

from pathlib import Path

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run-network",
        action="store_true",
        default=False,
        help="run tests that require network access")


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "network: test requires network access")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-network"):
        return
    skip_network = pytest.mark.skip(reason="requires --run-network")
    for item in items:
        if "network" in item.keywords:
            item.add_marker(skip_network)


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
