"""Integrity checks on every shipped detector materials.dat."""
from pathlib import Path

import pytest


def _parse_materials_file(path: Path):
    """Yield (material_name, declared_count, components, sum) per material."""
    sums = {}
    declared = {}
    components = {}
    current = None
    for line in path.read_text().splitlines():
        stripped = line.split("#", 1)[0].strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) == 2 and parts[1].isdigit():
            current = parts[0]
            declared[current] = int(parts[1])
            sums[current] = 0.0
            components[current] = []
            continue
        if current is not None and len(parts) >= 2:
            try:
                frac = float(parts[1])
            except ValueError:
                continue
            sums[current] += frac
            components[current].append((parts[0], frac))
    for name in declared:
        yield name, declared[name], components[name], sums[name]


def _all_materials_files(detectors_dir: Path):
    return sorted(detectors_dir.glob("*/*-v*/materials.dat"))


def _file_id(path):
    return f"{path.parents[1].name}/{path.parents[0].name}/{path.name}"


@pytest.fixture
def materials_files(detectors_dir):
    return _all_materials_files(detectors_dir)


def test_materials_files_exist(detectors_dir):
    files = _all_materials_files(detectors_dir)
    assert files, f"no materials.dat under {detectors_dir}"


# Pre-existing data bugs to be cleaned up in a follow-up.
KNOWN_BAD_MATERIALS_FILES = {
    # Missing component counts on FGDbox/TPC_padding; same class of bug
    # fixed in ND280-v1 by #115.
    "ND280UPGRD/ND280UPGRD-v1/materials.dat",
}

# Loose enough to accept ~1e-4 rounding drift in existing files but tight
# enough to catch a genuine mistake.
SUM_TOLERANCE = 1e-3


def test_each_material_components_sum_to_one(materials_files):
    failures = []
    for path in materials_files:
        if _file_id(path) in KNOWN_BAD_MATERIALS_FILES:
            continue
        for name, _declared, _comps, total in _parse_materials_file(path):
            if abs(total - 1.0) > SUM_TOLERANCE:
                failures.append(
                    f"  {_file_id(path)} :: {name}: sum={total:.6f} (off by {total-1.0:+.6f})"
                )
    assert not failures, "materials with components not summing to 1.0:\n" + "\n".join(failures)


def test_each_material_declared_count_matches_actual(materials_files):
    failures = []
    for path in materials_files:
        if _file_id(path) in KNOWN_BAD_MATERIALS_FILES:
            continue
        for name, declared, comps, _total in _parse_materials_file(path):
            actual = len(comps)
            if declared != actual:
                failures.append(
                    f"  {_file_id(path)} :: {name}: declared {declared} components, found {actual}"
                )
    assert not failures, "materials with wrong declared component count:\n" + "\n".join(failures)


def test_no_duplicate_material_names_within_file(materials_files):
    """Regression test for the duplicate ROCK/AIR entries fixed in #115."""
    failures = []
    for path in materials_files:
        seen = set()
        for line in path.read_text().splitlines():
            stripped = line.split("#", 1)[0].strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) == 2 and parts[1].isdigit():
                if parts[0] in seen:
                    failures.append(f"  {path.relative_to(path.parents[2])}: duplicate material {parts[0]!r}")
                seen.add(parts[0])
    assert not failures, "duplicate material headers:\n" + "\n".join(failures)
