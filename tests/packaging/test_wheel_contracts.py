"""Packaging regressions; run with CMake, pip, and scikit-build-core installed."""

import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace
import zipfile

import pytest


REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "wheel_smoke", REPO / "tools/wheels/test_hepmc3_wheel.py"
)
smoke = importlib.util.module_from_spec(spec)
spec.loader.exec_module(smoke)


def test_linux_maps_preserves_spaces_and_escaped_newline():
    maps = (
        "1000-2000 r-xp 00000000 08:01 123 /tmp/wheel env/siren.libs/libSIREN.so\n"
        r"2000-3000 r--p 00001000 08:01 124 /tmp/new\012line/libSIREN.so"
        "\n3000-4000 rw-p 00000000 00:00 0 [heap]\n"
        "4000-5000 rw-p 00000000 00:00 0\n"
    )
    assert smoke.linux_mapped_paths(maps) == [
        Path("/tmp/wheel env/siren.libs/libSIREN.so"),
        Path("/tmp/new\nline/libSIREN.so"),
    ]


def test_deleted_core_cannot_match_a_current_wheel_file():
    paths = smoke.linux_mapped_paths(
        "1000-2000 r-xp 00000000 08:01 123 /tmp/libSIREN.so (deleted)"
    )
    assert smoke.is_core_library(paths[0])
    assert paths[0] != Path("/tmp/libSIREN.so")


@pytest.mark.parametrize(
    "name",
    [
        "libSIREN.dylib",
        "libSIREN_python-abcd.so",
        "libSIREN.so.1",
        "SIREN_python-abcd.dll",
    ],
)
def test_repaired_core_names(name):
    assert smoke.is_core_library(Path(name))


def test_rejects_second_packaged_core(monkeypatch):
    wheel = SimpleNamespace(
        files=[Path("siren.libs/libSIREN.dylib"), Path("siren/.dylibs/libSIREN.dylib")],
        read_text=lambda name: (
            "Root-Is-Purelib: false\nTag: " + str(next(smoke.sys_tags())) + "\n"
        ),
    )
    monkeypatch.setattr(smoke.metadata, "distribution", lambda name: wheel)
    with pytest.raises(AssertionError, match="Expected one core library"):
        smoke.main()


@pytest.mark.skipif(
    sys.platform not in ("darwin", "linux"), reason="POSIX native fixture"
)
def test_wheel_tags_relocation_and_incremental_contents(tmp_path):
    pytest.importorskip("scikit_build_core")
    if not shutil.which("cmake"):
        pytest.skip("CMake is required")
    source = tmp_path / "native package"
    build = tmp_path / "native build"
    source.mkdir()
    for name in ("cmake", "package", "python", "resources"):
        (source / name).mkdir()
    for name in (
        "pyproject.toml",
        "README.md",
        "LICENSE",
        "cmake/siren_python_package.cmake",
        "cmake/siren_wheel_install.cmake",
        "package/CMakeLists.txt",
    ):
        shutil.copy2(REPO / name, source / name)
    (source / "python/__init__.py").write_text("")
    (source / "resources/probe.txt").write_text("resource\n")
    (source / "spline.cpp").write_text('extern "C" int spline_probe() { return 40; }\n')
    (source / "core.cpp").write_text(
        'extern "C" int spline_probe();\n'
        'extern "C" int native_probe() { return spline_probe() + 2; }\n'
    )
    (source / "CMakeLists.txt").write_text("""cmake_minimum_required(VERSION 3.20)
project(siren LANGUAGES CXX)
set(SIREN_PYTHON_PACKAGE ON)
if(APPLE)
    set(SIREN_RPATH_ORIGIN "@loader_path")
else()
    set(SIREN_RPATH_ORIGIN "$ORIGIN")
endif()
add_library(photospline SHARED spline.cpp)
set_target_properties(photospline PROPERTIES VERSION 2.4.1 SOVERSION 2)
add_library(SIREN_python SHARED core.cpp)
target_link_libraries(SIREN_python PRIVATE photospline)
set(SIREN_WHEEL_LIBRARIES SIREN_python photospline)
include(cmake/siren_wheel_install.cmake)
siren_install_wheel_libraries(${SIREN_WHEEL_LIBRARIES})
install(DIRECTORY python/ DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL)
install(DIRECTORY resources DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL)
include(cmake/siren_python_package.cmake)
""")
    env = {
        key: value
        for key, value in os.environ.items()
        if key
        not in (
            "MACOSX_DEPLOYMENT_TARGET",
            "ARCHFLAGS",
            "CMAKE_ARGS",
            "SKBUILD_CMAKE_ARGS",
            "PYTHONPATH",
            "DYLD_LIBRARY_PATH",
            "DYLD_FALLBACK_LIBRARY_PATH",
            "LD_LIBRARY_PATH",
        )
    }
    env["PIP_NO_INDEX"] = "1"
    commands = []

    def run(command, run_env=None):
        result = subprocess.run(
            command, env=run_env or env, text=True, capture_output=True, cwd=tmp_path
        )
        commands.append(command)
        (tmp_path / f"command-{len(commands)}.log").write_text(
            repr(command) + "\n" + result.stdout + result.stderr
        )
        assert result.returncode == 0, result.stdout + result.stderr
        return result

    configure = [
        "cmake",
        "-S",
        str(source),
        "-B",
        str(build),
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DPython_EXECUTABLE={sys.executable}",
        "-DSIREN_WHEEL_BUILD_ISOLATION=OFF",
    ]
    if sys.platform == "darwin":
        configure += [
            "-DCMAKE_CXX_COMPILER=/usr/bin/c++",
            "-DCMAKE_OSX_DEPLOYMENT_TARGET=11.0",
            "-DCMAKE_OSX_ARCHITECTURES=arm64;x86_64",
        ]
    run(configure)
    build_command = [
        "cmake",
        "--build",
        str(build),
        "--target",
        "python_package",
        "-j",
        "2",
    ]
    run(build_command)
    cmake_wheel = next((build / "dist_wheels").glob("*.whl"))
    original_stamp = cmake_wheel.stat().st_mtime_ns
    run(build_command)
    assert cmake_wheel.stat().st_mtime_ns == original_stamp

    probe = source / "python/added.py"
    for content in ("first = 1\n", "changed = 2\n", None):
        # Apple's bundled Make compares timestamps at whole-second precision.
        time.sleep(1.05)
        if content is None:
            probe.unlink()
        else:
            probe.write_text(content)
        run(build_command)
        with zipfile.ZipFile(cmake_wheel) as archive:
            if content is None:
                assert "siren/added.py" not in archive.namelist()
            else:
                assert archive.read("siren/added.py").decode() == content

    source_env = dict(env)
    if sys.platform == "darwin":
        source_env["MACOSX_DEPLOYMENT_TARGET"] = "11.0"
        source_env["ARCHFLAGS"] = "-arch arm64 -arch x86_64"
    source_wheels = tmp_path / "source wheels"
    run(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            "--no-deps",
            "--no-build-isolation",
            "-w",
            str(source_wheels),
            str(source),
        ],
        source_env,
    )
    for index, wheel in enumerate((cmake_wheel, next(source_wheels.glob("*.whl")))):
        installed = tmp_path / f"installed-{index}"
        with zipfile.ZipFile(wheel) as archive:
            libraries = [
                name for name in archive.namelist() if "/libphotospline" in name
            ]
            assert len(libraries) == 1, libraries
            assert libraries[0].endswith(
                ".2.dylib" if sys.platform == "darwin" else ".so.2"
            )
            metadata = archive.read("siren-0.1.0.dist-info/WHEEL").decode()
            assert "Root-Is-Purelib: false" in metadata
            if sys.platform == "darwin":
                assert "macosx_11_0_universal2" in metadata
            archive.extractall(installed)
        if sys.platform == "darwin":
            archs = run(
                ["lipo", "-archs", str(installed / "siren.libs/libSIREN_python.dylib")]
            )
            assert set(archs.stdout.split()) == {"arm64", "x86_64"}
        hidden = []
        try:
            for path in (source, build):
                destination = path.with_name(path.name + " hidden")
                path.rename(destination)
                hidden.append((path, destination))
            core = (
                installed
                / "siren.libs"
                / (
                    "libSIREN_python.dylib"
                    if sys.platform == "darwin"
                    else "libSIREN_python.so"
                )
            )
            run(
                [
                    sys.executable,
                    "-c",
                    "import ctypes,sys; lib=ctypes.CDLL(sys.argv[1]); assert lib.native_probe()==42",
                    str(core),
                ]
            )
        finally:
            for original, destination in reversed(hidden):
                destination.rename(original)
