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
spec = importlib.util.spec_from_file_location("native_loader", REPO / "python/_native.py")
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)


def test_linux_maps_preserves_spaces_and_escaped_newline():
    maps = (
        "1000-2000 r-xp 00000000 08:01 123 /tmp/wheel env/siren.libs/libSIREN.so\n"
        r"2000-3000 r--p 00001000 08:01 124 /tmp/new\012line/libSIREN.so"
        "\n3000-4000 rw-p 00000000 00:00 0 [heap]\n"
        "4000-5000 rw-p 00000000 00:00 0\n"
    )
    assert native.linux_mapped_paths(maps) == [
        Path("/tmp/wheel env/siren.libs/libSIREN.so"),
        Path("/tmp/new\nline/libSIREN.so"),
    ]


def test_deleted_core_cannot_match_a_current_wheel_file():
    paths = native.linux_mapped_paths(
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


@pytest.mark.parametrize("name", ["libphotospline.2.4.1.dylib", "libphotospline.so.2.4.1"])
def test_rejects_foreign_versioned_dependency(tmp_path, name):
    soname = "libphotospline.2.dylib" if name.endswith("dylib") else "libphotospline.so.2"
    relative = Path("siren.libs") / soname
    wheel = SimpleNamespace(files=[relative], locate_file=lambda path: tmp_path / path)
    smoke.verify_bundled_libraries(wheel, [tmp_path / relative])
    with pytest.raises(AssertionError, match="outside this wheel"):
        smoke.verify_bundled_libraries(wheel, [tmp_path / "native" / name])


def test_stripped_loader_override_fails(monkeypatch):
    monkeypatch.delenv("DYLD_LIBRARY_PATH", raising=False)
    with pytest.raises(AssertionError, match="did not reach"):
        smoke.check_loader_override("DYLD_LIBRARY_PATH", "/native/lib")


def test_extension_names_can_match_the_standard_library(tmp_path):
    relative = Path("siren/math.cpython-313-darwin.so")
    wheel = SimpleNamespace(files=[relative], locate_file=lambda path: tmp_path / path)
    smoke.verify_bundled_libraries(wheel, [tmp_path / "lib-dynload" / relative.name])


def test_distinct_dependency_abis_are_not_aliases():
    assert smoke.library_identity(Path("libHepMC3.4.dylib")) != smoke.library_identity(
        Path("libHepMC3.5.dylib"))


def test_other_wheel_copy_is_not_a_substitution(tmp_path, monkeypatch):
    own = Path("siren/.dylibs/libgfortran.5.dylib")
    scipy = Path("scipy/.dylibs/libgfortran.5.dylib")
    wheel = SimpleNamespace(files=[own], locate_file=lambda path: tmp_path / path)
    other = SimpleNamespace(files=[scipy], locate_file=wheel.locate_file,
                            metadata={"Name": "scipy"})
    monkeypatch.setattr(smoke.metadata, "distributions", lambda: [other])
    smoke.verify_bundled_libraries(wheel, [tmp_path / own, tmp_path / scipy])
    # Owning a file is not sufficient: SIREN's own copy must also be loaded.
    with pytest.raises(AssertionError, match="outside this wheel"):
        smoke.verify_bundled_libraries(wheel, [tmp_path / scipy])
    with pytest.raises(AssertionError, match="outside this wheel"):
        smoke.verify_bundled_libraries(wheel, [tmp_path / own, tmp_path / "libgfortran.5.dylib"])
    other.metadata = {"Name": "siren"}
    with pytest.raises(AssertionError, match="outside this wheel"):
        smoke.verify_bundled_libraries(wheel, [tmp_path / own, tmp_path / scipy])


def test_repair_hash_is_part_of_loader_name(tmp_path):
    own = Path("siren.libs/libcfitsio-a1b2c3d4.so.4")
    wheel = SimpleNamespace(files=[own], locate_file=lambda path: tmp_path / path)
    # The loader requests the hashed name: an unhashed library cannot replace it.
    smoke.verify_bundled_libraries(wheel, [tmp_path / own, tmp_path / "native/libcfitsio.so.4"])
    with pytest.raises(AssertionError, match="outside this wheel"):
        smoke.verify_bundled_libraries(wheel, [tmp_path / "native" / own.name])


@pytest.mark.parametrize("platform", ["freebsd14", "linux"])
def test_import_guard_allows_unavailable_probe_but_acceptance_does_not(monkeypatch, platform):
    monkeypatch.setattr(native.sys, "platform", platform)
    def no_procfs(path):
        raise FileNotFoundError("procfs is not mounted")
    monkeypatch.setattr(native.Path, "read_text", no_procfs)
    native.reject_standalone_runtime()
    with pytest.raises((native.NativeInspectionUnavailable, FileNotFoundError)):
        native.loaded_libraries()


@pytest.mark.parametrize("missing", ["platform", "backend", "cmake", "repair"])
def test_required_packaging_prerequisites_cannot_skip(tmp_path, monkeypatch, missing):
    monkeypatch.setenv("SIREN_TEST_REQUIRE_WHEEL_REPAIR", "1")
    monkeypatch.setattr(sys, "platform", "freebsd14" if missing == "platform" else "linux")
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None if missing == "backend" else True)
    missing_tool = {"cmake": "cmake", "repair": "auditwheel"}.get(missing)
    monkeypatch.setattr(shutil, "which", lambda name: None if name == missing_tool else name)
    with pytest.raises(pytest.fail.Exception, match="prerequisites missing"):
        test_wheel_tags_relocation_and_incremental_contents(tmp_path, "siren.libs")


@pytest.mark.parametrize("name", ["libSIREN.dylib", "libSIREN.so", "SIREN.dll"])
def test_standalone_core_rejected_before_extensions(monkeypatch, name):
    monkeypatch.setattr(native, "loaded_libraries", lambda: [Path("/native") / name])
    with pytest.raises(ImportError, match="separate processes"):
        native.reject_standalone_runtime()


@pytest.mark.parametrize("library_dir", ["siren.libs", "siren.native"])
def test_wheel_tags_relocation_and_incremental_contents(tmp_path, library_dir):
    required = os.environ.get("SIREN_TEST_REQUIRE_WHEEL_REPAIR") == "1"
    repair = "delocate-wheel" if sys.platform == "darwin" else "auditwheel"
    missing = []
    if sys.platform not in ("darwin", "linux"):
        missing.append("POSIX native fixture platform")
    if importlib.util.find_spec("scikit_build_core") is None:
        missing.append("scikit-build-core")
    for tool in ("cmake", repair):
        if not shutil.which(tool) and (tool == "cmake" or required):
            missing.append(tool)
    if missing:
        message = "Required packaging prerequisites missing: " + ", ".join(missing)
        if required:
            pytest.fail(message)
        pytest.skip(message)
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
        "cmake/build_wheel.py",
        "cmake/wheel_rpath.py",
        "package/CMakeLists.txt",
    ):
        shutil.copy2(REPO / name, source / name)
    # Change the single directory definition: modules, staging, and the nested
    # backend must agree without editing any of their rules.
    helper = source / "cmake/siren_wheel_install.cmake"
    helper.write_text(helper.read_text().replace('"siren.libs"', f'"{library_dir}"'))
    (source / "python/__init__.py").write_text("")
    (source / "resources/probe.txt").write_text("resource\n")
    (source / "spline.cpp").write_text(
        'extern "C" int spglam_probe();\n'
        'extern "C" int spline_probe() { return spglam_probe(); }\n'
    )
    (source / "spglam.cpp").write_text(
        'extern "C" int external_probe();\n'
        'extern "C" int spglam_probe() { return external_probe(); }\n'
    )
    (source / "core.cpp").write_text(
        'extern "C" int spline_probe();\n'
        'extern "C" int native_probe() { return spline_probe() + 2; }\n'
    )
    (source / "module.cpp").write_text(
        '#include <Python.h>\nextern "C" int native_probe();\n'
        'static PyObject *probe(PyObject *, PyObject *) { return PyLong_FromLong(native_probe()); }\n'
        'static PyMethodDef methods[] = {{"probe", probe, METH_NOARGS, nullptr}, {nullptr}};\n'
        'static PyModuleDef module = {PyModuleDef_HEAD_INIT, "probe", nullptr, -1, methods};\n'
        'PyMODINIT_FUNC PyInit_probe() { return PyModule_Create(&module); }\n'
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
add_library(external SHARED IMPORTED)
set_target_properties(external PROPERTIES IMPORTED_LOCATION "@external@")
add_library(spglam SHARED spglam.cpp)
set_target_properties(spglam PROPERTIES VERSION 2.4.1 SOVERSION 2)
target_link_libraries(spglam PRIVATE external)
target_link_libraries(photospline PRIVATE spglam external)
get_filename_component(external_directory "@external@" DIRECTORY)
set_target_properties(spglam PROPERTIES
    INSTALL_RPATH "${external_directory}" INSTALL_RPATH_USE_LINK_PATH TRUE)
set_target_properties(photospline PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS photospline spglam LIBRARY DESTINATION lib)
add_library(SIREN_python SHARED core.cpp)
target_link_libraries(SIREN_python PRIVATE photospline)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
Python_add_library(probe MODULE WITH_SOABI module.cpp)
target_link_libraries(probe PRIVATE SIREN_python)
set(SIREN_PYTHON_MODULES probe)
set(SIREN_WHEEL_LIBRARIES SIREN_python photospline spglam)
include(cmake/siren_wheel_install.cmake)
siren_install_wheel_libraries(${SIREN_WHEEL_LIBRARIES})
siren_install_wheel_modules(${SIREN_PYTHON_MODULES})
install(DIRECTORY python/ DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL)
install(DIRECTORY resources DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL)
include(cmake/siren_python_package.cmake)
""".replace("@external@", str(tmp_path / "external install/lib" / (
        "libexternal.dylib" if sys.platform == "darwin" else "libexternal.so"))))
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

    # Model a vendor dependency with an @rpath ID outside the native prefix.
    external_source = tmp_path / "external source"
    external_source.mkdir()
    (external_source / "external.cpp").write_text('extern "C" int external_probe() { return 40; }')
    (external_source / "CMakeLists.txt").write_text(
        'cmake_minimum_required(VERSION 3.20)\nproject(external LANGUAGES CXX)\n'
        'add_library(external SHARED external.cpp)\n'
        'set_target_properties(external PROPERTIES INSTALL_NAME_DIR "@rpath")\n'
        'install(TARGETS external LIBRARY DESTINATION lib)\n'
    )
    external_build = tmp_path / "external build"
    external_configure = ["cmake", "-S", str(external_source), "-B", str(external_build),
                          f"-DCMAKE_INSTALL_PREFIX={tmp_path / 'external install'}"]
    if sys.platform == "darwin":
        external_configure += ["-DCMAKE_CXX_COMPILER=/usr/bin/c++",
                               "-DCMAKE_OSX_DEPLOYMENT_TARGET=11.0",
                               "-DCMAKE_OSX_ARCHITECTURES=arm64;x86_64"]
    run(external_configure)
    run(["cmake", "--build", str(external_build)])
    run(["cmake", "--install", str(external_build)])

    configure = [
        "cmake",
        "-S",
        str(source),
        "-B",
        str(build),
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DCMAKE_INSTALL_PREFIX={tmp_path / 'native install'}",
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
    # A native component must never invoke pip in the active interpreter.
    # Block it before it can uninstall a developer's existing SIREN package.
    pip_blocker = tmp_path / "pip blocker"
    pip_blocker.mkdir()
    (pip_blocker / "pip.py").write_text(
        'raise RuntimeError("Native component install must not invoke pip")\n')
    run(["cmake", "--install", str(build), "--component", "Unspecified"],
        dict(env, PYTHONPATH=str(pip_blocker)))
    # Model an external prefix that also contains native photospline/spglam.
    # Its link-derived RPATH must not win during wheel staging or repair.
    for binary in (tmp_path / "native install/lib").iterdir():
        if binary.is_file():
            shutil.copy2(binary, tmp_path / "external install/lib" / binary.name)
    native_spline = tmp_path / "native install/lib" / (
        "libphotospline.2.dylib" if sys.platform == "darwin" else "libphotospline.so.2")
    run([sys.executable, "-c", "import ctypes,sys; "
         "assert ctypes.CDLL(sys.argv[1]).spline_probe() == 40", str(native_spline)])
    cmake_wheel = next((build / "dist_wheels").glob("*.whl"))
    original_stamp = cmake_wheel.stat().st_mtime_ns
    unchanged = run(build_command)
    assert cmake_wheel.stat().st_mtime_ns == original_stamp, unchanged.stdout + unchanged.stderr
    time.sleep(1.05)
    run(configure)
    unchanged = run(build_command)
    assert cmake_wheel.stat().st_mtime_ns == original_stamp, unchanged.stdout + unchanged.stderr

    # An install-only change must invalidate the wheel even with identical
    # target binaries and package source files.
    time.sleep(1.05)
    with (source / "CMakeLists.txt").open("a") as output:
        output.write('\ninstall(FILES resources/probe.txt DESTINATION siren/extra '
                     'COMPONENT PythonWheel EXCLUDE_FROM_ALL)\n')
    run(configure)
    run(build_command)
    with zipfile.ZipFile(cmake_wheel) as archive:
        assert archive.read("siren/extra/probe.txt") == b"resource\n"

    # Detect an equal-size restore even when its timestamp is preserved.
    resource = source / "resources/probe.txt"
    resource_stat = resource.stat()
    resource.write_text("replaced\n")
    os.utime(resource, ns=(resource_stat.st_atime_ns, resource_stat.st_mtime_ns))
    run(build_command)
    with zipfile.ZipFile(cmake_wheel) as archive:
        assert archive.read("siren/resources/probe.txt") == b"replaced\n"

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
    wheels = [cmake_wheel, next(source_wheels.glob("*.whl"))]
    repair = "delocate-wheel" if sys.platform == "darwin" else "auditwheel"
    if shutil.which(repair):
        for index, wheel in enumerate(list(wheels)):
            repaired = tmp_path / f"repaired-{index}"
            command = [repair] + (["repair"] if repair == "auditwheel" else [])
            run(command + ["-w", str(repaired), str(wheel)])
            wheels.append(next(repaired.glob("*.whl")))
    elif os.environ.get("SIREN_TEST_REQUIRE_WHEEL_REPAIR") == "1":
        pytest.fail(f"Required repair tool missing: {repair}")
    for index, wheel in enumerate(wheels):
        installed = tmp_path / f"installed-{index}"
        with zipfile.ZipFile(wheel) as archive:
            cores = [name for name in archive.namelist() if smoke.is_core_library(Path(name))]
            assert len(cores) == 1, cores
            core_path = installed / cores[0]
            libraries = [
                name for name in archive.namelist() if "/libphotospline" in name
            ]
            assert len(libraries) == 1, libraries
            spglam = [name for name in archive.namelist() if "/libspglam" in name]
            assert len(spglam) == 1, spglam
            assert libraries[0].endswith(
                ".2.dylib" if sys.platform == "darwin" else ".so.2"
            )
            metadata = archive.read("siren-0.1.0.dist-info/WHEEL").decode()
            assert "Root-Is-Purelib: false" in metadata
            if sys.platform == "darwin":
                assert "macosx_11_0_universal2" in metadata
            archive.extractall(installed)
        if sys.platform == "darwin":
            if index < 2:
                for relative in [cores[0], libraries[0], spglam[0]]:
                    load_commands = run(["otool", "-l", str(installed / relative)]).stdout
                    rpaths = smoke.re.findall(
                        r"cmd LC_RPATH\n\s+cmdsize \d+\n\s+path (.*?) \(offset \d+\)",
                        load_commands)
                    assert rpaths[0] == "@loader_path", rpaths
            archs = run(
                ["lipo", "-archs", str(core_path)]
            )
            assert set(archs.stdout.split()) == {"arm64", "x86_64"}
        hidden = []
        try:
            hidden_paths = [source, build]
            if index >= 2:  # Repaired wheels must also bundle the external dependency.
                hidden_paths.append(tmp_path / "external install")
            for path in hidden_paths:
                destination = path.with_name(path.name + " hidden")
                path.rename(destination)
                hidden.append((path, destination))
            run(
                [
                    sys.executable,
                    "-c",
                    "import ctypes,sys; lib=ctypes.CDLL(sys.argv[1]); assert lib.native_probe()==42",
                    str(core_path),
                ]
            )
            run([sys.executable, "-c",
                 "import sys; sys.path.insert(0, sys.argv[1]); "
                 "from siren.probe import probe; assert probe() == 42", str(installed)])
        finally:
            for original, destination in reversed(hidden):
                destination.rename(original)
