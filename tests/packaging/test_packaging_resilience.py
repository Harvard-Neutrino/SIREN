"""Failure recovery and installation boundaries for the native wheel helpers."""

import importlib.util
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from types import SimpleNamespace
import zipfile

import pytest


REPO = Path(__file__).resolve().parents[2]


def helper(name, relative):
    spec = importlib.util.spec_from_file_location(name, REPO / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


driver = helper("build_wheel", "cmake/build_wheel.py")
rpath = helper("wheel_rpath", "cmake/wheel_rpath.py")
native = helper("native_paths", "python/_native.py")
smoke = helper("smoke_check", "tools/wheels/test_hepmc3_wheel.py")


def run(command, **kwargs):
    result = subprocess.run(command, text=True, capture_output=True, **kwargs)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def tiny_wheel(directory, version, extra=None):
    """A real pip-installable SIREN sentinel with no backend or native deps."""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"siren-{version}-py3-none-any.whl"
    metadata = f"siren-{version}.dist-info"
    files = {
        "siren/__init__.py": f"marker = {version!r}\n",
        f"{metadata}/METADATA": f"Metadata-Version: 2.1\nName: siren\nVersion: {version}\n",
        f"{metadata}/WHEEL": "Wheel-Version: 1.0\nRoot-Is-Purelib: true\nTag: py3-none-any\n",
    }
    files.update(extra or {})
    files[f"{metadata}/RECORD"] = "".join(f"{name},,\n" for name in files) + f"{metadata}/RECORD,,\n"
    with zipfile.ZipFile(path, "w") as archive:
        for name, value in files.items():
            archive.writestr(name, value)
    return path


def fake_run(handler):
    """Stand in for subprocess.run in the driver: cmake succeeds, pip is handled."""
    def run(command, **kwargs):
        if "pip" in command:
            return handler(command)
        return SimpleNamespace(returncode=0, stdout="", stderr="")
    return run


def driver_project(tmp_path):
    source, build = tmp_path / "source", tmp_path / "build"
    (source / "package").mkdir(parents=True)
    for name in ("pyproject.toml", "README.md", "LICENSE", "package/CMakeLists.txt"):
        (source / name).write_text("new input\n")
    build.mkdir()
    argv = ["build_wheel", "--source", str(source), "--build", str(build),
            "--config", "Release", "--cmake", "cmake", "--library-dir", "siren.libs"]
    return source, build, argv


@pytest.mark.parametrize("component", [None, "Unspecified"])
def test_cmake_install_never_invokes_pip_or_installs_python(tmp_path, component):
    # The wheel is installed by the interpreter the user selects, never by
    # CMake: no install component may reach for pip or write Python files.
    source, build = tmp_path / "source", tmp_path / "build"
    source.mkdir()
    (source / "CMakeLists.txt").write_text(
        'cmake_minimum_required(VERSION 3.20)\nproject(install_probe LANGUAGES NONE)\n'
        'set(SIREN_PYTHON_PACKAGE ON)\n'
        f'include("{REPO / "cmake/siren_python_package.cmake"}")\n')
    pip_blocker = tmp_path / "pip blocker"
    pip_blocker.mkdir()
    (pip_blocker / "pip.py").write_text('raise RuntimeError("CMake installation must not invoke pip")\n')
    env = dict(os.environ, PYTHONPATH=str(pip_blocker))
    run(["cmake", "-S", str(source), "-B", str(build), f"-DPython_EXECUTABLE={sys.executable}"], env=env)
    tiny_wheel(build / "dist_wheels", "0.0.2")
    prefix = tmp_path / "separate install"
    command = ["cmake", "--install", str(build), "--prefix", str(prefix)]
    if component:
        command += ["--component", component]
    run(command, env=env)
    assert not prefix.exists() or not any(prefix.rglob("*"))
    assert list((build / "dist_wheels").glob("*.whl"))


def test_staged_tree_digest_hashes_what_the_backend_packages(tmp_path):
    # The wheel backend reads through links, so the cache key must too:
    # changed bytes behind an unchanged link are a content change.
    staging = tmp_path / "staging"
    (staging / "siren/resources").mkdir(parents=True)
    module = staging / "siren/__init__.py"
    module.write_text("first\n")
    external = tmp_path / "external"
    external.mkdir()
    (external / "data.txt").write_text("first payload\n")
    (staging / "siren/data.txt").symlink_to(external / "data.txt")
    (staging / "siren/resources/linked").symlink_to(external, target_is_directory=True)
    (staging / "siren/self").symlink_to(".")
    (staging / "siren/missing").symlink_to(tmp_path / "absent")
    first = driver.staged_digests(staging)
    assert set(first) == {"siren/__init__.py", "siren/data.txt", "siren/resources/linked/data.txt",
                          "siren/self", "siren/missing"}
    assert first["siren/data.txt"] == first["siren/resources/linked/data.txt"]
    assert first["siren/self"] == "cycle:." and first["siren/missing"].startswith("dangling:")
    (external / "data.txt").write_text("changed bytes\n")
    second = driver.staged_digests(staging)
    assert second["siren/data.txt"] != first["siren/data.txt"]
    assert second["siren/resources/linked/data.txt"] != first["siren/resources/linked/data.txt"]
    module.unlink()
    assert "siren/__init__.py" not in driver.staged_digests(staging)


def test_staging_is_rebuilt_from_scratch(tmp_path, monkeypatch):
    # A removed source file must not survive in the staging tree the wheel is
    # built from: staging starts empty on every invocation.
    source, build, argv = driver_project(tmp_path)
    staging = build / "python_staging"
    (staging / "siren").mkdir(parents=True)
    (staging / "siren/removed.py").write_text("removed upstream\n")
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(driver.subprocess, "run", fake_run(
        lambda command: tiny_wheel(Path(command[command.index("--wheel-dir") + 1]), "0.0.2")))
    driver.main()
    assert not (staging / "siren/removed.py").exists()
    assert (staging / "CMakeLists.txt").read_text() == "new input\n"
    assert sorted(driver.read_state(build / ".wheel_inputs.json")["staged"]) == [
        "CMakeLists.txt", "LICENSE", "README.md", "pyproject.toml"]


def test_failed_staging_preserves_previous_wheel(tmp_path, monkeypatch):
    source, build, argv = driver_project(tmp_path)
    state = build / ".wheel_inputs.json"
    state.write_text('{"staged": {"old": "signature"}}')
    wheel = tiny_wheel(build / "dist_wheels", "0.0.1")
    old_state, old_wheel = state.read_bytes(), wheel.read_bytes()
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(driver.subprocess, "run", lambda command, **kwargs: SimpleNamespace(
        returncode=3, stdout="", stderr="CMake Error: install failed\n"))
    with pytest.raises(subprocess.CalledProcessError):
        driver.main()
    assert wheel.read_bytes() == old_wheel and state.read_bytes() == old_state


@pytest.mark.parametrize("contents", ['{"payload":', '[]', 'null', b'\xff'])
def test_corrupt_cache_is_discarded(tmp_path, capsys, contents):
    state = tmp_path / ".wheel_inputs.json"
    state.write_bytes(contents if isinstance(contents, bytes) else contents.encode())
    assert driver.read_state(state) == {}
    assert "ignoring invalid cache" in capsys.readouterr().out


def test_failed_wheel_build_preserves_previous_output(tmp_path, monkeypatch):
    source, build, argv = driver_project(tmp_path)
    state = build / ".wheel_inputs.json"
    state.write_text('{"staged": {"old": "signature"}}')
    wheel = tiny_wheel(build / "dist_wheels", "0.0.1")
    old_state, old_wheel = state.read_bytes(), wheel.read_bytes()
    monkeypatch.setattr(sys, "argv", argv)
    def fail_pip(command):
        raise subprocess.CalledProcessError(1, command)
    monkeypatch.setattr(driver.subprocess, "run", fake_run(fail_pip))
    with pytest.raises(subprocess.CalledProcessError):
        driver.main()
    assert wheel.read_bytes() == old_wheel
    assert state.read_bytes() == old_state
    assert not list(build.glob(".wheel-build-*"))


@pytest.mark.parametrize("platform", ["linux", "darwin"])
def test_undecodable_image_paths_preserve_native_guard(tmp_path, monkeypatch, platform):
    raw = b"/tmp/non-utf8-\xff/libSIREN.so"
    monkeypatch.setattr(native.sys, "platform", platform)
    if platform == "linux":
        maps = tmp_path / "maps"
        maps.write_bytes(b"1000-2000 r-xp 0000 08:01 123 " + raw + b"\n")
        monkeypatch.setattr(native, "Path", lambda path:
                            maps if path == "/proc/self/maps" else Path(path))
    else:
        monkeypatch.setattr(native.ctypes, "CDLL", lambda _: SimpleNamespace(
            _dyld_image_count=lambda: 1, _dyld_get_image_name=lambda index: raw))
    assert os.fsencode(native.loaded_libraries()[0]) == raw
    with pytest.raises(ImportError, match="separate processes"):
        native.reject_standalone_runtime()


def macho(tmp_path, paths_by_arch):
    source = tmp_path / "probe.c"
    source.write_text("int probe(void) { return 42; }\n")
    slices = []
    for architecture, paths in paths_by_arch.items():
        path = tmp_path / (architecture + ".dylib")
        flags = [f"-Wl,-rpath,{entry}" for entry in paths]
        run(["/usr/bin/clang", "-arch", architecture, "-dynamiclib", "-Wl,-headerpad_max_install_names",
             "-mmacosx-version-min=11.0", *flags, str(source), "-o", str(path)])
        slices.append(path)
    library = tmp_path / "libprobe.dylib"
    if len(slices) == 1:
        shutil.copy2(slices[0], library)
    else:
        run(["/usr/bin/lipo", "-create", *map(str, slices), "-output", str(library)])
    return library


@pytest.mark.skipif(sys.platform != "darwin", reason="Mach-O tooling")
@pytest.mark.parametrize("paths_by_arch", [
    {"arm64": []},
    {"arm64": [], "x86_64": []},
    {"arm64": [], "x86_64": ["/external with spaces/lib"]},
    {"arm64": ["/arm/lib", "@loader_path"], "x86_64": ["/intel/lib", "@loader_path"]},
])
def test_rpath_normalization_handles_empty_and_differing_slices(tmp_path, paths_by_arch):
    library = macho(tmp_path, paths_by_arch)
    mode = library.stat().st_mode
    rpath.prefer_siblings(library)
    assert library.stat().st_mode == mode
    assert set(run(["/usr/bin/lipo", "-archs", str(library)]).split()) == set(paths_by_arch)
    for architecture, paths in paths_by_arch.items():
        output = run(["/usr/bin/otool", "-arch", architecture, "-l", str(library)])
        actual = re.findall(r"cmd LC_RPATH\n\s+cmdsize \d+\n\s+path (.*?) \(offset \d+\)", output)
        assert actual == ["@loader_path"] + [p for p in paths if p != "@loader_path"]
    run(["/usr/bin/codesign", "--verify", str(library)])
    before = library.read_bytes()
    rpath.prefer_siblings(library)
    assert library.read_bytes() == before


@pytest.mark.skipif(sys.platform != "darwin", reason="Mach-O tooling")
def test_failed_rpath_edit_preserves_original(tmp_path, monkeypatch):
    library = macho(tmp_path, {"arm64": ["/external/lib"], "x86_64": []})
    before = library.read_bytes()
    real_run = rpath.subprocess.run
    def failed_sign(command, **kwargs):
        if command[0] == "/usr/bin/codesign":
            raise subprocess.CalledProcessError(1, command)
        return real_run(command, **kwargs)
    monkeypatch.setattr(rpath.subprocess, "run", failed_sign)
    with pytest.raises(subprocess.CalledProcessError):
        rpath.prefer_siblings(library)
    assert library.read_bytes() == before
    assert not list(tmp_path.glob(".wheel-rpath-*"))


def test_unreadable_cache_is_a_cache_miss(tmp_path, monkeypatch, capsys):
    state = tmp_path / ".wheel_inputs.json"
    state.write_text("{}")
    def denied(path):
        raise PermissionError("cache is unreadable")
    monkeypatch.setattr(driver.Path, "read_text", denied)
    assert driver.read_state(state) == {}
    assert "cache is unreadable" in capsys.readouterr().out


def test_directory_at_cache_path_does_not_prevent_build(tmp_path, monkeypatch, capsys):
    source, build, argv = driver_project(tmp_path)
    state = build / ".wheel_inputs.json"
    state.mkdir()
    sentinel = state / "keep"
    sentinel.write_text("do not remove unrelated directory contents")
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(driver.subprocess, "run", fake_run(
        lambda command: tiny_wheel(Path(command[command.index("--wheel-dir") + 1]), "0.0.2")))
    driver.main()
    assert list((build / "dist_wheels").glob("siren-0.0.2-*.whl"))
    assert sentinel.read_text() == "do not remove unrelated directory contents"
    assert "input cache could not be saved" in capsys.readouterr().out
    assert not list(build.glob(".wheel-build-*"))
