"""Failure recovery and install isolation for the native wheel helpers."""

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


@pytest.mark.parametrize("component", [None, "PythonPackage"])
def test_prefix_install_preserves_interpreter_package(tmp_path, component):
    # Only ever run real pip installs in a new, private venv. An active SIREN
    # install is the witness; --prefix must not uninstall it from this venv.
    environment = tmp_path / "build interpreter"
    run([sys.executable, "-m", "venv", str(environment)])
    python = environment / "bin/python"
    env = {key: value for key, value in os.environ.items()
           if key not in ("PYTHONPATH", "PIP_TARGET", "PIP_PREFIX", "PIP_USER")}
    env.update(PIP_CONFIG_FILE=os.devnull, PIP_NO_INDEX="1", PIP_DISABLE_PIP_VERSION_CHECK="1")
    old = tiny_wheel(tmp_path / "old", "0.0.1")
    run([str(python), "-m", "pip", "install", "--no-deps", str(old)], env=env)
    site = Path(run([str(python), "-c", "import sysconfig; print(sysconfig.get_path('purelib'))"], env=env).strip())
    before = {str(path.relative_to(site)): path.read_bytes()
              for path in site.rglob("*") if path.is_file() and "siren" in str(path.relative_to(site))}
    source, build = tmp_path / "source", tmp_path / "build"
    source.mkdir()
    (source / "CMakeLists.txt").write_text(
        'cmake_minimum_required(VERSION 3.20)\nproject(install_probe LANGUAGES NONE)\n'
        'set(SIREN_PYTHON_PACKAGE ON)\n'
        f'include("{REPO / "cmake/siren_python_package.cmake"}")\n')
    run(["cmake", "-S", str(source), "-B", str(build), f"-DPython_EXECUTABLE={python}"], env=env)
    tiny_wheel(build / "dist_wheels", "0.0.2")
    prefix = tmp_path / "separate install"
    command = ["cmake", "--install", str(build), "--prefix", str(prefix)]
    if component:
        command += ["--component", component]
    run(command, env=env)
    assert "0.0.1" in run([str(python), "-c", "import siren; print(siren.marker)"], env=env)
    assert all((site / name).read_bytes() == contents for name, contents in before.items())
    installed = list(prefix.rglob("siren/__init__.py"))
    assert len(installed) == 1 and "0.0.2" in installed[0].read_text()


@pytest.mark.parametrize("contents", ['{"payload":', '[]', 'null', b'\xff'])
def test_corrupt_cache_is_discarded(tmp_path, capsys, contents):
    state = tmp_path / ".wheel_inputs.json"
    state.write_bytes(contents if isinstance(contents, bytes) else contents.encode())
    assert driver.read_state(state) == {}
    assert "ignoring invalid cache" in capsys.readouterr().out


def test_failed_wheel_build_preserves_previous_output(tmp_path, monkeypatch):
    source, build = tmp_path / "source", tmp_path / "build"
    (source / "package").mkdir(parents=True)
    for name in ("pyproject.toml", "README.md", "LICENSE", "package/CMakeLists.txt"):
        (source / name).write_text("new input\n")
    build.mkdir()
    inputs = build / "inputs.txt"
    inputs.write_text(str(source / "README.md") + "\n")
    (build / "cmake_install.cmake").write_text("install rules")
    state = build / ".wheel_inputs.json"
    state.write_text('{"payload": {"old": "signature"}}')
    wheel = tiny_wheel(build / "dist_wheels", "0.0.1")
    old_state, old_wheel = state.read_bytes(), wheel.read_bytes()
    monkeypatch.setattr(sys, "argv", ["build_wheel", "--source", str(source), "--build", str(build),
                                      "--inputs", str(inputs), "--config", "Release", "--cmake", "cmake",
                                      "--library-dir", "siren.libs"])
    def fail_pip(command, **kwargs):
        if "pip" in command:
            raise subprocess.CalledProcessError(1, command)
    monkeypatch.setattr(driver.subprocess, "run", fail_pip)
    with pytest.raises(subprocess.CalledProcessError):
        driver.main()
    assert wheel.read_bytes() == old_wheel
    assert state.read_bytes() == old_state
    assert not list(build.glob(".wheel-build-*"))


def test_missing_payload_is_distinct_from_missing_build_input(tmp_path):
    assert driver.payload_digest(tmp_path / "python/deleted.py", tmp_path) is None
    assert driver.payload_digest(tmp_path / "resources/deleted.txt", tmp_path) is None
    with pytest.raises(FileNotFoundError, match="Required wheel input missing"):
        driver.payload_digest(tmp_path / "build/libSIREN_python.so", tmp_path)
    with pytest.raises(FileNotFoundError, match="Required wheel input missing"):
        driver.payload_digest(tmp_path / "pyproject.toml", tmp_path)


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


@pytest.mark.parametrize("dist_name", [None, ""])
def test_broken_distribution_metadata_does_not_mask_substitution(tmp_path, monkeypatch, dist_name):
    own = Path("siren.libs/libgfortran.so.5")
    foreign = Path("other/libgfortran.so.5")
    wheel = SimpleNamespace(files=[own], locate_file=lambda path: tmp_path / path)
    broken = SimpleNamespace(files=[foreign], locate_file=wheel.locate_file,
                             metadata={"Name": dist_name})
    monkeypatch.setattr(smoke.metadata, "distributions", lambda: [broken])
    with pytest.raises(AssertionError, match="outside this wheel"):
        smoke.verify_bundled_libraries(wheel, [tmp_path / own, tmp_path / foreign])


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


def test_preloaded_non_pip_copy_requires_own_copy(tmp_path, monkeypatch):
    own = tmp_path / "siren.libs/libexample.so.1"
    other = tmp_path / "conda/lib/libexample.so.1"
    wheel = SimpleNamespace(files=[own], locate_file=lambda path: path)
    monkeypatch.setattr(smoke.metadata, "distributions", lambda: [])
    smoke.verify_bundled_libraries(wheel, [own, other], preexisting=[other])
    with pytest.raises(AssertionError, match="outside this wheel"):
        smoke.verify_bundled_libraries(wheel, [other], preexisting=[other])
    # A competing dependency loaded by this import is still a substitution.
    with pytest.raises(AssertionError, match="outside this wheel"):
        smoke.verify_bundled_libraries(wheel, [own, other])


def test_identical_dependency_reuse_checks_bytes_not_name(tmp_path, monkeypatch):
    name = "libcfitsio-a1b2c3d4.so.4"
    own, other = tmp_path / "siren.libs" / name, tmp_path / "healpy.libs" / name
    for path in (own, other):
        path.parent.mkdir()
        path.write_bytes(b"identical library contents")
    wheel = SimpleNamespace(files=[own], locate_file=lambda path: path)
    monkeypatch.setattr(smoke.metadata, "distributions", lambda: [])
    smoke.verify_bundled_libraries(wheel, [other])
    other.write_bytes(b"substituted library data!!")
    with pytest.raises(AssertionError, match="outside this wheel"):
        smoke.verify_bundled_libraries(wheel, [other])
    other.unlink()
    with pytest.raises(AssertionError, match="outside this wheel"):
        smoke.verify_bundled_libraries(wheel, [other])


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="ELF SONAME reuse")
def test_elf_loader_reuses_identical_dependency(tmp_path):
    name = "libshared-a1b2c3d4.so.1"
    own = tmp_path / "siren.libs" / name
    foreign = tmp_path / "other.libs" / name
    own.parent.mkdir()
    foreign.parent.mkdir()
    dependency = tmp_path / "dep.c"
    dependency.write_text("int dependency(void) { return 42; }\n")
    run(["cc", "-shared", "-fPIC", f"-Wl,-soname,{name}", str(dependency), "-o", str(own)])
    shutil.copy2(own, foreign)
    consumer = tmp_path / "consumer.c"
    consumer.write_text("extern int dependency(void); int probe(void) { return dependency(); }\n")
    library = tmp_path / "consumer.so"
    run(["cc", "-shared", "-fPIC", str(consumer), str(own), "-Wl,-rpath,$ORIGIN/siren.libs",
         "-o", str(library)])
    # A fresh child supplies actual link-map evidence, independent of our path
    # comparison. The expected copy remains on disk but is never mapped.
    program = '''
import ctypes, importlib.util, sys
from pathlib import Path
from types import SimpleNamespace
root, repo = map(Path, sys.argv[1:])
def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, repo / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
native = load("native_probe", "python/_native.py")
smoke = load("wheel_smoke", "tools/wheels/test_hepmc3_wheel.py")
own = next((root / "siren.libs").iterdir())
foreign = root / "other.libs" / own.name
ctypes.CDLL(str(foreign), mode=ctypes.RTLD_GLOBAL)
consumer = ctypes.CDLL(str(root / "consumer.so"))
assert consumer.probe() == 42
mapped = set(native.loaded_libraries())
assert foreign in mapped and own not in mapped, mapped
wheel = SimpleNamespace(files=[own], locate_file=lambda path: path)
smoke.verify_bundled_libraries(wheel, mapped)
'''
    run([sys.executable, "-c", program, str(tmp_path), str(REPO)])


@pytest.fixture
def prefix_installation(tmp_path):
    """Exercise actual installs without exposing the test runner's packages."""
    environment = tmp_path / "interpreter"
    run([sys.executable, "-m", "venv", str(environment)])
    python = environment / "bin/python"
    env = {key: value for key, value in os.environ.items()
           if key not in ("PYTHONPATH", "PIP_TARGET", "PIP_PREFIX", "PIP_USER", "PIP_ROOT", "DESTDIR")}
    env.update(PIP_CONFIG_FILE=os.devnull, PIP_NO_INDEX="1", PIP_DISABLE_PIP_VERSION_CHECK="1")
    old = tiny_wheel(tmp_path / "old", "0.0.1", {"siren/obsolete.py": "obsolete = True\n"})
    run([str(python), "-m", "pip", "install", "--no-deps", str(old)], env=env)
    site = Path(run([str(python), "-c", "import sysconfig; print(sysconfig.get_path('purelib'))"], env=env).strip())
    source, build = tmp_path / "source", tmp_path / "build"
    source.mkdir()
    (source / "CMakeLists.txt").write_text(
        'cmake_minimum_required(VERSION 3.20)\nproject(install_probe LANGUAGES NONE)\n'
        'set(SIREN_PYTHON_PACKAGE ON)\n'
        f'include("{REPO / "cmake/siren_python_package.cmake"}")\n')
    run(["cmake", "-S", str(source), "-B", str(build), f"-DPython_EXECUTABLE={python}"], env=env)
    return SimpleNamespace(environment=environment, python=python, site=site, build=build, env=env)


@pytest.mark.parametrize("version", ["0.0.1", "0.0.2"])
@pytest.mark.parametrize("spelling", ["direct", "symlink", "dotdot"])
def test_same_prefix_install_cleanly_replaces_previous_wheel(tmp_path, prefix_installation, version, spelling):
    install = prefix_installation
    prefix = install.environment
    if spelling == "symlink":
        prefix = tmp_path / "interpreter alias"
        prefix.symlink_to(install.environment, target_is_directory=True)
    elif spelling == "dotdot":
        prefix = install.environment / "bin/.."
    tiny_wheel(install.build / "dist_wheels", version)
    run(["cmake", "--install", str(install.build), "--prefix", str(prefix)], env=install.env)
    assert not (install.site / "siren/obsolete.py").exists()
    assert [path.name for path in install.site.glob("siren-*.dist-info")] == [f"siren-{version}.dist-info"]
    assert version in run([str(install.python), "-c", "import siren; print(siren.marker)"], env=install.env)
    run([str(install.python), "-m", "pip", "uninstall", "-y", "siren"], env=install.env)
    assert not (install.site / "siren").exists()
    assert not list(install.site.glob("siren-*.dist-info"))


@pytest.mark.parametrize("same_prefix", [True, False])
def test_destdir_stages_python_without_changing_live_prefix(tmp_path, prefix_installation, same_prefix):
    install = prefix_installation
    prefix = install.environment if same_prefix else tmp_path / "new live prefix"
    root = tmp_path / "package root"
    before = {str(path.relative_to(install.site)): path.read_bytes()
              for path in install.site.rglob("*") if path.is_file() and "siren" in str(path.relative_to(install.site))}
    tiny_wheel(install.build / "dist_wheels", "0.0.2")
    run(["cmake", "--install", str(install.build), "--prefix", str(prefix)],
        env=dict(install.env, DESTDIR=str(root)))
    staged = root / prefix.relative_to(prefix.anchor)
    modules = list(staged.rglob("siren/__init__.py"))
    assert len(modules) == 1 and "0.0.2" in modules[0].read_text()
    assert all((install.site / name).read_bytes() == value for name, value in before.items())
    if not same_prefix:
        assert not prefix.exists()
    assert "0.0.1" in run([str(install.python), "-c", "import siren; print(siren.marker)"], env=install.env)


@pytest.mark.parametrize("name", ["libphotospline.2.dylib", "libspglam.so.2", "libphotospline-abcd.so.2"])
@pytest.mark.parametrize("ownership", ["preloaded", "distribution", "identical"])
def test_vendored_library_cannot_be_exempted(tmp_path, monkeypatch, name, ownership):
    own, foreign = tmp_path / "siren.libs" / name, tmp_path / "native" / name
    for path in (own, foreign):
        path.parent.mkdir()
        path.write_bytes(b"same library" if ownership == "identical" else str(path).encode())
    wheel = SimpleNamespace(files=[own], locate_file=lambda path: path)
    other = SimpleNamespace(files=[foreign], locate_file=wheel.locate_file, metadata={"Name": "other"})
    monkeypatch.setattr(smoke.metadata, "distributions", lambda: [other] if ownership == "distribution" else [])
    with pytest.raises(AssertionError, match="outside this wheel"):
        smoke.verify_bundled_libraries(wheel, [own, foreign],
                                       preexisting=[foreign] if ownership == "preloaded" else [])


@pytest.mark.parametrize("spelling", ["direct", "relative", "dotdot", "source-symlink", "input-symlink"])
def test_removed_resource_path_is_canonicalized(tmp_path, monkeypatch, spelling):
    source = tmp_path / "source"
    (source / "resources").mkdir(parents=True)
    alias = tmp_path / "alias"
    alias.symlink_to(source, target_is_directory=True)
    removed = source / "resources/deleted.txt"
    monkeypatch.chdir(tmp_path)
    if spelling == "relative":
        source = Path("source")
    elif spelling == "dotdot":
        source = source / "resources/.."
    elif spelling == "source-symlink":
        source = alias
    elif spelling == "input-symlink":
        removed = alias / "resources/deleted.txt"
    assert driver.payload_digest(removed, source) is None


def test_unreadable_cache_is_a_cache_miss(tmp_path, monkeypatch, capsys):
    state = tmp_path / ".wheel_inputs.json"
    state.write_text("{}")
    def denied(path):
        raise PermissionError("cache is unreadable")
    monkeypatch.setattr(driver.Path, "read_text", denied)
    assert driver.read_state(state) == {}
    assert "cache is unreadable" in capsys.readouterr().out


def test_directory_at_cache_path_does_not_prevent_build(tmp_path, monkeypatch, capsys):
    source, build = tmp_path / "source", tmp_path / "build"
    (source / "package").mkdir(parents=True)
    for name in ("pyproject.toml", "README.md", "LICENSE", "package/CMakeLists.txt"):
        (source / name).write_text("new input\n")
    build.mkdir()
    inputs = build / "inputs.txt"
    inputs.write_text(str(source / "README.md") + "\n")
    state = build / ".wheel_inputs.json"
    state.mkdir()
    sentinel = state / "keep"
    sentinel.write_text("do not remove unrelated directory contents")
    monkeypatch.setattr(sys, "argv", ["build_wheel", "--source", str(source), "--build", str(build),
                                      "--inputs", str(inputs), "--config", "Release", "--cmake", "cmake",
                                      "--library-dir", "siren.libs"])
    def build_wheel(command, **kwargs):
        if "pip" in command:
            tiny_wheel(Path(command[command.index("--wheel-dir") + 1]), "0.0.2")
    monkeypatch.setattr(driver.subprocess, "run", build_wheel)
    driver.main()
    assert list((build / "dist_wheels").glob("siren-0.0.2-*.whl"))
    assert sentinel.read_text() == "do not remove unrelated directory contents"
    assert "input cache could not be saved" in capsys.readouterr().out
    assert not list(build.glob(".wheel-build-*"))
