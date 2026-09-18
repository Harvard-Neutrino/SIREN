"""Real worker completion, failure and cancellation without a display."""
from pathlib import Path
import importlib.util
import os
import shutil
import subprocess
import sys
import textwrap
import time

import pytest

from siren._visualization_job import MeshJob
from test_visualization_scene import fixture_gdml


def test_worker_from_installed_package_without_source_overlay(tmp_path):
    # Source overlays put the worker beside Python sources, without math.so.
    # Reproduce a wheel layout, where siren.math must not shadow stdlib math.
    import siren
    package_root = Path(siren.__file__).parent
    installed = tmp_path / 'installed'
    package = installed / 'siren'
    package.mkdir(parents=True)
    for entry in package_root.iterdir():
        if entry.name != '__pycache__':
            (package / entry.name).symlink_to(entry, target_is_directory=entry.is_dir())
    libraries = package_root.parent / 'siren.libs'
    if libraries.exists():
        (installed / 'siren.libs').symlink_to(libraries, target_is_directory=True)
    for name in ('visualization', '_visualization_job', '_visualization_worker',
                 '_visualization_scene', '_visualization_render', '_visualization_view'):
        destination = package / (name + '.py')
        destination.unlink(missing_ok=True)
        shutil.copyfile(importlib.util.find_spec('siren.' + name).origin, destination)
    environment = os.environ.copy()
    environment['PYTHONPATH'] = str(installed)
    environment.pop('SIREN_VIEWER_SOURCE', None)
    result = subprocess.run(
        [sys.executable, '-c', textwrap.dedent('''
            import sys
            import time
            from siren._visualization_job import MeshJob
            with MeshJob(sys.argv[1], dict(cache=False), preview_slices=12) as job:
                events = []
                deadline = time.monotonic() + 30
                while not job.finished:
                    assert time.monotonic() < deadline, 'worker timed out'
                    events.extend(job.poll())
                    time.sleep(.01)
                assert events[-1]['kind'] == 'done', events
                assert [e['stage'] for e in events if e['kind'] == 'stage_done'] == ['preview', 'detail']
        '''), str(fixture_gdml(tmp_path / 'a.gdml'))],
        cwd=tmp_path, env=environment, capture_output=True, text=True, timeout=40)
    assert result.returncode == 0, result.stdout + result.stderr


def drain(job):
    events = []
    deadline = time.monotonic() + 20
    while not job.finished:
        if time.monotonic() > deadline:
            pytest.fail('meshing worker did not finish')
        events.extend(job.poll())
        time.sleep(.01)
    return events


def test_worker_streams_preview_then_detail(tmp_path):
    path = fixture_gdml(tmp_path / 'a.gdml')
    with MeshJob(path, dict(mesh_slices=24, cache_dir=str(tmp_path / 'cache')), preview_slices=12) as job:
        directory = job.directory.name
        events = drain(job)
        assert [e['stage'] for e in events if e['kind'] == 'stage_done'] == ['preview', 'detail']
        assert events[-1]['kind'] == 'done'
        assert all(Path(e['path']).exists() for e in events if e['kind'] == 'mesh')
    assert not Path(directory).exists()
    with MeshJob(path, dict(mesh_slices=24, cache_dir=str(tmp_path / 'cache')), preview_slices=12) as job:
        events = drain(job)
        assert [e['stage'] for e in events if e['kind'] == 'stage_done'] == ['detail']
        assert next(e for e in events if e['kind'] == 'stage_done')['cache_hits'] == 1


def test_worker_failure_is_reported(tmp_path):
    with MeshJob(tmp_path / 'missing.gdml', dict(cache=False)) as job:
        events = drain(job)
        assert events[-1]['kind'] == 'failed'
        assert 'missing.gdml' in events[-1]['message']


def test_worker_cancel_removes_child_and_temporary_files(tmp_path):
    job = MeshJob(fixture_gdml(tmp_path / 'a.gdml'), dict(cache=False))
    directory = job.directory.name
    job.close()
    assert job.process.poll() is not None
    assert not Path(directory).exists()
