"""Actual Cocoa close-button acceptance in disposable Python processes."""
import ctypes
import gc
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest


@pytest.mark.skipif(sys.platform != 'darwin' or os.environ.get('SIREN_TEST_VTK_RENDER') != '1',
                    reason='requires the native macOS display')
@pytest.mark.parametrize('phase', ['early', 'ready', 'queued', 'exit'])
def test_native_window_close_and_reopen(tmp_path, phase):
    process = subprocess.Popen([sys.executable, __file__, phase, str(tmp_path)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, start_new_session=True)
    try:
        output, _ = process.communicate(timeout=40)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        output, _ = process.communicate()
        pytest.fail('Closing the viewer did not return to Python:\n' + output)
    assert process.returncode == 0, output
    assert 'REOPEN_OK' in output


def exercise_close(root, phase):
    from siren._visualization_view import ViewSession
    from test_visualization_scene import fixture_gdml

    objc = ctypes.CDLL('/usr/lib/libobjc.A.dylib')
    objc.sel_registerName.argtypes = [ctypes.c_char_p]
    objc.sel_registerName.restype = ctypes.c_void_p
    objc.objc_msgSend.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
    objc.objc_msgSend.restype = None

    for attempt in range(2):
        session = ViewSession(fixture_gdml(root / 'a.gdml'), cache_dir=root / 'cache')
        closed_at, rendered_after_close, jobs = [], [], []

        def rendered(*args):
            if closed_at:
                rendered_after_close.append(True)

        def close(*args):
            if closed_at or (phase == 'ready' and not session.detail_ready):
                return
            if phase == 'queued':
                events = Path(session.job.directory.name) / 'events.jsonl'
                if not events.exists() or '"kind": "mesh"' not in events.read_text():
                    return
            if session.job:
                jobs.append((session.job.process, session.job.directory.name))
            closed_at.append(time.monotonic())
            if phase == 'exit':
                session.viewer.iren.ExitCallback()
            else:
                # Invoke the action used by the native red close button. This
                # path calls TerminateApp but does NOT emit VTK's ExitEvent.
                window = session.viewer.renWin.GetRootWindow()
                pointer = int(window.split('_')[1], 16)
                objc.objc_msgSend(pointer, objc.sel_registerName(b'performClose:'), None)

        session.viewer.renWin.AddObserver('StartEvent', rendered)
        session.viewer.iren.AddObserver('TimerEvent', close, 1.)
        if phase == 'queued':
            tick = session.tick
            # Hold results until close, then dispatch the already-queued timer
            # callback. The regression must not draw or accept those results.
            session.tick = lambda *args: tick() if closed_at else None
        session.run()
        assert closed_at and time.monotonic() - closed_at[0] < 5
        assert session.closed and session.job is None
        assert not rendered_after_close
        assert session.timer is None
        assert not session.viewer.iren.HasObserver('TimerEvent')
        assert session.viewer.iren.GetRenderWindow() is None
        assert all(p.poll() is not None and not Path(d).exists() for p, d in jobs)
        assert session.metrics['cancelled'] == (phase != 'ready')
        # Late callbacks after returning must remain harmless, too.
        session.tick()
        session.request_geometry('all')
        session.receive(dict(kind='failed', message='late worker event'))
        session._render()
        assert session.error is None and session.job is None
        assert not rendered_after_close
        gc.collect()
    print('REOPEN_OK', flush=True)


if __name__ == '__main__':
    exercise_close(Path(sys.argv[2]), sys.argv[1])
