"""Nonblocking worker transport and deterministic cancellation/cleanup."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile


class MeshJob:
    def __init__(self, path, options, preview_slices=None):
        self.directory = tempfile.TemporaryDirectory(prefix="siren-view-")
        folder = Path(self.directory.name)
        self.offset, self.pending, self.finished = 0, b"", False
        self.log = (folder / "worker.log").open("wb")
        try:
            request = folder / "request.json"
            request.write_text(json.dumps(dict(path=str(path), options=options,
                                               preview_slices=preview_slices)))
            # File execution puts siren/ first on sys.path, so the installed
            # siren.math extension shadows stdlib math and registers types twice.
            # Module execution preserves normal package-qualified imports.
            self.process = subprocess.Popen(
                [sys.executable, "-m", "siren._visualization_worker", str(request)],
                stdin=subprocess.DEVNULL, stdout=self.log, stderr=self.log)
        except BaseException:
            self.log.close()
            self.directory.cleanup()
            raise

    def poll(self):
        if self.finished:
            return []
        path = Path(self.directory.name) / "events.jsonl"
        try:
            with path.open("rb") as stream:
                stream.seek(self.offset)
                data = stream.read()
                self.offset += len(data)
        except FileNotFoundError:
            data = b""
        lines = (self.pending + data).split(b"\n")
        self.pending = lines.pop()
        events = [json.loads(line) for line in lines if line]
        if any(e["kind"] in ("done", "failed") for e in events):
            self.finished = True
        elif self.process.poll() is not None:
            # Read once more after exit; the worker may have appended between
            # the read and poll. EOF must not erase a just-written result.
            with path.open("rb") if path.exists() else open(self.log.name, "rb") as stream:
                if path.exists():
                    stream.seek(self.offset)
                    tail = stream.read()
                else:
                    tail = b""
            if tail:
                self.offset += len(tail)
                lines = (self.pending + tail).split(b"\n")
                self.pending = lines.pop()
                events.extend(json.loads(line) for line in lines if line)
            self.finished = True
            if not any(e["kind"] in ("done", "failed") for e in events):
                message = Path(self.log.name).read_text(errors="replace")[-8000:]
                events.append(dict(kind="failed", message="mesh worker exited (%s): %s" %
                                   (self.process.returncode, message)))
        return events

    def close(self):
        if self.process.poll() is None:
            self.process.terminate()
            try:
                self.process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait(timeout=2)
        self.log.close()
        self.directory.cleanup()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
