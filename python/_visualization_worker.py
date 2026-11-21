"""Isolated meshing worker; events are JSON lines and meshes are numeric NPZ."""
import json
from pathlib import Path
import signal
import sys
import traceback


def exit_on_terminate():
    """Turn the parent's terminate() into SystemExit.

    An ordinary signal death skips Python cleanup, leaving partial cache files
    behind. Raising lets ``write_mesh``'s finally block remove its temporary.
    """
    def handler(signum, frame):
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, handler)


def run(request_path):
    from siren._visualization_scene import prepare_scene, cached_scene_available
    exit_on_terminate()
    request = json.loads(Path(request_path).read_text())
    directory = Path(request_path).parent
    with (directory / "events.jsonl").open("w") as stream:
        def emit(event):
            stream.write(json.dumps(event) + "\n")
            stream.flush()
        try:
            options = request["options"]
            detail = options.pop("mesh_slices", 48)
            preview = request.get("preview_slices")
            effective = detail
            if not effective:
                # 0/None keeps the mesher defaults; compare against those.
                import pyg4ometry as pg
                effective = max([getattr(o, "nslice", 0) or 0
                                 for o in vars(pg.config.SolidDefaults).values()] + [0])
            if (preview and effective and preview < effective and
                    not cached_scene_available(request['path'], mesh_slices=detail, **options)):
                prepare_scene(request["path"], directory, mesh_slices=preview,
                              emit=emit, stage="preview", **options)
            prepare_scene(request["path"], directory, mesh_slices=detail,
                          emit=emit, stage="detail", **options)
            emit(dict(kind="done"))
        except Exception:
            emit(dict(kind="failed", message=traceback.format_exc()))
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(run(sys.argv[1]))
