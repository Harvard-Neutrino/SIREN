"""Isolated meshing worker; events are JSON lines and meshes are numeric NPZ."""
import json
from pathlib import Path
import sys
import traceback


def run(request_path):
    from siren._visualization_scene import prepare_scene, cached_scene_available
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
            if (preview and detail and preview < detail and
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
