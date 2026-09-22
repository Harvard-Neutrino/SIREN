# Interactive geometry viewing

`siren.visualization.view` opens a GDML file or a SIREN `DetectorModel` in a
VTK window. It needs the `visualization-3d` extras (pyg4ometry, VTK) and a
display; pyvista is only required for models with `TriangularMesh` sectors.

```python
from siren import visualization

visualization.view("detector.gdml")
visualization.view(model, screenshot="detector.png")   # off-screen PNG, full detail
visualization.view("detector.gdml", timings=True)      # print phase times on return
```

The window opens immediately with a loading status. A worker process meshes a
coarse preview, then the requested detail, while the window stays interactive;
camera position, hidden materials and opacity survive refinement. Closing the
window cancels outstanding work, removes its temporary files and returns to
Python, where another `view()` can be opened at once. Each worker process pays
the pyg4ometry import, so on-demand loads take a few seconds to start.

## Options

| Option | Default | Effect |
| --- | --- | --- |
| `mesh_slices` | 48 | Full-detail slices/stacks for curved solids; `0`/`None` keeps the mesher defaults. |
| `preview_slices` | 12 | Preview resolution; `None` disables the preview. Skipped when full detail is cached. |
| `progressive` | `True` | `False` loads synchronously. Screenshots and `interactive=False` always wait for full detail. |
| `cache`, `cache_dir` | `True`, `~/.cache/siren/meshes` | Persistent mesh cache (or under `XDG_CACHE_HOME`); `cache=False` neither reads nor writes it. |
| `show_gas` | `False` | Volumes with density <= 0.05 g/cm3 start hidden; **g** loads and shows them. |
| `hidden_volumes` | `()` | Logical-volume name globs to defer; **v** loads and restores them. Their children keep their placements. |
| `instancing` | `True` | Repeated placements share GPU prototypes; reflections and shear use exact matrix actors. |
| `display` | `"full"` | `"exterior"` shows only outward faces of annotated region shells (below). |
| `cutter`, `clipper` | `False` | Section outlines / clipping widget (**k**), with `clip_origin` in mm and `clip_normal`; the negative side is displayed. Both force `display="full"`. |
| `coloured`, `axes`, `legend`, `bounding_box`, `picker`, `near_frac` | | Usual display controls. |

Cache files are disposable numeric NPZ files (no pickles), keyed by input and
included-file contents, resolution, backend and library versions. Corrupt
files are rebuilt, and new files are published atomically only when the inputs
were stable throughout the load. Volumes sharing a solid share one cached mesh.
GDML files with `<!DOCTYPE>`/`<!ENTITY>` declarations still load through
pyg4ometry's reader but are never cached (a warning says so). Files using
replicas, divisions or parameterised placements use pyg4ometry's eager
meshing, which limits deferral and cache savings. `<file>` references resolve
against the working directory, as in pyg4ometry.

## Exterior display of region shells

`display="exterior"` keeps only the faces whose normals point radially outward
from the local z axis for volumes annotated as shells of revolution, removing
inner walls and the artificial end caps between regions. The mapping is
supplied by the caller; SIREN carries no detector-specific region names.

```python
regions = {
    "auxtype": "pmt_region",           # <auxiliary auxtype=...> key on logical volumes
    "styles": {                        # region -> legend label, colour, opacity
        "external_tpb": {"label": "PMT TPB", "colour": [0.95, 0.94, 0.82], "alpha": 1.0},
        "bare_transparent_glass": {"label": "PMT bare glass", "colour": [0.72, 0.88, 0.94], "alpha": 0.25},
    },
    "hidden": ["vacuum"],              # regions not displayed at all
}
visualization.view("detector-cad.gdml", display="exterior", regions=regions,
                   hidden_volumes=["CADCryostat*"])
```

Unannotated volumes are rendered normally. Colours are display choices, not
optical properties. These open surfaces never alter the GDML or simulation
solids; use `display="full"` for volume inspection.

## Controls

- Left drag orbits, scroll zooms; arrow keys pan and fly, shift-up/down moves vertically.
- Right-click identifies the placed volume (name, material, density; a supplied
  SIREN model adds its sector query); shift-right-click hides that material.
- **g** gas volumes, **v** restore all, **n/m** opacity, **c** section outlines,
  **k** clip widget, **l** legend, **b** bounding box, **o** orientation cube,
  **h** help, **q**/**e** quit.

## Timings

Pass a dictionary as `timings` to receive phase times (import, export, parse,
mesh/cache, cache write, scene build, render), first-frame times
(`window_frame_seconds`, `first_geometry_frame_seconds`, `preview_frame_seconds`,
`detail_frame_seconds`), cache hits/misses per stage, shared-mesh reuses,
face and instance counts, and whether the load was cancelled. Times are
measured from the `view()` call. `direct_mesh_error` records why a model's
mesh sectors were not drawn (typically a missing pyvista), and
`request_errors` collects failures of on-demand loads, which leave the viewer
open with a warning.

The `backend="pyvista"` path (`view_pv`) keeps its own renderer and controls;
the options above apply to the pyg4ometry backend.
