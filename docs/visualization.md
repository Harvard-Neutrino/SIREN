# Interactive geometry viewing

`siren.visualization.view` opens GDML files or SIREN detector models. Install
SIREN's `visualization-3d` extras to use the pyg4ometry/VTK backend. Native
window rendering requires a working display; worker meshing does not.

```python
from siren import visualization

timing = {}
visualization.view("detector-cad.gdml", timings=timing)
print(timing)
```

The window appears with loading status. A worker process prepares a coarse
preview and then the requested detail, while the window remains interactive.
Camera movements, material hiding and opacity changes survive refinement.
Closing the window returns control to Python, cancels outstanding work and
removes its temporary files. Queued mesh updates are discarded after close;
the Python process stays alive and another `view()` can be opened immediately.
When full-detail meshes are already cached, the preview pass is skipped.

## Cache and resolution

Mesh files default to `~/.cache/siren/meshes`, or the corresponding directory
under `XDG_CACHE_HOME`. They are disposable numeric NPZ files, with no Python
pickle data. Input contents, referenced GDML files, effective resolution,
backend/library versions and cache format determine reuse. Invalid files are
rebuilt; new files are published atomically after checking input stability.
Logical volumes referencing the same solid share a cached mesh and a worker
transfer file. Different materials, placements and volume names are retained.
The shared-solid cache format creates a new namespace, so the first load after
this update rebuilds meshes once; older cache entries can be discarded.

```python
visualization.view("detector.gdml", cache_dir="/tmp/my-viewer-cache")
visualization.view("detector.gdml", cache=False)  # bypass persistent cache
visualization.view("detector.gdml", mesh_slices=64, preview_slices=12)
visualization.view("detector.gdml", progressive=False)  # synchronous loading
```

The default full resolution is 48 slices/stacks, with a 12-slice preview.
`preview_slices=None` disables the preview. `mesh_slices=0` or `None` retains
the mesher defaults. Temporary changes to pyg4ometry settings are restored on
success and failure. Interactive preparation uses an isolated process.

The preview adds work on a cold load: it improves time to visible geometry,
not time to the final high-resolution scene. Choose `progressive=False` when
only full-detail completion matters. Screenshots and `interactive=False`
always wait for full detail.

## Deferred volumes and display surfaces

Low-density volumes (density <= 0.05 g/cm3) start hidden and are prepared on
demand. Press **g** to load/show them, or start with `show_gas=True`. Name globs
can defer other logical volumes; **v** restores them. Children of a hidden
volume keep their original world placements and can still be displayed.
Necessary operands of visible Boolean solids must still be meshed.

```python
visualization.view(
    "detector-cad.gdml",
    hidden_volumes=["CADCryostat*", "CADLiquidArgon*"],
    display="exterior",
)
```

`display="exterior"` recognizes the `ccm_pmt_region` annotations on CCM's
PMT shells of revolution about local z. It removes the inner shell faces and
artificial end caps between regions, and suppresses PMT vacuum. TPB,
photocathode, reflector and bare glass receive distinct display colours;
these colours are not measured optical coefficients. The ordinary geometry is
retained for other volumes. Opaque coating regions cover the underlying glass
where both are present; the bare neck is translucent.

These open display surfaces do not replace simulation solids. Use
`display="full"` for volume inspection. `cutter=True` and `clipper=True`
automatically use full closed meshes and an expanded pipeline for accurate
section/clip calculations. `clip_origin` is in GDML/VTK millimetres and
`clip_normal` gives the normal; the displayed half-space is the negative side.

Files with replicas, divisions or parameterised placements use pyg4ometry's
eager-meshing compatibility path. Their hidden-volume and warm-cache savings
are consequently limited. GDML external entities are not supported by the
scene loader. Referenced `<file>` paths follow pyg4ometry's working-directory
resolution; include their contents when sharing a geometry.

## Rendering and controls

Repeated placements and logical volumes referencing the same solid share
prototype meshes through GPU glyphs. Actor groups also require matching material
and surface role; picking retains each original volume name and placement.
Reflections or
shear use exact matrix actors sharing a mapper. Set `instancing=False` to use
matrix actors throughout. Right-click picking intersects the shared prototype
meshes without expanding every copy, including translucent objects. A supplied
SIREN model additionally supplies the material/sector query.

Mesh notifications are consumed in bounded batches so large cached scenes yield
to native window input. Actor statistics and control updates run once per frame,
instead of rescanning the scene after every incoming mesh.

- Left drag orbits; scroll zooms. Arrow keys pan or move forwards/backwards;
  shift-up/down moves vertically.
- Right-click identifies a placed volume; shift-right-click hides its material.
- **g** toggles gas/vacuum, **v** restores materials and opacity.
- **n/m** decrease/increase opacity, **c** toggles requested section outlines.
- **k** toggles the clipping widget, **l** the legend, **b** the bounding box,
  **o** the orientation widget, and **h** the help overlay.

```python
visualization.view("detector.gdml", clipper=True, clip_normal=(1, 0, 0))
visualization.view("detector.gdml", screenshot="detector.png")
visualization.view("detector.gdml", interactive=False, timings=True)
```

`timings={}` collects import/export, per-stage parsing, meshing/cache and cache
write times, scene-building/render times, cache hits/misses, face/instance
counts and cancellation. `window_frame_seconds` measures the first actual
render, which can contain loading text; `first_geometry_frame_seconds` measures
the first rendered geometry; `preview_frame_seconds` and `detail_frame_seconds`
measure completed stages. Times are elapsed from the view call, excluding
earlier imports by the caller. `timings=True` prints the record on return.
Stage cache hits/misses count unique mesh keys; `shared_mesh_reuses` counts
additional logical volumes served by those meshes. `prototype_faces` retains
its per-logical-volume count and is not a count of unique shared GPU faces.
`control_update_seconds` records time applying visibility and opacity controls.

The separate `backend="pyvista"` path retains its existing renderer and
controls; the new cache, worker and instancing options apply to pyg4ometry.
