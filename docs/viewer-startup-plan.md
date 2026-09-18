# Viewer startup and rendering improvements

Work is on `fix/viewer-startup`, based on main `92d89d70`. The detector GDML
and simulation geometry remain authoritative. Display approximations are
explicit options and never get written back into simulation inputs.

The first change caches shared Boolean operands during a GDML read and removes
unrequested section pipelines. On the detailed CCM detector this reduces
loading/meshing from 191 s to 50 s without changing any of the 78 meshes.
Profiling the remaining work attributes 47 s to Boolean mesh construction.

## Implementation sequence and acceptance

1. **Persistent mesh cache.** Save reusable local meshes using input/dependency
   content hashes, effective meshing settings, backend/version and a format
   version. Use validated numeric data, atomic writes, and cache-miss recovery
   for corrupt entries. Check warm/cold mesh equality and invalidation after
   input, dependency or resolution changes. Expose cache location and disable.
2. **Scoped configuration and phase timings.** Restore pyg4ometry settings on
   success and failure. Report preparation, parsing, meshing/cache, scene
   construction, first preview frame and full-detail frame separately. Measure
   actual rendered frames before making window-startup claims.
3. **Deferred geometry.** Select display volumes before meshing, retain
   placements through hidden parents, and load newly requested volumes on
   demand. Handle replicas/divisions/parameterised volumes through an explicit
   compatible fallback. Verify hidden geometry does not get meshed unless it
   is an operand of a requested solid.
4. **Prototype instancing.** Triangulate each prototype once, share it across
   placement transforms, and retain per-instance identity for picking. Use
   a compatible expanded representation for clipping/sections or transforms
   that cannot be represented faithfully by the instancing mapper. Verify
   nested rotation/translation/scale, bounds, picking and cutaway behavior.
5. **Exterior surface display.** Add an explicit display mode for annotated PMT
   regions that suppresses vacuum and internal shell faces while preserving
   the outer surface and TPB/photocathode/reflector/bare-glass regions. Preserve
   full-detail display for geometry inspection and cutaways. Unknown volumes
   retain their ordinary geometry; do not infer optical regions from density.
6. **Progressive opening.** Open a responsive window with loading status and a
   coarse preview, then refine. Run meshing in a separate process, keep VTK
   operations on the UI thread, report worker failures, support cancellation
   on close, and preserve camera/visibility state during refinement. Check the
   synchronous path for scripts/screenshots and the interactive event loop.

## Validation and completion record

- Preserve the existing exact-mesh comparison and test suite.
- Add focused correctness regressions for cache corruption/invalidation,
  deferred traversal, transforms, surface selection and process lifecycle.
- Benchmark cold and warm runs on the same detailed detector and resolution;
  record prototype/placed face counts and native/Python runtime provenance.
- Exercise real VTK rendering, camera controls, visibility, picking, cutaways,
  progressive refinement and early-close cleanup where supported.
- Keep implementation, validation, commit/publication and installation status
  separate. Evidence lives in
  `artifacts/siren-viewer-startup-20260917` in the SBN_PRISM workspace.

## Implemented result

All six steps are implemented. Usage and compatibility boundaries are in
[visualization.md](visualization.md). The full Python suite passes 410 tests
with 49 optional dependency/data/network skips, including real native-window
tests. Full-detail cold/warm mesh comparisons retain all 78 prototype meshes
exactly. Shared settings, cache corruption/invalidation, hidden-parent
transforms, replicas/divisions, surface selection, picking, cutaways, camera
preservation, deferred loading, early close and SIREN model input are covered.

On the detailed detector, full mesh preparation is 56.19 s cold and 0.97 s warm.
The exterior representation reduces the complete scene's placed faces from
7,270,134 to 3,349,534 (including the world). A native exterior inspection view
with the cryostat and liquid volume hidden opens its loading window in 3.30 s,
renders its first geometry in 7.72 s, finishes its cold preview in 25.51 s and
refines in 68.41 s. Its cached reopen reaches full detail in 8.17 s. These are
single runs on macOS with pyg4ometry 1.4.2/CGAL and VTK; warm runs skip previews.
Cold refinement includes the extra preview computation. No frame-rate claim
is made. Screenshots and provenance are in the workspace evidence directory.

An installed-package regression exposed a worker import collision after the
initial source-overlay validation: launching `_visualization_worker.py` as a
file made `siren.math` shadow Python's standard-library `math`. Workers now
launch with `python -m siren._visualization_worker`. A wheel-layout test without
the overlay reproduces the original `Vector3D` registration failure and passes
with the fix. All 39 focused tests pass against the rebuilt wheel with native
rendering enabled. Evidence is under `worker-import-fix/` in the task artifacts.

Native window-close handling also watches VTK's `Done` state, since Cocoa's
close button does not emit `ExitEvent`. Closing stops queued rendering and
deferred loading, cancels the worker, and detaches window callbacks, widgets
and interactor references. Native tests cover early close, completed geometry,
queued mesh results, the VTK exit callback, and reopening in one Python process.
The old queued-results witness rendered after close; the repaired witness does
not. Evidence is under `window-close-fix/`.

All 43 focused viewer tests pass against the rebuilt wheel, including native
close-button tests. The installed public `view()` also returns to IPython after
closing the actual detector during loading and after completion; a second view
opens in the same process.

Status: implemented and validated; this branch includes the source, regression
tests and usage documentation. The corrected wheel is installed in CCM/local.
The import and
shutdown repairs change only Python viewer code and package metadata; native
libraries and geometry resources remain byte-for-byte unchanged.

## Large SBN scene scaling repair

SBND and ICARUS exposed a regression in the new renderer: each incoming mesh
rescanned all actors for statistics and controls, and flattened logical volumes
lost sharing in cache files, worker transfers and actor construction. Both
pre-repair native runs exceeded 180 s. Processing now yields between bounded
event batches, statistics/control scans run once per frame, and shared solids
reuse arrays, transfer files, VTK sources and glyph actors. Material/surface-role
boundaries, logical names, exact placements, deferred selection and picking are
preserved. Eager special-placement fallbacks retain separate logical mesh keys.

Native full-detail frame times at the same 48-slice resolution are:

| Standard SBN export | Original main | Repaired cold | Repaired cached |
| --- | ---: | ---: | ---: |
| ICARUS | 89.87 s | 12.33 s | 9.29 s |
| SBND | 64.04 s | 28.68 s | 15.78 s |

These single runs include the complete exported SBN site; dependency imports
precede timing. Cold first geometry appears at 8.42 s and 12.17 s respectively.
Independent direct-mesher comparisons match all 1,987 unique visible surfaces;
all 24,865 visible placement rows are retained, with reconstructed glyph
rotations agreeing within 6.2e-16. All 47 focused viewer tests pass against the
rebuilt wheel with native rendering enabled. The corrected wheel is installed
in CCM/local; native libraries and geometry resources remain unchanged.
Reproduction, runtime identity and native renders are under `sbn-performance/`
in the task evidence.
