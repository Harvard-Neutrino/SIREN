"""Progressive VTK window controller. All VTK calls run on the caller's thread."""
from collections import deque
import os
from pathlib import Path
import tempfile
import time

from ._visualization_job import MeshJob
from ._visualization_render import SceneRenderer
from ._visualization_scene import normalize_regions, prepare_scene, read_mesh
from .visualization import _LOW_DENSITY


class ViewSession:
    def __init__(self, model, gdml_path=None, *, screenshot=None, coloured=True,
                 axes=True, cutter=False, clipper=False, clip_origin=(0., 0., 0.),
                 clip_normal=(1., 0., 0.), legend=True, bounding_box=False,
                 picker=True, interactive=True, mesh_slices=48, near_frac=1e-3,
                 cache=True, cache_dir=None, progressive=True, preview_slices=12,
                 instancing=True, display="full", regions=None, show_gas=False,
                 hidden_volumes=(), timings=None, progress=True):
        self.started = time.perf_counter()
        self.metrics = {"stages": []}
        self.timings = timings
        self.job, self.controls, self.error = None, None, None
        self.pending = set()
        self.loaded = set()
        self.detail_ready = False
        self.camera_dirty = False
        self.camera_fitted = False
        self.closed = False
        self.timer = None
        self._timer_observer = None
        self._events = deque()
        self._mesh_files = {}
        self.interactive = interactive and not screenshot
        self.progressive = bool(progressive and self.interactive)
        self.preview_slices = preview_slices
        self.screenshot = screenshot
        self.model = None if isinstance(model, (str, os.PathLike)) else model
        self.axes, self.legend, self.bounding_box, self.picker = axes, legend, bounding_box, picker
        self.axes_added = False
        self.progress = bool(progress)
        self.near_frac = near_frac
        self.show_gas = show_gas
        self.direct_meshes_added = False
        if display not in ("full", "exterior"):
            raise ValueError("display must be 'full' or 'exterior'")
        regions = normalize_regions(regions)
        if display == "exterior" and regions is None:
            raise ValueError("display='exterior' requires a regions mapping")
        if preview_slices is not None and (isinstance(preview_slices, bool) or
                int(preview_slices) != preview_slices or preview_slices < 4):
            raise ValueError("preview_slices must be an integer >= 4 or None")
        if isinstance(hidden_volumes, str):
            hidden_volumes = [hidden_volumes]
        self.options = dict(mesh_slices=mesh_slices, cache=bool(cache),
                            cache_dir=None if cache_dir is None else str(cache_dir),
                            display="full" if cutter or clipper else display,
                            regions=regions, show_gas=show_gas,
                            hidden_volumes=list(hidden_volumes))
        self.metrics["display"] = self.options["display"]
        self.metrics["cutaway_full_detail"] = bool((cutter or clipper) and display != "full")
        import vtk
        from pyg4ometry.visualisation import VtkViewerNew
        from . import visualization as vis
        self.vtk = vtk
        self.vis = vis
        self.metrics["imports_seconds"] = time.perf_counter() - self.started
        self.viewer = VtkViewerNew(defaultCutters=bool(cutter), axisCubeWidget=bool(axes))
        # pyg4ometry's default style has a crash-prone right-click picker; a
        # plain trackball style guards the loading phase until _install_controls
        # replaces it. Each style carries its own camera observer while live.
        style = vtk.vtkInteractorStyleTrackballCamera()
        style.AddObserver('StartInteractionEvent', self._camera_changed)
        self.viewer.iren.SetInteractorStyle(style)
        if clipper:
            self.viewer.addClipper(list(clip_origin), list(clip_normal), True)
        self.coloured = bool(coloured)
        self.renderer = SceneRenderer(self.viewer, coloured=coloured, instancing=instancing,
                                      display=self.options["display"], regions=regions)
        if not picker:
            # picker=False keeps the point-only right-click of the plain viewer.
            self.viewer.pick_scene = None
        self.status = vis._text_actor(vtk, "Loading geometry...", .02, .96, size=18, anchor="tl")
        self.status.SetVisibility(bool(progress))
        self.viewer.ren.AddViewProp(self.status)
        self.viewer.request_geometry = self.request_geometry
        self.viewer.camera_moved = self._camera_changed
        self.viewer.iren.AddObserver("ExitEvent", self._exited)
        self.viewer.iren.Initialize()
        if self.progressive:
            self.viewer.renWin.Render()
            self.metrics["window_frame_seconds"] = time.perf_counter() - self.started
        before = time.perf_counter()
        try:
            if self.model is None:
                self.path = str(Path(model).resolve())
            else:
                if gdml_path is None:
                    with tempfile.NamedTemporaryFile(suffix=".gdml", delete=False) as out:
                        gdml_path = out.name
                self.path = str(Path(gdml_path).resolve())
                vis.to_gdml(model, self.path, skip_geo_types=("TriangularMesh",))
        except BaseException:
            # run() never starts, so its window teardown must happen here.
            self.closed = True
            self._release_window()
            raise
        self.metrics["export_seconds"] = time.perf_counter() - before

    def _camera_changed(self, *args):
        self.camera_dirty = True

    def _exited(self, *args):
        self.closed = True
        self.viewer.iren.TerminateApp()

    def _closing(self):
        # Native window-close (notably Cocoa) can call TerminateApp without
        # emitting ExitEvent. A TimerEvent already being dispatched must not
        # consume worker results or render into the closed window afterward.
        if self.viewer.iren.GetDone():
            self.closed = True
        return self.closed

    def _release_window(self):
        # These objects and callbacks belong to this one viewer window.
        # Detach widgets/styles before releasing the native window so a
        # later GC pass or another view() cannot dispatch stale callbacks.
        v = self.viewer
        v.iren.Disable()
        v.iren.EnableRenderOff()
        for name in ("axesWidget", "clipperPlaneWidget"):
            widget = getattr(v, name, None)
            if widget is not None:
                widget.SetEnabled(0)
                widget.SetInteractor(None)
        style = v.iren.GetInteractorStyle()
        if style is not None:
            style.RemoveAllObservers()
        v.iren.RemoveAllObservers()
        v.iren.SetInteractorStyle(None)
        v.request_geometry = v.camera_moved = None
        v.renWin.Finalize()
        v.iren.SetRenderWindow(None)
        v.renWin.SetInteractor(None)

    def _publish_metrics(self):
        self.metrics.update(self.renderer.stats)
        if isinstance(self.timings, dict):
            self.timings.update(self.metrics)

    def _render(self, stage=None):
        if self._closing():
            return
        if self.controls is not None:
            before = time.perf_counter()
            self.controls["refresh"]()
            self.metrics["control_update_seconds"] = self.metrics.get("control_update_seconds", 0.) + time.perf_counter() - before
        before = time.perf_counter()
        self.viewer.renWin.Render()
        now = time.perf_counter()
        self.metrics["render_seconds"] = self.metrics.get("render_seconds", 0.) + now - before
        self.metrics.setdefault("window_frame_seconds", now - self.started)
        if self.viewer.actors:
            self.metrics.setdefault("first_geometry_frame_seconds", now - self.started)
        if stage:
            self.metrics.setdefault(stage + "_frame_seconds", now - self.started)
        self._publish_metrics()

    def _setup_controls(self):
        from types import SimpleNamespace
        from pyg4ometry.visualisation import VisualisationOptions
        v = self.viewer
        if not self.direct_meshes_added and self.model is not None:
            self.direct_meshes_added = True
            try:
                self.vis._add_mesh_actors(self.vtk, v.ren, v.actors, self.model,
                                          self.renderer.options if self.coloured else None,
                                          self.renderer.registry)
            except Exception as exc:
                # TriangularMesh sectors need the optional pyvista dependency. The
                # GDML scene is complete without their direct actors, so report
                # the omission instead of failing the whole viewer.
                import warnings
                self.metrics["direct_mesh_error"] = str(exc)
                warnings.warn("direct mesh sectors were not rendered: %s" % exc, RuntimeWarning)
        for key, actor in v.actors.items():
            if key.startswith('mesh__'):
                self.renderer.external_actors[key] = actor
                material = actor._siren_material
                self.renderer.registry.materialDict[material] = SimpleNamespace(density=actor._siren_density)
                self.renderer.options[material] = VisualisationOptions(
                    colour=list(actor.GetProperty().GetColor()), alpha=actor.GetProperty().GetOpacity())
                if actor not in self.renderer.material_actors[material]:
                    self.renderer.material_actors[material].append(actor)
        bounds = self.vis._scene_bounds(v)
        if self.axes and not self.axes_added and bounds is not None:
            span = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
            v.addAxes(length=.12 * span, origin=[(bounds[i] + bounds[i+1]) / 2 for i in (0, 2, 4)])
            for getter in (v.axes[-1].GetXAxisCaptionActor2D,
                           v.axes[-1].GetYAxisCaptionActor2D,
                           v.axes[-1].GetZAxisCaptionActor2D):
                caption = getter()
                caption.GetTextActor().SetTextScaleModeToNone()
                caption.GetCaptionTextProperty().SetFontSize(12)
                caption.GetCaptionTextProperty().SetColor(.15, .15, .15)
            self.axes_added = True
        if self.controls is None:
            # coloured=False is the plain single-colour viewer: no material legend.
            self.controls = self.vis._install_controls(
                self.vtk, v, self.renderer.registry,
                self.renderer.options if self.coloured else None,
                self.model if self.picker else None, bounds, legend=self.legend,
                bounding_box=self.bounding_box, near_frac=self.near_frac,
                gas_visible=self.show_gas)
            v.interactorStyle.AddObserver('StartInteractionEvent', self._camera_changed)
        if v.bClipper and v.clipperPlaneWidget is None and v.clippers:
            v.addClipperWidget()
            v.clipperPlaneWidget.On()
            # The controls' state and key handler share this mutable handle.
            self.controls["clipper_widget"] = v.clipperPlaneWidget

    def receive(self, event):
        if self._closing():
            return
        kind = event["kind"]
        if kind == "structure":
            self._mesh_files.clear()
            self.renderer.set_structure(event["scene"])
            self.status.SetInput("Loading %s geometry..." % event["scene"]["stage"])
        elif kind == "mesh":
            before = time.perf_counter()
            path = event['path']
            if path not in self._mesh_files:
                self._mesh_files[path] = read_mesh(path)
            self.renderer.add_mesh(event["prototype"], *self._mesh_files[path])
            self.metrics["scene_build_seconds"] = self.metrics.get("scene_build_seconds", 0.) + time.perf_counter() - before
            if event["stage"] == "detail":
                self.loaded.add(event["prototype"])
            self.status.SetInput("%s: %d / %d shapes" % (event["stage"].capitalize(), event["completed"], event["total"]))
            if not self.camera_fitted and self.viewer.actors:
                # Input on the loading window already placed the camera; keep it.
                if not self.camera_dirty:
                    self.viewer.ren.ResetCamera()
                self.camera_fitted = True
        elif kind == "stage_done":
            before = time.perf_counter()
            self.renderer.finish_stage()
            self.metrics["scene_build_seconds"] = self.metrics.get("scene_build_seconds", 0.) + time.perf_counter() - before
            self.metrics["stages"].append(event)
            self._setup_controls()
            if not self.camera_dirty and not self.detail_ready:
                self.viewer.ren.ResetCamera()
            if event["stage"] == "detail":
                self.detail_ready = True
                self.status.SetInput("Ready")
                self.status.SetVisibility(False)
            self._render(event["stage"])
        elif kind == "warning":
            # Worker stderr goes to its log; cache/input warnings reach the caller here.
            import warnings
            self.metrics.setdefault("warnings", []).append(event["message"])
            warnings.warn(event["message"], RuntimeWarning)
        elif kind == "failed":
            if self.detail_ready:
                # An on-demand ('g'/'v') job failed after the scene was ready:
                # keep the loaded viewer open and report it instead of exiting.
                import warnings
                self.metrics.setdefault("request_errors", []).append(event["message"])
                warnings.warn("requested volumes failed to load: %s" % event["message"],
                              RuntimeWarning)
                self.status.SetInput("Requested volumes failed to load; see the Python warning.")
                self.status.SetVisibility(True)
                self._render()
                return
            self.error = event["message"]
            self.status.SetInput("Geometry loading failed; see the Python error.")
            self.status.SetVisibility(True)
            self._render()
            self.viewer.iren.TerminateApp()
        # A coarse-resolution Boolean failure is reported in stage statistics;
        # that prototype is still required to succeed at the final resolution.

    def _start_job(self, only=None):
        if self._closing():
            return
        options = dict(self.options)
        if only is not None:
            options.update(only=sorted(only), show_gas=True, hidden_volumes=[])
        self.job = MeshJob(self.path, options,
                           self.preview_slices if only is None and self.progressive else None)

    def request_geometry(self, group):
        if self._closing() or self.renderer.scene is None:
            return
        prototypes = self.renderer.scene["prototypes"]
        hidden_roles = (set(self.options["regions"]["hidden"])
                        if self.options["display"] == "exterior" and self.options["regions"] else set())
        needed = {name for name, p in prototypes.items()
                  if (group == "all" or p["density"] <= _LOW_DENSITY) and
                  p["role"] not in hidden_roles}
        self.pending.update(needed - self.loaded)
        if not self.pending and self.job is None:
            self.status.SetVisibility(False)
            return
        self.status.SetVisibility(self.progress)
        self.status.SetInput("Loading requested volumes...")
        if self.job is None and self.pending:
            only, self.pending = self.pending, set()
            self._start_job(only)

    def tick(self, *args):
        if self._closing() or self.job is None:
            return
        try:
            self._events.extend(self.job.poll())
            deadline = time.perf_counter() + .05
            processed = False
            # Leave time for native close/input events even when a fast worker
            # has queued thousands of meshes. Keep its files until consumed.
            # A trailing 'done' carries no work; consume it now so a close in
            # the next timer slot does not report a finished load as cancelled.
            while self._events and (not processed or time.perf_counter() < deadline or
                                    self._events[0]["kind"] == "done"):
                if self._closing():
                    break
                self.receive(self._events.popleft())
                processed = True
            if processed and not self.error and not self._closing():
                self._render()
            if self.job.finished and not self._events:
                self.job.close()
                self.job = None
                self.pending.difference_update(self.loaded)
                if self.pending and not self.error and not self._closing():
                    only, self.pending = self.pending, set()
                    self._start_job(only)
        except Exception as exc:
            self.receive(dict(kind="failed", message=str(exc)))

    def run(self):
        try:
            if self.progressive:
                self._start_job()
            else:
                with tempfile.TemporaryDirectory(prefix="siren-scene-") as directory:
                    prepare_scene(self.path, directory, emit=self.receive, **self.options)
            if self.interactive:
                self._timer_observer = self.viewer.iren.AddObserver("TimerEvent", self.tick)
                self.timer = self.viewer.iren.CreateRepeatingTimer(100)
                self.viewer.iren.Start()
            if self.error:
                raise RuntimeError(self.error)
            if self.screenshot:
                self._render()
                capture = self.vtk.vtkWindowToImageFilter()
                capture.SetInput(self.viewer.renWin)
                capture.ReadFrontBufferOff()
                capture.Update()
                writer = self.vtk.vtkPNGWriter()
                writer.SetFileName(str(self.screenshot))
                writer.SetInputConnection(capture.GetOutputPort())
                writer.Write()
            return self.screenshot or self.path
        except Exception as exc:
            self.error = str(exc)
            self.metrics['error'] = self.error
            raise
        finally:
            self.closed = True
            unfinished = bool(self.pending or any(e["kind"] != "done" for e in self._events) or
                              (self.job is not None and not self.job.finished))
            if self.timer is not None:
                self.viewer.iren.DestroyTimer(self.timer)
                self.timer = None
            if self._timer_observer is not None:
                self.viewer.iren.RemoveObserver(self._timer_observer)
                self._timer_observer = None
            if self.job is not None:
                self.job.close()
                self.job = None
            self.pending.clear()
            self._events.clear()
            self._mesh_files.clear()
            self.metrics["cancelled"] = bool((not self.detail_ready or unfinished) and not self.error)
            self._publish_metrics()
            if self.timings is True:
                import json
                print("[siren.view] " + json.dumps(self.metrics, sort_keys=True))
            self._release_window()
