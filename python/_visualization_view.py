"""Progressive VTK window controller. All VTK calls run on the caller's thread."""
from collections import deque
import os
from pathlib import Path
import re
import tempfile
import time

from ._visualization_job import MeshJob
from ._visualization_panel import Panel
from ._visualization_presets import BACKGROUND, OTHER, load_presets, normalize_section
from ._visualization_render import SceneRenderer
from ._visualization_scene import normalize_regions, prepare_scene, read_mesh
from .visualization import _GDML_MM_PER_M, _LOW_DENSITY


class ViewSession:
    def __init__(self, model, gdml_path=None, *, screenshot=None, coloured=True,
                 axes=True, cutter=False, clipper=False, clip_origin=(0., 0., 0.),
                 clip_normal=(1., 0., 0.), legend=True, bounding_box=False,
                 picker=True, interactive=True, mesh_slices=48, near_frac=1e-3,
                 cache=True, cache_dir=None, progressive=True, preview_slices=12,
                 instancing=True, display="full", regions=None, show_gas=False,
                 hidden_volumes=(), timings=None, progress=True, presets=None,
                 section=None, panel=True, edges=True):
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
        self.presets = load_presets(presets)
        hidden_volumes = list(hidden_volumes) + self.presets["hidden_volumes"]
        if section is not None and not isinstance(section, str):
            section = normalize_section(section)
            self.presets["sections"] = [r for r in self.presets["sections"]
                                        if r["label"] != section["label"]] + [section]
            section = section["label"]
        if section is not None and section not in [r["label"] for r in self.presets["sections"]]:
            raise ValueError("no preset section called %r" % section)
        if clipper and section is not None:
            raise ValueError("give either section or clipper, not both")
        self.first_section = section
        self.first_plane = (list(clip_origin), list(clip_normal)) if clipper else None
        self.current_section, self.slider_axis, self.slider_value = None, "x", 0.
        self.plane_widget = None
        cutaway = bool(cutter or clipper or section)
        self.options = dict(mesh_slices=mesh_slices, cache=bool(cache),
                            cache_dir=None if cache_dir is None else str(cache_dir),
                            display="full" if cutaway else display,
                            regions=regions, show_gas=show_gas,
                            hidden_volumes=list(hidden_volumes),
                            overlapping=self.model is not None)
        self.metrics["display"] = self.options["display"]
        self.metrics["cutaway_full_detail"] = bool(cutaway and display != "full")
        import vtk
        from pyg4ometry.visualisation import VtkViewerNew
        from . import visualization as vis
        self.vtk = vtk
        self.vis = vis
        self.metrics["imports_seconds"] = time.perf_counter() - self.started
        self.viewer = VtkViewerNew(defaultCutters=bool(cutter), axisCubeWidget=bool(axes))
        self.viewer.ren.SetBackground(*BACKGROUND)
        self.viewer.renWin.SetMultiSamples(8)
        self._set_lights()
        # pyg4ometry's default style has a crash-prone right-click picker; a
        # plain trackball style guards the loading phase until _install_controls
        # replaces it. Each style carries its own camera observer while live.
        style = vtk.vtkInteractorStyleTrackballCamera()
        style.AddObserver('StartInteractionEvent', self._camera_changed)
        self.viewer.iren.SetInteractorStyle(style)
        self.coloured = bool(coloured)
        self.renderer = SceneRenderer(self.viewer, coloured=coloured, instancing=instancing,
                                      display=self.options["display"], regions=regions,
                                      presets=self.presets, edges=edges)
        self.presets = self.renderer.presets
        self.viewer.layer_visible = lambda actor: self.renderer.layers.get(
            self.renderer.actor_layer.get(actor), True)
        self.viewer.exact_shapes = self._exact_shapes
        self.panel = None
        if panel and self.interactive and coloured and self.options["display"] == "full":
            self.viewer.renWin.SetSize(1400, 900)
            if getattr(self.viewer, "axesWidget", None) is not None:
                self.viewer.axesWidget.SetViewport(.88, 0., 1., .18)
            self.panel = Panel(vtk, self.viewer.renWin, self)
            self.panel.place(self.viewer.ren)
            self.viewer.renWin.AddObserver("ModifiedEvent", self._resized)
            self._size = None
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
                stats = vis.to_gdml(model, self.path, skip_geo_types=("TriangularMesh",))
                self.boxed = set(stats["bbox_sectors"])
        except BaseException:
            # run() never starts, so its window teardown must happen here.
            self.closed = True
            self._release_window()
            raise
        self.metrics["export_seconds"] = time.perf_counter() - before

    boxed = frozenset()

    def _sector_index(self, row):
        """Index of the sector behind an instance, from the name to_gdml gave it."""
        found = re.search(r"__(\d+)$", row["prototype"])
        return int(found.group(1)) if found else None

    def _exact_shapes(self, instances):
        """Shapes for cutting the sectors that the GDML holds only as boxes."""
        from functools import partial
        return {k: partial(self._exact_shape, i) for k, row in enumerate(instances)
                for i in [self._sector_index(row)] if i in self.boxed}

    def _exact_shape(self, index):
        shapes = self.__dict__.setdefault("_shapes", {})
        if index not in shapes:
            sector = self.vis._sectors(self.model)[index]
            shapes[index] = self.vis._exact_shape(sector.geo)
            if shapes[index] is None:
                import warnings
                warnings.warn("sector %s has no exact shape to cut and is left out of sections "
                              "(is pyvista installed?)" % sector.name, RuntimeWarning)
        return shapes[index]

    def _set_levels(self, scene):
        """Sector levels decide which of a model's overlapping sectors is cut."""
        if self.model is None:
            return
        levels = [s.level for s in self.vis._sectors(self.model)]
        for row in scene["instances"]:
            index = self._sector_index(row)
            if index is not None and index < len(levels):
                row["level"] = levels[index]

    # ---- views, sections and the panel ----------------------------------
    def _set_lights(self):
        """Soft lights that move with the camera, so every face is readable."""
        ren = self.viewer.ren
        ren.AutomaticLightCreationOff()
        ren.RemoveAllLights()
        for position, intensity in (((.5, 1., 1.), .75), ((-1., .2, .6), .35)):
            light = self.vtk.vtkLight()
            light.SetLightTypeToCameraLight()
            light.SetPosition(*position)
            light.SetFocalPoint(0, 0, 0)
            light.SetIntensity(intensity)
            ren.AddLight(light)

    def _resized(self, *args):
        size = self.viewer.renWin.GetSize()
        if self.panel is None or size == self._size or self._closing():
            return
        self._size = size
        self.panel.place(self.viewer.ren)
        self.panel.rebuild()

    def panel_changed(self):
        if self.panel is not None:
            self.panel.rebuild()

    def toggle_panel(self):
        if self.panel is not None:
            self.panel.visible = not self.panel.visible
            self.panel.place(self.viewer.ren)
            self.panel.rebuild()
            self.viewer.renWin.Render()

    def _aim(self, centre, direction, up, width, height, angle):
        """Look at *centre* along -direction so that width x height (mm) fits."""
        import math
        ren = self.viewer.ren
        # The renderer's own size is not known before its first frame.
        size, port = self.viewer.renWin.GetSize(), ren.GetViewport()
        aspect = size[0] * (port[2] - port[0]) / max(size[1] * (port[3] - port[1]), 1)
        half = math.tan(math.radians(angle) / 2)
        distance = max(height, width / max(aspect, 1e-6)) / 2 / half
        length = math.sqrt(sum(c * c for c in direction))
        camera = ren.GetActiveCamera()
        camera.SetViewAngle(angle)
        camera.SetFocalPoint(*centre)
        camera.SetPosition(*(c + distance * d / length for c, d in zip(centre, direction)))
        camera.SetViewUp(*up)
        ren.ResetCameraClippingRange()

    def go_home(self):
        home = self.presets["home"]
        if home is not None:
            self._aim([c * _GDML_MM_PER_M for c in home["centre"]], home["direction"], self.presets["up"],
                      home["width"] * _GDML_MM_PER_M, home["height"] * _GDML_MM_PER_M, 32.)
            return
        self.viewer.ren.GetActiveCamera().SetViewAngle(30.)
        self.viewer.ren.ResetCamera()

    def show_section(self, label, aim=True):
        """Cut along a preset section, or show the whole model for None."""
        self.current_section = label
        self._drop_plane_widget()
        row = next((r for r in self.presets["sections"] if r["label"] == label), None)
        if row is None:
            self.renderer.clear_section()
            if aim:
                self.go_home()
        else:
            origin = [c * _GDML_MM_PER_M for c in row["origin"]]
            up = row["up"] if row["up"] is not None else self.presets["up"]
            self.renderer.set_section(origin, row["normal"], up)
            if aim:
                self._aim_section(row, origin)
        self._after_cut()

    def _aim_section(self, row, origin):
        import numpy as np
        section = self.renderer.section
        centre = [c * _GDML_MM_PER_M for c in row["centre"]] if row["centre"] is not None else None
        width, height = ((row[k] * _GDML_MM_PER_M if row[k] is not None else None)
                         for k in ("width", "height"))
        pieces = list(section.pieces.values()) if section is not None else []
        if pieces and None in (centre, width, height):
            points = np.concatenate([q.reshape(-1, 2) for q in pieces])
            low, high = points.min(axis=0), points.max(axis=0)
            if centre is None:
                centre = section.to_world((low + high) / 2).tolist()
            width, height = width or 1.1 * (high[0] - low[0]), height or 1.1 * (high[1] - low[1])
        if section is None or None in (width, height):
            return
        self._aim(centre or origin, section.normal.tolist(), section.w.tolist(), width, height, 9.)

    def _after_cut(self):
        if self.controls is not None:
            self.controls["refresh"]()
        self.panel_changed()
        if not self._closing():
            self.viewer.renWin.Render()

    def slider_range(self):
        """Reach of the slider: the part of the model in view, so zooming refines it."""
        import math
        bounds = self.vis._scene_bounds(self.viewer)
        if bounds is None:
            return -1., 1.
        i = "xyz".index(self.slider_axis)
        camera = self.viewer.ren.GetActiveCamera()
        reach = 2 * camera.GetDistance() * math.tan(math.radians(camera.GetViewAngle()) / 2)
        focus = camera.GetFocalPoint()[i]
        low, high = max(bounds[2 * i], focus - reach), min(bounds[2 * i + 1], focus + reach)
        if low >= high:
            low, high = bounds[2 * i], bounds[2 * i + 1]
        return low / _GDML_MM_PER_M, high / _GDML_MM_PER_M

    def view_moved(self):
        self.panel_changed()

    def cut_axis(self, axis):
        """Cut across one axis, removing the side the camera is on."""
        self.slider_axis = axis
        i = "xyz".index(axis)
        low, high = self.slider_range()
        self.slider_value = min(max(self.slider_value, low), high)
        camera = self.viewer.ren.GetActiveCamera()
        focus, distance = camera.GetFocalPoint(), camera.GetDistance()
        view = [(p - f) / distance for p, f in zip(camera.GetPosition(), focus)]
        self.slider_sign = 1. if view[i] >= 0 else -1.
        if abs(view[i]) < .35:
            # Nearly edge-on: turn until the cut face can be seen.
            view = [.35 * v for v in view]
            view[i] = .94 * self.slider_sign
            length = sum(v * v for v in view) ** .5
            up = self.presets["up"]
            if abs(sum(u * v for u, v in zip(up, view))) > .9 * length * sum(u * u for u in up) ** .5:
                up = camera.GetViewUp()
            camera.SetPosition(*(f + distance * v / length for f, v in zip(focus, view)))
            camera.SetViewUp(*up)
            self.viewer.ren.ResetCameraClippingRange()
        self.slide(self.slider_value, True)

    slider_sign = 1.

    def slide(self, value, done=False):
        self.current_section, self.slider_value = "", value
        self._drop_plane_widget()
        i = "xyz".index(self.slider_axis)
        origin, normal = [0.] * 3, [0.] * 3
        origin[i], normal[i] = value * _GDML_MM_PER_M, self.slider_sign
        self.renderer.set_section(origin, normal, self.presets["up"])
        self._after_cut()

    def flip_section(self):
        """Keep the other half of the model at the current plane."""
        if self.renderer.plane is None:
            return
        self.slider_sign = -self.slider_sign
        self.renderer.set_section(self.renderer.section_origin,
                                  [-c for c in self.renderer.section_normal], self.renderer._section_up)
        self._after_cut()

    def used_layers(self):
        return {row.get("layer", OTHER) for row in self.renderer.scene["instances"]
                if row["prototype"] in self.renderer._drawn} if self.renderer.scene else set()

    def toggle_layer(self, label):
        self.renderer.layers[label] = not self.renderer.layers.get(label, True)
        self._after_cut()

    def show_all_layers(self):
        for label in self.renderer.layers:
            self.renderer.layers[label] = True
        self.panel_changed()

    def material_rows(self):
        """(label, colour, shown) of the drawn materials, densest first."""
        renderer, rows = self.renderer, {}
        if not renderer.scene or self.controls is None:
            return []
        hidden = self.controls["hidden"]
        for material, actors in renderer.material_actors.items():
            if not actors or material not in renderer.options:
                continue
            density = float(renderer.registry.materialDict[material].density)
            if density <= _LOW_DENSITY:
                continue
            row = rows.setdefault(renderer.labels.get(material, material),
                                  [density, renderer.options[material].colour, False])
            row[0] = max(row[0], density)
            row[2] = row[2] or material not in hidden
        return [(label, tuple(colour), on) for label, (_, colour, on)
                in sorted(rows.items(), key=lambda item: -item[1][0])]

    def toggle_material(self, label):
        hidden = self.controls["hidden"]
        names = {m for m in self.renderer.material_actors
                 if self.renderer.labels.get(m, m) == label}
        if names <= hidden:
            hidden -= names
        else:
            hidden |= names
        self._after_cut()

    def toggle_plane_widget(self):
        """Show a plane to drag and turn; its cut follows it."""
        if self.plane_widget is not None:
            self._drop_plane_widget()
            return
        bounds = self.vis._scene_bounds(self.viewer)
        if bounds is None:
            return
        if self.renderer.plane is None:
            camera = self.viewer.ren.GetActiveCamera()
            middle = [(bounds[2 * i] + bounds[2 * i + 1]) / 2 for i in range(3)]
            normal = [p - f for p, f in zip(camera.GetPosition(), camera.GetFocalPoint())]
            self._widget_cut(middle, normal, True)
        vtk = self.vtk
        representation = vtk.vtkImplicitPlaneRepresentation()
        representation.SetPlaceFactor(1.)
        representation.PlaceWidget(bounds)
        representation.SetOrigin(*self.renderer.section_origin)
        representation.SetNormal(*self.renderer.section_normal)
        representation.DrawPlaneOff()
        representation.OutlineTranslationOff()
        representation.ScaleEnabledOff()
        widget = vtk.vtkImplicitPlaneWidget2()
        widget.SetInteractor(self.viewer.iren)
        widget.SetDefaultRenderer(self.viewer.ren)
        widget.SetRepresentation(representation)
        def moved(*args):
            self._widget_cut(representation.GetOrigin(), representation.GetNormal(), False)
        def settled(*args):
            self._widget_cut(representation.GetOrigin(), representation.GetNormal(), True)
            self.viewer.renWin.Render()
        widget.AddObserver("InteractionEvent", moved)
        widget.AddObserver("EndInteractionEvent", settled)
        widget.On()
        self.plane_widget = widget
        if self.controls is not None:
            self.controls["clipper_widget"] = widget

    def _widget_cut(self, origin, normal, caps):
        self.current_section = ""
        self.renderer.set_section(list(origin), list(normal), self.presets["up"], caps=caps)
        if caps:
            self._after_cut()

    def _drop_plane_widget(self):
        widget, self.plane_widget = self.plane_widget, None
        if widget is not None:
            widget.Off()
            widget.RemoveAllObservers()
            widget.SetInteractor(None)
        if self.controls is not None:
            self.controls["clipper_widget"] = None

    def _open_view(self):
        """The first view, once there is geometry to look at."""
        if self.first_plane is not None:
            (origin, normal), self.first_plane = self.first_plane, None
            self._widget_cut(origin, normal, True)
            if self.interactive:
                self.toggle_plane_widget()
        if self.first_section is not None:
            label, self.first_section = self.first_section, None
            self.show_section(label, aim=not self.camera_dirty)
        elif self.current_section is None and not self.camera_dirty:
            self.go_home()

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
        self._drop_plane_widget()
        v.renWin.RemoveAllObservers()
        for name in ("axesWidget",):
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
        if self.axes and self.panel is None and not self.axes_added and bounds is not None:
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
                self.model if self.picker else None, bounds,
                legend=self.legend and self.panel is None,
                bounding_box=self.bounding_box, near_frac=self.near_frac,
                gas_visible=self.show_gas, session=self)
            v.interactorStyle.AddObserver('StartInteractionEvent', self._camera_changed)

    def receive(self, event):
        if self._closing():
            return
        kind = event["kind"]
        if kind == "structure":
            self._mesh_files.clear()
            self._set_levels(event["scene"])
            self.renderer.set_structure(event["scene"])
            self.status.SetInput("Loading %s geometry..." % event["scene"]["stage"])
        elif kind == "mesh":
            before = time.perf_counter()
            path = event['path']
            if path not in self._mesh_files:
                self._mesh_files[path] = read_mesh(path)
            self.renderer.add_mesh(event["prototype"], *self._mesh_files[path])
            self.metrics["scene_build_seconds"] = self.metrics.get("scene_build_seconds", 0.) + time.perf_counter() - before
            if event["stage"] == "detail" and event["prototype"] in self.renderer._drawn:
                self.loaded.add(event["prototype"])
            self.status.SetInput("%s: %d / %d shapes" % (event["stage"].capitalize(), event["completed"], event["total"]))
            if not self.camera_fitted and self.viewer.actors:
                # Input on the loading window already placed the camera; keep it.
                if not self.camera_dirty:
                    self.go_home()
                self.camera_fitted = True
        elif kind == "stage_done":
            before = time.perf_counter()
            self.renderer.finish_stage()
            self.metrics["scene_build_seconds"] = self.metrics.get("scene_build_seconds", 0.) + time.perf_counter() - before
            self.metrics["stages"].append(event)
            self._setup_controls()
            if not self.detail_ready:
                self._open_view()
            self.panel_changed()
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
