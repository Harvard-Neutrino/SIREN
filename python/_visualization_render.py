"""VTK prototype rendering and placement-aware picking (UI thread only)."""
from collections import defaultdict

from ._visualization_presets import (BACKGROUND, INK, assign_layers, assign_volumes, clean_name,
                                     load_presets, material_row)


def polydata(vertices, faces):
    import numpy as np
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(vertices, deep=True))
    cells = vtk.vtkCellArray()
    flat = np.column_stack((np.full(len(faces), 3, dtype=np.int64), faces)).ravel()
    cells.ImportLegacyFormat(numpy_to_vtkIdTypeArray(flat, deep=True))
    result = vtk.vtkPolyData()
    result.SetPoints(points)
    result.SetPolys(cells)
    return result


def glyph_components(matrix):
    """Return quaternion/scale for orthogonal positive transforms, else None.

    Reflections and shear retain their exact full matrix through a shared-mapper
    actor. Nested nonuniform scaling can produce shear even from GDML scales.
    """
    import numpy as np
    import vtk
    linear = np.asarray(matrix, dtype=float)[:3, :3]
    scale = np.linalg.norm(linear, axis=0)
    if not np.all(scale > 0) or not np.isfinite(linear).all():
        return None
    rotation = linear / scale
    if np.linalg.det(rotation) <= 0 or not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-12, rtol=0):
        return None
    quaternion = [0., 0., 0., 0.]
    vtk.vtkMath.Matrix3x3ToQuaternion(rotation.tolist(), quaternion)
    return quaternion, scale.tolist()


def feature_lines(vertices, faces, angle=28.):
    """Point pairs along edges where the surface turns by more than *angle* degrees."""
    import numpy as np
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    merge = vtk.vtkCleanPolyData()
    merge.SetInputData(polydata(vertices, faces))
    edges = vtk.vtkFeatureEdges()
    edges.SetInputConnection(merge.GetOutputPort())
    edges.BoundaryEdgesOff()
    edges.NonManifoldEdgesOff()
    edges.ManifoldEdgesOff()
    edges.FeatureEdgesOn()
    edges.ColoringOff()
    edges.SetFeatureAngle(angle)
    edges.Update()
    result = edges.GetOutput()
    if not result.GetNumberOfLines():
        return np.empty((0, 2, 3))
    ids = vtk_to_numpy(result.GetLines().GetConnectivityArray())
    return vtk_to_numpy(result.GetPoints().GetData()).astype(float)[ids].reshape(-1, 2, 3)


def line_polydata(pairs):
    import numpy as np
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    pairs = np.asarray(pairs, dtype=float)
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(pairs.reshape(-1, 3), deep=True))
    cells = vtk.vtkCellArray()
    ids = np.arange(2 * len(pairs), dtype=np.int64).reshape(-1, 2)
    flat = np.column_stack((np.full(len(ids), 2, dtype=np.int64), ids)).ravel()
    cells.ImportLegacyFormat(numpy_to_vtkIdTypeArray(flat, deep=True))
    result = vtk.vtkPolyData()
    result.SetPoints(points)
    result.SetLines(cells)
    return result


# Edge lines are skipped for a volume whose copies would draw more than this.
_MAX_EDGE_PAIRS = 20000
# Volumes smaller than this part of the scene get none; they would be a dark blot.
_MIN_EDGE_SIZE = .01


class SceneRenderer:
    def __init__(self, viewer, *, colored=True, instancing=True, display="full", regions=None,
                 presets=None, edges=True):
        self.viewer = viewer
        self.colored = colored
        self.instancing = instancing
        self.display = display
        # Normalized surface-region mapping (see _visualization_scene.normalize_regions).
        self.regions = regions
        self.presets = load_presets(presets)
        self.edges = edges
        self.scene = None
        self.layers = {row["label"]: row["shown"] for row in self.presets["layers"]}
        # Color, opacity ("solid" or "see-through") and shown, by kind and key:
        # ("material", label), ("layer", label), ("volume", rule label or "lv:" + logical name).
        self.styles = {}
        for kind in ("materials", "layers", "volumes"):
            for row in self.presets[kind]:
                self.styles[kind[:-1], row["label"]] = {k: row[k] for k in ("color", "opacity", "shown")}
        self.color_mode, self.surface = self.presets["color_mode"], self.presets["surfaces"]
        self.opacity_override = None
        self._separate = set()  # logical volumes styled on their own
        self.actor_layer, self.actor_material, self.actor_volume, self.actor_prototype = {}, {}, {}, {}
        self.line_actors, self.cap_actors = [], []
        self._edge_pairs = {}
        self.plane, self.section, self._sections = None, None, None
        self._section_up = None
        self.meshes = {}
        self.placements = defaultdict(list)
        self.prototype_actors = {}
        self.actor_instances = {}
        self._actor_is_glyph = {}
        self.external_actors = {}
        self._drawn = set()
        self.locators, self.inverses = {}, {}
        self.options = {}
        self._groups = defaultdict(list)
        self._group_representatives = {}
        self._group_loaded = {}
        self._source_data = {}
        self._source_users = defaultdict(int)  # groups drawing each polydata
        self.material_actors = defaultdict(list)
        self.viewer.material_actors = self.material_actors
        self.viewer.pick_scene = self.pick
        self._stats = dict(prototype_faces=0, placed_faces=0, glyph_instances=0, matrix_instances=0)
        self._stats_dirty = False

    @property
    def stats(self):
        # Reporting once per frame is linear in scene size. Recounting every
        # actor after every mesh made large flattened detector exports quadratic.
        if self._stats_dirty:
            self.update_stats()
        return self._stats

    def set_structure(self, scene):
        import numpy as np
        from types import SimpleNamespace
        from pyg4ometry.visualisation import VisualisationOptions
        from .visualization import _material_vis_options
        self.scene = scene
        self._stats_dirty = True
        self._groups.clear()
        self._source_data.clear()
        self._source_users.clear()
        selected = set(scene.get('selected', scene['prototypes']))
        active = selected | self._drawn | set(scene.get('holes', ()))
        self._drawn |= selected
        for name, p in scene['prototypes'].items():
            if name in active:
                self._groups[(p.get('mesh_key', name), p['material'], p['role'])].append(name)
        self.placements.clear()
        materials = {}
        self._sections = None
        for instance, layer, volume in zip(scene["instances"], assign_layers(scene, self.presets["layers"]),
                                           assign_volumes(scene, self.presets["volumes"])):
            instance["layer"] = layer
            instance["volume"] = ("lv:" + instance["prototype"]
                                  if instance["prototype"] in self._separate else volume)
            self.placements[instance["prototype"]].append(instance)
            self.inverses[instance["name"]] = np.linalg.inv(instance["matrix"])
        for p in scene["prototypes"].values():
            materials[p["material"]] = SimpleNamespace(density=p["density"])
        self.registry = SimpleNamespace(materialDict=materials)
        # Density colors; the caller's colors come through self.styles.
        self.options = _material_vis_options(self.registry, VisualisationOptions,
                                             (), self.presets["opacity"])
        rows = {name: material_row(name, self.presets["materials"]) for name in materials}
        self.labels = {name: row["label"] if row else clean_name(name) for name, row in rows.items()}
        legend_options, legend_materials = {}, {}
        for name in {p["material"] for k, p in scene["prototypes"].items() if k in self._drawn}:
            label, density = self.labels[name], materials[name].density
            if label not in legend_materials or density > legend_materials[label].density:
                legend_options[label] = VisualisationOptions(colour=list(self.material_color(name)),
                                                             alpha=self.options[name].alpha)
                legend_materials[label] = materials[name]
        styles = self.regions['styles'] if self.display == 'exterior' and self.regions else {}
        if styles:
            for material in materials:
                roles = {p['role'] for p in scene['prototypes'].values() if p['material'] == material}
                if roles and roles <= styles.keys():
                    legend_options.pop(self.labels[material], None)
            for name, p in scene['prototypes'].items():
                if p['role'] in styles:
                    label = styles[p['role']]['label']
                    legend_options[label] = self.vis_option(name)
                    legend_materials[label] = SimpleNamespace(density=p['density'])
        self.viewer.legend_options = legend_options
        self.viewer.legend_registry = SimpleNamespace(materialDict=legend_materials)

    def vis_option(self, name):
        import copy
        from pyg4ometry.visualisation import VisualisationOptions
        p = self.scene["prototypes"][name]
        vo = copy.copy(self.options[p["material"]]) if self.colored else VisualisationOptions(colour=[.6, .6, .6])
        if self.display == "exterior" and self.regions and p["role"] in self.regions["styles"]:
            # Caller-supplied display colors, not measured optical coefficients.
            style = self.regions["styles"][p["role"]]
            if self.colored:
                vo.colour = list(style["color"])
            vo.alpha = style["alpha"]
        return vo

    # ---- colors and surfaces -------------------------------------------
    def material_color(self, material):
        """A material's color under the current color choice."""
        own = self.styles.get(("material", self.labels.get(material, material)))
        if self.color_mode == "custom" and own and own.get("color") is not None:
            return tuple(own["color"])
        return tuple(self.options[material].colour)

    def style(self, material, layer=None, volume=None, prototype=None):
        """(color, opacity): a volume's own choice, else its layer's, else its material's."""
        from .visualization import _LOW_DENSITY
        if not self.colored or (prototype is not None and self.display == "exterior" and self.regions
                                 and self.scene["prototypes"][prototype]["role"] in self.regions["styles"]):
            vo = self.vis_option(prototype)
            return tuple(vo.colour), (vo.alpha if self.opacity_override is None else self.opacity_override)
        chain = [s for s in (self.styles.get(("volume", volume)) if volume else None,
                             self.styles.get(("layer", layer)),
                             self.styles.get(("material", self.labels.get(material, material)))) if s]
        color = self.material_color(material)
        if self.color_mode == "custom":
            color = next((tuple(s["color"]) for s in chain if s.get("color") is not None), color)
        if self.opacity_override is not None:
            return color, self.opacity_override
        if self.registry.materialDict[material].density <= _LOW_DENSITY:
            return color, self.options[material].alpha
        surface = next((s["opacity"] for s in chain if s.get("opacity")), self.surface)
        return color, self.presets["opacity"] if surface == "solid" else self.presets["see_through"]

    def actor_style(self, actor):
        return self.style(self.actor_material[actor], self.actor_layer.get(actor),
                          self.actor_volume.get(actor), self.actor_prototype.get(actor))

    def apply_styles(self):
        """Color and opacity of every surface, cut face and edge line."""
        keep = 1. if self.plane is None else .5  # surfaces behind a cut are paler
        for actor in self.actor_instances:
            color, alpha = self.actor_style(actor)
            color = tuple(keep * c + (1 - keep) * b for c, b in zip(color, BACKGROUND))
            prop = actor.GetProperty()
            if max(abs(a - b) for a, b in zip(prop.GetColor(), color)) > 1e-9:
                prop.SetColor(*color)
                if self._actor_is_glyph[actor]:
                    actor.GetMapper().Modified()  # glyph colors are cached per instance
            prop.SetOpacity(alpha)
        for actor in self.cap_actors:
            color, alpha = self.actor_style(actor)
            actor.GetProperty().SetColor(*color)
            actor.GetProperty().SetOpacity(alpha)
        for actor in self.line_actors:
            actor.GetProperty().SetOpacity(.38 if self.actor_style(actor)[1] >= .5 else .12)

    def set_style(self, kind, key, **change):
        """Change the color, opacity or shown of a material, layer or logical volume."""
        if kind == "volume":
            return self.set_volume(key, **change)
        self.styles.setdefault((kind, key), {}).update(change)
        self.apply_styles()

    def set_volume(self, prototype, **change):
        """Style one logical volume (all its copies) apart from its neighbors."""
        key = "lv:" + prototype
        if prototype not in self._separate:
            rows = self.placements[prototype]
            first = rows[0].get("volume") if rows else None
            self.styles["volume", key] = dict(self.styles.get(("volume", first), {}))
            self._separate.add(prototype)
            for row in rows:
                row["volume"] = key
            p = self.scene["prototypes"][prototype]
            group = (p.get("mesh_key", prototype), p["material"], p["role"])
            name = self._group_representatives.get(group)
            if name is not None and name in self.meshes:
                self.add_mesh(name, *self.meshes[name], force=True)
            self._build_edges()
            if self.plane is not None:
                self.build_caps()
        self.styles["volume", key].update(change)
        self.apply_styles()
        self.refresh_overlays()

    def volume_style(self, prototype):
        """(color, see-through, shown) of a logical volume as drawn."""
        rows = self.placements.get(prototype)
        if not rows:
            return None, False, True
        row = rows[0]
        p = self.scene["prototypes"][prototype]
        color, alpha = self.style(p["material"], row["layer"], row.get("volume"))
        own = self.styles.get(("volume", row.get("volume")), {})
        return color, alpha < self.presets["opacity"] - 1e-6, own.get("shown") is not False

    def _remember(self, actor, name, instances, glyph=False):
        prop = actor.GetProperty()
        if self.colored:
            prop.SetInterpolationToFlat()
            prop.SetAmbient(.5)
            prop.SetDiffuse(.55)
        if self.plane is not None:
            actor.GetMapper().AddClippingPlane(self.plane)
        self.actor_layer[actor] = instances[0].get("layer")
        self.actor_volume[actor] = instances[0].get("volume")
        self.actor_prototype[actor] = name
        self.actor_material[actor] = self.scene["prototypes"][name]["material"]
        color, alpha = self.actor_style(actor)
        keep = 1. if self.plane is None else .5
        prop.SetColor(*(keep * c + (1 - keep) * b for c, b in zip(color, BACKGROUND)))
        prop.SetOpacity(alpha)
        actor.SetVisibility(self.shown(actor))
        key = name + "#" + str(len(self.prototype_actors[name]))
        self.prototype_actors[name].append((key, actor))
        self.viewer.actors[key] = actor
        self.viewer.ren.AddActor(actor)
        self.material_actors[self.scene["prototypes"][name]["material"]].append(actor)
        self.actor_instances[actor] = instances
        self._actor_is_glyph[actor] = glyph

    def add_mesh(self, name, vertices, faces, force=False):
        import numpy as np
        import vtk
        from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
        p = self.scene['prototypes'][name]
        group = (p.get('mesh_key', name), p['material'], p['role'])
        aliases = self._groups[group]
        previous = self._group_loaded.get(group)
        if (not force and previous is not None and previous[0] is vertices and
                previous[1] is faces and previous[2] is aliases):
            return
        # Shared geometry becomes available for every selected alias at once.
        # Logical names and per-instance transforms remain in the picking rows.
        arrays = (vertices, faces)
        for alias in aliases:
            self.meshes[alias] = arrays
        self._group_loaded[group] = (vertices, faces, aliases)
        self._stats_dirty = True
        locator_key = (p.get('mesh_key', name), p['role'] if self.display == 'exterior' else '')
        self.locators.pop(locator_key, None)
        if not any(alias in self._drawn for alias in aliases):
            return  # a section opening: its shape is needed, not its surface
        name = self._group_representatives.setdefault(group, name)
        material = self.scene["prototypes"][name]["material"]
        old_actors = self.prototype_actors.get(name, [])
        for key, actor in old_actors:
            self.viewer.ren.RemoveActor(actor)
            self.viewer.actors.pop(key, None)
            self.actor_instances.pop(actor, None)
            self._actor_is_glyph.pop(actor, None)
            for table in (self.actor_layer, self.actor_material, self.actor_volume, self.actor_prototype):
                table.pop(actor, None)
            self.material_actors[material].remove(actor)
        self.prototype_actors[name] = []
        source_key = (id(vertices), id(faces))
        if source_key not in self._source_data:
            self._source_data[source_key] = (vertices, faces, polydata(vertices, faces))
        source = self._source_data[source_key][2]
        if previous is not None:
            # Release the preview polydata once no group draws it any more.
            old_key = (id(previous[0]), id(previous[1]))
            if old_key != source_key:
                self._source_users[old_key] -= 1
                if self._source_users[old_key] <= 0:
                    self._source_users.pop(old_key, None)
                    self._source_data.pop(old_key, None)
                self._source_users[source_key] += 1
        else:
            self._source_users[source_key] += 1
        by_layer, matrices = defaultdict(list), []
        for instance in (row for alias in aliases for row in self.placements[alias]):
            components = glyph_components(instance["matrix"]) if self.instancing else None
            if components is None:
                matrices.append(instance)
            else:
                by_layer[instance["layer"], instance.get("volume")].append((instance, components))
        for glyphs in by_layer.values():
            points = vtk.vtkPoints()
            points.SetData(numpy_to_vtk(np.array([np.asarray(i["matrix"])[:3, 3] for i, _ in glyphs]), deep=True))
            data = vtk.vtkPolyData()
            data.SetPoints(points)
            for key, values in (("rotation", [c[0] for _, c in glyphs]),
                                ("scale", [c[1] for _, c in glyphs])):
                array = numpy_to_vtk(np.asarray(values, dtype=float), deep=True)
                array.SetName(key)
                data.GetPointData().AddArray(array)
            ids = numpy_to_vtkIdTypeArray(np.arange(len(glyphs), dtype=np.int64), deep=True)
            ids.SetName("instance_id")
            data.GetPointData().AddArray(ids)
            mapper = vtk.vtkGlyph3DMapper()
            mapper.SetSourceData(source)
            mapper.SetInputData(data)
            mapper.SetOrientationArray("rotation")
            mapper.SetOrientationModeToQuaternion()
            mapper.OrientOn()
            mapper.SetScaleArray("scale")
            mapper.SetScaleModeToScaleByVectorComponents()
            mapper.ScalingOn()
            mapper.SetSelectionIdArray("instance_id")
            mapper.SetUseSelectionIds(True)
            mapper.ScalarVisibilityOff()
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            self._remember(actor, name, [i for i, _ in glyphs], glyph=True)
        if matrices:
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(source)
            mapper.ScalarVisibilityOff()
            for instance in matrices:
                matrix = vtk.vtkMatrix4x4()
                matrix.DeepCopy(np.asarray(instance["matrix"]).ravel())
                actor = vtk.vtkActor()
                actor.SetMapper(mapper)
                actor.SetUserMatrix(matrix)
                self._remember(actor, name, [instance])
        self.viewer.bBuiltPipelines = True

    def finish_stage(self):
        """Rebuild what is drawn once per stage: edge lines, outlines and caps."""
        self._sections = None
        self._build_edges()
        self._build_outlines()
        if self.plane is not None:
            self.build_caps()
        self.update_stats()

    # ---- edge lines and section outlines --------------------------------
    def _add_lines(self, pairs, color, width, opacity, layer=None, material=None, clip=True,
                   volume=None):
        import vtk
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(line_polydata(pairs))
        mapper.ScalarVisibilityOff()
        # VTK's default line offset (-1, -4) is meters in a building-size view.
        mapper.SetRelativeCoincidentTopologyLineOffsetParameters(0., 3.8)
        mapper.SetResolveCoincidentTopologyToPolygonOffset()
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        prop = actor.GetProperty()
        prop.LightingOff()
        prop.SetColor(*color)
        prop.SetLineWidth(width)
        prop.SetOpacity(opacity)
        actor.PickableOff()
        if clip and self.plane is not None:
            mapper.AddClippingPlane(self.plane)
        self.actor_layer[actor], self.actor_material[actor], self.actor_volume[actor] = layer, material, volume
        self.viewer.ren.AddActor(actor)
        return actor

    def _forget(self, actors):
        for actor in actors:
            self.viewer.ren.RemoveActor(actor)
            for table in (self.actor_layer, self.actor_material, self.actor_volume, self.actor_prototype):
                table.pop(actor, None)
        actors.clear()

    def _build_edges(self):
        import numpy as np
        self._forget(self.line_actors)
        if not self.edges or not self.colored:
            return
        from .visualization import _LOW_DENSITY
        from .visualization import _scene_bounds
        bounds = _scene_bounds(self.viewer)
        smallest = (_MIN_EDGE_SIZE * np.linalg.norm(np.subtract(bounds[1::2], bounds[::2]))
                    if bounds is not None else 0.)
        plain = {row["label"] for row in self.presets["layers"] if not row["edges"]}
        batches, kept = defaultdict(list), {}
        for name in self._drawn & set(self.meshes):
            p, arrays = self.scene["prototypes"][name], self.meshes[name]
            placements = [row for row in self.placements[name] if row["layer"] not in plain]
            if (p["density"] <= _LOW_DENSITY or not placements or
                    np.linalg.norm(np.ptp(arrays[0], axis=0)) < smallest):
                continue
            key = (id(arrays[0]), id(arrays[1]))
            if key not in kept:
                kept[key] = self._edge_pairs.get(key)
                if kept[key] is None:
                    kept[key] = (arrays, feature_lines(*arrays))
            pairs = kept[key][1]
            if not len(pairs) or len(pairs) * len(placements) > _MAX_EDGE_PAIRS:
                continue
            for instance in placements:
                m = np.asarray(instance["matrix"])
                batches[instance["layer"], p["material"], instance.get("volume")].append(
                    pairs @ m[:3, :3].T + m[:3, 3])
        self._edge_pairs = kept
        for (layer, material, volume), rows in batches.items():
            self.line_actors.append(self._add_lines(np.concatenate(rows), INK, 1., .38, layer, material,
                                                    volume=volume))
        self.apply_styles()

    def _build_outlines(self):
        """Colored outlines in the viewer's fixed x, y and z planes (``cutter=True``)."""
        import numpy as np
        self._forget([self.viewer.actors.pop(k) for k in list(self.viewer.actors)
                      if k.startswith("section_")])
        for name, origin in self.viewer.cutterOrigins.items():
            cut = self.sections().cut(origin, self.viewer.cutterNormals[name], drawn=self._drawn)
            if cut.outlines:
                pairs = cut.to_world(np.concatenate(list(cut.outlines.values())))
                self.viewer.actors["section_" + name] = self._add_lines(
                    pairs, self.viewer.cutterColors[name], 2., 1., clip=False)

    # ---- plane section --------------------------------------------------
    def sections(self):
        if self._sections is None:
            from ._visualization_section import SceneSections
            exact = getattr(self.viewer, "exact_shapes", None)
            self._sections = SceneSections(self.scene["instances"], self.meshes,
                                           exact(self.scene["instances"]) if exact else None)
        return self._sections

    def _clipped_actors(self):
        # Caps lie in the plane and fixed outlines have their own planes.
        return list(self.actor_instances) + self.line_actors + list(self.external_actors.values())

    def set_section(self, origin, normal, up=None, caps=True):
        """Remove everything on the positive side of a plane (mm) and fill the cut."""
        import numpy as np
        import vtk
        fresh = self.plane is None
        if fresh:
            self.plane = vtk.vtkPlane()
        normal = np.asarray(normal, dtype=float)
        normal = normal / np.linalg.norm(normal)
        self.plane.SetOrigin(*origin)
        self.plane.SetNormal(*(-normal))
        self.section_origin, self.section_normal = list(origin), list(normal)
        self._section_up = up
        if fresh:
            for actor in self._clipped_actors():
                actor.GetMapper().AddClippingPlane(self.plane)
            self.apply_styles()
        if caps:
            self.build_caps()
        else:
            self._forget(self.cap_actors)
            self.section = None

    def clear_section(self):
        if self.plane is None:
            return
        for actor in self._clipped_actors():
            actor.GetMapper().RemoveClippingPlane(self.plane)
        self.plane = self.section = None
        self._forget(self.cap_actors)
        self.apply_styles()

    def build_caps(self):
        import numpy as np
        import vtk
        from .visualization import _LOW_DENSITY
        self._forget(self.cap_actors)
        if not self.scene or not self.meshes:
            return
        solid = {name for name in self._drawn if self.scene["prototypes"][name]["density"] > _LOW_DENSITY}
        self.section = self.sections().cut(self.section_origin, self.section_normal,
                                           self._section_up, drawn=solid)
        instances = self.scene["instances"]
        groups = defaultdict(list)
        for index in self.section.pieces:
            row = instances[index]
            groups[row["layer"], self.scene["prototypes"][row["prototype"]]["material"],
                   row.get("volume")].append(index)
        for (layer, material, volume), rows in groups.items():
            quads = np.concatenate([self.section.pieces[i] for i in rows])
            points = self.section.to_world(quads).reshape(-1, 3)
            base = 4 * np.arange(len(quads))[:, None]
            faces = np.concatenate((base + [0, 1, 2], base + [0, 2, 3]))
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(polydata(points, faces.astype(np.int64)))
            mapper.ScalarVisibilityOff()
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().LightingOff()
            actor.PickableOff()
            self.actor_layer[actor], self.actor_material[actor], self.actor_volume[actor] = layer, material, volume
            self.actor_prototype[actor] = instances[rows[0]]["prototype"]
            self.viewer.ren.AddActor(actor)
            self.cap_actors.append(actor)
        self.apply_styles()
        self.refresh_overlays()

    # ---- visibility -----------------------------------------------------
    hidden_materials = frozenset()

    def shown(self, actor):
        """Shown by the layer and volume choices (materials and gases are the controls')."""
        own = self.styles.get(("volume", self.actor_volume.get(actor)))
        return bool(self.layers.get(self.actor_layer.get(actor), True) and
                    (own is None or own.get("shown") is not False))

    def visible(self, actor):
        return self.shown(actor) and self.actor_material.get(actor) not in self.hidden_materials

    def refresh_overlays(self, hidden=None):
        """Edge lines and caps follow the layer and material choices."""
        if hidden is not None:
            self.hidden_materials = hidden
        for actor in self.line_actors + self.cap_actors:
            actor.SetVisibility(self.visible(actor))

    def update_stats(self):
        drawn = {n: len(f) for n, (_, f) in self.meshes.items() if n in self._drawn}
        self._stats["prototype_faces"] = sum(drawn.values())
        self._stats["placed_faces"] = sum(k * len(self.placements[n]) for n, k in drawn.items())
        self._stats["glyph_instances"] = sum(len(rows) for actor, rows in self.actor_instances.items()
                                            if self._actor_is_glyph[actor])
        self._stats["matrix_instances"] = sum(len(rows) for actor, rows in self.actor_instances.items()
                                             if not self._actor_is_glyph[actor])
        self._stats_dirty = False

    def pick(self, x, y):
        """Ray-pick shared prototype meshes; works with translucent glyphs too."""
        import numpy as np
        endpoints = []
        for depth in (0., 1.):
            self.viewer.ren.SetDisplayPoint(x, y, depth)
            self.viewer.ren.DisplayToWorld()
            point = np.asarray(self.viewer.ren.GetWorldPoint())
            endpoints.append(point[:3] / point[3])
        return self.pick_ray(*endpoints)

    def pick_ray(self, start, end):
        """Intersect local BVHs with transformed rays, without copying meshes."""
        import numpy as np
        import vtk
        start, end = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
        if self.plane is not None:
            normal = np.asarray(self.section_normal)
            a, b = ((p - self.section_origin) @ normal for p in (start, end))
            if a > 0 and b > 0:
                return None
            if (a > 0) != (b > 0):
                crossing = start + (end - start) * (a / (a - b))
                if a > 0:
                    start = crossing
                    index = self.section.find(crossing) if self.section is not None else None
                    if index is not None:
                        row = self.scene["instances"][index]
                        p = self.scene["prototypes"][row["prototype"]]
                        own = self.styles.get(("volume", row.get("volume")), {})
                        if (self.layers.get(row["layer"], True) and own.get("shown") is not False and
                                p["material"] not in self.hidden_materials):
                            return dict(position=crossing.tolist(), name=row["name"],
                                        prototype=row["prototype"], material=p["material"],
                                        density=p["density"], section=True)
                else:
                    end = crossing
        # See-through surfaces are picked only when nothing solid is behind them.
        candidates = []
        for actor, rows in self.actor_instances.items():
            opacity = actor.GetProperty().GetOpacity()
            if actor.GetVisibility() and opacity > 0:
                candidates.extend((row, opacity < .5) for row in rows)
        best, result = [float('inf')] * 2, [None, None]
        for row, faint in candidates:
            name = row['prototype']
            p = self.scene['prototypes'][name]
            key = (p.get('mesh_key', name), p['role'] if self.display == 'exterior' else '')
            if key not in self.locators:
                locator = vtk.vtkStaticCellLocator()
                locator.SetDataSet(polydata(*self.meshes[name]))
                locator.BuildLocator()
                self.locators[key] = locator
            matrix = np.asarray(row['matrix'])
            inverse = self.inverses[row['name']]
            p0 = (inverse @ np.r_[start, 1.])[:3]
            p1 = (inverse @ np.r_[end, 1.])[:3]
            t, sub_id, cell_id = vtk.mutable(0.), vtk.mutable(0), vtk.mutable(0)
            hit, coords = [0., 0., 0.], [0., 0., 0.]
            found = self.locators[key].IntersectWithLine(
                p0, p1, 1e-8, t, hit, coords, sub_id, cell_id, vtk.vtkGenericCell())
            if not found:
                continue
            # Rank every hit by world-space distance from the ray origin. (The
            # local parameter t is affine-invariant too, but this is explicit.)
            world = (matrix @ np.r_[hit, 1.])[:3]
            distance = float(np.linalg.norm(world - start))
            if distance < best[faint]:
                best[faint] = distance
                result[faint] = dict(position=world.tolist(), name=row['name'], prototype=name,
                                     material=p['material'], density=p['density'])
        best, result = (best[0], result[0]) if result[0] is not None else (best[1], result[1])
        for name, actor in self.external_actors.items():
            if not actor.GetVisibility() or actor.GetProperty().GetOpacity() <= 0:
                continue
            key = 'external:' + name
            if key not in self.locators:
                locator = vtk.vtkStaticCellLocator()
                locator.SetDataSet(actor.GetMapper().GetInput())
                locator.BuildLocator()
                self.locators[key] = locator
            t, sub_id, cell_id = vtk.mutable(0.), vtk.mutable(0), vtk.mutable(0)
            hit, coords = [0., 0., 0.], [0., 0., 0.]
            found = self.locators[key].IntersectWithLine(
                start, end, 1e-8, t, hit, coords, sub_id, cell_id, vtk.vtkGenericCell())
            if found and float(np.linalg.norm(np.asarray(hit) - start)) < best:
                best = float(np.linalg.norm(np.asarray(hit) - start))
                result = dict(position=list(hit))
        return result
