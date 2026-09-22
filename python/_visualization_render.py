"""VTK prototype rendering and placement-aware picking (UI thread only)."""
from collections import defaultdict


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


class ArrayMesh:
    """The minimal mesh interface consumed by pyg4ometry's append pipeline."""
    def __init__(self, vertices, faces):
        self.vertices, self.faces = vertices, faces

    def toVerticesAndPolygons(self):
        return self.vertices.tolist(), self.faces.tolist(), len(self.faces)


class SceneRenderer:
    def __init__(self, viewer, *, coloured=True, instancing=True, display="full", regions=None):
        self.viewer = viewer
        self.coloured = coloured
        self.instancing = instancing
        self.display = display
        # Normalized surface-region mapping (see _visualization_scene.normalize_regions).
        self.regions = regions
        self.scene = None
        self.meshes = {}
        self.placements = defaultdict(list)
        self.prototype_actors = {}
        self.actor_instances = {}
        self._actor_is_glyph = {}
        self.external_actors = {}
        self.locators, self.inverses = {}, {}
        self.options = {}
        self._groups = defaultdict(list)
        self._group_representatives = {}
        self._group_loaded = {}
        self._source_data = {}
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
        active = set(scene.get('selected', scene['prototypes'])) | set(self.meshes)
        for name, p in scene['prototypes'].items():
            if name in active:
                self._groups[(p.get('mesh_key', name), p['material'], p['role'])].append(name)
        self.placements.clear()
        materials = {}
        for instance in scene["instances"]:
            self.placements[instance["prototype"]].append(instance)
            self.inverses[instance["name"]] = np.linalg.inv(instance["matrix"])
        for p in scene["prototypes"].values():
            materials[p["material"]] = SimpleNamespace(density=p["density"])
        self.registry = SimpleNamespace(materialDict=materials)
        self.options = _material_vis_options(self.registry, VisualisationOptions)
        legend_options, legend_materials = dict(self.options), dict(materials)
        styles = self.regions['styles'] if self.display == 'exterior' and self.regions else {}
        if styles:
            for material in materials:
                roles = {p['role'] for p in scene['prototypes'].values() if p['material'] == material}
                if roles and roles <= styles.keys():
                    legend_options.pop(material, None)
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
        vo = copy.copy(self.options[p["material"]]) if self.coloured else VisualisationOptions(colour=[.6, .6, .6])
        if self.display == "exterior" and self.regions and p["role"] in self.regions["styles"]:
            # Caller-supplied display colours, not measured optical coefficients.
            style = self.regions["styles"][p["role"]]
            if self.coloured:
                vo.colour = list(style["colour"])
            vo.alpha = style["alpha"]
        return vo

    def _remember(self, actor, name, instances, glyph=False):
        vo = self.vis_option(name)
        actor.GetProperty().SetColor(*vo.colour)
        actor.GetProperty().SetOpacity(vo.alpha)
        key = name + "#" + str(len(self.prototype_actors[name]))
        self.prototype_actors[name].append((key, actor))
        self.viewer.actors[key] = actor
        self.viewer.ren.AddActor(actor)
        self.material_actors[self.scene["prototypes"][name]["material"]].append(actor)
        self.actor_instances[actor] = instances
        self._actor_is_glyph[actor] = glyph

    def add_mesh(self, name, vertices, faces):
        import numpy as np
        import vtk
        from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
        p = self.scene['prototypes'][name]
        group = (p.get('mesh_key', name), p['material'], p['role'])
        aliases = self._groups[group]
        previous = self._group_loaded.get(group)
        if (previous is not None and previous[0] is vertices and
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
        if self.viewer.cutterOrigins or self.viewer.bClipper:
            return  # append once at stage completion for accurate cut surfaces
        name = self._group_representatives.setdefault(group, name)
        material = self.scene["prototypes"][name]["material"]
        old_actors = self.prototype_actors.get(name, [])
        for key, actor in old_actors:
            self.viewer.ren.RemoveActor(actor)
            self.viewer.actors.pop(key, None)
            self.actor_instances.pop(actor, None)
            self._actor_is_glyph.pop(actor, None)
            self.material_actors[material].remove(actor)
        self.prototype_actors[name] = []
        source_key = (id(vertices), id(faces))
        if source_key not in self._source_data:
            self._source_data[source_key] = (vertices, faces, polydata(vertices, faces))
        source = self._source_data[source_key][2]
        glyphs, matrices = [], []
        for instance in (row for alias in aliases for row in self.placements[alias]):
            components = glyph_components(instance["matrix"]) if self.instancing else None
            if components is None:
                matrices.append(instance)
            else:
                glyphs.append((instance, components))
        if glyphs:
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
        """Use pyg4ometry's closed-volume clipping/section path when requested."""
        import numpy as np
        if not (self.viewer.cutterOrigins or self.viewer.bClipper):
            return
        v = self.viewer
        for actor in list(v.actors.values()):
            v.ren.RemoveActor(actor)
        v.actors.clear()
        if v.clippers:
            # The widget callback moves the live clip functions only; carry the
            # current plane into the attributes buildPipelinesAppend rebuilds from.
            plane = v.clippers[0].GetClipFunction()
            v.clipperOrigin, v.clipperNormal = list(plane.GetOrigin()), list(plane.GetNormal())
        v.clippers.clear()
        v.cutters.clear()
        self.material_actors.clear()
        v.localmeshes = {name: ArrayMesh(*mesh) for name, mesh in self.meshes.items()}
        v.instancePlacements, v.instanceVisOptions = {}, {}
        for name in self.meshes:
            v.instancePlacements[name] = [dict(name=i["name"],
                transformation=np.asarray(i["matrix"])[:3, :3],
                translation=np.asarray(i["matrix"])[:3, 3]) for i in self.placements[name]]
            v.instanceVisOptions[name] = [self.vis_option(name)] * len(self.placements[name])
        v.buildPipelinesAppend()
        for key, actor in self.external_actors.items():
            v.actors[key] = actor
            v.ren.AddActor(actor)
        for name in self.meshes:
            actor = v.actors.get(str(self.vis_option(name)))
            material = self.scene["prototypes"][name]["material"]
            if actor is not None and actor not in self.material_actors[material]:
                self.material_actors[material].append(actor)
        self.update_stats()

    def update_stats(self):
        self._stats["prototype_faces"] = sum(len(f) for _, f in self.meshes.values())
        self._stats["placed_faces"] = sum(len(f) * len(self.placements[n]) for n, (_, f) in self.meshes.items())
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
        if self.viewer.bClipper and self.viewer.clippers:
            plane = self.viewer.clippers[0].GetClipFunction()
            a, b = plane.EvaluateFunction(start), plane.EvaluateFunction(end)
            # The append viewer draws vtkClipPolyData's negative half-space.
            if a > 0 and b > 0:
                return None
            if (a > 0) != (b > 0):
                crossing = start + (end - start) * (a / (a - b))
                if a > 0:
                    start = crossing
                else:
                    end = crossing
        candidates = []
        if self.actor_instances:
            for actor, rows in self.actor_instances.items():
                if actor.GetVisibility() and actor.GetProperty().GetOpacity() > 0:
                    candidates.extend(rows)
        else:
            for name in self.meshes:
                material = self.scene['prototypes'][name]['material']
                if any(a.GetVisibility() and a.GetProperty().GetOpacity() > 0
                       for a in self.material_actors[material]):
                    candidates.extend(self.placements[name])
        best, result = float('inf'), None
        for row in candidates:
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
            if distance < best:
                best = distance
                result = dict(position=world.tolist(), name=row['name'],
                              material=p['material'], density=p['density'])
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
