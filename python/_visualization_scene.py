"""Display-only mesh preparation. Optional geometry imports stay inside functions.

Numeric mesh files are shared by synchronous loading and the isolated worker.
Nothing in this module changes the source GDML or a SIREN DetectorModel.
"""
from contextlib import contextmanager
import fnmatch
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import tempfile
import threading
import time
import warnings
import xml.etree.ElementTree as ET
import zipfile

_FORMAT = 2
_CONFIG_LOCK = threading.RLock()
_PMT_SURFACES = {"external_tpb", "photocathode_inner_surface",
                 "reflector_inner_surface", "bare_transparent_glass"}


@contextmanager
def mesh_settings(pg, slices=None, meshing=None):
    """Scope pyg4ometry's global defaults, including exceptional exits.

    SIREN loads are serialized within a process. Interactive meshing uses a
    separate process so unrelated clients cannot observe temporary settings.
    """
    with _CONFIG_LOCK:
        defaults = pg.config.SolidDefaults
        old = [(obj, attr, getattr(obj, attr)) for obj in vars(defaults).values()
               for attr in ("nslice", "nstack") if hasattr(obj, attr)]
        old_meshing = pg.config.doMeshing
        try:
            if slices is not None and slices != 0:
                if isinstance(slices, bool) or int(slices) != slices or slices < 4:
                    raise ValueError("mesh_slices must be an integer >= 4, 0 or None")
                pg.config.setGlobalMeshSliceAndStack(int(slices))
            if meshing is not None:
                pg.config.doMeshing = meshing
            yield
        finally:
            for obj, attr, value in old:
                setattr(obj, attr, value)
            pg.config.doMeshing = old_meshing


def input_manifest(path, seen=None):
    """Hash GDML plus recursively included files; reject dependency cycles."""
    path = Path(path).resolve()
    seen = set() if seen is None else seen
    if path in seen:
        raise ValueError("cyclic GDML file dependency: %s" % path)
    data = path.read_bytes()
    root = ET.fromstring(data)
    # External entity resolution belongs to the GDML reader. It cannot safely
    # participate in a content cache whose full dependency set is unknown.
    if b"<!ENTITY" in data.upper() or b"<!DOCTYPE" in data.upper():
        raise ValueError("GDML with external entities is not supported by the scene loader")
    dependencies = []
    for ref in root.iter("file"):
        name = ref.get("name")
        if name:
            child = Path(name)
            # pyg4ometry resolves includes from the process working directory.
            # Use the same path, rather than hashing a different relative file.
            dependencies.append(input_manifest(child, seen | {path}))
    return {"path": str(path), "sha256": hashlib.sha256(data).hexdigest(),
            "dependencies": dependencies,
            "eager": any(n.tag in ("replicavol", "divisionvol", "paramvol")
                         for n in root.iter())}


def cache_key(manifest, pg):
    def content(m):
        return {"sha256": m["sha256"], "eager": m["eager"],
                "dependencies": [content(d) for d in m["dependencies"]]}
    defaults = {name: {attr: getattr(obj, attr) for attr in ("nslice", "nstack")
                       if hasattr(obj, attr)}
                for name, obj in vars(pg.config.SolidDefaults).items()
                if hasattr(obj, "nslice") or hasattr(obj, "nstack")}
    spec = {"format": _FORMAT, "input": content(manifest), "defaults": defaults,
            "pyg4ometry": importlib.metadata.version("pyg4ometry"),
            "vtk": importlib.metadata.version("vtk"), "backend": pg.config.backendName(),
            "two_pi_tolerance": pg.config.twoPiComparisonTolerance}
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()


def default_cache_dir():
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "siren" / "meshes"


def _digest(vertices, faces):
    return hashlib.sha256(vertices.tobytes() + faces.tobytes()).hexdigest()


def read_mesh(path):
    """Read checked float64 vertices/int64 triangles, never pickle objects."""
    import numpy as np
    with np.load(path, allow_pickle=False) as data:
        v, f = data["vertices"], data["faces"]
        digest = str(data["sha256"].item())
    if (v.dtype != np.dtype("<f8") or f.dtype != np.dtype("<i8") or
            v.ndim != 2 or v.shape[1] != 3 or f.ndim != 2 or f.shape[1] != 3 or
            not len(v) or not len(f) or not np.isfinite(v).all() or
            f.min() < 0 or f.max() >= len(v) or _digest(v, f) != digest):
        raise ValueError("invalid mesh cache: %s" % path)
    return v, f


def write_mesh(path, vertices, faces):
    import numpy as np
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    v, f = np.asarray(vertices, dtype="<f8"), np.asarray(faces, dtype="<i8")
    tmp = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".npz", delete=False) as out:
            tmp = Path(out.name)
            np.savez_compressed(out, vertices=v, faces=f, sha256=_digest(v, f))
        os.replace(tmp, path)
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)


def mesh_arrays(mesh):
    """Triangulate a prototype once; keep input coordinates at double precision."""
    import numpy as np
    vertices, polygons, _ = mesh.toVerticesAndPolygons()
    v = np.asarray(vertices, dtype="<f8")
    if not len(v) or not polygons:
        raise ValueError("null display mesh")
    if all(len(face) == 3 for face in polygons):
        return v, np.asarray(polygons, dtype="<i8")
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(v, deep=True))
    cells = vtk.vtkCellArray()
    flat = np.concatenate([np.asarray([len(face), *face], dtype=np.int64) for face in polygons])
    cells.ImportLegacyFormat(numpy_to_vtkIdTypeArray(flat, deep=True))
    pd = vtk.vtkPolyData()
    pd.SetPoints(points)
    pd.SetPolys(cells)
    triangle = vtk.vtkTriangleFilter()
    triangle.SetInputData(pd)
    triangle.Update()
    result = triangle.GetOutput()
    return (vtk_to_numpy(result.GetPoints().GetData()).astype("<f8"),
            vtk_to_numpy(result.GetPolys().GetData()).reshape(-1, 4)[:, 1:].astype("<i8"))


def exterior_faces(vertices, faces, role):
    """Outer faces of annotated CCM PMT shells of revolution about local z.

    Positive radial normal selects the outer wall; negative radial normal is
    the glass/vacuum interface, and axial faces close artificial region cuts.
    Unknown roles keep all faces. These open surfaces are display-only.
    """
    import numpy as np
    if role not in _PMT_SURFACES:
        return faces
    tri = vertices[faces]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    centres = tri.mean(axis=1)
    radial = (normals[:, :2] * centres[:, :2]).sum(axis=1)
    tolerance = np.linalg.norm(normals, axis=1) * np.linalg.norm(centres[:, :2], axis=1) * 1e-10
    result = faces[radial > tolerance]
    if not len(result):
        raise ValueError("annotated PMT shell has no outward radial faces")
    return result


def _auxiliary(lv):
    return {str(a.auxtype): str(a.auxvalue) for a in getattr(lv, "auxiliary", [])}


def scene_structure(reg, pg):
    """Traverse hidden parents too. Mesh-free for ordinary placements."""
    import numpy as np
    prototypes, instances, sources = {}, [], {}
    solid_keys = {}

    def prototype(key, lv, mesh=None):
        if key not in prototypes:
            try:
                density = float(lv.material.density)
            except (AttributeError, TypeError, ValueError):
                density = 1.0
            aux = _auxiliary(lv)
            if mesh is None and id(lv.solid) not in solid_keys:
                # Stable traversal ordinals distinguish different objects even
                # if separately included GDML registries reuse a solid name.
                solid_keys[id(lv.solid)] = 'solid:%d:%s' % (len(solid_keys), lv.solid.name)
            prototypes[key] = dict(name=lv.name, material=lv.material.name, density=density,
                                   role=aux.get("ccm_pmt_region", ""), auxiliary=aux,
                                   mesh_key=("logical:" + key if mesh is not None
                                             else solid_keys[id(lv.solid)]))
            sources[key] = (lv, mesh)
        return key

    def walk(lv, matrix, path, ancestors):
        if id(lv) in ancestors:
            raise ValueError("cyclic logical volume hierarchy")
        if lv.type == "logical":
            key = prototype(lv.name, lv)
            instances.append(dict(prototype=key, name=path, matrix=matrix.tolist()))
        for pv in lv.daughterVolumes:
            if pv.type == "placement":
                local = np.eye(4)
                local[:3, :3] = np.linalg.inv(pg.transformation.tbxyz2matrix(pv.rotation.eval()))
                if pv.scale:
                    local[:3, :3] = local[:3, :3] @ np.diag(pv.scale.eval())
                local[:3, 3] = pv.position.eval()
                walk(pv.logicalVolume, matrix @ local, path + "/" + pv.name, ancestors | {id(lv)})
            elif pv.type in ("replica", "division", "parametrised"):
                for i, (mesh, transform) in enumerate(zip(pv.meshes, pv.transforms)):
                    rotation, position = transform
                    if pv.type == "parametrised":
                        rotation, position = rotation.eval(), position.eval()
                    local = np.eye(4)
                    local[:3, :3] = pg.transformation.tbxyz2matrix(rotation)
                    local[:3, 3] = position
                    key = prototype(pv.name + "#" + str(i), pv.logicalVolume, mesh.localmesh)
                    instances.append(dict(prototype=key, name=path + "/" + pv.name + "#" + str(i),
                                          matrix=(matrix @ local).tolist()))
            else:
                raise ValueError("unsupported display placement: %s" % pv.type)

    world = reg.getWorldVolume()
    walk(world, np.eye(4), world.name, set())
    return dict(prototypes=prototypes, instances=instances, world=world.name), sources


def selected_prototypes(scene, show_gas=False, hidden_volumes=(), display="full", only=None):
    if display not in ("full", "exterior"):
        raise ValueError("display must be 'full' or 'exterior'")
    selected = []
    for key, p in scene["prototypes"].items():
        if only is not None and key not in only:
            continue
        if any(fnmatch.fnmatchcase(p["name"], pattern) for pattern in hidden_volumes):
            continue
        if not show_gas and p["density"] <= 0.05:
            continue
        if display == "exterior" and p["role"] == "vacuum":
            continue
        selected.append(key)
    return selected


def cached_scene_available(path, *, mesh_slices=48, cache=True, cache_dir=None,
                           show_gas=False, hidden_volumes=(), display="full", only=None):
    """Skip preview work when all requested full-detail prototypes are cached."""
    if not cache:
        return False
    import pyg4ometry as pg
    manifest = input_manifest(path)
    def eager(m):
        return m['eager'] or any(eager(d) for d in m['dependencies'])
    if eager(manifest):
        return False
    with mesh_settings(pg, mesh_slices, meshing=False):
        key = cache_key(manifest, pg)
        registry = pg.gdml.Reader(str(path)).getRegistry()
        metadata, _ = scene_structure(registry, pg)
        directory = (Path(cache_dir) if cache_dir is not None else default_cache_dir()) / key
        keys = {metadata['prototypes'][name]['mesh_key'] for name in
                selected_prototypes(metadata, show_gas, hidden_volumes, display, only)}
        for mesh_key in keys:
            token = hashlib.sha256(mesh_key.encode()).hexdigest()
            try:
                read_mesh(directory / (token + '.npz'))
            except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile):
                return False
    return input_manifest(path) == manifest


def prepare_scene(path, output_dir, *, mesh_slices=48, cache=True, cache_dir=None,
                  show_gas=False, hidden_volumes=(), display="full", only=None,
                  emit=None, stage="detail"):
    """Prepare selected meshes and emit small JSON events plus numeric files.

    Ordinary hidden volumes are never meshed merely to traverse their children.
    Required CSG operands still have to be evaluated. Complex placement files
    use eager meshing for compatibility, reported explicitly in the metadata.
    """
    import pyg4ometry as pg
    from .visualization import _cached_solid_meshes, _read_pyg4ometry_registry
    emit = (lambda event: None) if emit is None else emit
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    timings = {}
    start = time.perf_counter()
    manifest = input_manifest(path)
    timings["input_hash_seconds"] = time.perf_counter() - start

    def eager_input(m):
        return m["eager"] or any(eager_input(d) for d in m["dependencies"])

    eager = eager_input(manifest)
    with mesh_settings(pg, mesh_slices, meshing=eager):
        key = cache_key(manifest, pg)
        start = time.perf_counter()
        reg = (_read_pyg4ometry_registry(str(path), pg) if eager
               else pg.gdml.Reader(str(path)).getRegistry())
        timings["parse_seconds"] = time.perf_counter() - start
        scene, sources = scene_structure(reg, pg)
        if eager:
            # Division parameters can mutate a shared solid; only their eager
            # per-volume mesh snapshots are safe to reuse.
            for name, prototype in scene['prototypes'].items():
                prototype['mesh_key'] = 'logical:' + name
        selected = selected_prototypes(scene, show_gas, hidden_volumes, display, only)
        scene.update(selected=selected, eager_fallback=eager, cache_key=key, stage=stage)
        emit(dict(kind="structure", scene=scene))
        hits, misses, failed = 0, 0, []
        prepared, outputs = {}, {}
        reuses = 0
        pending_cache = []
        start = time.perf_counter()
        with _cached_solid_meshes(reg):
            for index, name in enumerate(selected):
                lv, explicit = sources[name]
                mesh_key = scene['prototypes'][name]['mesh_key']
                token = hashlib.sha256(mesh_key.encode()).hexdigest()
                # Meshes mutated by divisions must retain their eager snapshots.
                disk = (Path(cache_dir) if cache_dir is not None else default_cache_dir()) / key / (token + ".npz")
                arrays = prepared.get(mesh_key)
                if arrays is not None:
                    reuses += 1
                elif cache:
                    try:
                        arrays = read_mesh(disk)
                        hits += 1
                    except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile):
                        pass
                try:
                    if arrays is None:
                        mesh = explicit if explicit is not None else (
                            lv.mesh.localmesh if eager and lv.mesh is not None else lv.solid.mesh())
                        arrays = mesh_arrays(mesh)
                        misses += 1
                        if cache:
                            pending_cache.append((disk, arrays))
                    prepared[mesh_key] = arrays
                    v, f = arrays
                    role = scene['prototypes'][name]['role'] if display == 'exterior' else ''
                    variant = (mesh_key, role)
                    if variant not in outputs:
                        if display == "exterior":
                            f = exterior_faces(v, f, role)
                        suffix = hashlib.sha256(json.dumps(variant).encode()).hexdigest()
                        dest = output_dir / (stage + "-" + suffix + ".npz")
                        write_mesh(dest, v, f)
                        outputs[variant] = (dest, len(f))
                    dest, face_count = outputs[variant]
                    emit(dict(kind="mesh", prototype=name, path=str(dest), stage=stage,
                              completed=index + 1, total=len(selected), faces=face_count, vertices=len(v)))
                except Exception as exc:
                    if stage != "preview":
                        raise RuntimeError("meshing %s failed: %s" % (name, exc)) from exc
                    failed.append(name)
                    emit(dict(kind="preview_skip", prototype=name, message=str(exc)))
        timings["mesh_seconds"] = time.perf_counter() - start
        if input_manifest(path) != manifest:
            raise RuntimeError("GDML inputs changed while loading; retry with stable inputs")
        start = time.perf_counter()
        for disk, arrays in pending_cache:
            try:
                write_mesh(disk, *arrays)
            except OSError as exc:
                warnings.warn("mesh cache unavailable: %s" % exc, RuntimeWarning)
        timings["cache_write_seconds"] = time.perf_counter() - start
    stats = dict(kind="stage_done", stage=stage, timings=timings, cache_hits=hits,
                 cache_misses=misses, shared_mesh_reuses=reuses,
                 preview_skipped=failed, eager_fallback=eager)
    emit(stats)
    return scene, stats
