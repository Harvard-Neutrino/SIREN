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


def normalize_regions(regions):
    """Validate the surface-region mapping used by ``display='exterior'``.

    ``regions`` is a mapping with ``auxtype`` (the GDML ``<auxiliary auxtype=...>``
    key whose value names a logical volume's region), ``styles`` (region name ->
    ``{"label", "colour", "alpha"}``; these regions must be shells of revolution
    about their local z axis, of which only outward faces are displayed) and an
    optional ``hidden`` list of region names dropped from the exterior display.
    Detector-specific names belong to the caller, not to this module. Returns a
    JSON-compatible dict, or None.
    """
    if regions is None:
        return None
    try:
        auxtype = str(regions["auxtype"])
        styles = dict(regions["styles"])
    except (KeyError, TypeError) as exc:
        raise ValueError("regions needs 'auxtype' and 'styles' entries") from exc
    out = {"auxtype": auxtype, "styles": {},
           "hidden": [str(role) for role in regions.get("hidden", ())]}
    for role, style in styles.items():
        colour = [float(c) for c in style["colour"]]
        if len(colour) != 3:
            raise ValueError("region %r colour must have three components" % role)
        out["styles"][str(role)] = {"label": str(style.get("label", role)),
                                    "colour": colour, "alpha": float(style.get("alpha", 1.0))}
    return out


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


def input_manifest(path, seen=None, strict=True):
    """Hash GDML plus recursively included files; reject dependency cycles.

    With ``strict=False`` an entity-based file yields an uncacheable manifest
    (``cacheable=False``) instead of raising, so it can still be displayed.
    """
    path = Path(path).resolve()
    seen = set() if seen is None else seen
    if path in seen:
        raise ValueError("cyclic GDML file dependency: %s" % path)
    data = path.read_bytes()
    # External entity resolution belongs to the GDML reader. It cannot safely
    # participate in a content cache whose full dependency set is unknown.
    upper = data.upper()
    if b"<!ENTITY" in upper or b"<!DOCTYPE" in upper:
        if strict:
            raise ValueError("GDML with external entities is not supported by the scene loader")
        return {"path": str(path), "sha256": hashlib.sha256(data).hexdigest(),
                "dependencies": [], "cacheable": False,
                "eager": any(tag in data for tag in (b"<replicavol", b"<divisionvol", b"<paramvol"))}
    root = ET.fromstring(data)
    dependencies = []
    for ref in root.iter("file"):
        name = ref.get("name")
        if name:
            child = Path(name)
            # pyg4ometry resolves includes from the process working directory.
            # Use the same path, rather than hashing a different relative file.
            dependencies.append(input_manifest(child, seen | {path}))
    return {"path": str(path), "sha256": hashlib.sha256(data).hexdigest(),
            "dependencies": dependencies, "cacheable": True,
            "eager": any(n.tag in ("replicavol", "divisionvol", "paramvol")
                         for n in root.iter())}


def inputs_unchanged(manifest):
    """Re-hash the manifest's files without re-parsing any XML."""
    def check(m):
        try:
            digest = hashlib.sha256(Path(m["path"]).read_bytes()).hexdigest()
        except OSError:
            return False
        return digest == m["sha256"] and all(check(d) for d in m["dependencies"])
    return check(manifest)


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


_PARTIAL_PREFIX = ".partial-"
_STALE_PARTIAL_SECONDS = 3600


def sweep_partial_files(directory, now=None):
    """Remove abandoned partial writes left in a cache directory.

    ``write_mesh`` publishes atomically, but a worker killed mid-write (window
    closed during the cache flush) cannot run its cleanup. Only files older than
    an hour are removed so a concurrent writer's live temporary is untouched.
    """
    cutoff = (time.time() if now is None else now) - _STALE_PARTIAL_SECONDS
    removed = 0
    try:
        entries = list(Path(directory).glob(_PARTIAL_PREFIX + "*"))
    except OSError:
        return 0
    for entry in entries:
        try:
            if entry.stat().st_mtime < cutoff:
                entry.unlink()
                removed += 1
        except OSError:
            pass
    return removed


def stage_mesh(path, vertices, faces):
    """Write a partial file beside ``path``; publish it with ``os.replace``."""
    import numpy as np
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    v, f = np.asarray(vertices, dtype="<f8"), np.asarray(faces, dtype="<i8")
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=_PARTIAL_PREFIX,
                                     suffix=".npz", delete=False) as out:
        tmp = Path(out.name)
        try:
            np.savez_compressed(out, vertices=v, faces=f, sha256=_digest(v, f))
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
    return tmp


def write_mesh(path, vertices, faces):
    tmp = None
    try:
        tmp = stage_mesh(path, vertices, faces)
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


def exterior_faces(vertices, faces, role, regions):
    """Outer faces of a region shell of revolution about its local z axis.

    Positive radial normal selects the outer wall; negative radial normal is
    the inner interface, and axial faces close artificial region cuts. Roles
    absent from ``regions['styles']`` keep all faces. Display-only surfaces.
    """
    import numpy as np
    if regions is None or role not in regions["styles"]:
        return faces
    tri = vertices[faces]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    centres = tri.mean(axis=1)
    radial = (normals[:, :2] * centres[:, :2]).sum(axis=1)
    tolerance = np.linalg.norm(normals, axis=1) * np.linalg.norm(centres[:, :2], axis=1) * 1e-10
    result = faces[radial > tolerance]
    if not len(result):
        raise ValueError("region shell %r has no outward radial faces" % role)
    return result


def _auxiliary(lv):
    return {str(a.auxtype): str(a.auxvalue) for a in getattr(lv, "auxiliary", [])}


def scene_structure(reg, pg, region_auxtype=None):
    """Traverse hidden parents too. Mesh-free for ordinary placements.

    ``region_auxtype`` names the auxiliary key that classifies a logical volume's
    surface region; without it every prototype's role is empty.
    """
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
                                   role=aux.get(region_auxtype, "") if region_auxtype else "",
                                   auxiliary=aux,
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


def selected_prototypes(scene, show_gas=False, hidden_volumes=(), display="full", only=None,
                        hidden_roles=()):
    from .visualization import _LOW_DENSITY
    if display not in ("full", "exterior"):
        raise ValueError("display must be 'full' or 'exterior'")
    selected = []
    for key, p in scene["prototypes"].items():
        if only is not None and key not in only:
            continue
        if any(fnmatch.fnmatchcase(p["name"], pattern) for pattern in hidden_volumes):
            continue
        if not show_gas and p["density"] <= _LOW_DENSITY:
            continue
        if display == "exterior" and p["role"] in hidden_roles:
            continue
        selected.append(key)
    return selected


def cached_scene_available(path, *, mesh_slices=48, cache=True, cache_dir=None,
                           show_gas=False, hidden_volumes=(), display="full", only=None,
                           regions=None):
    """Skip preview work when all requested full-detail prototypes are cached."""
    if not cache:
        return False
    regions = normalize_regions(regions)
    import pyg4ometry as pg
    manifest = input_manifest(path, strict=False)
    def eager(m):
        return m['eager'] or any(eager(d) for d in m['dependencies'])
    if eager(manifest) or not manifest["cacheable"]:
        return False
    with mesh_settings(pg, mesh_slices, meshing=False):
        key = cache_key(manifest, pg)
        registry = pg.gdml.Reader(str(path)).getRegistry()
        metadata, _ = scene_structure(registry, pg, regions and regions["auxtype"])
        directory = (Path(cache_dir) if cache_dir is not None else default_cache_dir()) / key
        keys = {metadata['prototypes'][name]['mesh_key'] for name in
                selected_prototypes(metadata, show_gas, hidden_volumes, display, only,
                                    regions["hidden"] if regions else ())}
        # Validate content: a truncated or corrupt entry is a miss, so the
        # preview still runs instead of waiting for full-detail recovery.
        for mesh_key in keys:
            token = hashlib.sha256(mesh_key.encode()).hexdigest()
            try:
                read_mesh(directory / (token + '.npz'))
            except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile):
                return False
    return inputs_unchanged(manifest)


def prepare_scene(path, output_dir, *, mesh_slices=48, cache=True, cache_dir=None,
                  show_gas=False, hidden_volumes=(), display="full", only=None,
                  regions=None, emit=None, stage="detail"):
    """Prepare selected meshes and emit small JSON events plus numeric files.

    Ordinary hidden volumes are never meshed merely to traverse their children.
    Required CSG operands still have to be evaluated. Complex placement files
    use eager meshing for compatibility, reported explicitly in the metadata.
    ``display='exterior'`` needs a ``regions`` mapping (see normalize_regions).
    """
    import pyg4ometry as pg
    regions = normalize_regions(regions)
    if display == "exterior" and regions is None:
        raise ValueError("display='exterior' requires a regions mapping")
    from .visualization import _cached_solid_meshes, _read_pyg4ometry_registry
    emit = (lambda event: None) if emit is None else emit
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    timings = {}
    start = time.perf_counter()
    manifest = input_manifest(path, strict=False)
    if cache and not manifest["cacheable"]:
        warnings.warn("mesh cache disabled: GDML uses DOCTYPE/entity declarations", RuntimeWarning)
        cache = False
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
        scene, sources = scene_structure(reg, pg, regions and regions["auxtype"])
        if eager:
            # Division parameters can mutate a shared solid; only their eager
            # per-volume mesh snapshots are safe to reuse.
            for name, prototype in scene['prototypes'].items():
                prototype['mesh_key'] = 'logical:' + name
        selected = selected_prototypes(scene, show_gas, hidden_volumes, display, only,
                                       regions["hidden"] if regions else ())
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
                            f = exterior_faces(v, f, role, regions)
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
        if not inputs_unchanged(manifest):
            raise RuntimeError("GDML inputs changed while loading; retry with stable inputs")
        start = time.perf_counter()
        for directory in {disk.parent for disk, _ in pending_cache}:
            sweep_partial_files(directory)
        # Stage every file, re-check the inputs, then publish atomically. An
        # input edited during a long write must not be reused under the old key.
        staged = []
        try:
            for disk, arrays in pending_cache:
                try:
                    staged.append((stage_mesh(disk, *arrays), disk))
                except OSError as exc:
                    warnings.warn("mesh cache unavailable: %s" % exc, RuntimeWarning)
                    break
            if staged and inputs_unchanged(manifest):
                for tmp, disk in staged:
                    os.replace(tmp, disk)
            elif staged:
                warnings.warn("GDML inputs changed while loading; cache not published",
                              RuntimeWarning)
        finally:
            for tmp, _ in staged:
                tmp.unlink(missing_ok=True)
        timings["cache_write_seconds"] = time.perf_counter() - start
    stats = dict(kind="stage_done", stage=stage, timings=timings, cache_hits=hits,
                 cache_misses=misses, shared_mesh_reuses=reuses,
                 preview_skipped=failed, eager_fallback=eager)
    emit(stats)
    return scene, stats
