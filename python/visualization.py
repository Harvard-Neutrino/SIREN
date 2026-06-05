"""Rudimentary geometry visualization for SIREN ``DetectorModel`` objects.

Entry points (all with deferred heavy imports so ``import siren`` stays light):

* :func:`plot_cross_sections` -- matplotlib material-coloured 2-D slices through
  the geometry (samples ``GetMassDensity``; needs only matplotlib).
* :func:`to_gdml` -- export the sectors to a single, self-contained, *standard*
  GDML file viewable in ROOT (``TGeoManager::Import``), pyg4ometry, Geant4, or a
  browser (jsroot). Identical solid shapes are de-duplicated; every volume is
  named after its sector and every material carries its name + density, so a
  picking-capable viewer shows what you click on.
* :func:`view` -- interactive 3-D (or off-screen PNG) render via pyg4ometry,
  coloured by material, with a SIREN-backed right-click picker (robust where
  pyg4ometry's own picker segfaults) and keyboard controls for the legend,
  bounding box, section outlines, per-material visibility, and opacity.
* :func:`to_web` -- export the geometry to glTF / a three.js HTML page for the
  browser (needs ``pip install pygltflib``); headless-friendly.
* :func:`describe` -- print a material summary (name, density, sector count) for
  headless inspection.

Units/conventions: SIREN geometry is in metres; ``Cylinder``/``Cone`` store the
*full* z-extent, so the exported GDML uses the standard convention
(``<tube>`` z = full length) and renders correctly in external tools -- unlike
SIREN's own hand-written sub-GDMLs, whose parser treats ``<tube>`` z as a half.
"""
import math

__all__ = ["plot_cross_sections", "to_gdml", "view", "to_web", "describe", "at",
           "list_volumes"]

_PLANES = {
    #  name : (horizontal axis, vertical axis, fixed axis)
    "yz": ("z", "y", "x"),
    "xy": ("x", "y", "z"),
    "xz": ("x", "z", "y"),
}
_AXIS = {"x": 0, "y": 1, "z": 2}
_AXIS_LABEL = {"x": "x  (west)", "y": "y  (up)", "z": "z  (beam)"}


def _sectors(model):
    return list(model.Sectors)


def _material_name(model, mid):
    try:
        return model.Materials.GetMaterialName(mid)
    except Exception:
        return "mat%d" % mid


def _sector_density(s):
    """Scalar density (g/cm^3) of a sector, evaluated at its geo centre."""
    try:
        return float(s.density.Evaluate(s.geo.placement.Position))
    except Exception:
        return 1.0


# ---------------------------------------------------------------------------
# 1. matplotlib cross-sections
# ---------------------------------------------------------------------------
def plot_cross_sections(model, center=(0.0, 0.0, 0.0), half=15.0,
                        planes=("yz", "xy"), n=240, path=None, show=False,
                        title=None):
    """Plot material-coloured 2-D cross-sections of *model*.

    Samples the detector mass density on an ``n x n`` grid in each requested
    plane and colours each cell by its material. The legend lists only the
    materials that actually appear in the slices.

    :param center: (x, y, z) in metres (detector-local) about which to slice.
    :param half: half-width of each slice in metres.
    :param planes: any of ``"yz"``, ``"xy"``, ``"xz"``.
    :param path: if given, save the figure (PNG).
    :param show: open an interactive window.
    :returns: the matplotlib Figure.
    """
    import numpy as np
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    from siren.detector import DetectorPosition
    from siren.math import Vector3D

    def dens(x, y, z):
        return model.GetMassDensity(DetectorPosition(Vector3D(float(x), float(y), float(z))))

    # density -> material name, from the sectors (used only for labelling)
    name_of = {}
    for s in _sectors(model):
        name_of.setdefault(round(_sector_density(s), 6), _material_name(model, s.material_id))

    # sample every requested plane first
    planes = [p for p in planes if p in _PLANES]
    coord = np.linspace(-half, half, n)
    grids = []
    for pl in planes:
        ha, va, fa = _PLANES[pl]
        grid = np.empty((n, n))
        for iv, vv in enumerate(coord):
            for ih, hh in enumerate(coord):
                p = list(center)
                p[_AXIS[ha]] = center[_AXIS[ha]] + hh
                p[_AXIS[va]] = center[_AXIS[va]] + vv
                grid[iv, ih] = dens(*p)
        grids.append(grid)

    # legend = only the materials that actually appear in the slices
    present = sorted({round(float(v), 6) for g in grids for v in np.unique(g)}) or [1.0]

    def label_for(d):
        if d in name_of:
            return name_of[d]
        keys = list(name_of)
        if not keys:
            return "%.3g" % d
        k = min(keys, key=lambda x: abs(math.log((x or 1e-9) / (d or 1e-9))))
        return name_of.get(k, "%.3g" % d)

    labels = [label_for(d) for d in present]
    colors = cm.turbo(np.linspace(0.04, 0.96, len(present)))
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(len(present) + 1) - 0.5, len(present))
    idx_of = {d: i for i, d in enumerate(present)}

    def to_idx(g):
        out = np.zeros(g.shape, dtype=int)
        R = np.round(g, 6)
        for d, i in idx_of.items():
            out[R == d] = i
        return out

    fig, axes = plt.subplots(1, len(planes), figsize=(6.4 * len(planes), 6.4),
                             squeeze=False)
    for ax, pl, grid in zip(axes[0], planes, grids):
        ha, va, fa = _PLANES[pl]
        ax.imshow(to_idx(grid), origin="lower", cmap=cmap, norm=norm,
                  extent=[coord[0], coord[-1], coord[0], coord[-1]],
                  aspect="equal", interpolation="nearest")
        ax.set_xlabel(_AXIS_LABEL[ha] + " [m]")
        ax.set_ylabel(_AXIS_LABEL[va] + " [m]")
        ax.set_title("%s slice  (%s = %.3g m)" % (pl, fa, center[_AXIS[fa]]))
        ax.axhline(0, color="w", lw=0.4, alpha=0.3)
        ax.axvline(0, color="w", lw=0.4, alpha=0.3)

    handles = [Patch(facecolor=colors[i], edgecolor="k",
                     label="%s  (%.3g g/cm^3)" % (labels[i], present[i]))
               for i in range(len(present))]
    fig.legend(handles=handles, loc="lower center",
               ncol=min(len(present), 6), fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    if title:
        fig.suptitle(title)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    if path:
        fig.savefig(path, dpi=120, bbox_inches="tight")
    if show:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# 2. standard-GDML export (from sectors)
# ---------------------------------------------------------------------------
def _fmt(x):
    return repr(float(x))


def _sanitize(name):
    """A GDML-safe identifier from an arbitrary sector/material name."""
    out = "".join(c if (c.isalnum() or c == "_") else "_" for c in str(name))
    return out or "x"


def _solid_xml(geo, name):
    """Return GDML <solid> XML for *geo*, or None to fall back to a bbox."""
    cls = type(geo).__name__
    if cls == "Box":
        return ('<box name="%s" lunit="m" x="%s" y="%s" z="%s"/>'
                % (name, _fmt(geo.X), _fmt(geo.Y), _fmt(geo.Z)))
    if cls == "Sphere":
        return ('<sphere name="%s" lunit="m" aunit="rad" rmin="%s" rmax="%s" '
                'startphi="%s" deltaphi="%s" starttheta="%s" deltatheta="%s"/>'
                % (name, _fmt(geo.InnerRadius), _fmt(geo.Radius),
                   _fmt(geo.StartPhi), _fmt(geo.DeltaPhi),
                   _fmt(geo.StartTheta), _fmt(geo.DeltaTheta)))
    if cls == "Cylinder":
        return ('<tube name="%s" lunit="m" aunit="rad" rmin="%s" rmax="%s" '
                'z="%s" startphi="0" deltaphi="%s"/>'
                % (name, _fmt(geo.InnerRadius), _fmt(geo.Radius),
                   _fmt(geo.Z), _fmt(2 * math.pi)))
    if cls == "Cone":
        return ('<cone name="%s" lunit="m" aunit="rad" rmin1="%s" rmax1="%s" '
                'rmin2="%s" rmax2="%s" z="%s" startphi="0" deltaphi="%s"/>'
                % (name, _fmt(geo.Rmin1), _fmt(geo.Rmax1), _fmt(geo.Rmin2),
                   _fmt(geo.Rmax2), _fmt(geo.Z), _fmt(2 * math.pi)))
    if cls == "Polycone":
        zs, rmin, rmax = geo.ZPlanes, geo.Rmin, geo.Rmax
        planes = "".join(
            '<zplane rmin="%s" rmax="%s" z="%s"/>' % (_fmt(rmin[i]), _fmt(rmax[i]), _fmt(zs[i]))
            for i in range(len(zs)))
        return ('<polycone name="%s" lunit="m" aunit="rad" startphi="0" '
                'deltaphi="%s">%s</polycone>' % (name, _fmt(2 * math.pi), planes))
    return None


def _bbox(geo):
    """(center, half-sizes) of the world-frame AABB (``min_corner``/``max_corner``)."""
    try:
        bb = geo.GetWorldBoundingBox()
        lo, hi = bb.min_corner, bb.max_corner
        lo = (lo.GetX(), lo.GetY(), lo.GetZ())
        hi = (hi.GetX(), hi.GetY(), hi.GetZ())
    except Exception:
        return None
    center = tuple((lo[i] + hi[i]) / 2 for i in range(3))
    half = tuple(max((hi[i] - lo[i]) / 2, 1e-6) for i in range(3))
    return center, half


def _placement(geo):
    """(position xyz, euler xyz) for the GDML physvol (passive rotation)."""
    pl = geo.placement
    pos = pl.Position
    pos = (pos.GetX(), pos.GetY(), pos.GetZ())
    # GDML <physvol> rotation is passive -> use the inverse rotation's angles.
    rx, ry, rz = (~pl.Quaternion).GetEulerAnglesXYZs()
    return pos, (rx, ry, rz)


def to_gdml(model, path, world_margin=2.0, region=None, region_center=None):
    """Export *model*'s sectors to a self-contained standard GDML file.

    Each sector becomes a ``<volume>`` named after the sector, referencing a
    material (named after the SIREN material, carrying the sector density) and a
    solid. **Identical solid shapes are emitted once and shared** -- the many
    replicated placements in the composite (decay-pipe segments, absorber
    plates, the NuMI beamline) would otherwise produce hundreds of duplicate
    solid definitions. Box/Sphere/Cylinder/Cone/Polycone are exact; any other
    solid is approximated by its bounding box; unbounded sectors (the infinite
    universe-boundary sphere) are skipped.

    :param region: if given, a half-size (m) of a cube about *region_center* --
        only sectors whose centre is inside are exported. Use this to drop the
        huge site-geology blocks / far beamline and keep just the detector
        region (much smaller, faster, and avoids a pyg4ometry picker crash on
        very large volumes).
    :param region_center: (x, y, z) in geometry coords for the crop centre;
        defaults to the model's ``DetectorOrigin``.
    :returns: dict with ``n_sectors``, ``n_exact``, ``n_bbox``, ``n_skipped``,
        ``n_cropped``, ``n_solids`` (distinct), ``n_materials``, ``path``.
    """
    import re
    sectors = _sectors(model)
    crop = None
    if region is not None:
        if region_center is None:
            try:
                o = model.DetectorOrigin.get()
                region_center = (o.GetX(), o.GetY(), o.GetZ())
            except Exception:
                region_center = (0.0, 0.0, 0.0)
        crop = (region_center, float(region))
    solid_cache, solid_defs = {}, []       # signature -> name ; emitted <solid>s
    mat_cache, mat_defs, used = {}, [], set()
    structs, physvols = [], []
    lo, hi = [1e30] * 3, [-1e30] * 3
    stats = dict(n_sectors=len(sectors), n_exact=0, n_bbox=0, n_skipped=0,
                 n_cropped=0)

    def share_solid(xml_tmp):
        sig = re.sub(r'\sname="[^"]*"', "", xml_tmp, count=1)
        if sig not in solid_cache:
            nm = "sol%d" % len(solid_cache)
            solid_cache[sig] = nm
            # replace the placeholder name with the canonical shared name
            solid_defs.append("    " + re.sub(r'name="[^"]*"', 'name="%s"' % nm,
                                              xml_tmp, count=1))
        return solid_cache[sig]

    for i, s in enumerate(sectors):
        geo = s.geo
        bb = _bbox(geo)
        if bb is None or not all(math.isfinite(c) for c in (bb[0] + bb[1])):
            stats["n_skipped"] += 1
            continue
        if crop is not None:
            (cx, cy, cz), _ = bb
            (rx, ry, rz), rh = crop
            if abs(cx - rx) > rh or abs(cy - ry) > rh or abs(cz - rz) > rh:
                stats["n_cropped"] += 1
                continue
        # material: dedup by id, named after the SIREN material + its density
        mid = int(s.material_id)
        if mid not in mat_cache:
            base = _sanitize(_material_name(model, mid))
            nm = base if base not in used else "%s_%d" % (base, mid)
            used.add(nm)
            mat_cache[mid] = nm
            mat_defs.append('    <material name="%s" Z="1"><D value="%s" unit="g/cm3"/>'
                            '<atom value="1.008"/></material>'
                            % (nm, _fmt(max(_sector_density(s), 1e-30))))
        # solid: exact where possible, else bbox; placement carried on the physvol
        try:
            xml = _solid_xml(geo, "TMP")
        except Exception:
            xml = None
        if xml is not None and ("inf" in xml or "nan" in xml):
            xml = None
        if xml is not None:
            pos, rot = _placement(geo)
            stats["n_exact"] += 1
        else:
            (cx, cy, cz), (hx, hy, hz) = bb
            xml = '<box name="TMP" lunit="m" x="%s" y="%s" z="%s"/>' % (
                _fmt(hx), _fmt(hy), _fmt(hz))
            pos, rot = (cx, cy, cz), None
            stats["n_bbox"] += 1
        sname = share_solid(xml)
        base = _sanitize(s.name)
        rxml = ('<rotation unit="rad" x="%s" y="%s" z="%s"/>'
                % (_fmt(rot[0]), _fmt(rot[1]), _fmt(rot[2]))) if rot else ""
        structs.append('    <volume name="%s__%d"><materialref ref="%s"/>'
                       '<solidref ref="%s"/></volume>' % (base, i, mat_cache[mid], sname))
        physvols.append('      <physvol name="%s__pv%d"><volumeref ref="%s__%d"/>'
                        '<position unit="m" x="%s" y="%s" z="%s"/>%s</physvol>'
                        % (base, i, base, i, _fmt(pos[0]), _fmt(pos[1]), _fmt(pos[2]), rxml))
        (cx, cy, cz), (hx, hy, hz) = bb
        lo = [min(lo[k], (cx, cy, cz)[k] - (hx, hy, hz)[k]) for k in range(3)]
        hi = [max(hi[k], (cx, cy, cz)[k] + (hx, hy, hz)[k]) for k in range(3)]

    stats["n_solids"] = len(solid_cache)
    stats["n_materials"] = len(mat_cache)

    if hi[0] < lo[0]:
        lo, hi = [-1, -1, -1], [1, 1, 1]
    # World box is centred at the origin, so its half-extent on each axis must
    # reach the farthest sector (plus a margin). GDML <box> x/y/z are FULL
    # lengths, so double it -- else the world is half the size it needs to be and
    # the outer sectors poke out of it.
    wh = [2.0 * (max(abs(lo[k]), abs(hi[k])) + world_margin) for k in range(3)]

    mat_defs.insert(0, '    <material name="WORLD_VAC" Z="1"><D value="1e-25" '
                       'unit="g/cm3"/><atom value="1.008"/></material>')
    gdml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n<gdml>\n  <define/>\n'
        '  <materials>\n%s\n  </materials>\n  <solids>\n%s\n'
        '    <box name="WORLD" lunit="m" x="%s" y="%s" z="%s"/>\n  </solids>\n'
        '  <structure>\n%s\n    <volume name="World">\n'
        '      <materialref ref="WORLD_VAC"/>\n      <solidref ref="WORLD"/>\n'
        '%s\n    </volume>\n  </structure>\n'
        '  <setup name="Default" version="1.0"><world ref="World"/></setup>\n</gdml>\n'
        % ("\n".join(mat_defs), "\n".join(solid_defs),
           _fmt(wh[0]), _fmt(wh[1]), _fmt(wh[2]),
           "\n".join(structs), "\n".join(physvols)))
    with open(path, "w") as f:
        f.write(gdml)
    stats["path"] = path
    return stats


# ---------------------------------------------------------------------------
# 3. headless inspection
# ---------------------------------------------------------------------------
def describe(model, printout=True):
    """Summarise the geometry by material: (name, density, sector count).

    A viewer-agnostic way to see what is in the model and at what density.
    :returns: list of (material_name, density_g_cm3, n_sectors), most-used first.
    """
    rows = {}
    for s in _sectors(model):
        mid = int(s.material_id)
        if mid not in rows:
            rows[mid] = [_material_name(model, mid), _sector_density(s), 0]
        rows[mid][2] += 1
    table = sorted((tuple(r) for r in rows.values()), key=lambda r: -r[2])
    if printout:
        print("%-30s %14s %9s" % ("material", "density[g/cm^3]", "sectors"))
        print("-" * 56)
        for nm, d, c in table:
            print("%-30s %14.5g %9d" % (nm, d, c))
        print("-" * 56)
        print("%-30s %14s %9d" % ("TOTAL (%d materials)" % len(table), "", sum(r[2] for r in table)))
    return table


def at(model, x, y, z, printout=True):
    """Identify the volume at the geometry-frame point (*x*, *y*, *z*) [m].

    Queries SIREN's own geometry (``GetContainingSector`` / ``GetMassDensity``),
    so it is exact and robust for *every* volume -- including the large outer
    ones the VTK picker crashes on. :returns: (name, material, density).
    """
    from siren.detector import GeometryPosition
    from siren.math import Vector3D
    dp = model.GeoPositionToDetPosition(
        GeometryPosition(Vector3D(float(x), float(y), float(z))))
    sec = model.GetContainingSector(dp)
    nm, mat, d = sec.name, _material_name(model, sec.material_id), model.GetMassDensity(dp)
    if printout:
        print("(%g, %g, %g) m  ->  %s   material=%s   rho=%.5g g/cm^3"
              % (x, y, z, nm, mat, d))
    return nm, mat, d


def list_volumes(model, sort="size", n=40, printout=True):
    """List the model's sectors: name, material, density, centre, size (m).

    Sorted by *sort* -- ``"size"`` (largest first, the default) is handy for
    verifying the big outer volumes (site geology, beamline enclosures) by their
    numbers instead of the crash-prone 3-D picker; ``"name"`` sorts
    alphabetically. Positions/sizes are in geometry-frame metres.
    :returns: list of dicts (name, material, density, center, size, maxdim).
    """
    rows = []
    for s in _sectors(model):
        bb = _bbox(s.geo)
        if bb is None:
            continue
        (cx, cy, cz), (hx, hy, hz) = bb
        if not all(math.isfinite(v) for v in (cx, cy, cz, hx, hy, hz)):
            continue
        rows.append(dict(name=s.name, material=_material_name(model, s.material_id),
                         density=_sector_density(s), center=(cx, cy, cz),
                         size=(2 * hx, 2 * hy, 2 * hz), maxdim=2 * max(hx, hy, hz)))
    if sort == "size":
        rows.sort(key=lambda r: -r["maxdim"])
    elif sort == "name":
        rows.sort(key=lambda r: r["name"])
    shown = rows if n is None else rows[:n]
    if printout:
        print("%-26s %-16s %9s  %-26s %-24s" %
              ("volume", "material", "rho", "centre (m)", "size  L x W x H (m)"))
        print("-" * 108)
        for r in shown:
            c, sz = r["center"], r["size"]
            print("%-26.26s %-16.16s %9.4g  (%8.1f %8.1f %8.1f)  (%7.2f %7.2f %7.2f)"
                  % (r["name"], r["material"], r["density"],
                     c[0], c[1], c[2], sz[0], sz[1], sz[2]))
        if n is not None and len(rows) > n:
            print("... %d more (pass n=None for all)" % (len(rows) - n))
    return rows


# ---------------------------------------------------------------------------
# 4. pyg4ometry 3-D view
# ---------------------------------------------------------------------------
def _material_vis_options(registry, VisOpt):
    """Map every material in *registry* to a colour so the viewer is not all
    grey (the predefined-material viewer only colours known names like G4_*).

    Colour runs by log-density (blue = low -> red = high) and gases (rho < 0.05)
    are made semi-transparent so you can see through air to the volumes inside.
    """
    import colorsys
    dens = {}
    for nm, mt in registry.materialDict.items():
        try:
            dens[nm] = max(float(mt.density), 1e-30)
        except Exception:
            dens[nm] = 1.0
    hi = math.log(max(dens.values()))
    lo = math.log(1e-4)  # ~below air; vacuum/gases clamp to the blue end
    out = {}
    for nm, d in dens.items():
        t = 0.5 if hi <= lo else min(max((math.log(d) - lo) / (hi - lo), 0.0), 1.0)
        r, g, b = colorsys.hsv_to_rgb(0.66 * (1.0 - t), 0.85, 0.95)
        out[nm] = VisOpt(colour=[round(r, 3), round(g, 3), round(b, 3)],
                         alpha=(0.12 if d < 0.05 else 0.9))
    return out


# pyg4ometry/VTK work in millimetres; SIREN geometry is in metres. The exported
# GDML carries lunit="m", which the reader scales by this factor when meshing, so
# a VTK world coordinate divided by this is a SIREN geometry-frame metre.
_GDML_MM_PER_M = 1000.0
# Materials at or below this density (g/cm^3) are treated as gases/vacuum: drawn
# semi-transparent and toggled together by the 'g' key (matches the air cut in
# _material_vis_options so you can see through air to the volumes inside).
_LOW_DENSITY = 0.05

_HELP_TEXT = (
    "Controls\n"
    "  right-click       identify volume (name / material / density)\n"
    "  shift+right-click hide the clicked material\n"
    "  v                 restore all materials and opacity\n"
    "  g                 toggle gas / low-density volumes\n"
    "  n / m             less / more opacity (see inside)\n"
    "  c                 toggle x/y/z section outlines\n"
    "  b                 toggle world bounding box\n"
    "  l                 toggle material legend\n"
    "  o                 toggle orientation cube\n"
    "  k                 toggle clip-plane widget\n"
    "  w / s             wireframe / surface\n"
    "  r                 reset camera\n"
    "  e / q             quit\n"
    "  h                 toggle this help")


def _load_registry(model, gdml_path):
    """(vis module, registry, world LV, gdml path) for *model* (or a GDML str)."""
    import pyg4ometry
    import tempfile
    vis = pyg4ometry.visualisation
    if isinstance(model, str):
        path = model
    else:
        path = gdml_path or tempfile.NamedTemporaryFile(
            suffix=".gdml", delete=False).name
        to_gdml(model, path)
    reg = pyg4ometry.gdml.Reader(path).getRegistry()
    return vis, reg, reg.getWorldVolume(), path


def _build_viewer(model, gdml_path=None, coloured=True):
    """Export+load *model* and return (vis, viewer, registry, path, mvo) with the
    world logical volume added. The append pipeline is *not* built yet."""
    vis, reg, world, path = _load_registry(model, gdml_path)
    mvo = _material_vis_options(reg, vis.VisualisationOptions) if coloured else None
    if coloured:
        v = vis.VtkViewerColouredNew(materialVisOptions=mvo,
                                     defaultColour=[0.6, 0.6, 0.6])
    else:
        v = vis.VtkViewerNew()
    v.addLogicalVolume(world)
    return vis, v, reg, path, mvo


def _scene_bounds(viewer):
    """Aggregate (xmin,xmax,ymin,ymax,zmin,zmax) over the geometry actors, or None.

    Call before any 2-D overlay actors are added -- those are attached to the
    renderer, not ``viewer.actors``, so they never enter this aggregate anyway.
    """
    lo, hi = [1e30] * 3, [-1e30] * 3
    for a in viewer.actors.values():
        try:
            b = a.GetBounds()
        except Exception:
            continue
        if b is None:
            continue
        for i in range(3):
            if b[2 * i] <= b[2 * i + 1] and math.isfinite(b[2 * i]) and math.isfinite(b[2 * i + 1]):
                lo[i] = min(lo[i], b[2 * i])
                hi[i] = max(hi[i], b[2 * i + 1])
    if hi[0] < lo[0]:
        return None
    return (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2])


def _text_actor(vtk, text, nx, ny, size=14, color=(0.0, 0.0, 0.0), anchor="bl"):
    """A monospaced 2-D overlay ``vtkTextActor`` anchored to a window corner.

    *nx*/*ny* are normalized-viewport coordinates; *anchor* picks the corner the
    text grows from ("tl"/"tr"/"bl"/"br") so it stays put on window resize.
    """
    a = vtk.vtkTextActor()
    a.SetInput(text)
    tp = a.GetTextProperty()
    tp.SetFontFamilyToCourier()
    tp.SetFontSize(size)
    tp.SetColor(*color)
    tp.SetVerticalJustificationToTop() if "t" in anchor else tp.SetVerticalJustificationToBottom()
    tp.SetJustificationToRight() if "r" in anchor else tp.SetJustificationToLeft()
    pc = a.GetPositionCoordinate()
    pc.SetCoordinateSystemToNormalizedViewport()
    pc.SetValue(nx, ny)
    return a


def _legend_actor(vtk, reg, mvo, max_entries=18):
    """A ``vtkLegendBoxActor`` of material -> colour swatch, densest first, or None."""
    items = []
    for name, vo in mvo.items():
        if name == "WORLD_VAC":
            continue
        try:
            d = float(reg.materialDict[name].density)
        except Exception:
            d = 0.0
        items.append((name, d, list(vo.getColour())))
    items.sort(key=lambda t: -t[1])
    items = items[:max_entries]
    if not items:
        return None
    cube = vtk.vtkCubeSource()
    cube.Update()
    sym = cube.GetOutput()
    leg = vtk.vtkLegendBoxActor()
    leg.SetNumberOfEntries(len(items))
    for i, (name, d, col) in enumerate(items):
        leg.SetEntry(i, sym, "%-.20s %.3g" % (name, d), col)
    leg.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    leg.GetPositionCoordinate().SetValue(0.8, 0.05)
    leg.GetPosition2Coordinate().SetCoordinateSystemToNormalizedViewport()
    leg.GetPosition2Coordinate().SetValue(0.2, 0.92)
    leg.UseBackgroundOn()
    leg.SetBackgroundColor(0.12, 0.12, 0.12)
    leg.SetBackgroundOpacity(0.6)
    leg.BorderOn()
    return leg


def _bbox_actor(vtk, bounds, color=(0.0, 0.0, 0.0)):
    """A wireframe ``vtkOutlineSource`` actor for the scene AABB *bounds*."""
    src = vtk.vtkOutlineSource()
    src.SetBounds(*bounds)
    m = vtk.vtkPolyDataMapper()
    m.SetInputConnection(src.GetOutputPort())
    a = vtk.vtkActor()
    a.SetMapper(m)
    a.GetProperty().SetColor(*color)
    a.GetProperty().SetLineWidth(1.0)
    a.PickableOff()
    return a


def _world_to_sector(model, wx_mm, wy_mm, wz_mm, inward=None):
    """Map a VTK world point (mm) to SIREN's (xyz [m], name, material, density).

    Converts mm->m and queries SIREN's own geometry -- robust for *every* volume,
    including the huge outer site-geology blocks the pyg4ometry VTK picker
    segfaults on.

    A cell-picker hit lands exactly *on* a surface. The point-only
    ``GetContainingSector`` (BVH + ``IsInside``) is ambiguous there -- it can
    return the clicked volume or the one just outside it. When *inward* (a
    geometry-frame direction pointing into the clicked volume, i.e. camera->point)
    is given we instead build the ray's intersection list and let SIREN's
    direction-aware ``GetContainingSector(intersections, p0)`` select the sector
    the ray is *entering* -- the front (clicked) volume, unambiguously and with no
    geometric fudge. Density is read straight off that sector
    (``sector.density.Evaluate`` at the point), NOT via
    ``GetMassDensity(intersections, p0)``: that overload tie-breaks the boundary
    the other way (it picks the sector the ray is *exiting* at p0, often the outer
    vacuum/air -- DetectorModel.cxx:845 vs :1372), so it would report a density
    inconsistent with the volume above. Without *inward* it falls back to the
    point-only query (as :func:`at` does).
    """
    from siren.detector import GeometryPosition, GeometryDirection
    from siren.math import Vector3D
    g = (wx_mm / _GDML_MM_PER_M, wy_mm / _GDML_MM_PER_M, wz_mm / _GDML_MM_PER_M)
    dp = model.GeoPositionToDetPosition(GeometryPosition(Vector3D(*g)))
    if inward is not None:
        n = (inward[0] ** 2 + inward[1] ** 2 + inward[2] ** 2) ** 0.5
        if n > 0:
            u = (inward[0] / n, inward[1] / n, inward[2] / n)
            ddir = model.GeoDirectionToDetDirection(GeometryDirection(Vector3D(*u)))
            inter = model.GetIntersections(dp, ddir)
            sec = model.GetContainingSector(inter, dp)
            try:
                density = sec.density.Evaluate(Vector3D(*g))
            except Exception:
                density = model.GetMassDensity(dp)
            return g, sec.name, _material_name(model, sec.material_id), density
    sec = model.GetContainingSector(dp)
    return g, sec.name, _material_name(model, sec.material_id), model.GetMassDensity(dp)


def _install_controls(vtk, viewer, reg, mvo, model, bounds,
                      legend=True, bounding_box=False, clipper_widget=None):
    """Attach overlays + a robust interactor style to an already-built *viewer*.

    Installs a custom ``vtkInteractorStyleTrackballCamera`` that (a) replaces
    pyg4ometry's segfault-prone right-click picker with a SIREN geometry query
    (see :func:`_world_to_sector`) and (b) adds keyboard toggles (see
    ``_HELP_TEXT``). *model* is the ``DetectorModel`` captured for the picker, or
    None to disable identification (right-click then just reports the 3-D point).
    Returns the state dict (mostly for testing).
    """
    ren, iren = viewer.ren, viewer.iren

    # ---- overlay actors -------------------------------------------------
    pick_txt = _text_actor(vtk, "right-click a volume to identify it",
                           0.01, 0.01, size=15, anchor="bl")
    help_txt = _text_actor(vtk, _HELP_TEXT, 0.01, 0.99, size=14, anchor="tl")
    help_txt.SetVisibility(False)
    hint_txt = _text_actor(vtk, "h: help   right-click: identify",
                           0.99, 0.99, size=13, color=(0.3, 0.3, 0.3), anchor="tr")
    for a in (pick_txt, help_txt, hint_txt):
        ren.AddViewProp(a)

    # ---- material -> body-actor map (coloured viewer only) --------------
    # buildPipelinesAppend keys each merged body actor by str(visOptions); since
    # we built mvo, mvo[name] is that exact options object, so str() matches.
    mat_actors, body_actors, gas_actors, orig_opacity = {}, [], [], {}
    if mvo is not None:
        seen = set()
        for name, vo in mvo.items():
            if name == "WORLD_VAC":
                continue
            act = viewer.actors.get(str(vo))
            if act is None:
                continue
            mat_actors[name] = act
            if id(act) not in seen:
                seen.add(id(act))
                body_actors.append(act)
                orig_opacity[act] = act.GetProperty().GetOpacity()
                try:
                    if float(reg.materialDict[name].density) <= _LOW_DENSITY:
                        gas_actors.append(act)
                except Exception:
                    pass

    # x/y/z section outlines pyg4ometry draws by default (keys end in _yz/_xz/_xy)
    cutter_actors = [a for k, a in viewer.actors.items()
                     if k.endswith(("_yz", "_xz", "_xy"))]

    legend_actor = None
    if legend and mvo is not None:
        legend_actor = _legend_actor(vtk, reg, mvo)
        if legend_actor is not None:
            ren.AddViewProp(legend_actor)

    bbox_actor = None
    if bounds is not None:
        bbox_actor = _bbox_actor(vtk, bounds)
        bbox_actor.SetVisibility(bool(bounding_box))
        ren.AddViewProp(bbox_actor)

    picker = vtk.vtkCellPicker()
    picker.SetTolerance(0.0005)

    state = dict(pick=pick_txt, help=help_txt, hint=hint_txt, legend=legend_actor,
                 bbox=bbox_actor, cutters=cutter_actors, mat_actors=mat_actors,
                 body_actors=body_actors, gas_actors=gas_actors,
                 orig_opacity=orig_opacity, hidden=set(), opacity=1.0,
                 gas_visible=True, cutters_visible=True, picker=picker,
                 clipper_widget=clipper_widget)

    def render():
        ren.GetRenderWindow().Render()

    def on_right(obj, event):
        x, y = iren.GetEventPosition()
        if not state["picker"].Pick(x, y, 0, ren) or state["picker"].GetActor() is None:
            pick_txt.SetInput("(nothing under cursor)")
            render()
            return
        wx, wy, wz = state["picker"].GetPickPosition()
        if model is None:
            pick_txt.SetInput("point  (%.2f, %.2f, %.2f) m\n(no model: identify unavailable)"
                              % (wx / _GDML_MM_PER_M, wy / _GDML_MM_PER_M, wz / _GDML_MM_PER_M))
            render()
            return
        # The pick lands exactly on a surface (a sector boundary). Pass the view
        # ray (camera->point, which continues into the front-most volume) so the
        # direction-aware GetContainingSector returns the clicked volume, not the
        # one just outside that surface.
        cam = ren.GetActiveCamera()
        cpx, cpy, cpz = cam.GetPosition()
        inward = (wx - cpx, wy - cpy, wz - cpz)
        try:
            (gx, gy, gz), name, material, density = _world_to_sector(
                model, wx, wy, wz, inward=inward)
        except Exception as e:
            pick_txt.SetInput("query failed: %s" % e)
            render()
            return
        print("[siren.view] (%.2f, %.2f, %.2f) m -> %s  [%s, rho=%.4g g/cm^3]"
              % (gx, gy, gz, name, material, density))
        pick_txt.SetInput("picked  (%.2f, %.2f, %.2f) m\nvolume    %s\nmaterial  %s\n"
                          "density   %.4g g/cm^3" % (gx, gy, gz, name, material, density))
        if iren.GetShiftKey() and material in mat_actors:
            mat_actors[material].SetVisibility(False)
            state["hidden"].add(material)
        render()

    def on_key(obj, event):
        key = (iren.GetKeySym() or "").lower()
        changed = True
        if key == "h":
            help_txt.SetVisibility(not help_txt.GetVisibility())
        elif key == "b" and bbox_actor is not None:
            bbox_actor.SetVisibility(not bbox_actor.GetVisibility())
        elif key == "l" and legend_actor is not None:
            legend_actor.SetVisibility(not legend_actor.GetVisibility())
        elif key == "c" and cutter_actors:
            state["cutters_visible"] = not state["cutters_visible"]
            for a in cutter_actors:
                a.SetVisibility(state["cutters_visible"])
        elif key == "g":
            state["gas_visible"] = not state["gas_visible"]
            for a in gas_actors:
                a.SetVisibility(state["gas_visible"])
        elif key in ("n", "m"):
            state["opacity"] = min(1.0, max(0.05, state["opacity"] + (0.15 if key == "m" else -0.15)))
            for a in body_actors:
                a.GetProperty().SetOpacity(state["opacity"])
        elif key == "v":
            for a in body_actors:
                a.SetVisibility(True)
                a.GetProperty().SetOpacity(orig_opacity.get(a, 1.0))
            state.update(opacity=1.0, gas_visible=True)
            state["hidden"].clear()
        elif key == "o":
            w = getattr(viewer, "axesWidget", None)
            if w is not None:
                w.EnabledOff() if w.GetEnabled() else w.EnabledOn()
        elif key == "k" and clipper_widget is not None:
            try:
                clipper_widget.Off() if clipper_widget.GetEnabled() else clipper_widget.On()
            except Exception:
                pass
        else:
            changed = False
        if changed:
            render()

    class _SirenInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
        # Observing RightButtonPressEvent (instead of pyg4ometry's broken
        # MouseInteractorNamePhysicalVolume) both fixes the picker and frees us
        # from its vtkImplicitPolyDataDistance loop; left-drag still orbits.
        def __init__(self):
            self.AddObserver("RightButtonPressEvent", on_right)
            self.AddObserver("KeyPressEvent", on_key)

    style = _SirenInteractorStyle()
    style.SetDefaultRenderer(ren)
    iren.SetInteractorStyle(style)
    viewer.interactorStyle = style  # keep a reference so VTK does not GC it
    state["style"] = style
    return state


def view(model, gdml_path=None, screenshot=None, coloured=True, axes=True,
         cutter=False, clipper=False, clip_origin=(0.0, 0.0, 0.0),
         clip_normal=(1.0, 0.0, 0.0), legend=True, bounding_box=False,
         picker=True, interactive=True):
    """Render *model* in 3-D with pyg4ometry, coloured by material.

    Exports to a temporary (or *gdml_path*) standard GDML, loads it, builds the
    append pipeline, and opens an interactive window with a SIREN-backed picker
    and keyboard controls (see :func:`_install_controls` / ``_HELP_TEXT``).

    Right-click identifies the volume under the cursor by querying SIREN's own
    geometry (segfault-free for every volume, unlike pyg4ometry's built-in
    picker); shift+right-click hides the clicked material. Keys toggle the
    legend, bounding box, section outlines, gas visibility, opacity, etc.

    :param model: a ``DetectorModel`` (or a str path to an existing GDML).
    :param coloured: colour by material; else a plain single-colour viewer.
    :param axes: add an orientation axis triad (auto-scaled to the scene).
    :param cutter: add an interactive cutting-plane (section) widget.
    :param clipper: add an interactive clipping-plane widget to slice the scene
        open; plane given by *clip_origin* / *clip_normal*. Toggle with 'k'.
    :param legend: show an in-window material colour legend (toggle with 'l').
    :param bounding_box: start with the world bounding box shown (toggle 'b').
    :param picker: enable the right-click SIREN identification (needs a model,
        not a bare GDML path); the custom interactor -- which also replaces
        pyg4ometry's crash-prone built-in picker -- is installed regardless. The
        clicked surface is disambiguated by the view direction (the front volume
        is reported, not the one just outside it).
    :param screenshot: write an off-screen PNG instead of opening a window
        (uses the legacy coloured viewer, which supports ``exportScreenShot``;
        needs a display-capable VTK). On a headless box prefer
        :func:`plot_cross_sections` or :func:`to_gdml` + an external viewer.
    """
    import vtk
    vis, reg, world, path = _load_registry(model, gdml_path)
    mvo = _material_vis_options(reg, vis.VisualisationOptions) if coloured else None
    model_obj = None if isinstance(model, str) else model

    # Off-screen PNG: the "New" pipeline has no exportScreenShot, so use the
    # legacy coloured viewer (which also honours our material->colour map).
    if screenshot:
        lv = vis.VtkViewerColoured(materialVisOptions=mvo) if coloured else vis.VtkViewer()
        lv.addLogicalVolume(world)
        if axes:
            try:
                lv.addAxes(20.0)
            except Exception:
                pass
        lv.exportScreenShot(screenshot)
        return screenshot

    # Supply our own material->colour map; the predefined-material viewer leaves
    # custom names (MINERAL_OIL, env_GlacialTill, ...) grey.
    if coloured:
        v = vis.VtkViewerColouredNew(materialVisOptions=mvo,
                                     defaultColour=[0.6, 0.6, 0.6])
    else:
        v = vis.VtkViewerNew()
    v.addLogicalVolume(world)
    if clipper:
        v.addClipper(list(clip_origin), list(clip_normal), True)  # before build
    v.buildPipelinesAppend()

    bounds = _scene_bounds(v)
    clipper_widget = None
    if clipper:
        try:
            v.addClipperWidget()
            clipper_widget = v.clipperPlaneWidget
        except Exception:
            pass
    if axes:
        try:
            if bounds is not None:
                span = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
                origin = ((bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2,
                          (bounds[4] + bounds[5]) / 2)
                v.addAxes(length=0.12 * span, origin=origin)
            else:
                v.addAxes(length=20.0)
        except Exception:
            pass
    if cutter:
        try:
            v.addCutterWidget()
        except Exception:
            pass

    _install_controls(vtk, v, reg, mvo, model_obj if picker else None, bounds,
                       legend=legend, bounding_box=bounding_box,
                       clipper_widget=clipper_widget)
    v.view(interactive=interactive)
    return path


def to_web(model, path, gdml_path=None, single_instance=False):
    """Export *model* for the browser via pyg4ometry's own exporters.

    The *path* extension selects the format:

    * ``.gltf`` / ``.glb`` -- a glTF scene (``exportGLTFScene``). Geometry and
      structure are exact; note pyg4ometry assigns *random* PBR colours here.
    * ``.html`` -- a self-contained three.js page (``exportThreeJSScene``,
      which also emits the ``.gltf`` and ``.css`` beside it).

    Both paths need ``pip install pygltflib`` (the HTML path also needs jinja2);
    a clear error is raised if a dependency is missing. Unlike the interactive
    :func:`view`, this needs no display, so it works headless.
    :returns: *path*.
    """
    vis, v, reg, gpath, mvo = _build_viewer(model, gdml_path, coloured=True)
    low = str(path).lower()
    if low.endswith((".gltf", ".glb")):
        try:
            import pygltflib  # noqa: F401
        except Exception:
            raise RuntimeError("glTF export needs: pip install pygltflib")
        v.exportGLTFScene(path, singleInstance=single_instance)
    elif low.endswith((".html", ".htm")):
        try:
            import pygltflib  # noqa: F401
            import jinja2  # noqa: F401
        except Exception:
            raise RuntimeError("HTML export needs: pip install pygltflib jinja2")
        v.exportThreeJSScene(str(path).rsplit(".", 1)[0])
    else:
        raise ValueError("unknown web-export extension (use .gltf/.glb/.html): %r" % path)
    return path
