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
  coloured by material.
* :func:`describe` -- print a material summary (name, density, sector count) for
  headless inspection.

Units/conventions: SIREN geometry is in metres; ``Cylinder``/``Cone`` store the
*full* z-extent, so the exported GDML uses the standard convention
(``<tube>`` z = full length) and renders correctly in external tools -- unlike
SIREN's own hand-written sub-GDMLs, whose parser treats ``<tube>`` z as a half.
"""
import math

__all__ = ["plot_cross_sections", "to_gdml", "view", "describe"]

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


def to_gdml(model, path, world_margin=2.0):
    """Export *model*'s sectors to a self-contained standard GDML file.

    Each sector becomes a ``<volume>`` named after the sector, referencing a
    material (named after the SIREN material, carrying the sector density) and a
    solid. **Identical solid shapes are emitted once and shared** -- the many
    replicated placements in the composite (decay-pipe segments, absorber
    plates, the NuMI beamline) would otherwise produce hundreds of duplicate
    solid definitions. Box/Sphere/Cylinder/Cone/Polycone are exact; any other
    solid is approximated by its bounding box; unbounded sectors (the infinite
    universe-boundary sphere) are skipped.

    :returns: dict with ``n_sectors``, ``n_exact``, ``n_bbox``, ``n_skipped``,
        ``n_solids`` (distinct), ``n_materials``, ``path``.
    """
    import re
    sectors = _sectors(model)
    solid_cache, solid_defs = {}, []       # signature -> name ; emitted <solid>s
    mat_cache, mat_defs, used = {}, [], set()
    structs, physvols = [], []
    lo, hi = [1e30] * 3, [-1e30] * 3
    stats = dict(n_sectors=len(sectors), n_exact=0, n_bbox=0, n_skipped=0)

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
    # World box is centred at the origin; size it to contain every sector.
    wh = [max(abs(lo[k]), abs(hi[k])) + world_margin for k in range(3)]

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


def view(model, gdml_path=None, screenshot=None, coloured=True, axes=True,
         cutter=False, clipper=False, clip_origin=(0.0, 0.0, 0.0),
         clip_normal=(1.0, 0.0, 0.0), interactive=True):
    """Render *model* in 3-D with pyg4ometry, coloured by material.

    Uses the fuller ``VtkViewerColouredMaterialNew`` pipeline (interactive
    section/clip widgets, three.js/GLTF export). Exports to a temporary (or
    *gdml_path*) standard GDML, loads it, and opens an interactive window.

    :param model: a ``DetectorModel`` (or a str path to an existing GDML).
    :param coloured: colour by material; else a plain single-colour viewer.
    :param axes: add an orientation axis triad.
    :param cutter: add an interactive cutting-plane (section) widget.
    :param clipper: add an interactive clipping-plane widget to slice the scene
        open; plane given by *clip_origin* / *clip_normal*.
    :param screenshot: write an off-screen PNG instead of opening a window
        (uses the legacy coloured viewer, which supports ``exportScreenShot``;
        needs a display-capable VTK). On a headless box prefer
        :func:`plot_cross_sections` or :func:`to_gdml` + an external viewer.
    """
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
    world = reg.getWorldVolume()
    mvo = _material_vis_options(reg, vis.VisualisationOptions) if coloured else None

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
    # Use the APPEND pipeline: the right-click "name physical volume" picker
    # walks a vtkAppendPolyData. buildPipelinesSeparate's deeper pipeline makes
    # that picker over-traverse and raise AttributeError on right-click.
    if clipper:
        v.addClipper(list(clip_origin), list(clip_normal), True)  # before build
    v.buildPipelinesAppend()
    if clipper:
        try:
            v.addClipperWidget()
        except Exception:
            pass
    if axes:
        try:
            v.addAxes(length=20.0)
        except Exception:
            pass
    if cutter:
        try:
            v.addCutterWidget()
        except Exception:
            pass
    v.view(interactive=interactive)
    return path
