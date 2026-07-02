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
  pyg4ometry's own picker segfaults) and keyboard controls for fly/pan
  navigation, the legend, bounding box, section outlines, per-material
  visibility, and opacity.
* :func:`view_pv` (or ``view(backend="pyvista")``) -- pure-pyvista render with
  no pyg4ometry dependency. Every SIREN solid type -- including hollow/phi-cut
  solids and (nested) BooleanGeometry -- is tessellated parametrically; see
  the builder section around :func:`_pv_parametric`.
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

__all__ = ["plot_cross_sections", "to_gdml", "view", "view_pv", "to_web",
           "describe", "at", "list_volumes"]

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
    # turbo only exists in matplotlib >= 3.3; degrade to viridis on older ones
    cmap_src = getattr(cm, "turbo", cm.viridis)
    colors = cmap_src(np.linspace(0.04, 0.96, len(present)))
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
                % (name, _fmt(geo.X / 2), _fmt(geo.Y / 2), _fmt(geo.Z / 2)))
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


def to_gdml(model, path, world_margin=2.0, region=None, region_center=None,
            skip_geo_types=()):
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
                 n_cropped=0, n_mesh=0)

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
        # Meshes are skipped here when requested: pyg4ometry would CSG them
        # (slow) or, lacking a <tessellated> emitter, fall back to a bounding box.
        # The viewer renders them straight to VTK instead (see _add_mesh_actors).
        if type(geo).__name__ in skip_geo_types:
            stats["n_mesh"] += 1
            continue
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
            # _bbox returns HALF-extents; GDML <box> x/y/z are FULL lengths, so
            # double them -- else bbox-approximated solids (booleans, etc.) render
            # at half size, shrunk toward their centre and visibly misaligned with
            # the full-size exact solids around them.
            (cx, cy, cz), (hx, hy, hz) = bb
            xml = '<box name="TMP" lunit="m" x="%s" y="%s" z="%s"/>' % (
                _fmt(2 * hx), _fmt(2 * hy), _fmt(2 * hz))
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
# pyg4ometry meshes curved solids (Sphere/Orb/Tubs/Cons) with this many angular
# slices/stacks. The library default for a <sphere> is only 10x10, which makes
# the MiniBooNE oil shells visibly faceted and their intersections hard to read;
# raise it so spheres/pipes render smoothly. Higher = smoother but more polygons.
_MESH_SLICES = 48
# Near/far clipping is driven by the camera's distance to its focal point
# (cam.GetDistance()), NOT by the scene's total extent. VTK's default ties the
# near plane to the far plane (near ~ tol x far); with an Earth-sized model far
# is thousands of km, so even a tiny tolerance leaves the near plane kilometres
# out and you cannot zoom into metre- (let alone mm-) scale features. We instead
# put the near plane at this fraction of the view distance, so zooming in (which
# shrinks that distance) shrinks the near plane in lockstep. Depth resolution at
# the focal point stays ~ view_distance / (frac * 2^depthbits) -- a fixed
# fraction of whatever you are looking at, hence crisp at every scale.
_NEAR_DIST_FRAC = 1e-3
# Absolute near-plane floor (mm) so it never reaches zero when zoomed to microns.
_NEAR_FLOOR_MM = 1e-3

_HELP_TEXT = (
    "Controls\n"
    "  left-drag         orbit;  scroll = zoom\n"
    "  arrow keys        pan L/R; fly in/out (up/down)\n"
    "  shift+up/down     pan up / down\n"
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


def _load_registry(model, gdml_path, skip_geo_types=("TriangularMesh",)):
    """(vis module, registry, world LV, gdml path) for *model* (or a GDML str).

    *skip_geo_types* names geometry classes left out of the GDML (default the
    TriangularMesh, which the viewer renders directly in VTK via _add_mesh_actors
    rather than through pyg4ometry's CSG).
    """
    import pyg4ometry
    import tempfile
    vis = pyg4ometry.visualisation
    if isinstance(model, str):
        path = model
    else:
        path = gdml_path or tempfile.NamedTemporaryFile(
            suffix=".gdml", delete=False).name
        to_gdml(model, path, skip_geo_types=skip_geo_types)
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


def _apply_placement(tris, geo):
    """Transform LOCAL-frame triangle vertices to WORLD frame via the geo's placement."""
    import numpy as np
    pl = geo.placement
    pos = pl.Position
    t = np.array([pos.GetX(), pos.GetY(), pos.GetZ()])
    q = pl.Quaternion
    rx, ry, rz = q.GetEulerAnglesXYZs()
    cx, sx = np.cos(rx), np.sin(rx)
    cy, sy = np.cos(ry), np.sin(ry)
    cz, sz = np.cos(rz), np.sin(rz)
    R = (np.array([[cz, -sz, 0.], [sz, cz, 0.], [0., 0., 1.]])
         @ np.array([[cy, 0., sy], [0., 1., 0.], [-sy, 0., cy]])
         @ np.array([[1., 0., 0.], [0., cx, -sx], [0., sx, cx]]))
    pts = tris.reshape(-1, 3)
    pts = (R @ pts.T).T + t
    return pts.reshape(tris.shape)


def _mesh_sector_polydata(model, top_only=False):
    """Yield (sector, material_name, pyvista.PolyData) for each TriangularMesh sector."""
    import numpy as np
    import pyvista as pv
    for s in _sectors(model):
        if type(s.geo).__name__ != "TriangularMesh":
            continue
        tris = np.asarray(s.geo.GetTriangles(), dtype=float)        # (N, 3, 3)
        if tris.size == 0:
            continue
        tris = _apply_placement(tris, s.geo)                        # local -> world frame
        if top_only:                                                # drop skirt+base (frame-agnostic)
            bb = s.geo.GetWorldBoundingBox()
            lo = min(bb.min_corner.GetX(), bb.min_corner.GetY(), bb.min_corner.GetZ())
            mask = tris.min(axis=(1, 2)) > lo + 500
            tris = tris[mask]
        n = len(tris)
        faces = np.empty((n, 4), dtype=np.int64)
        faces[:, 0] = 3
        faces[:, 1:] = np.arange(3 * n).reshape(n, 3)
        pd = pv.PolyData(tris.reshape(-1, 3), faces.ravel())
        yield s, _material_name(model, int(s.material_id)), pd


def _add_mesh_actors(vtk, renderer, actors, model, mvo, reg):
    """Render TriangularMesh sectors straight to *renderer*, bypassing pyg4ometry.

    Mesh sectors are skipped in the GDML (see _load_registry); here their triangles
    (TriangularMesh.GetTriangles(), via pyvista PolyData) become VTK actors coloured
    on the same log-density scale as _material_vis_options.
    """
    import colorsys
    items = list(_mesh_sector_polydata(model))
    if not items:
        return 0
    dens = []
    for mt in reg.materialDict.values():
        try:
            dens.append(max(float(mt.density), 1e-30))
        except Exception:
            pass
    dens += [max(_sector_density(s), 1e-30) for s, _, _ in items]
    hi = math.log(max(dens)) if dens else 0.0
    lo = math.log(1e-4)

    def colour_for(name, d):
        if mvo is not None and name in mvo:
            try:
                return tuple(mvo[name].getColour())
            except Exception:
                pass
        t = 0.5 if hi <= lo else min(max((math.log(d) - lo) / (hi - lo), 0.0), 1.0)
        return colorsys.hsv_to_rgb(0.66 * (1.0 - t), 0.85, 0.95)

    for s, mat, pd in items:
        # pyg4ometry's VTK scene is in mm; our mesh polydata is in m.
        pd = pd.copy()
        pd.points = pd.points * _GDML_MM_PER_M
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(pd)
        mapper.ScalarVisibilityOff()
        a = vtk.vtkActor()
        a.SetMapper(mapper)
        a.GetProperty().SetColor(*colour_for(mat, max(_sector_density(s), 1e-30)))
        renderer.AddActor(a)
        actors["mesh__%s" % _sanitize(s.name)] = a
    return len(items)


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
                      legend=True, bounding_box=False, clipper_widget=None,
                      near_frac=_NEAR_DIST_FRAC):
    """Attach overlays + a robust interactor style to an already-built *viewer*.

    Installs a custom ``vtkInteractorStyleTrackballCamera`` that (a) replaces
    pyg4ometry's segfault-prone right-click picker with a SIREN geometry query
    (see :func:`_world_to_sector`) and (b) adds keyboard toggles (see
    ``_HELP_TEXT``). *model* is the ``DetectorModel`` captured for the picker, or
    None to disable identification (right-click then just reports the 3-D point).
    Returns the state dict (mostly for testing).
    """
    ren, iren = viewer.ren, viewer.iren

    # ---- distance-based near/far clipping -------------------------------
    # Drive the clip planes off the camera's distance to its focal point (see
    # _NEAR_DIST_FRAC) so small features stay visible as you zoom in, even inside
    # an Earth-sized model where the scene-extent-based default would clip
    # everything near the camera. Re-applied on every camera move by a
    # ModifiedEvent observer (added with the style below); the style's
    # AutoAdjustCameraClippingRange is turned off so it cannot fight this.
    _span = (max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
             if bounds is not None else 1000.0)
    _scene_R = 0.5 * _span
    _clip_busy = [False]

    def _update_clip(obj=None, event=None):
        if _clip_busy[0] or not near_frac or near_frac <= 0:
            return
        cam = ren.GetActiveCamera()
        d = cam.GetDistance()
        if not math.isfinite(d) or d <= 0:
            d = max(_scene_R, 1.0)
        near = max(d * near_frac, _NEAR_FLOOR_MM)
        far = d + 2.0 * _scene_R + near    # always reaches the far side of the scene
        cur = cam.GetClippingRange()
        if abs(cur[0] - near) <= 1e-6 * near and abs(cur[1] - far) <= 1e-6 * far:
            return                          # unchanged -> skip (avoids set storms)
        _clip_busy[0] = True                # guard: SetClippingRange re-fires ModifiedEvent
        try:
            cam.SetClippingRange(near, far)
        finally:
            _clip_busy[0] = False

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

    # ---- fly / pan keyboard navigation ----------------------------------
    # (_span / _scene_R are computed above for the clipping helper.)
    def _norm(v):
        n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) ** 0.5
        return (v[0] / n, v[1] / n, v[2] / n) if n > 1e-12 else (0.0, 0.0, 0.0)

    def _cross(a, b):
        return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0])

    def navigate(key, shift):
        """Translate the camera: arrows pan left/right and fly forward/back
        (up/down); shift+up/down pans vertically. The step scales with the view
        distance so it feels the same zoomed in or out. Both the camera and its
        focal point move together (true fly/pan, not a dolly toward the centre),
        so orientation is preserved. Returns True if the camera moved."""
        cam = ren.GetActiveCamera()
        pos, fp = cam.GetPosition(), cam.GetFocalPoint()
        fwd = _norm(cam.GetDirectionOfProjection())   # camera -> focal point
        right = _norm(_cross(fwd, cam.GetViewUp()))
        vup = _norm(_cross(right, fwd))
        dist = cam.GetDistance()
        if not math.isfinite(dist) or dist <= 0:
            dist = 0.5 * _span
        step = 0.12 * dist
        if key == "up":
            d = vup if shift else fwd
        elif key == "down":
            d = tuple(-c for c in (vup if shift else fwd))
        elif key == "left":
            d = tuple(-c for c in right)
        elif key == "right":
            d = right
        else:
            return False
        cam.SetPosition(*(pos[i] + d[i] * step for i in range(3)))
        cam.SetFocalPoint(*(fp[i] + d[i] * step for i in range(3)))
        _update_clip()
        return True

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
        elif key in ("up", "down", "left", "right"):
            changed = navigate(key, bool(iren.GetShiftKey()))
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
    # Manage clipping ourselves: stop the style's scene-extent reset from
    # overriding the distance-based range, then re-apply it on every camera
    # change (orbit, zoom, fly/pan, reset) and once now for the initial view.
    if near_frac and near_frac > 0:
        style.AutoAdjustCameraClippingRangeOff()
        ren.GetActiveCamera().AddObserver("ModifiedEvent", _update_clip)
        _update_clip()
    state["style"] = style
    return state


def view(model, gdml_path=None, screenshot=None, coloured=True, axes=True,
         cutter=False, clipper=False, clip_origin=(0.0, 0.0, 0.0),
         clip_normal=(1.0, 0.0, 0.0), legend=True, bounding_box=False,
         picker=True, interactive=True, mesh_slices=_MESH_SLICES,
         near_frac=_NEAR_DIST_FRAC, backend="pyg4ometry"):
    """Render *model* in 3-D with pyg4ometry, coloured by material.

    Exports to a temporary (or *gdml_path*) standard GDML, loads it, builds the
    append pipeline, and opens an interactive window with a SIREN-backed picker
    and keyboard controls (see :func:`_install_controls` / ``_HELP_TEXT``).

    Right-click identifies the volume under the cursor by querying SIREN's own
    geometry (segfault-free for every volume, unlike pyg4ometry's built-in
    picker); shift+right-click hides the clicked material. Keys toggle the
    legend, bounding box, section outlines, gas visibility, and opacity; the
    arrow keys pan and fly through the scene (press 'h' for the full list).

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
    :param mesh_slices: angular slices/stacks for curved solids (spheres, tubes).
        pyg4ometry's default for a sphere is only 10, which looks faceted; higher
        is smoother but adds polygons. 0/None keeps the library default.
    :param near_frac: near clip plane as a fraction of the camera's distance to
        its focal point (default 1e-3), set dynamically on every camera move.
        Decoupled from the scene's total size, so you can zoom into metre- and
        mm-scale features even in an Earth-sized model (where VTK's scene-extent
        default leaves the near plane kilometres out). Smaller lets the camera
        approach closer before clipping; larger gives more depth precision. Set
        0/None to keep VTK's default scene-extent clipping.
    :param screenshot: write an off-screen PNG instead of opening a window
        (uses the legacy coloured viewer, which supports ``exportScreenShot``;
        needs a display-capable VTK). On a headless box prefer
        :func:`plot_cross_sections` or :func:`to_gdml` + an external viewer.
    """
    if backend == "pyvista":                         # pure-pyvista path (no pyg4ometry)
        return view_pv(model, screenshot=screenshot, axes=axes, legend=legend,
                       picker=picker)
    import vtk
    import pyg4ometry
    # Smooth the curved solids before anything is meshed: pyg4ometry meshes a
    # <sphere> with only 10x10 slices/stacks by default, so the MiniBooNE oil
    # shells look faceted and their intersections are hard to read. Must precede
    # buildPipelinesAppend (meshing reads the global config).
    if mesh_slices:
        try:
            pyg4ometry.config.setGlobalMeshSliceAndStack(int(mesh_slices))
        except Exception:
            pass
    vis, reg, world, path = _load_registry(model, gdml_path)
    mvo = _material_vis_options(reg, vis.VisualisationOptions) if coloured else None
    model_obj = None if isinstance(model, str) else model

    # Off-screen PNG: the "New" pipeline has no exportScreenShot, so use the
    # legacy coloured viewer (which also honours our material->colour map).
    if screenshot:
        lv = vis.VtkViewerColoured(materialVisOptions=mvo) if coloured else vis.VtkViewer()
        lv.addLogicalVolume(world)
        if model_obj is not None:
            try:
                _add_mesh_actors(vtk, lv.ren, getattr(lv, "actors", {}), model_obj, mvo, reg)
            except Exception:
                pass
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

    if model_obj is not None:                        # mesh sectors -> direct VTK
        try:
            _add_mesh_actors(vtk, v.ren, v.actors, model_obj, mvo, reg)
        except Exception:
            pass

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
                       clipper_widget=clipper_widget, near_frac=near_frac)
    v.view(interactive=interactive)
    return path


# ---------------------------------------------------------------------------
# Pure-pyvista renderer (no pyg4ometry). Every SIREN solid type has a crisp
# parametric tessellation (see _PV_BUILDERS), including hollow/phi-cut solids
# and BooleanGeometry (operands are meshed recursively and combined by surface
# CSG). Marching cubes on SIREN's own IsInside remains only as a last-resort
# fallback for unknown types or failed builds; every fallback is counted in
# _MC_FALLBACKS and reported by view_pv, so regressions are visible.
# ---------------------------------------------------------------------------
# Angular segments per full circle for curved parametric solids (revolutions,
# tubes, spheres); partial arcs get a proportional share. Polyhedra ignore this
# (their side count is exact).
_PV_REV_RES = 48
def _euler_transform(euler, pos):
    """4x4 active transform: rotate by XYZ *euler* (rad), then translate to *pos*."""
    import numpy as np
    rx, ry, rz = euler
    cx, sx = math.cos(rx), math.sin(rx)
    cy, sy = math.cos(ry), math.sin(ry)
    cz, sz = math.cos(rz), math.sin(rz)
    R = (np.array([[cz, -sz, 0.], [sz, cz, 0.], [0., 0., 1.]])
         @ np.array([[cy, 0., sy], [0., 1., 0.], [-sy, 0., cy]])
         @ np.array([[1., 0., 0.], [0., cx, -sx], [0., sx, cx]]))
    T = np.eye(4); T[:3, :3] = R; T[:3, 3] = pos
    return T


def _pv_assemble(pts, faces):
    """PolyData from a vertex list + per-face index lists, triangulated.

    Faces may be polygons of any length >= 3 (vtkTriangleFilter handles
    non-convex polygons, e.g. the phi cap of a hollow polycone). Bit-identical
    duplicate points (revolution seams, cap rims) are merged so the result is
    a closed manifold wherever the construction provides one.
    """
    import numpy as np, pyvista as pv
    flat = []
    for f in faces:
        flat.append(len(f))
        flat.extend(f)
    pd = pv.PolyData(np.asarray(pts, float), faces=np.asarray(flat, np.int64))
    pd = pd.triangulate().clean()
    if pd.n_lines or pd.n_verts:        # degenerate cells demoted by clean()
        pd = pv.PolyData(pd.points, faces=pd.faces)
    return pd


def _dedup_face(f):
    """Drop consecutive repeated indices (collapsed quad corners on the axis)."""
    g = [f[0]]
    for v in f[1:]:
        if v != g[-1]:
            g.append(v)
    if len(g) > 1 and g[-1] == g[0]:
        g.pop()
    return g


def _rz_area(rz):
    """Twice the signed (shoelace) area of a closed (r, z) polygon."""
    a = 0.0
    for i in range(len(rz)):
        r0, z0 = rz[i]
        r1, z1 = rz[(i + 1) % len(rz)]
        a += r0 * z1 - r1 * z0
    return a


def _pv_revolve_rz(rz, start_phi=0.0, delta_phi=None, n=None, caps=True,
                   orient=True):
    """Solid of revolution of a closed polygon in the (r, z) half-plane.

    The contour (last vertex implicitly connects to the first; r >= 0) is
    normalised to counter-clockwise so every swept quad gets an outward
    normal -- the signed-distance CSG in :func:`_pv_csg` relies on that.
    Segments lying on the axis sweep to nothing and are skipped; on-axis
    vertices collapse to a single point so cones/poles stay watertight.
    *delta_phi* = None (or ~2*pi) revolves fully, wrapping indices at the seam
    (exact closure, no reliance on float coincidence); otherwise the two flat
    phi faces are added as polygon caps (*caps*) -- valid whenever the contour
    is a simple polygon. *n* = angular segments (default scales with the arc;
    pass the side count for prisms/polyhedra).
    """
    import numpy as np
    full = delta_phi is None or abs(delta_phi - 2.0 * math.pi) < 1e-9
    dphi = 2.0 * math.pi if full else float(delta_phi)
    if n is None:
        n = max(8, int(round(_PV_REV_RES * dphi / (2.0 * math.pi))))
    clean = []
    for r, z in rz:
        if not clean or abs(r - clean[-1][0]) > 1e-12 or abs(z - clean[-1][1]) > 1e-12:
            clean.append((float(r), float(z)))
    while len(clean) > 1 and abs(clean[0][0] - clean[-1][0]) < 1e-12 \
            and abs(clean[0][1] - clean[-1][1]) < 1e-12:
        clean.pop()
    rz = clean
    if len(rz) < 3:
        return None
    if orient and _rz_area(rz) < 0:
        rz.reverse()
    m = len(rz)
    ang = np.linspace(start_phi, start_phi + dphi, n + 1)
    ca, sa = np.cos(ang), np.sin(ang)
    pts, faces, index = [], [], {}

    def P(i, j):
        r, z = rz[i]
        if r < 1e-12:
            key = (-1, i)                       # on-axis: one vertex per entry
        else:
            if full and j == n:
                j = 0                           # wrap the seam exactly
            key = (i, j)
        k = index.get(key)
        if k is None:
            k = len(pts)
            pts.append((0.0, 0.0, z) if r < 1e-12 else (r * ca[j], r * sa[j], z))
            index[key] = k
        return k

    for i in range(m):
        i2 = (i + 1) % m
        if rz[i][0] < 1e-12 and rz[i2][0] < 1e-12:
            continue                            # axis segment sweeps to nothing
        for j in range(n):
            f = _dedup_face([P(i, j), P(i, j + 1), P(i2, j + 1), P(i2, j)])
            if len(f) >= 3:
                faces.append(f)
    if caps and not full:
        # CCW contour mapped into the (r_hat, z) plane has normal r_hat x z_hat
        # = -phi_hat: outward at the start face as-is, reversed at the far face.
        cap0 = _dedup_face([P(i, 0) for i in range(m)])
        cap1 = _dedup_face([P(i, n) for i in range(m)])
        if len(cap0) >= 3:
            faces.append(cap0)
        if len(cap1) >= 3:
            faces.append(cap1[::-1])
    return _pv_assemble(pts, faces)


def _pv_hexa(v):
    """Closed hexahedron from 8 corners ordered (-z: --, +-, ++, -+; then +z)."""
    return _pv_assemble(v, [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
                            (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)])


def _phi_arg(geo):
    """(start_phi, delta_phi-or-None) of a solid; None means a full circle."""
    dphi = geo.DeltaPhi
    return geo.StartPhi, (dphi if dphi < 2.0 * math.pi - 1e-9 else None)


def _pv_box(geo):
    import pyvista as pv
    return pv.Cube(center=(0, 0, 0), x_length=geo.X, y_length=geo.Y,
                   z_length=geo.Z)


def _pv_cylinder(geo):
    h = 0.5 * geo.Z                 # SIREN Cylinder stores the FULL z extent
    phi0, dphi = _phi_arg(geo)
    rz = [(geo.Radius, -h), (geo.Radius, h), (geo.InnerRadius, h),
          (geo.InnerRadius, -h)]
    return _pv_revolve_rz(rz, phi0, dphi)


def _pv_cone(geo):
    h = 0.5 * geo.Z                 # full z extent; Rm*1 at -h, Rm*2 at +h
    phi0, dphi = _phi_arg(geo)
    rz = [(geo.Rmax1, -h), (geo.Rmax2, h), (geo.Rmin2, h), (geo.Rmin1, -h)]
    return _pv_revolve_rz(rz, phi0, dphi)


def _pv_polycone(geo):
    zs, rmin, rmax = list(geo.ZPlanes), list(geo.Rmin), list(geo.Rmax)
    phi0, dphi = _phi_arg(geo)
    rz = [(rmax[i], zs[i]) for i in range(len(zs))]
    rz += [(rmin[i], zs[i]) for i in reversed(range(len(zs)))]
    return _pv_revolve_rz(rz, phi0, dphi)


def _pv_generic_polycone(geo):
    phi0, dphi = _phi_arg(geo)
    return _pv_revolve_rz(list(zip(geo.R, geo.Z)), phi0, dphi)


def _pv_polyhedra(geo):
    ns = int(geo.NumSides)
    if ns < 3:
        return None
    phi0, dphi = _phi_arg(geo)
    # rmin/rmax are APOTHEMS (distance to the flat side); the cross-section
    # vertices sit at the circumradius, at angles start_phi + k*face_dphi
    # (matching Polyhedra::precompute_trig).
    face_dphi = (dphi if dphi is not None else 2.0 * math.pi) / ns
    scale = 1.0 / math.cos(0.5 * face_dphi)
    zs, rmin, rmax = list(geo.ZPlanes), list(geo.Rmin), list(geo.Rmax)
    rz = [(rmax[i] * scale, zs[i]) for i in range(len(zs))]
    rz += [(rmin[i] * scale, zs[i]) for i in reversed(range(len(zs)))]
    return _pv_revolve_rz(rz, phi0, dphi, n=ns)


def _pv_sphere(geo):
    import numpy as np, pyvista as pv
    rin, rout = geo.InnerRadius, geo.Radius
    phi0, dphi = _phi_arg(geo)
    th0 = min(max(geo.StartTheta, 0.0), math.pi)
    th1 = min(max(th0 + geo.DeltaTheta, 0.0), math.pi)
    theta_cut = th0 > 1e-9 or th1 < math.pi - 1e-9
    if rin <= 1e-12 and dphi is None and not theta_cut:
        return pv.Sphere(radius=rout, theta_resolution=_PV_REV_RES,
                         phi_resolution=_PV_REV_RES)
    nt = max(8, int(round(_PV_REV_RES * (th1 - th0) / math.pi)))
    thetas = np.linspace(th0, th1, nt + 1)

    def arc(rad, t):                # snap sin(pi) ~ 1e-16 to a true pole point
        s = math.sin(t)
        return (0.0 if abs(s) < 1e-12 else rad * s, rad * math.cos(t))

    rz = [arc(rout, t) for t in thetas]
    if rin > 1e-12:
        rz += [arc(rin, t) for t in reversed(thetas)]
    elif theta_cut:
        rz.append((0.0, 0.0))       # close the sector through the apex
    return _pv_revolve_rz(rz, phi0, dphi)


def _pv_torus(geo):
    import numpy as np, pyvista as pv
    R, rout, rin = geo.MajorRadius, geo.MinorRadius, geo.InnerRadius
    phi0, dphi = _phi_arg(geo)
    ts = np.linspace(0.0, 2.0 * math.pi, _PV_REV_RES + 1)[:-1]
    outer = [(R + rout * math.cos(t), rout * math.sin(t)) for t in ts]
    if rin <= 1e-12:
        return _pv_revolve_rz(outer, phi0, dphi)
    # Hollow tube: the cross-section is an annulus that does not touch the
    # axis, so it cannot be a single simple contour. Revolve the two circles
    # separately (inner pre-reversed so its normals face the tube axis) and,
    # for a partial torus, close the ends with quad strips pairing the rings.
    inner = [(R + rin * math.cos(t), rin * math.sin(t)) for t in ts]
    parts = [_pv_revolve_rz(outer, phi0, dphi, caps=False),
             _pv_revolve_rz(inner[::-1], phi0, dphi, caps=False, orient=False)]
    if dphi is not None:
        for phi, flip in ((phi0, False), (phi0 + dphi, True)):
            c, s = math.cos(phi), math.sin(phi)
            pts, faces = [], []
            for r, z in outer + inner:
                pts.append((r * c, r * s, z))
            m = len(outer)
            for k in range(m):
                k1 = (k + 1) % m
                q = [k, k1, m + k1, m + k]      # outward -phi_hat at start
                faces.append(q[::-1] if flip else q)
            parts.append(_pv_assemble(pts, faces))
    parts = [p for p in parts if p is not None]
    if not parts:
        return None
    return pv.merge(parts).clean()


def _pv_cut_tube(geo):
    import numpy as np
    rmin, rmax, dz = geo.Rmin, geo.Rmax, geo.Dz     # dz is the HALF height
    ln, hn = geo.LowNorm, geo.HighNorm
    lnx, lny, lnz = ln.GetX(), ln.GetY(), ln.GetZ()
    hnx, hny, hnz = hn.GetX(), hn.GetY(), hn.GetZ()
    if abs(lnz) < 1e-12 or abs(hnz) < 1e-12:
        return None
    phi0, dphi = _phi_arg(geo)
    full = dphi is None
    span = 2.0 * math.pi if full else dphi
    n = max(8, int(round(_PV_REV_RES * span / (2.0 * math.pi))))
    ang = np.linspace(phi0, phi0 + span, n + 1)

    def zlo(r, a):                  # cut plane through (0,0,-dz), normal LowNorm
        return -dz - r * (lnx * math.cos(a) + lny * math.sin(a)) / lnz

    def zhi(r, a):                  # cut plane through (0,0,+dz), normal HighNorm
        return dz - r * (hnx * math.cos(a) + hny * math.sin(a)) / hnz

    pts, faces, index = [], [], {}

    def P(r, j, low):
        if r < 1e-12:
            key = ("axis", low)
        else:
            jj = 0 if (full and j == n) else j
            key = (r, jj, low)
        k = index.get(key)
        if k is None:
            a = ang[0 if (full and j == n) else j]
            z = zlo(r, a) if low else zhi(r, a)
            k = len(pts)
            pts.append((r * math.cos(a), r * math.sin(a), z))
            index[key] = k
        return k

    for j in range(n):
        ob0, ob1 = P(rmax, j, True), P(rmax, j + 1, True)
        ot0, ot1 = P(rmax, j, False), P(rmax, j + 1, False)
        ib0, ib1 = P(rmin, j, True), P(rmin, j + 1, True)
        it0, it1 = P(rmin, j, False), P(rmin, j + 1, False)
        faces.append([ob0, ob1, ot1, ot0])              # outer barrel, outward
        if rmin > 1e-12:
            faces.append([it0, it1, ib1, ib0])          # inner barrel, inward
        for f in ([ob0, ib0, ib1, ob1], [ot0, ot1, it1, it0]):  # low/high caps
            f = _dedup_face(f)
            if len(f) >= 3:
                faces.append(f)
    if not full:
        for j, flip in ((0, False), (n, True)):
            f = [P(rmax, j, True), P(rmax, j, False), P(rmin, j, False),
                 P(rmin, j, True)]                      # CCW in (r_hat, z)
            f = _dedup_face(f[::-1] if flip else f)
            if len(f) >= 3:
                faces.append(f)
    return _pv_assemble(pts, faces)


def _pv_ellipsoid(geo):
    import numpy as np
    ax, by, cz = geo.Ax, geo.By, geo.Cz
    z1, z2 = geo.Zcut1, geo.Zcut2
    th_lo = math.acos(min(max(z2 / cz, -1.0), 1.0))     # top cut -> small theta
    th_hi = math.acos(min(max(z1 / cz, -1.0), 1.0))
    if th_hi - th_lo < 1e-12:
        return None
    nt = max(8, int(round(_PV_REV_RES * (th_hi - th_lo) / math.pi)))
    n = _PV_REV_RES
    thetas = np.linspace(th_lo, th_hi, nt + 1)
    ang = np.linspace(0.0, 2.0 * math.pi, n + 1)
    pts, faces, index = [], [], {}

    def P(i, j):
        j = j % n
        key = (i, j)
        k = index.get(key)
        if k is None:
            st, ct = math.sin(thetas[i]), math.cos(thetas[i])
            k = len(pts)
            pts.append((ax * st * math.cos(ang[j]), by * st * math.sin(ang[j]),
                        cz * ct))
            index[key] = k
        return k

    for i in range(nt):
        for j in range(n):
            f = _dedup_face([P(i, j), P(i + 1, j), P(i + 1, j + 1), P(i, j + 1)])
            if len(f) >= 3:
                faces.append(f)
    if z2 < cz - 1e-12:                                 # flat top cap at z2
        c = len(pts)
        pts.append((0.0, 0.0, cz * math.cos(th_lo)))
        for j in range(n):
            faces.append([c, P(0, j), P(0, j + 1)])
    if z1 > -cz + 1e-12:                                # flat bottom cap at z1
        c = len(pts)
        pts.append((0.0, 0.0, cz * math.cos(th_hi)))
        for j in range(n):
            faces.append([c, P(nt, j + 1), P(nt, j)])
    return _pv_assemble(pts, faces)


def _pv_elliptical_tube(geo):
    import numpy as np
    dx, dy, dz = geo.Dx, geo.Dy, geo.Dz                 # semi-axes + HALF height
    n = _PV_REV_RES
    ang = np.linspace(0.0, 2.0 * math.pi, n + 1)
    pts = []
    for z in (-dz, dz):
        for j in range(n):
            pts.append((dx * math.cos(ang[j]), dy * math.sin(ang[j]), z))
    lo, hi = 0, n
    faces = []
    for j in range(n):
        j1 = (j + 1) % n
        faces.append([lo + j, lo + j1, hi + j1, hi + j])
    cb, ct = len(pts), len(pts) + 1
    pts += [(0.0, 0.0, -dz), (0.0, 0.0, dz)]
    for j in range(n):
        j1 = (j + 1) % n
        faces.append([cb, lo + j1, lo + j])
        faces.append([ct, hi + j, hi + j1])
    return _pv_assemble(pts, faces)


def _pv_extrpoly(geo):
    poly = [(float(p[0]), float(p[1])) for p in geo.polygon]
    zsec = [(float(s.zpos), float(s.scale), (float(s.offset[0]), float(s.offset[1])))
            for s in geo.zsections]
    if len(poly) < 3 or len(zsec) < 2:
        return None
    area = 0.0
    for i in range(len(poly)):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % len(poly)]
        area += x0 * y1 - x1 * y0
    if area < 0:
        poly.reverse()                                  # normalise to CCW
    m = len(poly)
    pts, faces = [], []
    for z, s, (ox, oy) in zsec:
        for x, y in poly:
            pts.append((x * s + ox, y * s + oy, z))
    for k in range(len(zsec) - 1):
        a, b = k * m, (k + 1) * m
        for i in range(m):
            i1 = (i + 1) % m
            faces.append([a + i, a + i1, b + i1, b + i])
    bottom = list(range(m))
    top = list(range((len(zsec) - 1) * m, len(zsec) * m))
    faces.append(bottom[::-1])                          # outward -z
    faces.append(top)                                   # outward +z
    return _pv_assemble(pts, faces)


def _pv_para(geo):
    dx, dy, dz = geo.Dx, geo.Dy, geo.Dz                 # HALF lengths
    ta, tt = math.tan(geo.Alpha), math.tan(geo.Theta)
    cp, sp = math.cos(geo.Phi), math.sin(geo.Phi)
    e1 = (dx, 0.0, 0.0)
    e2 = (dy * ta, dy, 0.0)
    e3 = (dz * tt * cp, dz * tt * sp, dz)               # matches Para::ComputeFaceData
    v = []
    for s3 in (-1, 1):
        for s1, s2 in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            v.append(tuple(s1 * e1[c] + s2 * e2[c] + s3 * e3[c] for c in range(3)))
    return _pv_hexa(v)


def _pv_trd(geo):
    x1, x2, y1, y2, dz = geo.Dx1, geo.Dx2, geo.Dy1, geo.Dy2, geo.Dz  # HALF lengths
    v = [(-x1, -y1, -dz), (x1, -y1, -dz), (x1, y1, -dz), (-x1, y1, -dz),
         (-x2, -y2, dz), (x2, -y2, dz), (x2, y2, dz), (-x2, y2, dz)]
    return _pv_hexa(v)


def _pv_trap(geo):
    import numpy as np
    dz, dy1, dy2 = geo.Dz, geo.Dy1, geo.Dy2             # HALF lengths
    x1, x2, x3, x4 = geo.Dx1, geo.Dx2, geo.Dx3, geo.Dx4
    tt = math.tan(geo.Theta)
    cp, sp = math.cos(geo.Phi), math.sin(geo.Phi)
    ta1, ta2 = math.tan(geo.Alpha1), math.tan(geo.Alpha2)
    v = []
    for s3, dy, xa, xb, ta in ((-1, dy1, x1, x2, ta1), (1, dy2, x3, x4, ta2)):
        xc, yc = s3 * dz * tt * cp, s3 * dz * tt * sp   # matches Trap::ComputeVertices
        v += [(xc - dy * ta - xa, yc - dy, s3 * dz),
              (xc - dy * ta + xa, yc - dy, s3 * dz),
              (xc + dy * ta + xb, yc + dy, s3 * dz),
              (xc + dy * ta - xb, yc + dy, s3 * dz)]
    v = np.asarray(v, float)
    # With alpha1 != alpha2 (or mismatched dx ratios) the side quads through
    # these corners are NOT planar; SIREN's solid is the intersection of six
    # half-spaces, each through the same 3-vertex triple Trap::ComputePlane
    # uses. Rebuild the corners as the triple intersections of those planes
    # (z-face, y-face, x-face) -- exact for the planar case too.
    triples = [(0, 1, 2), (7, 6, 5), (0, 1, 5), (3, 2, 6), (0, 3, 7), (1, 2, 6)]
    centroid = v.mean(axis=0)
    planes = []
    for i, j, k in triples:
        nrm = np.cross(v[j] - v[i], v[k] - v[i])
        nrm /= np.linalg.norm(nrm)
        d = -nrm.dot(v[i])
        if nrm.dot(centroid) + d > 0:                   # outward (centroid inside)
            nrm, d = -nrm, -d
        planes.append((nrm, d))
    corners = []
    for iz, iy, ix in ((0, 2, 4), (0, 2, 5), (0, 3, 5), (0, 3, 4),
                       (1, 2, 4), (1, 2, 5), (1, 3, 5), (1, 3, 4)):
        A = np.array([planes[iz][0], planes[iy][0], planes[ix][0]])
        rhs = -np.array([planes[iz][1], planes[iy][1], planes[ix][1]])
        corners.append(np.linalg.solve(A, rhs))
    return _pv_hexa(corners)


def _pv_distance_fn(surface):
    """Batched signed-distance callable to a closed surface (< 0 inside).

    vtkImplicitPolyDataDistance signs by the angle-weighted pseudonormal of
    the nearest cell, so the surface winding must be outward -- which every
    builder here guarantees. The locator is built once and reused, and the
    whole batch is evaluated in C++ (FunctionValue) -- a python-level
    EvaluateFunction loop is 10-50x slower and dominates big booleans.
    """
    import numpy as np, vtk
    from pyvista.core.utilities.arrays import convert_array
    f = vtk.vtkImplicitPolyDataDistance()
    f.SetInput(surface)

    def dist(pts):
        pts = np.ascontiguousarray(pts, float)
        if not len(pts):
            return np.empty(0)
        out = vtk.vtkDoubleArray()
        f.FunctionValue(convert_array(pts), out)
        return convert_array(out)
    return dist


def _pv_intersection_curve(a, b, step):
    """(N, 3) points densely sampling the a/b surface intersection curve.

    Uses vtkIntersectionPolyDataFilter for the curve ONLY -- its triangle-
    triangle intersection tests are solid, while its mesh-splitting outputs
    can be corrupt (duplicated/unsplit faces), so those are never used.
    Polyline segments are resampled to ~*step* spacing so a point-to-curve
    distance is a faithful curve distance. Returns None when there is no
    transversal intersection (disjoint, contained, or only coincident faces).
    """
    import numpy as np, pyvista as pv, vtk
    warn = vtk.vtkObject.GetGlobalWarningDisplay()
    vtk.vtkObject.GlobalWarningDisplayOff()
    try:
        f = vtk.vtkIntersectionPolyDataFilter()
        f.SetInputData(0, a)
        f.SetInputData(1, b)
        f.SplitFirstOutputOff()
        f.SplitSecondOutputOff()
        f.Update()
        curve = pv.wrap(f.GetOutput(0))
    except Exception:
        return None
    finally:
        vtk.vtkObject.SetGlobalWarningDisplay(warn)
    pts = np.asarray(curve.points, float)
    if not len(pts):
        return None
    out = [pts]
    lines = np.asarray(curve.lines)
    i = 0
    while i < len(lines):                       # walk [n, i0, ..., in-1] cells
        n = lines[i]
        seg = pts[lines[i + 1:i + 1 + n]]
        i += n + 1
        for p0, p1 in zip(seg[:-1], seg[1:]):
            length = np.linalg.norm(p1 - p0)
            k = int(length / step)
            if k:
                t = np.linspace(0.0, 1.0, k + 2)[1:-1, None]
                out.append(p0 + t * (p1 - p0))
    pts = np.concatenate(out, axis=0)
    return pts[:200000]                          # hard cap, plenty for refinement


def _tri_subdivide(T):
    """Midpoint 4-split of an (n, 3, 3) triangle array (winding preserved)."""
    import numpy as np
    m01 = 0.5 * (T[:, 0] + T[:, 1])
    m12 = 0.5 * (T[:, 1] + T[:, 2])
    m20 = 0.5 * (T[:, 2] + T[:, 0])
    return np.concatenate([np.stack([T[:, 0], m01, m20], 1),
                           np.stack([m01, T[:, 1], m12], 1),
                           np.stack([m20, m12, T[:, 2]], 1),
                           np.stack([m01, m12, m20], 1)], axis=0)


def _tri_cut(t, d):
    """Cut one triangle along the linearly interpolated zero of vertex values
    *d* (mixed signs). Returns 3 sub-triangles with the original winding;
    zero-valued vertices yield degenerate slivers that are filtered later."""
    import numpy as np
    pos = d > 0
    i = int(np.flatnonzero(pos)[0]) if pos.sum() == 1 else int(np.flatnonzero(~pos)[0])
    j, k = (i + 1) % 3, (i + 2) % 3
    pij = t[i] + (d[i] / (d[i] - d[j])) * (t[j] - t[i])
    pik = t[i] + (d[i] / (d[i] - d[k])) * (t[k] - t[i])
    return [np.stack([t[i], pij, pik]), np.stack([pij, t[j], t[k]]),
            np.stack([pij, t[k], pik])]


def _pv_csg_split(mesh, dist, curve, eps_len, max_depth=6):
    """Triangle soup of *mesh*, refined and cut so every output triangle lies
    on a single side of the surface behind *dist* (within ~*eps_len*).

    Triangles near the intersection *curve* (or whose vertex distances change
    sign -- the curve can pierce a big face without crossing its edges, so the
    sign test alone is insufficient, and vice versa) are midpoint-subdivided
    until smaller than *eps_len*; remaining sign-straddlers are cut along the
    interpolated zero of *dist*. The midpoint splits leave zero-width
    T-junctions against unrefined neighbours -- invisible and harmless for
    rendering, picking, and reuse as a nested CSG operand. Returns (n, 3, 3).
    """
    import numpy as np
    T = mesh.points[mesh.faces.reshape(-1, 4)[:, 1:]]
    done = []
    for depth in range(max_depth + 1):
        if not len(T):
            break
        d = dist(T.reshape(-1, 3)).reshape(-1, 3)
        mixed = (d.max(axis=1) > 0) & (d.min(axis=1) < 0)
        diam = np.linalg.norm(T - np.roll(T, 1, axis=1), axis=2).max(axis=1)
        near = mixed
        if curve is not None:
            cen = T.mean(axis=1)
            dc = np.empty(len(cen))
            for i in range(0, len(cen), 256):   # blockwise point-to-curve dist
                blk = cen[i:i + 256]
                dc[i:i + 256] = np.sqrt(
                    ((blk[:, None, :] - curve[None]) ** 2).sum(2)).min(1)
            near = near | (dc < diam)
        refine = near & (diam > eps_len) & (depth < max_depth)
        final = ~refine
        done.append(T[final & ~mixed])
        for t, dv in zip(T[final & mixed], d[final & mixed]):
            done.extend(_tri_cut(t, dv))
        if not refine.any():
            break
        T = _tri_subdivide(T[refine])
    done = [np.asarray(x, float).reshape(-1, 3, 3) for x in done if len(x)]
    if not done:
        return np.empty((0, 3, 3))
    return np.concatenate(done, axis=0)


def _pv_csg(op, a, b):
    """Combine two closed operand meshes under a BooleanOperation.

    Deliberately does NOT use VTK's exact mesh booleans:
    vtkBooleanOperationPolyDataFilter dies with a native bus error on some
    SBND building booleans (a signal -- uncatchable from python), takes
    minutes on others, and rarely produces a watertight result around the
    coincident faces ubiquitous in GDML (flush cut-outs, touching union
    halves); vtkLoopBooleanPolyDataFilter is no better. Instead the surfaces
    are trimmed: each operand is refined and cut along the other's surface
    (:func:`_pv_csg_split`), and whole cells are then kept or dropped by
    probing the other operand a small step along the cell's own outward
    normal:

      union:        keep iff just-outside stays vacuum  (outside the other)
      intersection: keep iff just-inside gains material (inside the other)
      difference:   keep A-cells iff just-inside keeps material (outside B);
                    keep B-cells iff just-outside-B is carved out of A.

    The probe direction makes coincident faces exact (no hidden internal
    walls, no cracks): the material configuration across the cell decides,
    not the sign of a ~0 centroid distance. Same-orientation coincident pairs
    would survive twice, so B-cells sitting ON dA are dropped for the
    symmetric ops; for a difference such cells are exactly the flipped notch
    boundary and must stay. The B piece of a difference is flipped so normals
    stay outward -- nested booleans reuse the result as an operand, where
    orientation matters.
    """
    import numpy as np, pyvista as pv
    name = str(op)
    kind = ("union" if "UNION" in name else
            "difference" if "SUBTRACTION" in name else "intersection")
    a = a.triangulate()
    b = b.triangulate()
    # Replicated placements (the 48 BNB target fins, repeated building
    # modules) re-pose the SAME boolean: operand placements are internal, the
    # differing placement is the boolean's own and is applied AFTER this call,
    # so the operand meshes are bit-identical -- memoise on their content.
    key = (kind, a.n_points, b.n_points,
           hash(a.points.tobytes()), hash(b.points.tobytes()))
    if key in _PV_CSG_CACHE:
        hit = _PV_CSG_CACHE[key]
        return hit.copy() if hit is not None else None
    diag = max(np.linalg.norm(np.asarray(a.bounds[1::2]) - np.asarray(a.bounds[::2])),
               np.linalg.norm(np.asarray(b.bounds[1::2]) - np.asarray(b.bounds[::2])))
    diag = max(diag, 1e-6)
    delta = 1e-6 * diag
    eps_len = 4e-3 * diag
    dist_a = _pv_distance_fn(a)
    dist_b = _pv_distance_fn(b)
    curve = _pv_intersection_curve(a, b, eps_len)
    TA = _pv_csg_split(a, dist_b, curve, eps_len)
    TB = _pv_csg_split(b, dist_a, curve, eps_len)

    def classify(T, dist_other, inward, want_inside):
        if not len(T):
            return T, np.zeros(0, bool)
        nrm = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
        area2 = np.linalg.norm(nrm, axis=1)
        ok = area2 > (1e-9 * diag) ** 2
        nrm[~ok] = 0.0
        nrm[ok] /= area2[ok, None]
        probe = T.mean(axis=1) + (-delta if inward else delta) * nrm
        dp = dist_other(probe)
        return T, ok & ((dp < 0) if want_inside else (dp > 0))

    if kind == "union":
        TA, ka = classify(TA, dist_b, inward=False, want_inside=False)
        TB, kb = classify(TB, dist_a, inward=False, want_inside=False)
        kb &= np.abs(dist_a(TB.mean(axis=1))) > delta if len(TB) else kb
    elif kind == "intersection":
        TA, ka = classify(TA, dist_b, inward=True, want_inside=True)
        TB, kb = classify(TB, dist_a, inward=True, want_inside=True)
        kb &= np.abs(dist_a(TB.mean(axis=1))) > delta if len(TB) else kb
    else:                                               # difference: A - B
        TA, ka = classify(TA, dist_b, inward=True, want_inside=False)
        TB, kb = classify(TB, dist_a, inward=False, want_inside=True)
        TB = TB[:, ::-1, :]                             # flip kept B outward

    def soup(T):
        n = len(T)
        fc = np.empty((n, 4), np.int64)
        fc[:, 0] = 3
        fc[:, 1:] = np.arange(3 * n).reshape(n, 3)
        return pv.PolyData(T.reshape(-1, 3), fc.ravel())

    kept = np.concatenate([TA[ka], TB[kb]], axis=0)
    out = soup(kept).clean() if len(kept) else None
    if out is not None and out.n_cells == 0:
        out = None
    _PV_CSG_CACHE[key] = out.copy() if out is not None else None
    return out


def _pv_boolean(geo):
    # Children's placements live in the boolean's local frame (the recursive
    # _pv_parametric applies them), so combine first, then place the result.
    a = _pv_parametric(geo.Left)
    b = _pv_parametric(geo.Right)
    if a is None or b is None:
        return None
    return _pv_csg(geo.Operation, a, b)


_PV_BUILDERS = {
    "Box": _pv_box,
    "Sphere": _pv_sphere,
    "Cylinder": _pv_cylinder,
    "Cone": _pv_cone,
    "Polycone": _pv_polycone,
    "GenericPolycone": _pv_generic_polycone,
    "Polyhedra": _pv_polyhedra,
    "ExtrPoly": _pv_extrpoly,
    "Para": _pv_para,
    "Trd": _pv_trd,
    "Trap": _pv_trap,
    "Ellipsoid": _pv_ellipsoid,
    "EllipticalTube": _pv_elliptical_tube,
    "CutTube": _pv_cut_tube,
    "Torus": _pv_torus,
    "BooleanGeometry": _pv_boolean,
}


def _pv_parametric(geo):
    """Crisp parametric pyvista mesh for any SIREN solid, or None (-> MC).

    Every solid is built at the origin in its local frame and then moved by
    its placement; a BooleanGeometry recurses into its operands (whose
    placements live in the boolean's local frame) and combines the surfaces.
    TriangularMesh triangles are already in the enclosing frame.
    """
    import numpy as np, pyvista as pv
    cls = type(geo).__name__
    if cls == "TriangularMesh":
        tris = np.asarray(geo.GetTriangles(), float)
        if tris.size == 0:
            return None
        tris = _apply_placement(tris, geo)                          # local -> world frame
        n = len(tris)
        f = np.empty((n, 4), np.int64); f[:, 0] = 3; f[:, 1:] = np.arange(3 * n).reshape(n, 3)
        return pv.PolyData(tris.reshape(-1, 3), f.ravel())
    builder = _PV_BUILDERS.get(cls)
    if builder is None:
        return None
    m = builder(geo)
    if m is None or m.n_cells == 0:
        return None
    pl = geo.placement; p = pl.Position
    T = _euler_transform(pl.Quaternion.GetEulerAnglesXYZs(),
                         (p.GetX(), p.GetY(), p.GetZ()))
    return m.transform(T, inplace=False)


# geometry class name -> number of marching-cubes fallbacks (diagnostic; the
# goal is for every solid to mesh parametrically -- view_pv reports the diff).
_MC_FALLBACKS = {}
# (kind, operand content hash) -> CSG result mesh (or None). Replicated
# placements of the same boolean (target fins, repeated building modules) are
# bit-identical at this level, so this collapses them to one computation.
_PV_CSG_CACHE = {}


def _pv_marching(geo, res=48, pad=0.04):
    """Surface of ANY geometry by marching cubes on IsInside over its world AABB.

    Last-resort fallback (blocky and slow); every call is counted in
    _MC_FALLBACKS so :func:`view_pv` can report which solids lack a
    parametric mesh.
    """
    import numpy as np, pyvista as pv
    from siren.math import Vector3D
    cls = type(geo).__name__
    _MC_FALLBACKS[cls] = _MC_FALLBACKS.get(cls, 0) + 1
    print("[siren.view_pv] marching-cubes fallback for %s" % cls)
    bb = _bbox(geo)
    if bb is None:
        return None
    (cx, cy, cz), (hx, hy, hz) = bb
    hx, hy, hz = hx * (1 + pad) + 1.0, hy * (1 + pad) + 1.0, hz * (1 + pad) + 1.0
    xs = np.linspace(cx - hx, cx + hx, res)
    ys = np.linspace(cy - hy, cy + hy, res)
    zs = np.linspace(cz - hz, cz + hz, res)
    inside = geo.IsInside
    vol = np.empty((res, res, res), np.float32)
    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                vol[i, j, k] = inside(Vector3D(float(x), float(y), float(z)))
    if vol.max() < 0.5:
        return None
    grid = pv.ImageData(dimensions=(res, res, res),
                        spacing=(xs[1] - xs[0], ys[1] - ys[0], zs[1] - zs[0]),
                        origin=(xs[0], ys[0], zs[0]))
    grid["v"] = vol.ravel(order="F")
    return grid.contour([0.5], scalars="v")


def _pv_sector_mesh(geo, mc_res=48):
    try:
        m = _pv_parametric(geo)
    except Exception as e:                  # never let a bad build drop a sector
        print("[siren.view_pv] parametric build failed for %s (%s)"
              % (type(geo).__name__, e))
        m = None
    return m if m is not None else _pv_marching(geo, mc_res)


def _pv_material_colours(model):
    """material name -> (rgb, alpha) on the same log-density scale as the legend."""
    import colorsys
    dens = {}
    for s in _sectors(model):
        nm = _material_name(model, int(s.material_id))
        dens[nm] = max(dens.get(nm, 0.0), max(_sector_density(s), 1e-30))
    if not dens:
        return {}
    hi, lo = math.log(max(dens.values())), math.log(1e-4)
    out = {}
    for nm, d in dens.items():
        t = 0.5 if hi <= lo else min(max((math.log(d) - lo) / (hi - lo), 0.0), 1.0)
        out[nm] = (colorsys.hsv_to_rgb(0.66 * (1.0 - t), 0.85, 0.95),
                   0.15 if d < _LOW_DENSITY else 0.92)
    return out


def view_pv(model, screenshot=None, region=None, region_center=None, max_extent=None,
            axes=True, legend=True, picker=True, mc_res=48, z_exag=1.0, dark=True,
            window_size=(1500, 950)):
    """Render *model* entirely in pyvista -- no pyg4ometry.

    Every solid type is meshed parametrically (crisp and fast), including
    hollow/phi-cut solids and booleans (operand surfaces combined by CSG);
    marching cubes on SIREN's IsInside (*mc_res* grid) remains only as a
    fallback for unknown types or failed builds, and a summary line reports
    how many sectors needed it. Coloured by material (log-density); gases
    semi-transparent. Right-click (press 'p' over a surface) identifies the
    volume via SIREN. *region*/*region_center* (half-width, centre) or
    *max_extent* drop huge background volumes (PREM shells, world boxes).
    *z_exag* applies a vertical exaggeration (terrain).
    """
    import pyvista as pv
    from siren.math import Vector3D
    from siren.detector import DetectorPosition
    mc0 = dict(_MC_FALLBACKS)
    cols = _pv_material_colours(model)
    crop = None
    if region is not None:
        if region_center is None:
            try:
                o = model.DetectorOrigin.get(); region_center = (o.GetX(), o.GetY(), o.GetZ())
            except Exception:
                region_center = (0.0, 0.0, 0.0)
        crop = (region_center, float(region))
    pl = pv.Plotter(off_screen=bool(screenshot), window_size=window_size)
    pl.set_background("#1a1d24" if dark else "white")
    drawn = 0
    for s in _sectors(model):
        bb = _bbox(s.geo)
        if bb is None:
            continue
        (cx, cy, cz), (hx, hy, hz) = bb
        if max_extent is not None and max(hx, hy, hz) > max_extent:
            continue
        if crop is not None:
            (rx, ry, rz), rh = crop
            if abs(cx - rx) > rh or abs(cy - ry) > rh or abs(cz - rz) > rh:
                continue
        try:
            mesh = _pv_sector_mesh(s.geo, mc_res)
        except Exception:
            mesh = None
        if mesh is None or mesh.n_cells == 0:
            continue
        nm = _material_name(model, int(s.material_id))
        rgb, alpha = cols.get(nm, ((0.6, 0.6, 0.6), 0.9))
        pl.add_mesh(mesh, color=rgb, opacity=alpha, smooth_shading=True)
        drawn += 1
    mc = {k: v - mc0.get(k, 0) for k, v in _MC_FALLBACKS.items() if v - mc0.get(k, 0)}
    n_mc = sum(mc.values())
    print("[siren.view_pv] %d sectors drawn: %d parametric, %d marching-cubes%s"
          % (drawn, drawn - n_mc, n_mc,
             " (%s)" % ", ".join("%s:%d" % kv for kv in sorted(mc.items())) if mc else ""))
    if z_exag != 1.0:
        pl.set_scale(zscale=z_exag)
    if axes:
        pl.add_axes()
    if legend and cols:
        pl.add_legend([[nm, c[0]] for nm, c in sorted(cols.items())],
                      bcolor="k" if dark else "w", size=(0.20, 0.28))
    if picker and not isinstance(model, str):
        fg = "w" if dark else "k"
        pl.add_text("", position="upper_left", font_size=10, color=fg, name="pick")
        # use_picker=True makes pyvista invoke callback(point, picker); the
        # point is already divided by plotter.scale (i.e. z_exag is undone,
        # see pyvista picking.py "HACK: handle scale"), so use it as-is.
        def _cb(pt, *_):
            try:
                sec = model.GetContainingSector(DetectorPosition(
                    Vector3D(float(pt[0]), float(pt[1]), float(pt[2]))))
                pl.add_text("%s [%s]  rho=%.3g g/cc" % (
                    sec.name, _material_name(model, int(sec.material_id)),
                    _sector_density(sec)), position="upper_left",
                    font_size=10, color=fg, name="pick")
            except Exception:
                pass
        try:
            pl.enable_point_picking(callback=_cb, show_point=True, use_picker=True)
        except Exception:
            pass
    if screenshot:
        pl.show(screenshot=screenshot)
        return screenshot
    pl.show()
    return drawn


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
