"""
Homestake/SURF overburden + PREM Earth model for the DUNE far detector.

Adds real topography (graded, anti-aliased DEM mesh), Homestake geological
formations (Poorman/Yates/rhyolite), and a full PREM Earth to a DetectorModel
that was loaded from the DUNE FD GDML. The GDML volumes always take priority
over these sectors wherever both exist (same pattern as SBN earth_model.py).

Site frame (shared with detector.py): +y = up, +z = LBNF beam downstream
(geographic azimuth ~277.15 deg), +x = right-handed. The overburden is built in
a local ENU frame (z = up) then rotated by the ENU->site placement so ENU-up
maps to +y and the ENU horizontal plane carries +z along the beam azimuth.

Call ``add_earth_model(model)`` after ``model.LoadGDML(path)``.
"""
from __future__ import annotations

import math
import os
import importlib.util
import sys
import numpy as np

_THIS_DIR = os.path.dirname(os.path.realpath(__file__))

# SURF/Homestake 4850L site
SITE_LAT = 44.352
SITE_LON = -103.751
DEPTH = 1478.0          # vertical depth (m) to 4850L
R_PREM = 6371000.0


def _load_sibling(name, filename):
    fqn = f"siren._dunefd.{name}"
    if fqn in sys.modules:
        return sys.modules[fqn]
    spec = importlib.util.spec_from_file_location(
        fqn, os.path.join(_THIS_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[fqn] = mod
    try:
        spec.loader.exec_module(mod)
    except Exception:
        del sys.modules[fqn]
        raise
    return mod


# ---- DEM access (Copernicus GLO-30) ----
_DEM_DIR = os.path.join(_THIS_DIR, "dem")
_TILES = {}
_SAT = {}

_phi = math.radians(SITE_LAT)
M_LAT = 111132.954 - 559.822 * math.cos(2 * _phi) + 1.175 * math.cos(4 * _phi)
M_LON = (math.pi / 180.0) * 6378137.0 * math.cos(_phi) / math.sqrt(
    1 - 0.00669437999014 * math.sin(_phi) ** 2)


def _ensure_dem():
    """Download Copernicus GLO-30 tiles if not cached."""
    if _TILES:
        return
    os.makedirs(_DEM_DIR, exist_ok=True)
    import tifffile
    for ul_lon, fn in [(-104.0, "cop30_N44_W104.tif"), (-105.0, "cop30_N44_W105.tif")]:
        path = os.path.join(_DEM_DIR, fn)
        if not os.path.isfile(path):
            import urllib.request
            url = (f"https://copernicus-dem-30m.s3.amazonaws.com/"
                   f"Copernicus_DSM_COG_10_N44_00_W{abs(int(ul_lon)):03d}_00_DEM/"
                   f"Copernicus_DSM_COG_10_N44_00_W{abs(int(ul_lon)):03d}_00_DEM.tif")
            urllib.request.urlretrieve(url, path)
        arr = tifffile.imread(path).astype(np.float64)
        _TILES[ul_lon] = (arr, 45.0, 1.0 / 3600.0)
        sat = np.pad(arr.cumsum(0).cumsum(1), ((1, 0), (1, 0)))
        _SAT[ul_lon] = (sat, 45.0, 1.0 / 3600.0)


def _sample_dem(lat, lon):
    _ensure_dem()
    lat = np.atleast_1d(np.asarray(lat, float))
    lon = np.atleast_1d(np.asarray(lon, float))
    out = np.full(lat.shape, np.nan)
    for ul_lon, (arr, ul_lat, dx) in _TILES.items():
        m = (lon >= ul_lon) & (lon < ul_lon + 1.0)
        if not m.any():
            continue
        nr, nc = arr.shape
        col = (lon[m] - ul_lon) / dx; row = (ul_lat - lat[m]) / dx
        c0 = np.clip(np.floor(col).astype(int), 0, nc - 2)
        r0 = np.clip(np.floor(row).astype(int), 0, nr - 2)
        fc = np.clip(col - c0, 0, 1); fr = np.clip(row - r0, 0, 1)
        out[m] = (arr[r0, c0] * (1-fc) * (1-fr) + arr[r0, c0+1] * fc * (1-fr)
                  + arr[r0+1, c0] * (1-fc) * fr + arr[r0+1, c0+1] * fc * fr)
    return out


Z_SURF_SITE = None
Z_DET_ASL = None

def _init_site():
    global Z_SURF_SITE, Z_DET_ASL
    if Z_SURF_SITE is not None:
        return
    _ensure_dem()
    Z_SURF_SITE = float(_sample_dem(SITE_LAT, SITE_LON)[0])
    Z_DET_ASL = Z_SURF_SITE - DEPTH


_M_PER_PIX_LAT = (1.0 / 3600.0) * M_LAT
_M_PER_PIX_LON = (1.0 / 3600.0) * M_LON


def _avg_elev_asl(x, y, window):
    _ensure_dem()
    x = np.atleast_1d(np.asarray(x, float)); y = np.atleast_1d(np.asarray(y, float))
    window = np.broadcast_to(np.atleast_1d(np.asarray(window, float)), x.shape)
    lat = SITE_LAT + y / M_LAT; lon = SITE_LON + x / M_LON
    out = np.full(x.shape, np.nan)
    for ul_lon, (sat, ul_lat, dx) in _SAT.items():
        m = (lon >= ul_lon) & (lon < ul_lon + 1.0)
        if not m.any():
            continue
        nr, nc = sat.shape[0] - 1, sat.shape[1] - 1
        col = (lon[m] - ul_lon) / dx; row = (ul_lat - lat[m]) / dx
        hr = np.maximum(window[m] / 2 / _M_PER_PIX_LAT, 0.5)
        hc = np.maximum(window[m] / 2 / _M_PER_PIX_LON, 0.5)
        r0 = np.clip(np.floor(row - hr).astype(int), 0, nr - 1)
        r1 = np.clip(np.ceil(row + hr).astype(int), 1, nr)
        c0 = np.clip(np.floor(col - hc).astype(int), 0, nc - 1)
        c1 = np.clip(np.ceil(col + hc).astype(int), 1, nc)
        S = sat[r1, c1] - sat[r0, c1] - sat[r1, c0] + sat[r0, c0]
        area = ((r1 - r0) * (c1 - c0)).astype(float)
        out[m] = S / area
    return out


# ---- graded mesh builder (ENU frame, z=up, origin at detector) ----
H_MIN = 80.0
L_HALF = 20480.0
BASE_Z = -3000.0
C_Z_ENU = -(R_PREM - DEPTH)  # Earth center in ENU
TAPER_START = 15000.0
_BANDS = [(2560, 80.0), (5120, 160.0), (10240, 320.0), (15360, 640.0), (1e9, 1280.0)]


def _target_cell(xc, yc):
    r = max(abs(xc), abs(yc))
    for rmax, c in _BANDS:
        if r < rmax:
            return c
    return _BANDS[-1][1]


def _surface_sphere_z(x, y):
    return C_Z_ENU + math.sqrt(max(R_PREM * R_PREM - x * x - y * y, 0.0))


def _elev(x, y, window):
    """Terrain surface elevation (ENU z above the detector), with curvature correction.

    The curvature drop d^2/(2R) accounts for the Earth curving away from the
    tangent plane. The surface is the real DEM terrain everywhere we have data;
    no blending to the PREM sphere (which would add fake rock where the terrain
    is lower than the sphere). The mesh BASE tapers to the PREM sphere instead
    (see _base_z), so leaving the mesh laterally falls through the bottom into
    PREM crust -- physically correct.
    """
    _init_site()
    z = float(_avg_elev_asl(x, y, window)[0]) - Z_DET_ASL
    z -= (x * x + y * y) / (2.0 * R_PREM)   # Earth-curvature drop
    return z


def _base_z(x, y):
    """Mesh base elevation (ENU z): flat at BASE_Z so the mesh is a filled rock
    volume from the terrain surface down to BASE_Z. (A graded spherical-cap base
    is the planned refinement.) The local_crust box continues the rock below."""
    return BASE_Z


_vcache = {}
U = H_MIN
ROOT = int(round(L_HALF / U))


def _v(ix, iy, V3):
    key = (ix, iy)
    if key not in _vcache:
        x, y = ix * U, iy * U
        _vcache[key] = V3(x, y, _elev(x, y, _target_cell(x, y)))
    return _vcache[key]


def _build_graded_mesh(V3):
    _vcache.clear()
    _init_site()
    leaves = set()

    def subdivide(x0, y0, s):
        cx, cy = (x0 + s / 2) * U, (y0 + s / 2) * U
        if s > 1 and s * U > _target_cell(cx, cy):
            h = s // 2
            for dx in (0, h):
                for dy in (0, h):
                    subdivide(x0 + dx, y0 + dy, h)
        else:
            leaves.add((x0, y0, s))

    subdivide(-ROOT, -ROOT, 2 * ROOT)
    # 2:1 balance
    changed = True
    while changed:
        changed = False
        for cell in list(leaves):
            x0, y0, s = cell
            if s == 1:
                continue
            pts = []
            for t in (0.25, 0.75):
                pts += [(x0 + s + 0.01, y0 + s * t), (x0 - 0.01, y0 + s * t),
                        (x0 + s * t, y0 + s + 0.01), (x0 + s * t, y0 - 0.01)]
            def _lsz(px, py):
                sz = 1
                while sz <= 2 * ROOT:
                    lx = int(np.floor(px / sz)) * sz
                    ly = int(np.floor(py / sz)) * sz
                    if (lx, ly, sz) in leaves:
                        return sz
                    sz *= 2
                return None
            if any((sz := _lsz(px, py)) is not None and sz < s / 2 for px, py in pts):
                leaves.discard(cell)
                h = s // 2
                for dx in (0, h):
                    for dy in (0, h):
                        leaves.add((x0 + dx, y0 + dy, h))
                changed = True

    def _has_mid(x0, y0, s, side):
        if side == "R":   px, py = x0 + s + 0.01, y0 + s * 0.5
        elif side == "L": px, py = x0 - 0.01,     y0 + s * 0.5
        elif side == "T": px, py = x0 + s * 0.5,  y0 + s + 0.01
        else:             px, py = x0 + s * 0.5,  y0 - 0.01
        sz = 1
        while sz <= 2 * ROOT:
            lx = int(np.floor(px / sz)) * sz
            ly = int(np.floor(py / sz)) * sz
            if (lx, ly, sz) in leaves:
                return sz == s // 2
            sz *= 2
        return False

    tris = []
    for (x0, y0, s) in leaves:
        h = s // 2
        cx, cy = (x0 + h) * U, (y0 + h) * U
        cz = _elev(cx, cy, s * U)
        C = V3(cx, cy, cz)
        corners = {"B": [(x0, y0), (x0 + s, y0)],
                   "R": [(x0 + s, y0), (x0 + s, y0 + s)],
                   "T": [(x0 + s, y0 + s), (x0, y0 + s)],
                   "L": [(x0, y0 + s), (x0, y0)]}
        mids = {"B": (x0 + h, y0), "R": (x0 + s, y0 + h),
                "T": (x0 + h, y0 + s), "L": (x0, y0 + h)}
        for side in ("B", "R", "T", "L"):
            (a, b) = corners[side]
            if _has_mid(x0, y0, s, side):
                mx, my = mids[side]
                tris.append([C, _v(*a, V3), _v(mx, my, V3)])
                tris.append([C, _v(mx, my, V3), _v(*b, V3)])
            else:
                tris.append([C, _v(*a, V3), _v(*b, V3)])

    step = int(round(_BANDS[-1][1] / U))
    n = (2 * ROOT) // step
    ring = ([(-ROOT + k * step, -ROOT) for k in range(n)]
            + [(ROOT, -ROOT + k * step) for k in range(n)]
            + [(ROOT - k * step, ROOT) for k in range(n)]
            + [(-ROOT, ROOT - k * step) for k in range(n)])
    def Bv(ix, iy):
        return V3(ix * U, iy * U, _base_z(ix * U, iy * U))
    for k in range(len(ring)):
        a = ring[k]; b = ring[(k + 1) % len(ring)]
        tris.append([_v(*a, V3), Bv(*a), Bv(*b)])
        tris.append([_v(*a, V3), Bv(*b), _v(*b, V3)])
    for k in range(1, len(ring) - 1):
        tris.append([Bv(*ring[0]), Bv(*ring[k + 1]), Bv(*ring[k])])
    return tris, leaves


# ---- PREM shells ----
_PREM = [
    ("innercore",    1221500, "INNERCORE", 13.0),
    ("outercore",    3480000, "OUTERCORE", 11.0),
    ("lower_mantle", 5701000, "MANTLE",     4.9),
    ("upper_mantle", 6336000, "MANTLE",     3.4),
    ("inner_crust",  6356000, "ROCK",       2.9),
    ("upper_crust",  6371000, "ROCK",       2.65),
]
# Atmosphere shells are generated at load time from a pluggable air-density model
# (default US Standard Atmosphere 1976; see atmosphere.py) as mass-conserving
# spherical shells, replacing the old 3 hand-set AIR shells. See add_earth_model.

# ---- schematic geology ----
CONTACT_Z = -150.0    # Poorman/Yates contact (ENU z rel. detector)
CONTACT_DIP = 55.0
CONTACT_DIPDIR = 45.0
DIKES = [(-400.0, 200.0, 315.0, 40.0, 30000.0),
         (1200.0, -800.0, 300.0, 25.0, 30000.0)]


# Beam downstream azimuth at the FD (from the FNAL->SURF chord). Must match
# detector.py BEAM_AZIMUTH_DEG: +z (geometry) points along this geographic bearing.
BEAM_AZIMUTH_DEG = 277.15


def add_earth_model(model, atmosphere=None, atmo_top_km=100.0):
    """Add the Homestake overburden + PREM Earth to a DUNE FD DetectorModel.

    Must be called after the detector GDML is loaded. The model is in the SITE
    FRAME: +y = up, +z = LBNF beam downstream (geographic azimuth ~277.15 deg),
    +x = right-handed. The overburden is built in a local ENU frame (z = up) and
    rotated into the site frame, so the real Black Hills topography is oriented
    correctly relative to the beam.

    The atmosphere is built as mass-conserving spherical shells from a pluggable
    air-density model (``atmosphere``; default ``atmosphere.USStandard1976``) out
    to ``atmo_top_km`` km ASL, and the near-surface ``local_air`` box is set from
    the same model. Pass a custom ``atmosphere.Atmosphere`` subclass (e.g. an
    NRLMSIS / GDAS / CORSIKA-Linsley profile) for atmospheric-neutrino production
    studies where a specific or seasonal air column matters.
    """
    from siren.math import Vector3D, Quaternion, Matrix3D
    from siren.geometry import Box, Sphere, Placement, TriangularMesh, BooleanGeometry, BooleanOperation
    from siren.detector import DetectorSector, ConstantDensityDistribution
    V3 = Vector3D

    # ---- load materials (supplements the GDML's own) ----
    geo_mod = _load_sibling("homestake_geology", "homestake_geology.py")
    mat_path = os.path.join(_THIS_DIR, "DUNEFD_homestake_materials.dat")
    if not os.path.isfile(mat_path):
        geo_mod.write_materials_dat(mat_path)
    model.LoadMaterialModel(mat_path)
    mats = model.GetMaterials()

    # ---- ENU -> site rotation -------------------------------------------
    # ENU: x=East, y=North, z=Up.  Site: x=right-handed, y=Up, z=beam azimuth.
    # R maps ENU vectors to site coords (rows = site axes expressed in ENU):
    #   site_z (beam) = (sin az, cos az, 0);  site_y (up) = (0,0,1);
    #   site_x = site_y x site_z = (-cos az, sin az, 0).
    az = math.radians(BEAM_AZIMUTH_DEG)
    ca, sa = math.cos(az), math.sin(az)
    R = (-ca, sa, 0.0,    0.0, 0.0, 1.0,    sa, ca, 0.0)
    q_enu_site = Quaternion()
    q_enu_site.SetMatrix(Matrix3D(*R))
    enu_placement = Placement(V3(0.0, 0.0, 0.0), q_enu_site)

    # Earth center is DEPTH below the surface; the detector (geometry origin,
    # y~0) is the 4850L, so the center sits at y = -(R_PREM - DEPTH) (straight
    # down). The PREM surface sphere (R_PREM) then passes through y = +DEPTH.
    earth_place = Placement(V3(0.0, -(R_PREM - DEPTH), 0.0))

    def sector(name, material, level, geo, rho):
        s = DetectorSector()
        s.name = name
        s.material_id = mats.GetMaterialId(material)
        s.level = level
        s.geo = geo
        s.density = ConstantDensityDistribution(float(rho))
        model.AddSector(s)

    # ---- level strategy -------------------------------------------------
    # The composite has the vacuum "World" box at level 0 and all detector
    # volumes above it (the DUSEL_Rock world boxes are unwrapped away by
    # as_assembly). We want:
    #   PREM (<0) < World(0) < local_crust/air < overburden < detector volumes.
    # Shift every detector volume up to open a gap for the overburden.
    sectors = model.Sectors
    dup = {"local_crust", "local_air", "Poorman_bulk", "Yates_lens", "dike_0", "dike_1"}
    dup |= set(n for n, *_ in _PREM)
    if any(s.name in dup or s.name.startswith("atmo_") for s in sectors):
        raise RuntimeError("add_earth_model already applied to this model")
    SHIFT = 100000
    for s in sectors:
        if s.name not in ("UNIVERSE", "World"):
            s.level += SHIFT
    model.Sectors = sorted(sectors, key=lambda s: s.level)

    # ---- atmosphere shells from the air-density model (mass-conserving) ---
    # Built as spheres above R_PREM; R_PREM <-> the site surface (Z_SURF_SITE
    # ASL), so a shell whose outer edge is at altitude h (ASL) has radius
    # R_PREM + (h - Z_SURF_SITE). Each density is the mass-conserving (column-
    # preserving) mean over its band -- the value that matters for slant depth.
    atm_mod = _load_sibling("atmosphere", "atmosphere.py")
    if atmosphere is None:
        atmosphere = atm_mod.USStandard1976()
    _init_site()                                   # sets Z_SURF_SITE, Z_DET_ASL
    bnds = atm_mod.default_boundaries(Z_SURF_SITE, atmo_top_km * 1000.0)
    atmo_shells = [(f"atmo_{k}", R_PREM + (asl_top - Z_SURF_SITE), "AIR", rho)
                   for k, (asl_top, rho) in enumerate(atmosphere.shells(bnds))]

    # ---- PREM + atmosphere (far negative levels; innermost = highest wins) -
    shells = _PREM + atmo_shells                   # innermost first (ascending R)
    n_sh = len(shells)
    for i, (name, Rsh, material, rho) in enumerate(shells):
        sector(name, material, -200 + (n_sh - 1 - i), Sphere(earth_place, float(Rsh), 0.0), rho)

    # ---- local crust bridge + air (built in ENU, placed into site frame) --
    deep = -15000.0
    sky = 5000.0
    sector("local_crust", "ROCK", 1,
           Box(Placement(V3(0, (BASE_Z + deep) / 2, 0), q_enu_site),
               2 * L_HALF, 2 * L_HALF, BASE_Z - deep), 2.75)
    # local_air density from the same atmosphere model: mass-conserving over the
    # box's air column (site surface up to the box top; ENU sky -> ASL via Z_DET_ASL).
    air_rho = atmosphere.shell_density_gcc(Z_SURF_SITE, Z_DET_ASL + sky)
    sector("local_air", "AIR", 2,
           Box(Placement(V3(0, (BASE_Z + sky) / 2, 0), q_enu_site),
               2 * L_HALF, 2 * L_HALF, sky - BASE_Z), air_rho)

    # ---- graded terrain mesh + formations -------------------------------
    tris_enu, leaves = _build_graded_mesh(V3)
    mesh = TriangularMesh(enu_placement, tris_enu)
    sector("Poorman_bulk", "Poorman_phyllite", 3, mesh, 2.80)

    dd = math.radians(CONTACT_DIP); cz = math.radians(CONTACT_DIPDIR)
    nx, ny, nz = -math.sin(dd)*math.sin(cz), -math.sin(dd)*math.cos(cz), math.cos(dd)
    Hz = 8000.0
    rot_cont = Quaternion.rotation_between(V3(0, 0, 1), V3(nx, ny, nz))
    cen_enu = V3(-nx*Hz/2, -ny*Hz/2, CONTACT_Z - nz*Hz/2)
    yates_slab = Box(Placement(cen_enu, rot_cont), 120000.0, 120000.0, Hz)
    yates_bool = BooleanGeometry(enu_placement, BooleanOperation.INTERSECTION,
                                 TriangularMesh(tris_enu), yates_slab)
    sector("Yates_lens", "Yates_amphibolite", 4, yates_bool, 2.95)

    ceiling = DEPTH + 700
    for i, (x0, y0, strike_az, t, L) in enumerate(DIKES):
        a = math.radians(strike_az)
        rotd = Quaternion.rotation_between(V3(1, 0, 0), V3(math.cos(a), -math.sin(a), 0.0))
        dike = Box(Placement(V3(x0, y0, 0.5*(BASE_Z + ceiling)), rotd),
                   t, L, ceiling - BASE_Z)
        dike_bool = BooleanGeometry(enu_placement, BooleanOperation.INTERSECTION,
                                    TriangularMesh(tris_enu), dike)
        sector(f"dike_{i}", "Rhyolite_dike", 5 + i, dike_bool, 2.62)
