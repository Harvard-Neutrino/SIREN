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

# Generated/downloaded data (DEM tiles, materials .dat) goes in the writable
# mirror of this resource directory: _THIS_DIR itself may be a read-only
# site-packages install. Sibling .py sources are still read from _THIS_DIR.
from siren.download import writable_data_dir
_ABS_DIR = writable_data_dir(_THIS_DIR)

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
_DEM_DIR = os.path.join(_ABS_DIR, "dem")
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
    try:
        import tifffile
    except ImportError as e:
        raise ImportError(
            "The DUNE FD terrain model reads Copernicus GLO-30 DEM tiles, "
            "which requires the optional 'tifffile' package. Install it with "
            "'pip install tifffile' to use earth_model=True.") from e
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


# The PREM shells and the atmosphere are the canonical, shared ones from
# siren.earth (see add_earth_model); this module supplies only the DUNE-specific
# Homestake overburden (DEM terrain, formations, dikes) built on top of them.

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
    from siren.geometry import Box, Placement, TriangularMesh, BooleanGeometry, BooleanOperation
    from siren.detector import DetectorSector, ConstantDensityDistribution
    from siren.earth import (
        CANONICAL_PREM_LAYERS,
        USStandard1976,
        AtmosphereForm,
        build_earth_sectors,
        insert_sectors_above,
    )
    V3 = Vector3D

    # ---- load materials (supplements the GDML's own) ----
    geo_mod = _load_sibling("homestake_geology", "homestake_geology.py")
    mat_path = os.path.join(_ABS_DIR, "DUNEFD_homestake_materials.dat")
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

    # ---- canonical PREM + exponential atmosphere via the shared builder ----
    # siren.earth places PREM/atmosphere below every existing volume using the
    # predict-count-and-shift-up level scheme (World and the detector volumes
    # shift up to make room). The site surface sits at radius R_PREM because the
    # Earth center is DEPTH below the detector origin (y ~ 0 == 4850L).
    if atmosphere is None:
        atmosphere = USStandard1976()
    _init_site()                                   # sets Z_SURF_SITE, Z_DET_ASL

    build_earth_sectors(
        model,
        surface_offset=DEPTH,
        prem_layers=CANONICAL_PREM_LAYERS,
        atmosphere=AtmosphereForm(model=atmosphere, mode="exponential",
                                  top_m=atmo_top_km * 1000.0),
        up_axis=(0.0, 1.0, 0.0),
        r_prem=R_PREM,
    )

    # ---- site overburden: built in ENU, inserted just above "World" --------
    # The local crust/air and Homestake formations must override the vacuum
    # World box but stay below the detector volumes. insert_sectors_above shifts
    # the detector volumes up and gives the overburden the freed levels
    # (ascending priority: local_crust lowest, dikes highest) -- no absolute or
    # negative level constants.
    def _geo(name, material, geo, rho):
        s = DetectorSector()
        s.name = name
        s.material_id = mats.GetMaterialId(material)
        s.geo = geo
        s.density = ConstantDensityDistribution(float(rho))
        return s

    overburden = []

    # ---- local crust bridge + air (built in ENU, placed into site frame) --
    deep = -15000.0
    sky = 5000.0
    overburden.append(_geo("local_crust", "ROCK",
           Box(Placement(V3(0, (BASE_Z + deep) / 2, 0), q_enu_site),
               2 * L_HALF, 2 * L_HALF, BASE_Z - deep), 2.75))
    # local_air density from the same atmosphere model: mass-conserving mean of
    # the air column from the site surface up to the box top (altitudes measured
    # above the site surface, ASL top = Z_DET_ASL + sky).
    air_rho = atmosphere.shell_mean_density(0.0, max(1.0, (Z_DET_ASL + sky) - Z_SURF_SITE))
    overburden.append(_geo("local_air", "AIR",
           Box(Placement(V3(0, (BASE_Z + sky) / 2, 0), q_enu_site),
               2 * L_HALF, 2 * L_HALF, sky - BASE_Z), air_rho))

    # ---- graded terrain mesh + formations -------------------------------
    # Build the DEM mesh (and its BVH) once in each frame it is needed:
    # `mesh` carries the ENU->site placement and is the Poorman_bulk sector
    # geometry; `terrain_enu` is the same triangles with identity placement,
    # shared by the boolean intersections below, whose components live in the
    # boolean's local (ENU) frame -- the enclosing BooleanGeometry applies the
    # ENU->site placement, so a placed mesh there would rotate twice.
    tris_enu, leaves = _build_graded_mesh(V3)
    mesh = TriangularMesh(enu_placement, tris_enu)
    terrain_enu = TriangularMesh(tris_enu)
    overburden.append(_geo("Poorman_bulk", "Poorman_phyllite", mesh, 2.80))

    dd = math.radians(CONTACT_DIP); cz = math.radians(CONTACT_DIPDIR)
    nx, ny, nz = -math.sin(dd)*math.sin(cz), -math.sin(dd)*math.cos(cz), math.cos(dd)
    Hz = 8000.0
    rot_cont = Quaternion.rotation_between(V3(0, 0, 1), V3(nx, ny, nz))
    cen_enu = V3(-nx*Hz/2, -ny*Hz/2, CONTACT_Z - nz*Hz/2)
    yates_slab = Box(Placement(cen_enu, rot_cont), 120000.0, 120000.0, Hz)
    yates_bool = BooleanGeometry(enu_placement, BooleanOperation.INTERSECTION,
                                 terrain_enu, yates_slab)
    overburden.append(_geo("Yates_lens", "Yates_amphibolite", yates_bool, 2.95))

    ceiling = DEPTH + 700
    for i, (x0, y0, strike_az, t, L) in enumerate(DIKES):
        a = math.radians(strike_az)
        rotd = Quaternion.rotation_between(V3(1, 0, 0), V3(math.cos(a), -math.sin(a), 0.0))
        dike = Box(Placement(V3(x0, y0, 0.5*(BASE_Z + ceiling)), rotd),
                   t, L, ceiling - BASE_Z)
        dike_bool = BooleanGeometry(enu_placement, BooleanOperation.INTERSECTION,
                                    terrain_enu, dike)
        overburden.append(_geo(f"dike_{i}", "Rhyolite_dike", dike_bool, 2.62))

    insert_sectors_above(model, "World", overburden)
