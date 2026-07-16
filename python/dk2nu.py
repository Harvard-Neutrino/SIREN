"""
Read beam-parent kinematics from dk2nu ROOT files (siren.dk2nu).

dk2nu is a standardized format for neutrino flux simulation output
(see https://github.com/NuSoftHEP/dk2nu).  Each entry describes a
neutrino produced in a particle decay, recording the parent's identity,
momentum, decay vertex, and importance weight, plus the ancestor chain
and the parent-of-parent kinematics.

This module extracts that information and provides it in forms suitable
for SIREN injection:

    - As numpy arrays for direct use in flux calculations, including the
      neutrino production time from the ancestor chain and the muon spin
      axis in the bsim::calcEnuWgt convention
    - As a PrimaryExternalDistribution that injects the recorded parents
      at their decay vertices, times, and momenta
    - As a TabulatedFluxDistribution or energy spectra binned by parent
      species

Reading the ROOT files requires uproot, imported when read_dk2nu is
called; the rest of the module has no optional dependencies.
"""

import math

import numpy as np

from . import distributions as _distributions
from .detector import GeometryPosition, GeometryDirection
from .errors import ConfigurationError
from .math import Vector3D


# PDG codes for parent mesons
PTYPE_PIPLUS = 211
PTYPE_PIMINUS = -211
PTYPE_KPLUS = 321
PTYPE_KMINUS = -321
PTYPE_K0L = 130

_PARENT_NAMES = {
    PTYPE_PIPLUS: "pi+",
    PTYPE_PIMINUS: "pi-",
    PTYPE_KPLUS: "K+",
    PTYPE_KMINUS: "K-",
    PTYPE_K0L: "K0L",
    13: "mu-",
    -13: "mu+",
}

# dk2nu branch paths (nested under dk2nu/ in NuMI files, flat in BNB).
# The muon production kinematics (ppdxdz, ppdydz, pppz with ppenergy) and
# the parent-of-muon momentum (muparpx/y/z, mupare) are part of the
# standard bsim::Decay block; they carry the information for the muon
# polarization axis.
_BRANCH_SETS = [
    {
        "ptype": "dk2nu/decay/decay.ptype",
        "pdpx": "dk2nu/decay/decay.pdpx",
        "pdpy": "dk2nu/decay/decay.pdpy",
        "pdpz": "dk2nu/decay/decay.pdpz",
        "ppenergy": "dk2nu/decay/decay.ppenergy",
        "ppdxdz": "dk2nu/decay/decay.ppdxdz",
        "ppdydz": "dk2nu/decay/decay.ppdydz",
        "pppz": "dk2nu/decay/decay.pppz",
        "mupare": "dk2nu/decay/decay.mupare",
        "muparpx": "dk2nu/decay/decay.muparpx",
        "muparpy": "dk2nu/decay/decay.muparpy",
        "muparpz": "dk2nu/decay/decay.muparpz",
        "vx": "dk2nu/decay/decay.vx",
        "vy": "dk2nu/decay/decay.vy",
        "vz": "dk2nu/decay/decay.vz",
        "nimpwt": "dk2nu/decay/decay.nimpwt",
        "ntype": "dk2nu/decay/decay.ntype",
        "ndecay": "dk2nu/decay/decay.ndecay",
    },
    {
        "ptype": "decay.ptype",
        "pdpx": "decay.pdpx",
        "pdpy": "decay.pdpy",
        "pdpz": "decay.pdpz",
        "ppenergy": "decay.ppenergy",
        "ppdxdz": "decay.ppdxdz",
        "ppdydz": "decay.ppdydz",
        "pppz": "decay.pppz",
        "mupare": "decay.mupare",
        "muparpx": "decay.muparpx",
        "muparpy": "decay.muparpy",
        "muparpz": "decay.muparpz",
        "vx": "decay.vx",
        "vy": "decay.vy",
        "vz": "decay.vz",
        "nimpwt": "decay.nimpwt",
        "ntype": "decay.ntype",
        "ndecay": "decay.ndecay",
    },
]

# Decay-mode codes for muon decay in flight (bsim::dkproc_t).
_NDECAY_MUPLUS = 11
_NDECAY_MUMINUS = 12

# The muon mass value used by bsim::calcEnuWgt, reproduced here so the
# polarization axis matches the dk2nu reference computation exactly.
_MUON_MASS = 0.1056583715


def _muon_polarization(data, branches):
    """Muon spin axis in the muon rest frame, following bsim::calcEnuWgt.

    A muon from meson two-body decay is fully polarized: in the meson rest
    frame its spin points opposite its momentum for mu+ and along it for
    mu-. Boosted to the muon rest frame, that axis is the direction of the
    meson momentum there, with the same orientation convention for both
    charges (+spin for mu+, -spin for mu-), which is why calcEnuWgt can
    measure all decay angles from the parent direction. This reproduces
    that computation: boost the parent-of-muon momentum (muparpx/y/z,
    mupare) into the muon rest frame using the muon production kinematics
    (ppdxdz, ppdydz, pppz, ppenergy), and return the unit vector in beam
    coordinates. Rows that are not muon decays in flight, or that carry no
    parent-of-muon information, get a zero vector, which downstream models
    read as unpolarized.
    """
    ptype = data[branches["ptype"]].astype(int)
    ndecay = data[branches["ndecay"]].astype(int)
    pp_energy = np.asarray(data[branches["ppenergy"]], dtype=float)
    pppz = np.asarray(data[branches["pppz"]], dtype=float)
    beta = np.stack([
        np.asarray(data[branches["ppdxdz"]], dtype=float) * pppz,
        np.asarray(data[branches["ppdydz"]], dtype=float) * pppz,
        pppz,
    ], axis=1)
    mupar = np.stack([
        np.asarray(data[branches["muparpx"]], dtype=float),
        np.asarray(data[branches["muparpy"]], dtype=float),
        np.asarray(data[branches["muparpz"]], dtype=float),
    ], axis=1)
    mupare = np.asarray(data[branches["mupare"]], dtype=float)

    axis = np.zeros_like(beta)
    rows = ((np.abs(ptype) == 13)
            & np.isin(ndecay, (_NDECAY_MUPLUS, _NDECAY_MUMINUS))
            & (pp_energy > _MUON_MASS)
            & (mupare > 0))
    if np.any(rows):
        gamma = pp_energy[rows] / _MUON_MASS
        b = beta[rows] / pp_energy[rows, None]
        partial = gamma * np.einsum("ij,ij->i", b, mupar[rows])
        partial = mupare[rows] - partial / (gamma + 1.0)
        p_pcm = mupar[rows] - b * (gamma * partial)[:, None]
        norm = np.linalg.norm(p_pcm, axis=1)
        unit = np.zeros_like(p_pcm)
        good = norm > 1e-12
        unit[good] = p_pcm[good] / norm[good, None]
        axis[rows] = unit
    return {"pol_x": axis[:, 0], "pol_y": axis[:, 1], "pol_z": axis[:, 2]}

# Ancestor branch paths, used for the neutrino production time. The last
# entry of the ancestor chain is the neutrino itself; its start time is the
# parent decay time in nanoseconds relative to the primary proton.
_ANCESTOR_SETS = [
    {"pdg": "dk2nu/ancestor/ancestor.pdg",
     "startt": "dk2nu/ancestor/ancestor.startt"},
    {"pdg": "ancestor.pdg", "startt": "ancestor.startt"},
]


def _detect_branches(tree):
    """Auto-detect which branch naming convention the file uses."""
    keys = set(tree.keys())
    for bset in _BRANCH_SETS:
        if all(v in keys for v in bset.values()):
            return bset
    raise ValueError(
        "Could not find dk2nu decay branches. "
        f"Available branches: {sorted(keys)[:20]}..."
    )


def read_dk2nu(
    filenames,
    parent_pdg=None,
    decay_modes=None,
    nu_pdg=None,
    read_time=True,
    read_polarization=True,
    entry_start=None,
    entry_stop=None,
):
    """
    Read parent meson kinematics from one or more dk2nu ROOT files.

    Parameters
    ----------
    filenames : str or list of str
        Path(s) to dk2nu ROOT files.  Files without a dk2nuTree (for
        example a truncated write) are skipped.
    parent_pdg : int or list of int, optional
        Filter to specific parent PDG code(s).  Default: all parents.
    decay_modes : int or list of int, optional
        Filter to specific dk2nu decay-mode code(s) (the ndecay branch),
        e.g. 13 for pi+ -> mu+ nu_mu or 5 for K+ -> mu+ nu_mu.  Use this
        when injecting parents whose SIREN decay model implements one
        channel, so every row belongs to that channel.
    nu_pdg : int or list of int, optional
        Filter on the recorded neutrino PDG code (the ntype branch).  A
        muon decay writes one row per neutrino; selecting one flavor keeps
        one row per decay.
    read_time : bool
        Read the neutrino production time from the ancestor chain into a
        "t0" key [ns relative to the primary proton].  Rows whose last
        ancestor is not the recorded neutrino are dropped.  Files without
        ancestor branches leave "t0" absent (with a notice).
    read_polarization : bool
        Compute the muon spin axis for muon-decay rows into pol_x, pol_y,
        pol_z keys (unit vector in beam coordinates, muon rest frame; see
        _muon_polarization).  Non-muon rows get a zero vector, read
        downstream as unpolarized.
    entry_start, entry_stop : int, optional
        Limit the number of entries read (per file).

    Returns
    -------
    dict with keys:
        ptype      : int array, parent PDG code
        E          : float array, parent energy [GeV]
        px, py, pz : float arrays, parent momentum components [GeV]
        vx, vy, vz : float arrays, decay vertex [cm]
        nimpwt     : float array, importance weight
        ntype      : int array, neutrino PDG code
        ndecay     : int array, dk2nu decay-mode code
        t0         : float array, production time [ns] (when read_time and
                     the files carry ancestor branches)
        pol_x, pol_y, pol_z : float arrays, muon spin axis (when
                     read_polarization)
        pot        : float, total simulated POT across all files
    """
    try:
        import uproot
    except ImportError:
        raise ImportError(
            "uproot is required to read dk2nu files: pip install uproot"
        )

    if isinstance(filenames, str):
        filenames = [filenames]

    def _as_list(value):
        if value is None:
            return None
        if not hasattr(value, "__iter__"):
            return [value]
        return list(value)

    parent_pdg = _as_list(parent_pdg)
    decay_modes = _as_list(decay_modes)
    nu_pdg = _as_list(nu_pdg)

    keys = ["ptype", "E", "px", "py", "pz", "vx", "vy", "vz",
            "nimpwt", "ntype", "ndecay"]
    if read_polarization:
        keys += ["pol_x", "pol_y", "pol_z"]
    all_data = {k: [] for k in keys + ["t0"]}
    have_time = read_time
    total_pot = 0.0

    for fname in filenames:
        f = uproot.open(fname)
        if "dk2nuTree" not in [k.split(";")[0] for k in f.keys()]:
            continue
        tree = f["dk2nuTree"]
        branches = _detect_branches(tree)

        kw = {}
        if entry_start is not None:
            kw["entry_start"] = entry_start
        if entry_stop is not None:
            kw["entry_stop"] = entry_stop

        data = tree.arrays(
            list(branches.values()),
            library="np",
            **kw,
        )

        columns = {
            "ptype": data[branches["ptype"]].astype(int),
            "E": data[branches["ppenergy"]],
            "px": data[branches["pdpx"]],
            "py": data[branches["pdpy"]],
            "pz": data[branches["pdpz"]],
            "vx": data[branches["vx"]],
            "vy": data[branches["vy"]],
            "vz": data[branches["vz"]],
            "nimpwt": data[branches["nimpwt"]],
            "ntype": data[branches["ntype"]].astype(int),
            "ndecay": data[branches["ndecay"]].astype(int),
        }
        if read_polarization:
            columns.update(_muon_polarization(data, branches))

        mask = np.ones(len(columns["ptype"]), dtype=bool)
        if parent_pdg is not None:
            mask &= np.isin(columns["ptype"], parent_pdg)
        if decay_modes is not None:
            mask &= np.isin(columns["ndecay"], decay_modes)
        if nu_pdg is not None:
            mask &= np.isin(columns["ntype"], nu_pdg)

        if have_time:
            anc = None
            tree_keys = set(tree.keys())
            for aset in _ANCESTOR_SETS:
                if all(v in tree_keys for v in aset.values()):
                    anc = aset
                    break
            if anc is None:
                print("  %s carries no ancestor branches; production "
                      "times unavailable" % fname)
                have_time = False
            else:
                arr = tree.arrays(list(anc.values()), library="ak", **kw)
                last_pdg = np.asarray(arr[anc["pdg"]][:, -1], dtype=int)
                last_t = np.asarray(arr[anc["startt"]][:, -1], dtype=float)
                mask &= last_pdg == columns["ntype"]
                columns["t0"] = last_t

        for k, v in columns.items():
            all_data[k].append(v[mask])

        # uproot method to get total POT from dkmetaTree (if available)
        if "dkmetaTree" in f:
            meta_tree = f["dkmetaTree"]
            if "dkmeta/pots" in meta_tree.keys():
                pots = meta_tree["dkmeta/pots"].array(library="np")
                if len(pots) > 0:
                    pots = pots[0]
                else:
                    pots = 0.0
                total_pot += pots

    if not have_time:
        all_data.pop("t0", None)
    result = {k: np.concatenate(v) if v else np.array([])
              for k, v in all_data.items()}
    result["pot"] = total_pot
    return result


def meson_energy_spectrum(
    dk2nu_data,
    parent_pdg,
    n_bins=100,
    E_range=None,
):
    """
    Build a weighted energy spectrum of a given parent meson from dk2nu data.

    Parameters
    ----------
    dk2nu_data : dict
        Output of read_dk2nu().
    parent_pdg : int
        PDG code of the parent meson to histogram.
    n_bins : int
        Number of energy bins.
    E_range : tuple of (float, float), optional
        Energy range in GeV.  Default: auto from data.

    Returns
    -------
    E_centers : array, bin centers [GeV]
    dN_dE     : array, weighted counts per GeV per POT
    """
    mask = dk2nu_data["ptype"] == parent_pdg
    E = dk2nu_data["E"][mask]
    w = dk2nu_data["nimpwt"][mask]
    pot = dk2nu_data["pot"]

    if E_range is None:
        E_range = (E.min() * 0.9, E.max() * 1.1) if len(E) > 0 else (0.0, 5.0)

    counts, bin_edges = np.histogram(E, bins=n_bins, range=E_range, weights=w)
    dE = bin_edges[1] - bin_edges[0]
    E_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    if pot > 0:
        dN_dE = counts / (dE * pot)
    else:
        dN_dE = counts / dE

    return E_centers, dN_dE


def dk2nu_to_tabulated_flux(
    dk2nu_data,
    parent_pdg,
    n_bins=100,
    E_range=None,
    physically_normalized=True,
):
    """
    Convert dk2nu parent meson data into a SIREN TabulatedFluxDistribution
    (energy spectrum of the parent meson, per POT).

    Parameters
    ----------
    dk2nu_data : dict
        Output of read_dk2nu().
    parent_pdg : int
        PDG code of the parent meson.
    n_bins : int
        Number of energy bins.
    E_range : tuple, optional
        Energy range in GeV.
    physically_normalized : bool
        Whether the flux is physically normalized.

    Returns
    -------
    siren.distributions.TabulatedFluxDistribution
    """
    E_centers, dN_dE = meson_energy_spectrum(
        dk2nu_data, parent_pdg, n_bins=n_bins, E_range=E_range
    )

    return _distributions.TabulatedFluxDistribution(
        float(E_centers[0]),
        float(E_centers[-1]),
        list(E_centers),
        list(dN_dE),
        physically_normalized,
    )


def dk2nu_to_primary_distribution(
    dk2nu_data,
    detector_model,
    parent_pdg=None,
    sampling_bias=None,
):
    """
    Build a PrimaryExternalDistribution directly from dk2nu data.

    Converts positions from BNB (geometry) coordinates to detector-local
    coordinates using the detector model's ToDet transform, and from cm
    to meters.  No intermediate CSV file is needed.

    Parameters
    ----------
    dk2nu_data : dict
        Output of read_dk2nu().
    detector_model : siren.detector.DetectorModel
        Detector model (provides the geometry-to-detector transform).
    parent_pdg : int or list of int, optional
        Filter to specific parent PDG code(s).
    sampling_bias : callable, optional
        Function f(E, px, py, pz, vx, vy, vz) -> weight that computes
        per-entry sampling weights from the pion kinematics (in geometry
        coordinates, before the detector transform).  Entries are selected
        with probability proportional to these weights.  The generation
        probability accounts for the bias so event weights remain correct.
        Arguments are numpy arrays; the return value should broadcast to
        the same length.  When None (default), uniform selection is used.

    Returns
    -------
    siren.distributions.PrimaryExternalDistribution
    """
    ptype = dk2nu_data["ptype"]
    if parent_pdg is not None:
        if not hasattr(parent_pdg, "__iter__"):
            parent_pdg = [parent_pdg]
        mask = np.isin(ptype, parent_pdg)
    else:
        mask = np.ones(len(ptype), dtype=bool)

    simulated_pot = dk2nu_data["pot"]
    # Per-POT weights are meaningless without a positive POT. read_dk2nu leaves
    # pot at 0.0 when a file has no dkmetaTree/pots branch; dividing by it would
    # emit inf/nan weights silently. Fail loud at the point the weights are
    # formed rather than propagate a corrupt distribution.
    if not (simulated_pot > 0):
        raise ConfigurationError(
            "dk2nu_data['pot'] is %r; a positive simulated POT is required to "
            "compute per-POT weights. The input file(s) carried no POT metadata "
            "(no dkmetaTree/pots branch)." % (simulated_pot,))

    E = dk2nu_data["E"][mask]
    px = dk2nu_data["px"][mask]
    py = dk2nu_data["py"][mask]
    pz = dk2nu_data["pz"][mask]
    vx = dk2nu_data["vx"][mask]
    vy = dk2nu_data["vy"][mask]
    vz = dk2nu_data["vz"][mask]
    nimpwt = dk2nu_data["nimpwt"][mask]
    pt = ptype[mask]
    t0 = dk2nu_data["t0"][mask] if "t0" in dk2nu_data else None

    # Muon spin axes ride along as pol_x/y/z columns (and from there into
    # each record's interaction parameters) when any selected row carries
    # one; a table of meson rows stays free of them.
    pol = None
    if all(k in dk2nu_data for k in ("pol_x", "pol_y", "pol_z")):
        pol = np.stack([dk2nu_data["pol_x"][mask],
                        dk2nu_data["pol_y"][mask],
                        dk2nu_data["pol_z"][mask]], axis=1)
        if not np.any(np.abs(pol) > 0):
            pol = None

    weight = nimpwt / simulated_pot

    mass_map = {
        211: 0.13957039, -211: 0.13957039,
        321: 0.49368,    -321: 0.49368,
        130: 0.49761,
        13: 0.10566,     -13: 0.10566,
    }

    # The decay time from the ancestor chain becomes the primary's initial
    # time (the t0 column of PrimaryExternalDistribution), so interaction
    # times downstream are physical.
    keys = ["E", "px", "py", "pz", "x", "y", "z", "m", "weight"]
    if t0 is not None:
        keys.append("t0")
    if pol is not None:
        keys += ["pol_x", "pol_y", "pol_z"]
    data = []
    for i in range(len(E)):
        # Convert position from geometry (BNB) to detector coordinates
        geo_pos = GeometryPosition(Vector3D(
            vx[i] * 0.01, vy[i] * 0.01, vz[i] * 0.01
        ))
        det_pos = detector_model.GeoPositionToDetPosition(geo_pos).get()

        # Convert momentum direction from geometry to detector coordinates.
        # Energy is a scalar and is unchanged; the 3-momentum direction
        # must be rotated if the detector axes differ from geometry axes.
        p_mag = math.sqrt(float(px[i])**2 + float(py[i])**2 + float(pz[i])**2)
        if p_mag > 0:
            geo_dir = GeometryDirection(Vector3D(
                float(px[i]) / p_mag, float(py[i]) / p_mag, float(pz[i]) / p_mag))
            det_dir = detector_model.GeoDirectionToDetDirection(geo_dir).get()
            px_det = det_dir.GetX() * p_mag
            py_det = det_dir.GetY() * p_mag
            pz_det = det_dir.GetZ() * p_mag
        else:
            px_det = py_det = pz_det = 0.0

        m = mass_map.get(int(pt[i]), 0.13957)
        # dk2nu stores momentum at decay (pdpx/pdpy/pdpz) but energy
        # at production (ppenergy). Compute on-shell energy from the
        # decay-point momentum and known mass.
        E_decay = math.sqrt(p_mag * p_mag + m * m)
        row = [
            E_decay, px_det, py_det, pz_det,
            det_pos.GetX(), det_pos.GetY(), det_pos.GetZ(),
            m, float(weight[i]),
        ]
        if t0 is not None:
            row.append(float(t0[i]))
        if pol is not None:
            # The spin axis is a rest-frame direction expressed in beam
            # coordinates; rotate it into detector coordinates like the
            # momentum (both frames are reached by pure boosts, so their
            # spatial bases are related by the same rotation).
            norm = float(np.linalg.norm(pol[i]))
            if norm > 0:
                geo_ax = GeometryDirection(Vector3D(*(pol[i] / norm)))
                det_ax = detector_model.GeoDirectionToDetDirection(geo_ax).get()
                row += [det_ax.GetX() * norm, det_ax.GetY() * norm,
                        det_ax.GetZ() * norm]
            else:
                row += [0.0, 0.0, 0.0]
        data.append(row)

    if sampling_bias is not None:
        sw = np.asarray(
            sampling_bias(E, px, py, pz, vx, vy, vz),
            dtype=float,
        )
        np.maximum(sw, 0.0, out=sw)
        return _distributions.PrimaryExternalDistribution(
            keys, data, sw.tolist()
        )

    return _distributions.PrimaryExternalDistribution(keys, data)


def dk2nu_to_csv(
    dk2nu_data,
    output_path,
    parent_pdg=None,
    position_transform=None,
    units_cm=True,
):
    """
    Write dk2nu parent meson kinematics to a CSV file suitable for
    SIREN's PrimaryExternalDistribution.

    The CSV has columns: E, px, py, pz, x0, y0, z0, m, nimpwt

    Parameters
    ----------
    dk2nu_data : dict
        Output of read_dk2nu().
    output_path : str
        Path to write the CSV file.
    parent_pdg : int or list of int, optional
        Filter to specific parent PDG code(s).  Default: use all entries
        in dk2nu_data (which may already be filtered).
    position_transform : callable, optional
        Function that takes (vx, vy, vz) arrays in dk2nu coordinates
        and returns (x0, y0, z0) arrays in detector coordinates.
        dk2nu positions are in cm.  If None, positions are used as-is.
    units_cm : bool
        If True (default), positions in the CSV are in cm.
        If False, positions are converted to meters.

    Returns
    -------
    int
        Number of rows written.
    """
    ptype = dk2nu_data["ptype"]
    if parent_pdg is not None:
        if not hasattr(parent_pdg, "__iter__"):
            parent_pdg = [parent_pdg]
        mask = np.isin(ptype, parent_pdg)
    else:
        mask = np.ones(len(ptype), dtype=bool)

    E = dk2nu_data["E"][mask]
    px = dk2nu_data["px"][mask]
    py = dk2nu_data["py"][mask]
    pz = dk2nu_data["pz"][mask]
    vx = dk2nu_data["vx"][mask]
    vy = dk2nu_data["vy"][mask]
    vz = dk2nu_data["vz"][mask]
    nimpwt = dk2nu_data["nimpwt"][mask]
    pt = ptype[mask]

    if position_transform is not None:
        vx, vy, vz = position_transform(vx, vy, vz)

    scale = 1.0 if units_cm else 0.01

    mass_map = {
        211: 0.13957039, -211: 0.13957039,
        321: 0.49368,    -321: 0.49368,
        130: 0.49761,
        13: 0.10566,     -13: 0.10566,
    }

    with open(output_path, "w") as f:
        f.write("E,px,py,pz,x,y,z,m,nimpwt\n")
        for i in range(len(E)):
            m = mass_map.get(int(pt[i]), 0.13957)
            f.write(
                f"{E[i]:.8e},{px[i]:.8e},{py[i]:.8e},{pz[i]:.8e},"
                f"{vx[i]*scale:.8e},{vy[i]*scale:.8e},{vz[i]*scale:.8e},"
                f"{m:.8e},{nimpwt[i]:.8e}\n"
            )

    return int(np.sum(mask))


def print_summary(dk2nu_data):
    """Print a summary of the dk2nu data."""
    ptypes = dk2nu_data["ptype"]
    unique, counts = np.unique(ptypes, return_counts=True)
    n_total = len(ptypes)
    pot = dk2nu_data["pot"]

    print(f"Total entries: {n_total}")
    print(f"Total POT: {pot:.3e}")
    print(f"Parent breakdown:")
    for pdg, count in sorted(zip(unique, counts), key=lambda x: -x[1]):
        name = _PARENT_NAMES.get(int(pdg), str(int(pdg)))
        frac = 100.0 * count / n_total if n_total > 0 else 0
        E = dk2nu_data["E"][ptypes == pdg]
        print(f"  {name:>5s} ({int(pdg):>4d}): {count:>7d} ({frac:5.1f}%)  "
              f"E = [{E.min():.3f}, {E.max():.3f}] GeV")
