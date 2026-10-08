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

Coordinate systems: a dk2nu file is written in its beam simulation's own
frame (G4BNB files in BNB coordinates, g4numi files in NuMI coordinates,
and so on), which is not in general the frame of the detector model the
events are injected into. dk2nu_to_primary_distribution accepts a
FrameTransform that maps the file's coordinates into the detector
model's geometry (or directly detector) coordinates; the default assumes
the file frame and the model's geometry frame coincide. Detector
resources that know their beamline surveys can hand any object carrying
a rotation and a translation straight to that argument.

Reading the ROOT files requires uproot, imported when read_dk2nu is
called; the rest of the module has no optional dependencies.
"""

import math
from collections.abc import Mapping
from numbers import Integral

import numpy as np

from . import distributions as _distributions
from . import particles
from .dataclasses import ParticleType
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


class FrameTransform:
    """Rigid map from a dk2nu file's coordinate system into a detector
    model's frames.

    Positions transform as r_out = rotation @ r_file + translation and
    directions (momenta, spin axes) as d_out = rotation @ d_file, all in
    meters and applied after the file's centimeters are converted. The
    ``target`` names the frame the transform lands in: ``"geometry"``
    (default) for the detector model's geometry coordinates, after which
    the model's geometry-to-detector conversion runs as usual, or
    ``"detector"`` for detector-local coordinates directly, in which case
    the detector model is not consulted at all.

    The rotation must be a proper rotation (orthonormal, determinant +1);
    a reflection would silently flip the geometry's handedness, so it is
    rejected. Beam-frame conventions published as GNuMIFlux beamdir and
    beampos elements, site survey matrices, or a detector resource's
    frame-graph transforms all provide exactly this (rotation,
    translation) pair.
    """

    def __init__(self, rotation=None, translation=None, target="geometry"):
        if target not in ("geometry", "detector"):
            raise ConfigurationError(
                "FrameTransform target must be 'geometry' or 'detector', "
                "got %r" % (target,))
        self.target = target
        if rotation is None:
            self.rotation = np.eye(3)
        else:
            try:
                self.rotation = np.array(rotation, dtype=float)
            except (TypeError, ValueError) as exc:
                raise ConfigurationError(
                    "FrameTransform rotation is not numeric: %s" % exc)
            if self.rotation.shape != (3, 3):
                raise ConfigurationError(
                    "FrameTransform rotation must be a 3x3 matrix, got "
                    "shape %r" % (self.rotation.shape,))
            if not np.allclose(self.rotation @ self.rotation.T, np.eye(3),
                               atol=1e-8):
                raise ConfigurationError(
                    "FrameTransform rotation is not orthonormal")
            if np.linalg.det(self.rotation) < 0:
                raise ConfigurationError(
                    "FrameTransform rotation is a reflection (determinant "
                    "-1); coordinate frames must keep their handedness")
        if translation is None:
            self.translation = np.zeros(3)
        else:
            try:
                self.translation = np.array(translation, dtype=float)
            except (TypeError, ValueError) as exc:
                raise ConfigurationError(
                    "FrameTransform translation is not numeric: %s" % exc)
            if self.translation.shape != (3,):
                raise ConfigurationError(
                    "FrameTransform translation must be a 3-vector, got "
                    "shape %r" % (self.translation.shape,))

    def position(self, r):
        """Transform a position [m] from the file frame."""
        return self.rotation @ np.asarray(r, dtype=float) + self.translation

    def direction(self, d):
        """Transform a direction or momentum from the file frame."""
        return self.rotation @ np.asarray(d, dtype=float)

    def __repr__(self):
        return ("FrameTransform(rotation=%s, translation=%s, target=%r)"
                % (self.rotation.tolist(), self.translation.tolist(),
                   self.target))


def _as_frame_transform(frame):
    """Normalize the ``frame`` argument of dk2nu_to_primary_distribution.

    Accepts None or "geometry" (the file frame is the detector model's
    geometry frame), "detector" (the file frame is the detector-local
    frame), a FrameTransform, a (rotation, translation) pair, or any
    object carrying the pair as ``rotation``/``translation`` or ``R``/
    ``t`` attributes -- the shape a detector resource's frame-graph
    transform naturally has. Duck-typed objects always target the
    geometry frame; wrap in a FrameTransform to target the detector
    frame explicitly.
    """
    if frame is None or (isinstance(frame, str) and frame == "geometry"):
        return FrameTransform()
    if isinstance(frame, str):
        if frame == "detector":
            return FrameTransform(target="detector")
        raise ConfigurationError(
            "Unknown frame %r; expected 'geometry', 'detector', a "
            "FrameTransform, a (rotation, translation) pair, or an object "
            "with rotation/translation (or R/t) attributes" % (frame,))
    if isinstance(frame, FrameTransform):
        return frame
    if isinstance(frame, (tuple, list)) and len(frame) == 2:
        return FrameTransform(frame[0], frame[1])

    def _attribute_pair(rotation_name, translation_name):
        # An attribute pair only counts when it is data: a class may carry
        # a method of the same name (a constructor helper, say), which is
        # not a transform component.
        rotation = getattr(frame, rotation_name, None)
        translation = getattr(frame, translation_name, None)
        if callable(rotation) or callable(translation):
            return None, None
        return rotation, translation

    rotation, translation = _attribute_pair("rotation", "translation")
    if rotation is None and translation is None:
        rotation, translation = _attribute_pair("R", "t")
    if rotation is not None or translation is not None:
        return FrameTransform(rotation, translation)
    raise ConfigurationError(
        "Cannot interpret frame %r; expected 'geometry', 'detector', a "
        "FrameTransform, a (rotation, translation) pair, or an object "
        "with rotation/translation (or R/t) attributes" % (frame,))


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
        E          : float array, parent energy at production [GeV]
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


# Retain the beam simulation masses used by existing dk2nu conversions.
_PARENT_MASSES = {
    111: 0.1349768,
    211: 0.13957039, -211: 0.13957039,
    321: 0.49368,    -321: 0.49368,
    130: 0.49761,
    13: 0.10566,     -13: 0.10566,
}


def _parent_masses(ptype, masses):
    """Resolve selected PDGs without changing conventional beam masses."""
    known = dict(_PARENT_MASSES)
    if masses is not None:
        if not isinstance(masses, Mapping):
            raise ConfigurationError("masses must map integer PDG codes to masses in GeV")
        for pdg, value in masses.items():
            if isinstance(pdg, bool) or not isinstance(pdg, Integral):
                raise ConfigurationError("masses keys must be integer PDG codes")
            try:
                value = float(value)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ConfigurationError(f"Invalid mass for PDG {pdg}") from exc
            if not math.isfinite(value) or value < 0:
                raise ConfigurationError(f"Mass for PDG {pdg} must be finite and non-negative")
            known[int(pdg)] = value
    result = np.empty(len(ptype), dtype=float)
    for value in np.unique(ptype):
        if not np.isfinite(value) or value != int(value):
            raise ConfigurationError(f"Parent PDG must be a finite integer, got {value!r}")
        pdg = int(value)
        if pdg not in known:
            try:
                known[pdg] = particles.mass(ParticleType(pdg))
            except (ValueError, TypeError, OverflowError, ConfigurationError) as exc:
                raise ConfigurationError(
                    f"No mass for parent PDG {pdg}; supply masses={{PDG: mass_in_GeV}} "
                    "or register it with particles.define().") from exc
        result[ptype == value] = known[pdg]
    return result


def _parent_mask(data, parent_pdg):
    """Apply the same species and importance-weight selection to both exports."""
    ptype = np.asarray(data["ptype"])
    mask = (np.ones(len(ptype), dtype=bool) if parent_pdg is None
            else np.isin(ptype, parent_pdg))
    weights = np.asarray(data["nimpwt"])
    bad = ~np.isfinite(weights) | (weights < 0)
    if np.any(bad & mask):
        print("  dropping %d row(s) with negative or non-finite nimpwt "
              "(importance-reweighting bookkeeping artifacts)"
              % int(np.sum(bad & mask)))
        mask &= ~bad
    return mask


def _extra_columns(columns, mask):
    """Validate metadata before the native table converts it to doubles."""
    if columns is None:
        return {}
    if not isinstance(columns, Mapping):
        raise ConfigurationError("extra_columns must map column names to full-length arrays")
    reserved = {"E", "px", "py", "pz", "x", "y", "z", "x0", "y0", "z0",
                "m", "weight", "nimpwt", "t0", "pol_x", "pol_y", "pol_z"}
    result = {}
    for name, values in columns.items():
        if (not isinstance(name, str) or not name or name.strip() != name
                or any(c in name for c in ",\r\n") or name in reserved
                or name.startswith("PrimaryExternalDistribution_")):
            raise ConfigurationError(f"Invalid or reserved extra column name {name!r}")
        supplied = values
        values = np.asarray(values)
        if values.shape != mask.shape or values.dtype.kind not in "biuf":
            raise ConfigurationError(f"extra_columns[{name!r}] must be a real numeric array "
                                     f"of shape {mask.shape}")
        if (np.ma.isMaskedArray(supplied)
                and np.any(np.ma.getmaskarray(supplied)[mask])):
            raise ConfigurationError(f"extra_columns[{name!r}] contains masked values "
                                     "in selected rows")
        values = values[mask]
        # Check integer arrays regardless of the original container. Compare
        # bounds directly so abs(int64_min) cannot overflow.
        invalid_integer = values.dtype.kind in "iu" and (
            np.any(values > 2**53)
            or (values.dtype.kind == "i" and np.any(values < -(2**53))))
        # A sequence mixing signed/unsigned integers or integers/floats can
        # become float64 in np.asarray before its integer IDs are checked.
        # Inspect the supplied elements without that lossy promotion, including
        # zero-dimensional arrays used as scalar elements.
        if not isinstance(supplied, np.ndarray):
            originals = np.asarray(supplied, dtype=object)[mask]
            for value in originals:
                if isinstance(value, np.ndarray) and value.ndim == 0:
                    value = value.item()
                if isinstance(value, Integral) and abs(int(value)) > 2**53:
                    invalid_integer = True
                    break
        if invalid_integer:
            raise ConfigurationError(f"Integer extra column {name!r} exceeds the exact "
                                     "double-precision ID range [-2**53, 2**53]")
        values = values.astype(float)
        if not np.all(np.isfinite(values)):
            raise ConfigurationError(f"extra_columns[{name!r}] must contain finite values")
        result[name] = values
    return result


def dk2nu_to_primary_distribution(
    dk2nu_data,
    detector_model,
    parent_pdg=None,
    sampling_bias=None,
    frame=None,
    *,
    masses=None,
    extra_columns=None,
):
    """
    Build a PrimaryExternalDistribution directly from dk2nu data.

    Converts positions from the file's coordinate system to detector-local
    coordinates and from centimeters to meters.  No intermediate CSV file
    is needed.  By default the file's frame is taken to be the detector
    model's geometry frame (true for G4BNB files with the SBN model, whose
    geometry frame is the BNB beam frame); files written in any other
    beamline's coordinates supply the rigid map through ``frame``.

    Each row carries its physical weight nimpwt/POT in the table's
    ``weight`` column.  The recorded parent list is treated as the
    intended injection ensemble -- a production's importance reweighting
    deliberately oversamples the phase space it cares about -- so rows
    are sampled UNIFORMLY unless an explicit ``sampling_bias`` reshapes
    the selection.  The importance weights encode how to return to the
    physical distribution given that sampling: the distribution reports
    the physical row density on the physical side of the weight ratio
    and the weight total (the recorded parents per POT) as its physical
    normalization, so event weights from the standard shared
    injection/physical assembly come out in events per POT, with the
    per-row nimpwt spread carried by the weights.  Importance-reweighted
    files (nimpwt spanning orders of magnitude) and unweighted files
    (nimpwt identically 1) then yield consistent rates and spectra.

    Parameters
    ----------
    dk2nu_data : dict
        Output of read_dk2nu(), or equivalent plain, unmasked NumPy arrays.
        Masked raw kinematics and weights are unsupported; only extra_columns
        receives explicit selected-row mask validation.
    detector_model : siren.detector.DetectorModel
        Detector model (provides the geometry-to-detector transform).
        May be None when ``frame`` targets detector coordinates directly.
    parent_pdg : int or list of int, optional
        Filter to specific parent PDG code(s).
    sampling_bias : callable or array-like, optional
        Additional biasing of the ROW SELECTION, independent of the
        importance weights.  A callable f(E, px, py, pz, vx, vy, vz)
        computes per-row selection weights from the parent kinematics in
        the file's own coordinate system (before any frame transform);
        arguments are numpy arrays and the return value should broadcast
        to the same length.  An array-like must align with the arrays in
        ``dk2nu_data`` (full length; it is filtered alongside them), so
        ``sampling_bias=dk2nu_data["nimpwt"]`` selects rows proportional
        to their importance weights (the near-equal-event-weight
        configuration).  Rows are selected with probability proportional
        to the result; the generation density accounts for the bias while
        the physical row density stays proportional to nimpwt, so event
        weights remain correct even when one instance is shared between
        the injection and physical sides.  When None (default), rows are
        selected uniformly. Biases must be finite and non-negative and
        retain positive support for every row with positive physical weight.
        Callable E is the stored production energy for compatibility.
    frame : optional
        Coordinate system of the dk2nu file.  None or "geometry": the
        file frame is the detector model's geometry frame.  "detector":
        the file frame is the detector-local frame (the model is not
        consulted).  A FrameTransform, a (rotation, translation) pair, or
        any object with rotation/translation (or R/t) attributes: the
        rigid map from the file frame into the geometry frame, in meters
        (a FrameTransform may instead target the detector frame).  For
        transforms that are not rigid, transform the arrays in
        ``dk2nu_data`` before calling.

    masses : mapping of int to float, optional
        PDG-to-mass overrides in GeV. Otherwise use the conventional beam
        masses, then particles.mass (including particles.define registrations).
        Unknown masses raise ConfigurationError.
    extra_columns : mapping of str to array-like, optional
        Numeric metadata aligned with the full input rows, filtered alongside
        them and copied by name into interaction_parameters. Reserved native
        columns cannot be overridden. Integer IDs must lie within [-2**53, 2**53].
        Inclusion probabilities are metadata only: nimpwt must already include
        any source-subsampling correction. No additional weight is applied.

    Returns
    -------
    siren.distributions.PrimaryExternalDistribution
    """
    xform = _as_frame_transform(frame)
    if xform.target == "geometry" and detector_model is None:
        raise ConfigurationError(
            "detector_model is required to convert geometry coordinates "
            "into detector coordinates. Pass the model, or a frame that "
            "targets 'detector' if the file coordinates are already "
            "detector-local.")

    ptype = np.asarray(dk2nu_data["ptype"])
    mask = _parent_mask(dk2nu_data, parent_pdg)

    simulated_pot = dk2nu_data["pot"]
    # Per-POT weights are meaningless without a positive POT. read_dk2nu leaves
    # pot at 0.0 when a file has no dkmetaTree/pots branch; dividing by it would
    # emit inf/nan weights silently. Fail loud at the point the weights are
    # formed rather than propagate a corrupt distribution.
    if not (math.isfinite(simulated_pot) and simulated_pot > 0):
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
    # Normalize in binary64 even when the beam input stores float32 weights.
    nimpwt = np.asarray(dk2nu_data["nimpwt"], dtype=float)[mask]
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

    row_masses = _parent_masses(pt, masses)
    metadata = _extra_columns(extra_columns, mask)

    # Map the file's coordinates into the frame the transform targets:
    # positions (after cm -> m) pick up the rotation and translation,
    # momenta and spin axes the rotation alone. The default transform is
    # the identity into the geometry frame.
    pos_m = np.stack([vx, vy, vz], axis=1) * 0.01
    mom = np.stack([px, py, pz], axis=1)
    pos_m = pos_m @ xform.rotation.T + xform.translation
    mom = mom @ xform.rotation.T
    if pol is not None:
        pol = pol @ xform.rotation.T
    to_detector = xform.target == "detector"

    # The decay time from the ancestor chain becomes the primary's initial
    # time (the t0 column of PrimaryExternalDistribution), so interaction
    # times downstream are physical.
    keys = ["E", "px", "py", "pz", "x", "y", "z", "m", "weight"]
    if t0 is not None:
        keys.append("t0")
    if pol is not None:
        keys += ["pol_x", "pol_y", "pol_z"]
    keys += list(metadata)
    data = []
    for i in range(len(E)):
        # Convert position from geometry to detector coordinates, unless
        # the frame transform already landed in detector coordinates.
        if to_detector:
            x_det = (float(pos_m[i][0]), float(pos_m[i][1]),
                     float(pos_m[i][2]))
        else:
            geo_pos = GeometryPosition(Vector3D(*pos_m[i]))
            det_pos = detector_model.GeoPositionToDetPosition(geo_pos).get()
            x_det = (det_pos.GetX(), det_pos.GetY(), det_pos.GetZ())

        # Convert momentum direction from geometry to detector coordinates.
        # Energy is a scalar and is unchanged; the 3-momentum direction
        # must be rotated if the detector axes differ from geometry axes.
        p_mag = math.sqrt(float(mom[i][0])**2 + float(mom[i][1])**2
                          + float(mom[i][2])**2)
        if to_detector or p_mag == 0.0:
            px_det, py_det, pz_det = (float(mom[i][0]), float(mom[i][1]),
                                      float(mom[i][2]))
        else:
            geo_dir = GeometryDirection(Vector3D(
                float(mom[i][0]) / p_mag, float(mom[i][1]) / p_mag,
                float(mom[i][2]) / p_mag))
            det_dir = detector_model.GeoDirectionToDetDirection(geo_dir).get()
            px_det = det_dir.GetX() * p_mag
            py_det = det_dir.GetY() * p_mag
            pz_det = det_dir.GetZ() * p_mag

        m = row_masses[i]
        # dk2nu stores momentum at decay (pdpx/pdpy/pdpz) but energy
        # at production (ppenergy). Compute on-shell energy from the
        # decay-point momentum and known mass.
        E_decay = math.sqrt(p_mag * p_mag + m * m)
        row = [
            E_decay, px_det, py_det, pz_det,
            x_det[0], x_det[1], x_det[2],
            m, float(weight[i]),
        ]
        if t0 is not None:
            row.append(float(t0[i]))
        if pol is not None:
            # The spin axis is a rest-frame direction expressed in the
            # file's spatial basis; rotate it like the momentum (all these
            # frames are related by rotations and pure boosts, which share
            # their spatial bases).
            norm = float(np.linalg.norm(pol[i]))
            if to_detector or norm == 0.0:
                row += [float(pol[i][0]), float(pol[i][1]), float(pol[i][2])]
            else:
                geo_ax = GeometryDirection(Vector3D(*(pol[i] / norm)))
                det_ax = detector_model.GeoDirectionToDetDirection(geo_ax).get()
                row += [det_ax.GetX() * norm, det_ax.GetY() * norm,
                        det_ax.GetZ() * norm]
        row += [float(values[i]) for values in metadata.values()]
        data.append(row)

    if sampling_bias is not None:
        # The bias shapes the row selection only; the physical row weights
        # (nimpwt/POT, already in the table) de-bias it through the
        # distribution's physical density, so no composition happens here.
        if callable(sampling_bias):
            sw = np.asarray(
                sampling_bias(E, px, py, pz, vx, vy, vz),
                dtype=float,
            )
            if sw.shape != E.shape:
                sw = np.broadcast_to(sw, E.shape).astype(float).copy()
        else:
            sw = np.asarray(sampling_bias, dtype=float)
            if sw.shape != np.shape(dk2nu_data["E"]):
                raise ConfigurationError(
                    "array-like sampling_bias must align with the arrays in "
                    "dk2nu_data (shape %r), got %r"
                    % (np.shape(dk2nu_data["E"]), sw.shape))
            sw = sw[mask]
        sw = np.array(sw, dtype=float)
        if not np.all(np.isfinite(sw)) or np.any(sw < 0):
            raise ConfigurationError("sampling_bias must be finite and non-negative")
        with np.errstate(over="ignore"):
            total_bias = float(np.sum(sw))
        if not math.isfinite(total_bias) or total_bias <= 0:
            raise ConfigurationError("sampling_bias must have a finite positive sum")
        if np.any((weight > 0) & (sw == 0)):
            raise ConfigurationError("sampling_bias has missing physical support")
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
    *,
    masses=None,
    extra_columns=None,
):
    """
    Write dk2nu parent meson kinematics to a CSV file suitable for
    SIREN's PrimaryExternalDistribution.

    Energies use the decay-point momentum and resolved mass. Columns include
    E, px, py, pz, x, y, z, m, nimpwt, and weight=nimpwt/POT, plus available
    time, polarization, and explicit metadata. Use units_cm=False when loading
    the file directly into PrimaryExternalDistribution, which expects metres.

    Parameters
    ----------
    dk2nu_data : dict
        Output of read_dk2nu(), or equivalent plain, unmasked NumPy arrays.
        Masked raw kinematics and weights are unsupported; only extra_columns
        receives explicit selected-row mask validation.
    output_path : str or os.PathLike
        Filesystem path to write plain UTF-8 CSV with LF newlines, regardless
        of the filename suffix. File-like objects are unsupported.
    parent_pdg : int or list of int, optional
        Filter to specific parent PDG code(s).  Default: use all entries
        in dk2nu_data (which may already be filtered).
    position_transform : callable, optional
        Function that maps (vx, vy, vz) position arrays, in centimeters.
        It does not rotate momenta or polarization; transform those input
        arrays as well when changing the spatial basis.
        dk2nu positions are in cm.  If None, positions are used as-is.
    units_cm : bool
        If True (default), positions in the CSV are in cm.
        If False, positions are converted to meters.

    masses : mapping of int to float, optional
        PDG-to-mass overrides in GeV. Otherwise use the conventional beam
        masses, then particles.mass (including particles.define registrations).
        Unknown masses raise ConfigurationError.
    extra_columns : mapping of str to array-like, optional
        Numeric metadata aligned with the full input rows, filtered alongside
        them and copied by name into interaction_parameters. Reserved native
        columns cannot be overridden. Integer IDs must lie within [-2**53, 2**53].
        Inclusion probabilities are metadata only: nimpwt must already include
        any source-subsampling correction. No additional weight is applied.

    Returns
    -------
    int
        Number of rows written.
    """
    ptype = np.asarray(dk2nu_data["ptype"])
    mask = _parent_mask(dk2nu_data, parent_pdg)
    simulated_pot = dk2nu_data["pot"]
    if not (math.isfinite(simulated_pot) and simulated_pot > 0):
        raise ConfigurationError("A finite positive simulated POT is required for CSV weights")

    px = np.asarray(dk2nu_data["px"], dtype=float)[mask]
    py = np.asarray(dk2nu_data["py"], dtype=float)[mask]
    pz = np.asarray(dk2nu_data["pz"], dtype=float)[mask]
    vx = dk2nu_data["vx"][mask]
    vy = dk2nu_data["vy"][mask]
    vz = dk2nu_data["vz"][mask]
    # Promote before POT division so small positive weights do not underflow
    # in the input array's precision.
    nimpwt = np.asarray(dk2nu_data["nimpwt"], dtype=float)[mask]
    pt = ptype[mask]

    if position_transform is not None:
        vx, vy, vz = position_transform(vx, vy, vz)

    scale = 1.0 if units_cm else 0.01

    row_masses = _parent_masses(pt, masses)
    metadata = _extra_columns(extra_columns, mask)
    # ppenergy describes production, while px/py/pz describe the decay point.
    energy = np.sqrt(px**2 + py**2 + pz**2 + row_masses**2)
    keys = ["E", "px", "py", "pz", "x", "y", "z", "m", "nimpwt", "weight"]
    columns = [energy, px, py, pz, vx * scale, vy * scale, vz * scale,
               row_masses, nimpwt, nimpwt / simulated_pot]
    for name in ("t0", "pol_x", "pol_y", "pol_z"):
        if name in dk2nu_data:
            keys.append(name)
            columns.append(dk2nu_data[name][mask])
    keys += list(metadata)
    columns += list(metadata.values())
    # Round-trip all binary64 values, including large but exactly representable IDs.
    # The native loader reads plain text. Passing a filename to savetxt would
    # implicitly compress .gz/.bz2 paths and break the native round trip.
    table = np.column_stack(columns)
    with open(output_path, "w", encoding="utf-8", newline="\n") as output:
        np.savetxt(output, table, delimiter=",",
                   header=",".join(keys), comments="", fmt="%.17g")

    return int(np.sum(mask))


def print_summary(dk2nu_data):
    """Print a summary of the dk2nu data."""
    ptypes = dk2nu_data["ptype"]
    unique, counts = np.unique(ptypes, return_counts=True)
    n_total = len(ptypes)
    pot = dk2nu_data["pot"]

    print(f"Total entries: {n_total}")
    print(f"Total POT: {pot:.3e}")
    print("Parent breakdown:")
    for pdg, count in sorted(zip(unique, counts), key=lambda x: -x[1]):
        name = _PARENT_NAMES.get(int(pdg), str(int(pdg)))
        frac = 100.0 * count / n_total if n_total > 0 else 0
        E = dk2nu_data["E"][ptypes == pdg]
        print(f"  {name:>5s} ({int(pdg):>4d}): {count:>7d} ({frac:5.1f}%)  "
              f"E = [{E.min():.3f}, {E.max():.3f}] GeV")
