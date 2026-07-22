"""Merge split decay-in-flight / decay-at-rest dk2nu productions.

Analysis-level helper for the example4 beam-timing scripts.

The standard g4numi production kills every particle other than a neutrino
once its kinetic energy falls below KillTrackingThreshold (0.05 GeV), so
its files contain no decays from parents below that energy.  A dedicated
decay-at-rest production tracks everything to rest and records decays at
all parent energies -- the two samples overlap above the threshold, and
they carry different protons-on-target (POT).  The merge implemented here
makes the split exact and the normalization consistent:

- decay-in-flight rows are kept where the parent kinetic energy is at or
  above the cut;
- decay-at-rest rows are kept where it is below the cut;
- decay-at-rest importance weights are rescaled by the POT ratio, so the
  merged rows behave exactly like a single sample whose POT is the
  decay-in-flight POT.

Downstream consumers that form per-POT weights as nimpwt / pot (for
example siren.dk2nu.dk2nu_to_primary_distribution) therefore need no
changes.
"""

import numpy as np

from siren import dk2nu

# Parent rest masses [GeV] for kinetic-energy classification.  Distinct
# from the mass map inside dk2nu_to_primary_distribution, which feeds
# injected kinematics.
PARENT_MASSES = {
    13: 0.1056583755,     # mu
    211: 0.13957039,      # pi+-
    321: 0.493677,        # K+-
    130: 0.497611,        # K0L
    310: 0.497611,        # K0S
    2112: 0.9395654205,   # n
    2212: 0.9382720882,   # p
}


def parent_kinetic_energy(dk2nu_data):
    """Parent kinetic energy at decay [GeV] for each row.

    Computed as E - m(ptype) from the recorded parent energy.  Raises
    KeyError for a parent species without a rest mass on record.
    """
    codes = np.abs(np.asarray(dk2nu_data["ptype"], dtype=np.int64))
    masses = np.empty(codes.shape, dtype=float)
    for code in np.unique(codes):
        if int(code) not in PARENT_MASSES:
            raise KeyError(
                "No rest mass on record for dk2nu parent pdg %d; extend "
                "PARENT_MASSES to classify it" % int(code))
        masses[codes == code] = PARENT_MASSES[int(code)]
    return dk2nu_data["E"] - masses


def combine_dif_dar(dif_data, dar_data, kinetic_energy_cut=0.05):
    """Merge a decay-in-flight and a decay-at-rest dk2nu sample.

    Parameters
    ----------
    dif_data, dar_data : dict
        Outputs of siren.dk2nu.read_dk2nu() for the decay-in-flight and
        decay-at-rest file sets, read with the same filters.
    kinetic_energy_cut : float
        Parent kinetic energy boundary between the samples [GeV].  Use
        the production's kill threshold (default 0.05).

    Returns
    -------
    dict with the intersection of the two samples' array keys, plus:
        dar     : bool array, True for rows from the decay-at-rest sample
        pot     : the decay-in-flight POT (the merged set's normalization)
        dif_pot, dar_pot : the raw per-sample POT
    """
    dif_pot = float(dif_data.get("pot", 0.0))
    dar_pot = float(dar_data.get("pot", 0.0))
    if not (dif_pot > 0.0) or not (dar_pot > 0.0):
        raise ValueError(
            "combine_dif_dar requires a positive POT in both samples "
            "(got dif_pot=%r, dar_pot=%r); the input files carried no "
            "usable dkmetaTree/pots metadata" % (dif_pot, dar_pot))

    dif_keep = parent_kinetic_energy(dif_data) >= kinetic_energy_cut
    dar_keep = parent_kinetic_energy(dar_data) < kinetic_energy_cut

    n_below = int(np.count_nonzero(~dif_keep))
    if n_below > 0:
        print(
            "combine_dif_dar: dropped %d decay-in-flight rows below the "
            "%.3g GeV cut (the decay-at-rest sample covers them)"
            % (n_below, kinetic_energy_cut))

    def _array_keys(data):
        return {k for k, v in data.items()
                if isinstance(v, np.ndarray) and k != "pot"}

    dif_keys = _array_keys(dif_data)
    dar_keys = _array_keys(dar_data)
    for key in sorted(dif_keys ^ dar_keys):
        print(
            "combine_dif_dar: column '%s' is present in only one sample "
            "and is dropped from the merge" % key)

    merged = {}
    for key in sorted(dif_keys & dar_keys):
        merged[key] = np.concatenate(
            [dif_data[key][dif_keep], dar_data[key][dar_keep]])
    merged["nimpwt"] = np.concatenate(
        [dif_data["nimpwt"][dif_keep],
         dar_data["nimpwt"][dar_keep] * (dif_pot / dar_pot)])
    merged["dar"] = np.concatenate(
        [np.zeros(int(np.count_nonzero(dif_keep)), dtype=bool),
         np.ones(int(np.count_nonzero(dar_keep)), dtype=bool)])
    merged["pot"] = dif_pot
    merged["dif_pot"] = dif_pot
    merged["dar_pot"] = dar_pot
    return merged


def read_dif_dar(
    dif_filenames,
    dar_filenames,
    kinetic_energy_cut=0.05,
    **read_kwargs,
):
    """Read matched decay-in-flight and decay-at-rest dk2nu productions.

    Both file sets are read with siren.dk2nu.read_dk2nu using the same
    filters, then merged with combine_dif_dar.
    """
    dif_data = dk2nu.read_dk2nu(dif_filenames, **read_kwargs)
    dar_data = dk2nu.read_dk2nu(dar_filenames, **read_kwargs)
    return combine_dif_dar(
        dif_data, dar_data, kinetic_energy_cut=kinetic_energy_cut)
