r"""G4BNB/G4NuMI dk2nu -> BeamDecay neutrinos -> DarkNews HNL timing in SBN.

This example keeps the dk2nu parent meson as the SIREN primary.  The first
vertex uses the BeamDecays resource to generate

    pi+ -> mu+ nu_mu,

directing the neutrino toward the selected detector with a physical fallback.
The neutrino then upscatters through a DarkNews dipole-portal model and the N4
decays to nu + gamma.  The complete event chain is therefore

    dk2nu pi+ -> mu+ nu_mu -> N4 + Ar -> nu + gamma.

No intermediate neutrino CSV or flux histogram is required.  The pion decay
time from the dk2nu ancestor chain becomes the initial event time, and SIREN
propagates the neutrino and HNL times of flight between later vertices.

Every SBN detector model uses the BNB/SAND frame for geometry. G4BNB files
therefore use the default ``--beam-frame BNB``. For G4NuMI input, pass
``--beam-frame NuMI`` to apply the surveyed NuMI-to-BNB transform carried by
the SBN geometry resource.

By default the DarkNews interpolation tables are filled before injection,
up to the forward Doppler bound of the loaded dk2nu rows, and saved next to
the DarkNewsTables resource. No DarkNews cross-section evaluation then
happens inside the injection loop, and later runs load the tables from disk
instead of recomputing them. Use ``--table-emax`` to override the fill
range, or ``--no-precompute-tables`` to build the tables lazily during
injection.

Example::

    python DarkNewsHNL_SBN_dk2nu_timing.py \
        '/data/g4bnb/*dk2nu*.root' --detector SBND --events 500 \
        --m4 0.10 --mu-tr-mu4 2.5e-6 \
        --output output/sbnd_darknews_hnl_timing.png

    python DarkNewsHNL_SBN_dk2nu_timing.py \
        '/data/g4numi/*dk2nu*.root' --beam-frame NuMI \
        --detector ICARUS --events 500 \
        --output output/icarus_numi_darknews_hnl_timing.png
"""

import argparse
import glob
import os

import numpy as np

import siren
from siren import _util, channels, dk2nu, expand
from siren.Injector import Injector
from siren.Weighter import Weighter


_FIDUCIALS = {
    "SBND": {
        "center": (0.0, 0.59, -0.415),
        "widths": (4.026, 4.074645, 5.01),
    },
    "ICARUS": {
        "center": (0.0, 0.0, 0.0),
        "widths": (7.20, 3.16, 17.95),
    },
}

# G4BNB geometry coordinates, in meters.  The SBN detector loader performs the
# conversion into each detector's local coordinates.
BNB_TARGET_CENTER = (0.0, 0.0, 0.3905)


def _prompt_origin_geometry(beam_frame, frame_transform):
    """Reference point for the beta=1 timing subtraction, in BNB geometry."""
    if beam_frame == "BNB":
        return np.asarray(BNB_TARGET_CENTER, dtype=float)
    # G4NuMI positions and ancestor times use MCZERO, the NuMI frame origin.
    return np.asarray(frame_transform.apply([0.0, 0.0, 0.0]), dtype=float)


def _expand_input_paths(values):
    paths = []
    for value in values:
        matches = sorted(glob.glob(value))
        if matches:
            paths.extend(matches)
        elif os.path.isfile(value):
            paths.append(value)
    return list(dict.fromkeys(os.path.abspath(path) for path in paths))


def _load_beam_decay_model():
    """Load the module-only BeamDecays resource."""
    module = siren.resources.processes.BeamDecays
    return module.MesonTwoBodyLeptonicDecay(211)


def _load_beam_samples():
    """Load the sibling decay-in-flight / decay-at-rest merging helper."""
    path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "_beam_samples.py")
    return _util.load_module("example4_beam_samples", path)


def _max_neutrino_energy(data):
    """Forward Doppler bound on the decay neutrino energy over all rows."""
    decay = _load_beam_decay_model()
    m_meson = decay.m_meson
    e_cm = (m_meson ** 2 - decay.m_lepton ** 2) / (2.0 * m_meson)
    energy = np.asarray(data["E"], dtype=float)
    momentum = np.sqrt(np.maximum(energy ** 2 - m_meson ** 2, 0.0))
    return float(np.max((energy + momentum) / m_meson) * e_cm)


def _make_threshold_bias(m4):
    """Sampling bias that skips rows unable to reach the N4 threshold.

    The forward Doppler bound (E + p) / m * E_cm is the largest neutrino
    energy a recorded parent can produce.  A row whose bound does not
    exceed the N4 mass can never drive the upscatter, so it gets zero
    sampling weight: its physical contribution is exactly zero, and the
    generation probability accounts for the remaining bias.  Without
    this, a merged decay-at-rest sample dominates the sampling with
    stopped pions whose 30 MeV neutrinos are far below threshold.
    """
    decay = _load_beam_decay_model()
    m_meson = decay.m_meson
    e_cm = (m_meson ** 2 - decay.m_lepton ** 2) / (2.0 * m_meson)

    def bias(E, px, py, pz, vx, vy, vz):
        energy = np.asarray(E, dtype=float)
        momentum = np.sqrt(np.maximum(energy ** 2 - m_meson ** 2, 0.0))
        bound = (energy + momentum) / m_meson * e_cm
        return (bound > m4).astype(float)

    return bias


def _darknews_bundle(detector_model, m4, mu_tr_mu4,
                     table_emax=None, save_tables=True):
    model_kwargs = {
        "m4": m4,
        "mu_tr_mu4": mu_tr_mu4,
        "UD4": 0,
        "Umu4": 0,
        "epsilon": 0.0,
        "gD": 0.0,
        "decay_product": "photon",
        "noHC": True,
        "HNLtype": "dirac",
    }
    table_name = "DarkNewsTables-v%s/" % siren.utilities.darknews_version()
    table_name += "Dipole_M%2.2e_mu%2.2e" % (m4, mu_tr_mu4)
    bundle = siren.load_processes(
        "DarkNewsTables",
        primary_type=siren.particles.NuMu,
        detector_model=detector_model,
        # The Ar40-only model gives zero upscatter density in other sectors.
        # Add selected upstream nuclei here when production outside the active
        # liquid argon is part of the study.
        nuclear_targets=["Ar40"],
        table_name=table_name,
        # Fill the interpolation tables over the full beam energy range up
        # front, so no DarkNews cross-section evaluation happens mid-injection.
        fill_tables_at_start=table_emax is not None,
        Emax=table_emax,
        **model_kwargs,
    )
    # The N4 total width is computed lazily on first use; trigger it here so
    # the width integral is also paid before injection starts.
    for particle_type, decays in bundle.secondary.items():
        for decay in decays:
            decay.TotalDecayWidthAllFinalStates(particle_type)
    if table_emax is not None and save_tables:
        # Persist the filled tables next to the resource; later runs with
        # the same model parameters load them from disk.
        darknews_tables = siren.resources.processes.DarkNewsTables
        table_dir = os.path.join(
            _util.resource_package_dir(), "processes", "DarkNewsTables",
            table_name)
        primary_ups_keys, secondary_dec_keys = bundle.metadata
        darknews_tables.SaveDarkNewsProcesses(
            table_dir,
            bundle.primary, primary_ups_keys,
            bundle.secondary, secondary_dec_keys)
    return bundle


def build_vertices(detector_model, external, fiducial, m4, mu_tr_mu4,
                   table_emax=None, save_tables=True):
    """Build pi -> nu, nu -> N4, and N4 -> nu gamma vertices."""
    pion_decay = _load_beam_decay_model()
    bundle = _darknews_bundle(detector_model, m4, mu_tr_mu4,
                              table_emax=table_emax, save_tables=save_tables)
    hnl_models = bundle.secondary[siren.particles.N4]
    hnl_interactions = siren.interactions.InteractionCollection(
        siren.particles.N4, hnl_models)

    pion = siren.Vertex(
        # The dk2nu parent is a pi+.
        siren.dataclasses.ParticleType(211),
        # The only interaction is the two-body leptonic decay to mu+ nu_mu.
        pion_decay,
        distributions=[external],
        physical=[external],
        # Fixed vertex weighting is appropriate because dk2nu already records
        # the pion decay position.
        weighting=siren.Fixed(),
        kinematics=(
            # Primarily bias the neutrino towards the fiducial volume
            0.99 * channels.toward("NuMu", fiducial)
            # Retain a small fraction of physically distributed decays.
            + 0.01 * channels.physical()
        ),
        # The neutrino is the only secondary we expand into another vertex.
        expand=(expand.child("NuMu"),),
    )

    neutrino = siren.Vertex(
        # This vertex is the neutrino upscatter to N4.
        "NuMu",
        # Load the DarkNews dipole-portal model from its precomputed tables.
        bundle.primary[siren.particles.NuMu],
        # Bias the neutrino upscatter by the chance that a collinear N4 with
        # approximately the neutrino energy subsequently decays in the active
        # volume. Both legs use interaction depth, so this remains valid for
        # nonuniform material profiles. The Weighter removes the full bias.
        position=siren.dist.DecayRangeVertex(
            fiducial, hnl_interactions, m4,
            daughter_energy_fraction=1.0,
            max_length=1000.0),
        expand=(expand.child("N4"),),
    )

    hnl = siren.Vertex(
        "N4",
        hnl_models,
        # Bias the N4 decay to the active volume.  The Weighter removes this
        # position bias using the physical DarkNews decay width.
        position=siren.dist.BoundedVertex(fiducial, np.inf),
        expand=(expand.depth_below(0),),
    )
    return pion, (neutrino, hnl)


def _detector_position(detector_model, geometry_position):
    point = siren.detector.GeometryPosition(
        siren.math.Vector3D(*geometry_position))
    point = detector_model.GeoPositionToDetPosition(point).get()
    return np.array([point.GetX(), point.GetY(), point.GetZ()], dtype=float)


def _record_for_primary(event, particle_type):
    wanted = int(particle_type)
    for datum in event.tree:
        if int(datum.record.signature.primary_type) == wanted:
            return datum.record
    return None


def collect_timing(results, detector_model,
                   prompt_origin_geometry=BNB_TARGET_CENTER):
    """Return plot-ready arrays for the pion, upscatter, and N4 decay."""
    target = _detector_position(detector_model, prompt_origin_geometry)
    c_m_per_ns = float(siren.utilities.Constants.c)
    columns = {
        "weight": [],
        "parent_energy_GeV": [],
        "neutrino_energy_GeV": [],
        "hnl_energy_GeV": [],
        "parent_decay_ns": [],
        "upscatter_ns": [],
        "hnl_decay_ns": [],
        "upscatter_delay_ns": [],
        "hnl_decay_delay_ns": [],
        "hnl_flight_ns": [],
    }

    for event, weight in results:
        root = event.tree[0].record
        upscatter = _record_for_primary(event, siren.particles.NuMu)
        hnl_decay = _record_for_primary(event, siren.particles.N4)

        nu_energy = np.nan
        hnl_energy = np.nan
        upscatter_time = np.nan
        hnl_decay_time = np.nan
        upscatter_delay = np.nan
        hnl_decay_delay = np.nan
        hnl_flight = np.nan

        if upscatter is not None:
            nu_energy = float(upscatter.primary_momentum[0])
            upscatter_time = float(upscatter.interaction_time)
            distance = np.linalg.norm(
                np.asarray(upscatter.interaction_vertex, dtype=float) - target)
            upscatter_delay = upscatter_time - distance / c_m_per_ns

            for i, ptype in enumerate(upscatter.signature.secondary_types):
                if ptype == siren.particles.N4:
                    hnl_energy = float(upscatter.secondary_momenta[i][0])
                    break

        if hnl_decay is not None:
            hnl_energy = float(hnl_decay.primary_momentum[0])
            hnl_decay_time = float(hnl_decay.interaction_time)
            distance = np.linalg.norm(
                np.asarray(hnl_decay.interaction_vertex, dtype=float) - target)
            hnl_decay_delay = hnl_decay_time - distance / c_m_per_ns

        if upscatter is not None and hnl_decay is not None:
            hnl_flight = hnl_decay_time - upscatter_time

        columns["weight"].append(float(weight))
        columns["parent_energy_GeV"].append(float(root.primary_momentum[0]))
        columns["neutrino_energy_GeV"].append(nu_energy)
        columns["hnl_energy_GeV"].append(hnl_energy)
        columns["parent_decay_ns"].append(float(root.interaction_time))
        columns["upscatter_ns"].append(upscatter_time)
        columns["hnl_decay_ns"].append(hnl_decay_time)
        columns["upscatter_delay_ns"].append(upscatter_delay)
        columns["hnl_decay_delay_ns"].append(hnl_decay_delay)
        columns["hnl_flight_ns"].append(hnl_flight)

    return {name: np.asarray(values, dtype=float)
            for name, values in columns.items()}


def _weighted_hist(ax, values, weights, *, label=None, bins=50):
    mask = np.isfinite(values) & np.isfinite(weights) & (weights >= 0.0)
    if not np.any(mask) or np.sum(weights[mask]) <= 0.0:
        ax.text(0.5, 0.5, "no completed vertices", ha="center", va="center",
                transform=ax.transAxes)
        return
    ax.hist(values[mask], bins=bins, weights=weights[mask], density=True,
            histtype="step", linewidth=1.8, label=label)


def plot_timing(columns, output_path, detector_name, beam_frame="BNB"):
    import matplotlib.pyplot as plt

    weights = columns["weight"]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), constrained_layout=True)

    _weighted_hist(axes[0, 0], columns["parent_decay_ns"], weights)
    axes[0, 0].set_xlabel("pion decay time after proton [ns]")
    axes[0, 0].set_ylabel("weighted density")

    _weighted_hist(axes[0, 1], columns["upscatter_ns"], weights,
                   label=r"$\nu_\mu\to N_4$")
    _weighted_hist(axes[0, 1], columns["hnl_decay_ns"], weights,
                   label=r"$N_4\to\nu\gamma$")
    axes[0, 1].set_xlabel("absolute vertex time after proton [ns]")
    axes[0, 1].legend()

    _weighted_hist(axes[1, 0], columns["upscatter_delay_ns"], weights,
                   label=r"$\nu_\mu\to N_4$")
    _weighted_hist(axes[1, 0], columns["hnl_decay_delay_ns"], weights,
                   label=r"$N_4\to\nu\gamma$")
    axes[1, 0].set_xlabel(r"delay relative to prompt $\beta=1$ [ns]")
    axes[1, 0].set_ylabel("weighted density")
    axes[1, 0].legend()

    _weighted_hist(axes[1, 1], columns["hnl_flight_ns"], weights)
    axes[1, 1].set_xlabel(r"$N_4$ decay time after upscatter [ns]")

    for ax in axes.flat:
        ax.grid(alpha=0.25)
    fig.suptitle("%s: %s dk2nu DarkNews HNL timing"
                 % (detector_name, beam_frame))

    output_path = os.path.abspath(output_path)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    return output_path


def save_timing_csv(columns, output_path):
    names = list(columns)
    table = np.column_stack([columns[name] for name in names])
    np.savetxt(output_path, table, delimiter=",", header=",".join(names),
               comments="")
    return output_path


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dk2nu_files", nargs="+",
                        help="G4BNB/G4NuMI dk2nu ROOT files or quoted glob(s)")
    parser.add_argument(
        "--beam-frame", choices=("BNB", "NuMI"), default="BNB",
        help="coordinate frame used by the dk2nu files (default: BNB)")
    parser.add_argument("--detector", choices=sorted(_FIDUCIALS),
                        default="SBND")
    parser.add_argument("--events", type=int, default=250)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--entry-stop", type=int,
                        help="read at most this many rows per ROOT file")
    parser.add_argument("--m4", type=float, default=0.10,
                        help="HNL mass in GeV (default: 0.10)")
    parser.add_argument("--mu-tr-mu4", type=float, default=2.5e-6,
                        help="transition dipole coupling in GeV^-1")
    parser.add_argument("--output",
                        help=("PNG path (default: output/<detector>_<beam>_"
                              "darknews_hnl_timing.png)"))
    parser.add_argument(
        "--table-emax", type=float,
        help=("fill the DarkNews interpolation tables up to this neutrino "
              "energy in GeV before injecting (default: the forward Doppler "
              "bound of the loaded dk2nu rows)"))
    parser.add_argument(
        "--no-precompute-tables", action="store_true",
        help=("skip the up-front table fill and save; build the DarkNews "
              "tables lazily during injection"))
    parser.add_argument(
        "--dar-files", nargs="+",
        help=("decay-at-rest dk2nu ROOT files or quoted glob(s), merged "
              "with the (decay-in-flight) positional files: each sample "
              "keeps its side of the kinetic-energy cut and is normalized "
              "by its own POT"))
    parser.add_argument(
        "--dar-ke-cut", type=float, default=0.05,
        help=("parent kinetic energy boundary in GeV between the "
              "decay-in-flight and decay-at-rest samples (default: 0.05, "
              "the g4numi kill threshold)"))
    args = parser.parse_args(argv)

    files = _expand_input_paths(args.dk2nu_files)
    if not files:
        parser.error("none of the dk2nu paths/globs matched a file")

    read_kwargs = dict(
        parent_pdg=dk2nu.PTYPE_PIPLUS,
        decay_modes=13,
        nu_pdg=14,
        read_time=True,
        entry_stop=args.entry_stop,
    )
    if args.dar_files:
        dar_files = _expand_input_paths(args.dar_files)
        if not dar_files:
            parser.error("none of the --dar-files paths/globs matched a file")
        data = _load_beam_samples().read_dif_dar(
            files, dar_files,
            kinetic_energy_cut=args.dar_ke_cut,
            **read_kwargs)
    else:
        data = dk2nu.read_dk2nu(files, **read_kwargs)
    dk2nu.print_summary(data)
    if "t0" not in data:
        raise RuntimeError(
            "These dk2nu files have no ancestor start-time branches. "
            "They can drive event generation, but not this timing example.")
    if len(data["E"]) == 0:
        raise RuntimeError("No pi+ -> mu+ nu_mu rows survived the filters")

    detector_model = siren.load_detector("SBN", detector=args.detector)
    frame = siren.resources.detectors.SBN.geo.transform(args.beam_frame, "BNB")
    bias = _make_threshold_bias(args.m4)
    open_rows = bias(data["E"], data["px"], data["py"], data["pz"],
                     data["vx"], data["vy"], data["vz"])
    n_open = int(np.count_nonzero(open_rows))
    if n_open == 0:
        raise RuntimeError(
            "No dk2nu row can reach the N4 threshold for m4 = %g GeV"
            % args.m4)
    print("Threshold bias keeps %d of %d rows for sampling"
          % (n_open, len(open_rows)))
    external = dk2nu.dk2nu_to_primary_distribution(
        data, detector_model, frame=frame, sampling_bias=bias)
    prompt_origin = _prompt_origin_geometry(args.beam_frame, frame)
    spec = _FIDUCIALS[args.detector]
    fiducial = siren.geometry.Box(
        widths=spec["widths"], center=spec["center"])
    if args.no_precompute_tables:
        table_emax = None
    else:
        table_emax = args.table_emax or _max_neutrino_energy(data)
        print("Precomputing DarkNews tables up to %.3f GeV" % table_emax)
    primary, secondaries = build_vertices(
        detector_model, external, fiducial, args.m4, args.mu_tr_mu4,
        table_emax=table_emax)

    injector = Injector(
        detector=detector_model,
        primary=primary,
        secondaries=secondaries,
        events=args.events,
        seed=args.seed,
    )
    weighter = Weighter(injector, primary_physical=primary.physical)
    results = siren.generate(
        injector, weighter, events=args.events, on_shortfall="warn")
    results.summary()
    print(injector.report())

    columns = collect_timing(
        results, detector_model, prompt_origin_geometry=prompt_origin)
    output = args.output or os.path.join(
        "output", "%s_%s_darknews_hnl_timing.png"
        % (args.detector.lower(), args.beam_frame.lower()))
    output = plot_timing(
        columns, output, args.detector, beam_frame=args.beam_frame)
    csv_path = os.path.splitext(output)[0] + ".csv"
    save_timing_csv(columns, csv_path)
    print("Wrote %s" % output)
    print("Wrote %s" % csv_path)


if __name__ == "__main__":
    main()
