r"""G4BNB/G4NuMI dk2nu -> SBN vector-portal events -> timing plots.

This is a small companion example to VectorPortal_SBND_dk2nu.py.
It uses the same off-shell Dutta-Kim chain, but accepts either the SBND or
ICARUS detector and writes timing plots plus a CSV table of the plotted
quantities.

The important data path is direct and does not require an intermediate CSV::

    dk2nu ROOT files
        -> siren.dk2nu.read_dk2nu(...)
        -> siren.dk2nu.dk2nu_to_primary_distribution(...)
        -> siren.distributions.PrimaryExternalDistribution in memory

Every SBN detector model uses BNB/SAND coordinates as its geometry frame.
G4BNB files therefore need no extra transform.  G4NuMI files use NuMI
coordinates; pass ``--beam-frame NuMI`` and this example loads the surveyed
NuMI-to-BNB transform from the SBN geometry resource.

Example::

    python VectorPortal_SBN_dk2nu_timing.py \
        /data/g4bnb/*dk2nu*.root --detector SBND --events 500 \
        --output output/sbnd_vector_portal_timing.png

    python VectorPortal_SBN_dk2nu_timing.py \
        /data/g4numi/*dk2nu*.root --beam-frame NuMI \
        --detector ICARUS --events 500 \
        --output output/icarus_numi_vector_portal_timing.png

Times are nanoseconds relative to the primary proton. The dk2nu ancestor
chain supplies the parent-decay time; SIREN propagates later vertices using
the sampled particles' time of flight.
"""

import argparse
import glob
import os

import numpy as np

import siren
from siren import dk2nu
from siren.Injector import Injector
from siren.Weighter import Weighter

from VectorPortal_SBND_dk2nu import build_models, build_vertices


# Detector-coordinate boxes. SBND uses the two active TPC drift volumes from
# the existing dk2nu chain example. The ICARUS box is the active-LAr envelope
# tabulated by the SBN detector resource; load_detector places detector (0,0,0)
# at its center.
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

# The BNB frame origin is the G4BNB SAND-world center. The physical Be target
# is about 39.05 cm downstream. This is the reference used for the prompt
# beta=1 comparison in the plots; moving it by 39 cm changes the reference by
# only about 1.3 ns, but spelling it out avoids hiding that convention.
_BNB_TARGET_GEOMETRY = (0.0, 0.0, 0.3905)


def _beam_frame_to_sbn_geometry(beam_frame):
    """Map the file frame into the SBN model's BNB geometry frame."""
    if beam_frame == "BNB":
        return None
    return siren.resources.detectors.SBN.geo.transform(beam_frame, "BNB")


def _prompt_origin_geometry(beam_frame, frame_transform):
    """Reference point for the beta=1 timing subtraction, in BNB geometry."""
    if beam_frame == "BNB":
        return np.asarray(_BNB_TARGET_GEOMETRY, dtype=float)
    # The NuMI frame origin is MCZERO (the Horn 1 upstream face), the time and
    # position reference used by the G4NuMI beam simulation.
    return np.asarray(frame_transform.apply([0.0, 0.0, 0.0]), dtype=float)


def _expand_input_paths(values):
    """Expand quoted globs while also accepting shell-expanded filenames."""
    paths = []
    for value in values:
        matches = sorted(glob.glob(value))
        if matches:
            paths.extend(matches)
        elif os.path.isfile(value):
            paths.append(value)
    # Preserve order while removing duplicates.
    return list(dict.fromkeys(os.path.abspath(path) for path in paths))


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


def collect_timing(results, detector_model, chi_type, visible_type,
                   prompt_origin_geometry=_BNB_TARGET_GEOMETRY):
    """Extract one timing row per successfully generated event.

    ``scatter_delay_ns`` and ``visible_delay_ns`` subtract the arrival time of
    a beta=1 particle traveling directly from the selected beam simulation's
    origin to the same vertex. They therefore retain the dk2nu
    parent-production time and every massive-particle/decay delay in the
    simulated chain.
    """
    target = _detector_position(detector_model, prompt_origin_geometry)
    c_m_per_ns = float(siren.utilities.Constants.c)

    columns = {
        "weight": [],
        "parent_energy_GeV": [],
        "parent_decay_ns": [],
        "scatter_ns": [],
        "visible_decay_ns": [],
        "scatter_delay_ns": [],
        "visible_delay_ns": [],
        "visible_after_scatter_ns": [],
    }

    for event, weight in results:
        root = event.tree[0].record
        scatter = _record_for_primary(event, chi_type)
        visible = _record_for_primary(event, visible_type)

        scatter_time = np.nan
        visible_time = np.nan
        scatter_delay = np.nan
        visible_delay = np.nan
        visible_after_scatter = np.nan

        if scatter is not None:
            scatter_time = float(scatter.interaction_time)
            distance = np.linalg.norm(
                np.asarray(scatter.interaction_vertex, dtype=float) - target)
            scatter_delay = scatter_time - distance / c_m_per_ns

        if visible is not None:
            visible_time = float(visible.interaction_time)
            distance = np.linalg.norm(
                np.asarray(visible.interaction_vertex, dtype=float) - target)
            visible_delay = visible_time - distance / c_m_per_ns

        if scatter is not None and visible is not None:
            visible_after_scatter = visible_time - scatter_time

        columns["weight"].append(float(weight))
        columns["parent_energy_GeV"].append(float(root.primary_momentum[0]))
        columns["parent_decay_ns"].append(float(root.interaction_time))
        columns["scatter_ns"].append(scatter_time)
        columns["visible_decay_ns"].append(visible_time)
        columns["scatter_delay_ns"].append(scatter_delay)
        columns["visible_delay_ns"].append(visible_delay)
        columns["visible_after_scatter_ns"].append(visible_after_scatter)

    return {name: np.asarray(values, dtype=float)
            for name, values in columns.items()}


def _weighted_hist(ax, values, weights, *, label=None, bins=50):
    mask = np.isfinite(values) & np.isfinite(weights) & (weights >= 0.0)
    if not np.any(mask):
        ax.text(0.5, 0.5, "no completed vertices", ha="center", va="center",
                transform=ax.transAxes)
        return
    ax.hist(values[mask], bins=bins, weights=weights[mask], density=True,
            histtype="step", linewidth=1.8, label=label)


def plot_timing(columns, output_path, detector_name, beam_frame="BNB"):
    """Write four weighted timing projections and return the output path."""
    import matplotlib.pyplot as plt

    weights = columns["weight"]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), constrained_layout=True)

    _weighted_hist(axes[0, 0], columns["parent_decay_ns"], weights)
    axes[0, 0].set_xlabel("parent decay time after proton [ns]")
    axes[0, 0].set_ylabel("weighted density")

    _weighted_hist(axes[0, 1], columns["scatter_ns"], weights,
                   label=r"$\chi$ scatter")
    _weighted_hist(axes[0, 1], columns["visible_decay_ns"], weights,
                   label=r"$V_1\to e^+e^-$")
    axes[0, 1].set_xlabel("absolute vertex time after proton [ns]")
    axes[0, 1].legend()

    _weighted_hist(axes[1, 0], columns["scatter_delay_ns"], weights,
                   label=r"$\chi$ scatter")
    _weighted_hist(axes[1, 0], columns["visible_delay_ns"], weights,
                   label=r"$V_1\to e^+e^-$")
    axes[1, 0].set_xlabel(r"delay relative to prompt $\beta=1$ [ns]")
    axes[1, 0].set_ylabel("weighted density")
    axes[1, 0].legend()

    _weighted_hist(axes[1, 1], columns["visible_after_scatter_ns"], weights)
    axes[1, 1].set_xlabel(r"$V_1\to e^+e^-$ time after $\chi$ scatter [ns]")

    for ax in axes.flat:
        ax.grid(alpha=0.25)
    fig.suptitle("%s: %s dk2nu vector-portal timing"
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
    parser.add_argument("--tune", action="store_true",
                        help="tune phase-space channel weights first")
    parser.add_argument("--output",
                        help=("PNG path (default: output/<detector>_<beam>_"
                              "vector_portal_timing.png)"))
    args = parser.parse_args(argv)

    files = _expand_input_paths(args.dk2nu_files)
    if not files:
        parser.error("none of the dk2nu paths/globs matched a file")

    # This MesonThreeBodySIRENDecay models pi+ -> mu+ nu_mu V1. Filtering the
    # dk2nu rows to the matching ordinary decay avoids silently mixing parent
    # rows generated for a different channel.
    data = dk2nu.read_dk2nu(
        files,
        parent_pdg=dk2nu.PTYPE_PIPLUS,
        decay_modes=13,
        nu_pdg=14,
        read_time=True,
        entry_stop=args.entry_stop,
    )
    dk2nu.print_summary(data)
    if "t0" not in data:
        raise RuntimeError(
            "These dk2nu files have no ancestor start-time branches. "
            "They can drive event generation, but not this timing example.")
    if len(data["E"]) == 0:
        raise RuntimeError(
            "No pi+ -> mu+ nu_mu rows survived the dk2nu filters")

    detector_model = siren.utilities.load_detector(
        "SBN", detector=args.detector)
    # All SBN models use BNB as their geometry frame.  G4BNB is already in
    # that frame; G4NuMI receives the surveyed NuMI -> BNB rigid transform.
    frame = _beam_frame_to_sbn_geometry(args.beam_frame)
    external = dk2nu.dk2nu_to_primary_distribution(
        data, detector_model, frame=frame)
    prompt_origin = _prompt_origin_geometry(args.beam_frame, frame)

    fid = _FIDUCIALS[args.detector]
    fiducial = siren.geometry.Box(
        widths=fid["widths"], center=fid["center"])
    models = build_models()
    # NuMI-to-ICARUS paths are about 800 m in the shared BNB geometry frame;
    # leave margin for the detector's full projected extent.
    primary, secondaries = build_vertices(
        models, external, fiducial, max_length=1000.0)

    injector = Injector(
        detector=detector_model,
        primary=primary,
        secondaries=secondaries,
        events=args.events,
        seed=args.seed,
    )
    weighter = Weighter(injector, primary_physical=primary.physical)
    if args.tune:
        print(siren.tune.tune(injector, weighter, events=200, rounds=3))
        injector.reset()

    results = siren.generate(
        injector, weighter, events=args.events, on_shortfall="warn")
    results.summary()
    print(injector.report())

    # build_models/build_vertices register these public particle names.
    chi_type = siren.particles.resolve("chi")
    visible_type = siren.particles.resolve("V1_sig")
    columns = collect_timing(
        results, detector_model, chi_type, visible_type,
        prompt_origin_geometry=prompt_origin)

    output = args.output or os.path.join(
        "output", "%s_%s_vector_portal_timing.png"
        % (args.detector.lower(), args.beam_frame.lower()))
    output = plot_timing(
        columns, output, args.detector, beam_frame=args.beam_frame)
    csv_path = os.path.splitext(output)[0] + ".csv"
    save_timing_csv(columns, csv_path)
    print("Wrote %s" % output)
    print("Wrote %s" % csv_path)


if __name__ == "__main__":
    main()
