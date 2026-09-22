"""Dark-neutrino e+e- events at an SBN detector from G4BNB dk2nu beam parents.

The single dark neutrino model of the MicroBooNE search for dark-sector e+e-
explanations of the MiniBooNE anomaly (Abdullahi et al., PRL 136, 121804
(2026)): a muon neutrino upscatters on argon to a heavy neutral lepton N4
through a light dark photon Z', and the N4 decays to nu e+e- through the
on-shell Z'. The benchmark is the paper's Fig. 3 example point, m4 = 106 MeV,
mZ' = 30 MeV, |U_mu4|^2 = 2e-10, epsilon = 8e-4, g_D = 2, |U_D4|^2 = 1, with
a Dirac N4; DarkNews defaults cover the rest.

Instead of a neutrino flux table, each run injects one species of recorded
beam parent (pi+, K+, or mu-) at its dk2nu decay vertex and decays it with the
BeamDecays models, so the chain is

    parent -> nu_mu + ...,   nu_mu N -> N4 N,   N4 -> nu e+ e-

The parent vertex is weighted in Fixed mode (the dk2nu row already records the
decay) and aims the neutrino at the active volume with a physical fallback.
The neutrino upscatters where the N4 can still decay inside the active volume
(the physical interaction density times the N4's decay probability in the
volume), on argon by default or, with --all-targets, on every nucleus of the
site so that long-lived N4s produced upstream count too. The weighter removes
every bias. Event weights are per POT, so the outputs of the species add.
MicroBooNE itself is not in the SBN detector resource; SBND (default) or
ICARUS stand in.

Requires G4BNB dk2nu ROOT files (not shipped with SIREN): pass --dk2nu or set
DK2NU_FILES to a path or glob.
"""
import argparse
import glob
import math
import os

import siren
from siren import channels, dk2nu, expand

parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--dk2nu", default=os.environ.get("DK2NU_FILES"),
                    help="G4BNB dk2nu ROOT file(s), a path or a glob")
parser.add_argument("--detector", choices=("SBND", "ICARUS"), default="SBND")
parser.add_argument("--events", type=int, default=100,
                    help="events per parent species")
parser.add_argument("--seed", type=int, default=None)
parser.add_argument("--output", default="output/DarkNeutrino")
parser.add_argument("--all-targets", action="store_true",
                    help="upscatter on every nucleus of the detector model "
                         "(soil, concrete, cryostat, argon) instead of argon "
                         "only; generates more DarkNews tables")
args = parser.parse_args()

files = sorted(glob.glob(args.dk2nu or ""))
if not files:
    raise SystemExit("no dk2nu files: pass --dk2nu or set DK2NU_FILES")

model_kwargs = {
    "m4": 0.106,
    "mzprime": 0.030,
    "Umu4": math.sqrt(2.0e-10),
    "UD4": 1.0,
    "epsilon": 8e-4,
    "gD": 2.0,
    "decay_product": "e+e-",
    "HNLtype": "dirac",
}

# Active liquid-argon volumes in detector coordinates (metres).
ACTIVE_VOLUME = {
    "SBND": dict(widths=(4.026, 4.074645, 5.01), center=(0.0, 0.59, -0.415)),
    "ICARUS": dict(widths=(7.20, 3.16, 17.95), center=(0.0, 0.0, 0.0)),
}

NuMu, N4 = siren.particles.NuMu, siren.particles.N4
BeamDecays = siren.resources.processes.BeamDecays

# Parent species: PDG code, dk2nu row selection, and the SIREN decay model.
# Row selection keeps the channel the model implements (K+ -> mu+ nu_mu is
# dk2nu decay mode 5, pi+ -> mu+ nu_mu mode 13) or, for muons, the nu_mu row
# of each three-body decay.
PARENTS = {
    "piplus": (211, {"decay_modes": 13}, BeamDecays.MesonTwoBodyLeptonicDecay(211)),
    "Kplus": (321, {"decay_modes": 5}, BeamDecays.MesonTwoBodyLeptonicDecay(321)),
    "muminus": (13, {"nu_pdg": 14}, BeamDecays.MuonThreeBodyDecay(13)),
}

detector_model = siren.load_detector("SBN", detector=args.detector)
fiducial = siren.geometry.Box(**ACTIVE_VOLUME[args.detector])

dn_version = siren.utilities.darknews_version()
table_name = f"DarkNewsTables-v{dn_version}/"
table_name += "DarkNeutrino_M%2.2e_Z%2.2e_U%2.2e_eps%2.2e" % (
    model_kwargs["m4"], model_kwargs["mzprime"],
    model_kwargs["Umu4"], model_kwargs["epsilon"])

# Upscattering targets: argon only by default; --all-targets uses every
# nucleus of the detector model, so the neutrino can also upscatter in the
# soil, the concrete and the cryostat upstream of the active volume.
bundle = siren.load_processes(
    "DarkNewsTables",
    primary_type=NuMu,
    detector_model=detector_model,
    nuclear_targets=None if args.all_targets else ["Ar40"],
    table_name=table_name,
    **model_kwargs,
)
hnl_decays = siren.interactions.InteractionCollection(N4, bundle.secondary[N4])

# Longest path a beam particle needs: from the BNB target, upstream of every
# parent decay, to the far end of the active volume.
reach = (detector_model.GetDetectorOrigin().get().GetZ()
         + ACTIVE_VOLUME[args.detector]["widths"][2])

neutrino = siren.Vertex(
    NuMu, bundle.primary[NuMu],
    # Upscatter where the N4, taken collinear at the neutrino's energy, can
    # still decay inside the active volume.
    position=siren.dist.DecayRangeVertex(
        fiducial, hnl_decays, model_kwargs["m4"], max_length=reach),
    expand=(expand.child("N4"),),
)

hnl = siren.Vertex(
    N4, bundle.secondary[N4],
    position=siren.dist.BoundedVertex(fiducial, reach),  # decay in the active volume
    expand=(expand.depth_below(0),),  # terminal: nu e+ e- are final
)


def parent_vertex(pdg, decay, parents):
    """A recorded parent decaying at its dk2nu vertex, aiming nu_mu at the TPC."""
    if abs(pdg) == 13:  # three-body: the electron is the spectator
        aim = channels.toward_3body("NuMu", fiducial, spectator="EMinus")
    else:
        aim = channels.toward("NuMu", fiducial)
    return siren.Vertex(
        siren.dataclasses.ParticleType(pdg), decay,
        distributions=[parents],
        physical=[parents],
        weighting=siren.Fixed(),  # the row already records the decay vertex
        kinematics=0.99 * aim + 0.01 * channels.physical(),
        expand=(expand.child("NuMu"),),
    )


for species, (pdg, rows, decay) in PARENTS.items():
    data = dk2nu.read_dk2nu(files, parent_pdg=pdg, **rows)
    if len(data["E"]) == 0:
        print(f"{species}: no rows in the dk2nu files, skipped")
        continue
    parents = dk2nu.dk2nu_to_primary_distribution(data, detector_model)
    print(f"{species}: {parents.GetPhysicalNumEvents()} recorded parents")

    sim = siren.Simulation(
        events=args.events, seed=args.seed, detector=detector_model,
        primary=parent_vertex(pdg, decay, parents), secondaries=[neutrino, hnl])
    results = sim.run()
    results.summary()
    results.save(f"{args.output}_{species}")
