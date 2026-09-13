# Beam timing from dk2nu in SBND or ICARUS

The short answer to the CSV question is: **do not make a CSV for the dk2nu
workflow.** `PrimaryExternalDistribution` can read a CSV, but it also accepts
an in-memory table. `siren.dk2nu.dk2nu_to_primary_distribution` builds that
table directly from the ROOT data and performs the unit, frame, and POT
bookkeeping.

There are two runnable examples.  The first,
[`VectorPortal_SBN_dk2nu_timing.py`](VectorPortal_SBN_dk2nu_timing.py),
models the off-shell Dutta-Kim chain

```text
pi+ -> mu+ nu V1
V1 -> chi chi
chi Ar -> chi V1 Ar
V1 -> e+ e-
```

using G4BNB or G4NuMI dk2nu pion rows as the source and either the SBND or
ICARUS SBN geometry. It makes event-weighted histograms of the pion decay
time, the dark matter scattering and visible-decay times, and their delays
relative to a prompt beta=1 particle from the selected beam reference.

The second,
[`DarkNewsHNL_SBN_dk2nu_timing.py`](DarkNewsHNL_SBN_dk2nu_timing.py), uses
the same dk2nu pion rows with the in-tree `BeamDecays` model and DarkNews:

```text
pi+ -> mu+ nu_mu
nu_mu Ar -> N4 Ar
N4 -> nu gamma
```

This is a single event chain rather than an energy-only neutrino flux
approximation, so it retains the pion decay time and propagates both the
neutrino and HNL times of flight.

## Running it

From an environment with SIREN installed:

```bash
python resources/examples/example4/VectorPortal_SBN_dk2nu_timing.py \
    '/path/to/g4bnb/*dk2nu*.root' \
    --detector SBND \
    --events 500 \
    --output output/sbnd_vector_portal_timing.png
```

For ICARUS, only the detector selection changes:

```bash
python resources/examples/example4/VectorPortal_SBN_dk2nu_timing.py \
    '/path/to/g4bnb/*dk2nu*.root' \
    --detector ICARUS \
    --events 500 \
    --output output/icarus_vector_portal_timing.png
```

The same vector-portal example accepts G4NuMI files. Their coordinates are in
the NuMI beam frame, so select the transform explicitly:

```bash
python resources/examples/example4/VectorPortal_SBN_dk2nu_timing.py \
    '/path/to/g4numi/*dk2nu*.root' \
    --beam-frame NuMI \
    --detector ICARUS \
    --events 500 \
    --output output/icarus_numi_vector_portal_timing.png
```

Productions that split decay-in-flight and decay-at-rest into separate file
sets are supported by both timing examples. The standard g4numi
configuration kills every particle other than a neutrino once its kinetic
energy drops below 0.05 GeV, so those files contain no decays from slower
parents; a dedicated decay-at-rest production tracks everything to rest and
records decays at all parent energies, with its own protons-on-target. Pass
the decay-at-rest files through `--dar-files` and the merging helper beside
the examples (`_beam_samples.py`, analysis-level code, not part of the
`siren.dk2nu` reader) keeps each sample on its side of the kinetic-energy
boundary (`--dar-ke-cut`, default 0.05 GeV, the g4numi kill threshold) and
rescales the decay-at-rest weights by the ratio of the two samples'
protons-on-target, so the merged sample is normalized consistently:

```bash
python resources/examples/example4/VectorPortal_SBN_dk2nu_timing.py \
    '/path/to/g4numi/g4numi_*_fhc_1*.root' \
    --dar-files '/path/to/g4numi/g4numi_*_fhc_dar_*.root' \
    --beam-frame NuMI \
    --detector ICARUS \
    --events 500 \
    --output output/icarus_numi_vector_portal_timing.png
```

For the DarkNews HNL example:

```bash
python resources/examples/example4/DarkNewsHNL_SBN_dk2nu_timing.py \
    '/path/to/g4bnb/*dk2nu*.root' \
    --detector SBND \
    --events 500 \
    --m4 0.10 \
    --mu-tr-mu4 2.5e-6 \
    --output output/sbnd_darknews_hnl_timing.png
```

For G4NuMI-driven HNL production, add `--beam-frame NuMI` in the same way:

```bash
python resources/examples/example4/DarkNewsHNL_SBN_dk2nu_timing.py \
    '/path/to/g4numi/*dk2nu*.root' \
    --beam-frame NuMI \
    --detector ICARUS \
    --events 500 \
    --m4 0.10 \
    --mu-tr-mu4 2.5e-6 \
    --output output/icarus_numi_darknews_hnl_timing.png
```

Use `--detector ICARUS` for ICARUS with either beam. The HNL example defaults
to a Dirac dipole-portal model with `N4 -> nu gamma`; edit
`_darknews_bundle()` when a different DarkNews model is desired. It
deliberately requests only the `Ar40` DarkNews target, so its nonzero
upscattering density is confined to argon-containing sectors. Loading every
material found in the complete detector-building GDML would construct many
irrelevant processes. Add the nuclei from selected upstream sectors when the
study should include HNL production before the active liquid argon.

The HNL example samples dk2nu rows through a threshold bias: a row whose
forward Doppler bound (the largest neutrino energy its parent can produce)
does not exceed the N4 mass gets zero sampling weight, since it cannot
drive the upscatter and its physical contribution is exactly zero. This
matters mostly with `--dar-files`: stopped pions make 30 MeV neutrinos,
far below the default N4 threshold, and would otherwise dominate the
sampling with rows that can never produce an event. The generation
probability accounts for the bias, so event weights are unchanged.

The HNL example also fills the DarkNews interpolation tables before injection
begins, up to the forward Doppler bound of the loaded dk2nu rows (the largest
neutrino energy any surviving row can produce), and then saves the filled
tables next to the DarkNewsTables resource. The first run therefore pays the
DarkNews evaluation cost up front instead of stalling inside the injection
loop, and every later run with the same model parameters loads the tables
from disk. `--table-emax` overrides the fill range in GeV, and
`--no-precompute-tables` builds the tables lazily during injection. The
vector-portal examples need no such step: their Dutta-Kim models are
analytic and compute their widths at construction.

The quoted glob lets the script expand the files itself. `--entry-stop 10000`
is useful for a quick first test. The script writes both the PNG and a CSV with
one row per generated event, so the plotted quantities can be studied without
rerunning the simulation.

The Python dependencies beyond SIREN are `uproot` for the ROOT input and
`matplotlib` for the plot. The HNL example also requires the DarkNews version
supported by that SIREN installation. Loading the SBN detector may download
its GDML data the first time, and the first DarkNews run may need to build
interpolation tables.

## What the loader consumes

`read_dk2nu()` extracts these quantities from each selected row:

| dk2nu value | Meaning | Unit in file |
|---|---|---|
| `ptype` | parent species | PDG code |
| `pdpx,pdpy,pdpz` | parent momentum at decay | GeV |
| `vx,vy,vz` | parent decay vertex | cm |
| `nimpwt` | beam-simulation importance weight | dimensionless |
| `ntype,ndecay` | neutrino species and decay-channel tags | integer codes |
| last ancestor `startt` | parent decay time after the primary proton | ns |
| `dkmeta/pots` | simulated exposure | POT |

The example deliberately requests `parent_pdg=211`, `decay_modes=13`, and
`nu_pdg=14`, matching the pi+ -> mu+ nu_mu production model. This is important:
filtering only on the parent can mix rows representing a different ordinary
decay channel.

`dk2nu_to_primary_distribution()` then constructs the in-memory external
distribution with:

- on-shell parent energy and momentum in GeV;
- the recorded decay vertex, converted from cm to m, then from the selected
  beam frame through BNB geometry into detector coordinates;
- the parent mass;
- `weight = nimpwt / simulated_POT`;
- `t0`, in ns, when the ancestor branches are present.

For this primary vertex the recorded decay location is fixed and the
distribution is weighted with `siren.Fixed()`: the dk2nu row already says that
the pion decayed at that location and into the selected ordinary channel.
The vector-portal example applies its BSM meson decay at that vertex. The HNL
example instead applies `BeamDecays.MesonTwoBodyLeptonicDecay(211)` and directs
the resulting `NuMu` daughter toward the active detector volume with a small
physical fallback. Later vertices receive physical time-of-flight propagation
from SIREN.

## G4BNB and G4NuMI coordinate frames

The `geometry` frame of every SBN detector model is the G4BNB frame: its origin
is the SAND world center and its axes follow the BNB beam. G4BNB dk2nu files
therefore use `frame=None`.

G4NuMI dk2nu files instead store decay positions and momenta in the NuMI frame,
whose origin is MCZERO at the Horn 1 upstream face. The examples implement
`--beam-frame NuMI` by loading the surveyed NuMI-to-BNB transform already
carried by the SBN detector resource:

```python
import siren

numi_to_bnb = siren.resources.detectors.SBN.geo.transform("NuMI", "BNB")

external = siren.dk2nu.dk2nu_to_primary_distribution(
    data,
    detector_model,
    frame=numi_to_bnb,
)
```

The helper rotates momenta and polarization axes, and rotates plus translates
decay positions after converting cm to m. The detector model then performs its
usual BNB-geometry-to-detector conversion. This transform is independent of
whether the selected detector is SBND or ICARUS. Do not translate directly to
the detector and then also pass the result through the detector model; that
would transform positions twice.

## Reading the plots

SIREN's time unit is the nanosecond. The plotted absolute times are relative to
the primary proton used by the beam simulation. The prompt-subtracted panels
use

```text
delay = simulated vertex time
        - distance(beam reference, vertex) / c
```

so they include the parent production/decay time and delays from massive BSM
particles and displaced decays. The histograms use `Results.weights`; an
unweighted histogram would describe the biased injection sample rather than
the physical prediction.

For G4BNB the beam reference is the physical Be target, 0.3905 m downstream of
the SAND-frame origin. For G4NuMI it is the NuMI MCZERO origin, transformed into
the SBN model's BNB geometry frame with the same survey transform used for the
dk2nu rows.

The dk2nu `t0` does **not** include the accelerator bunch position within a
spill. To compare to detector trigger time, sample or supply the bunch time and
add it as a final offset. Keep that beam-clock smearing separate from the
particle transport delay so each contribution can be inspected independently.

## If a CSV is needed anyway

The CSV constructor remains useful for a table produced by another program.
For a fixed parent-decay vertex its header would be

```text
E,px,py,pz,x,y,z,m,t0,weight
```

with GeV, meters in **detector coordinates**, ns, and a physical per-POT weight.
For a propagating neutrino ray, use `x0,y0,z0` for its initial/reference point
instead of `x,y,z`; the latter fixes the interaction vertex and is therefore
wrong when SIREN is supposed to sample where the neutrino upscatters.

## How the DarkNews chain uses BeamDecay

`dk2nu_to_primary_distribution()` correctly returns parent mesons at their
recorded decay vertices; those rows must not be relabeled as neutrinos. The HNL
example makes the neutrino as an ordinary secondary vertex instead:

```python
import siren

beam_decays = siren.resources.processes.BeamDecays

pion = siren.Vertex(
    siren.dataclasses.ParticleType(211),
    beam_decays.MesonTwoBodyLeptonicDecay(211),
    distributions=[external],
    physical=[external],
    weighting=siren.Fixed(),
    kinematics=(
        0.99 * siren.channels.toward("NuMu", fiducial)
        + 0.01 * siren.channels.physical()
    ),
    expand=(siren.expand.child("NuMu"),),
)

neutrino = siren.Vertex(
    "NuMu",
    darknews_bundle.primary[siren.particles.NuMu],
    position=siren.dist.DecayRangeVertex(
        fiducial,
        siren.interactions.InteractionCollection(
            siren.particles.N4,
            darknews_bundle.secondary[siren.particles.N4],
        ),
        daughter_mass=m4,
        daughter_energy_fraction=1.0,
        max_length=1000.0,
    ),
    expand=(siren.expand.child("N4"),),
)
```

The directed two-body channel samples the exact pion-decay kinematics and its
generation density is removed by the `Weighter`; it does not substitute the
pion energy for the neutrino energy. The neutrino interaction vertex then uses
the DarkNews cross sections. Its `DecayRangeVertex` density is proportional to
the neutrino first-interaction density times the conditional probability for a
collinear, equal-energy proxy `N4` to decay inside the fiducial volume. The
parent survival depth, local interaction rate, and downstream decay depth are
all evaluated through the detector model. This reduces to the familiar
piecewise two-exponential density in a uniform one-dimensional medium, while
also handling layered or continuously varying density profiles. A final
bounded `N4` vertex samples the actual decay location; the Weighter removes
both position biases.

The proxy is used only to improve variance: the upscatter model still samples
the actual `N4` energy and direction, and the reported numerical generation
density exactly matches the sampled interaction-depth distribution. Change
`daughter_energy_fraction` when a model has a reliably different typical
energy transfer. To permit upstream production, load DarkNews cross sections
for the nuclei in those upstream detector sectors as well as `Ar40`; a zero
cross section for a material correctly gives zero upscatter density there.
