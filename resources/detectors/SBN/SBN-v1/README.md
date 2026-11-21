# SBN detector and beamline geometry

The SBN loader places the selected detector, the BNB beamline, the NuMI beamline
and site geology in a common BNB coordinate frame.

```python
import siren

model = siren.utilities.load_detector("SBN", detector="ICARUS", numi_config="ME")
```

`numi_config` selects NuMI mechanical geometry. It defaults to `"ME"`
(medium energy / NOvA), and is case-insensitive. This is currently the only
supported configuration; values such as `"LE"`, `"HE"` or `"me000z200i"` raise
`ValueError` before downloading or composing geometry. It does not select
neutrino flux files, beam polarity, currents or a particular run period.
`earth_model=True` and `lbnf=True` remain independent options.

The ME asset is `numi_ME_g4export_2026-09-17.gdml`, verified by SHA256
`730466f287196d65a7fee074203014471faee6be0fbfa3da4769046d92355ed7`.
It uses explicit ME target, horn and baffle positions, and the production
template's alternate horn-1 geometry. The export macro, source provenance and
validation are maintained with the asset in
[SIREN-data](https://github.com/SIREN-Generator/SIREN-data/tree/328413f691b33052a77b1955c7de0820a17cb7a4/detectors/SBN/v1/NuMI).

The older `numi_g4export_2026-05-19.gdml` combined the tall ME target with legacy
positions and contained real target–horn overlaps. It is retained in SIREN-data
for provenance, but is not an available beam configuration. Existing calls
that omit `numi_config` now use the corrected ME asset. Its new filename also
gives the composed GDML a distinct cache identity; existing cached files remain
available for reproducing older work.

The data provenance documents residual horn-envelope and downstream shielding/
containment findings. Correcting the ME placement does not establish that the
complete beamline is overlap-free.

## Detectors

`detector` selects `"ICARUS"`, `"SBND"`, `"MicroBooNE"`, `"MiniBooNE"` or
`"DUNE_ND"`. ICARUS, SBND and DUNE_ND load their GDML exports from SIREN-data;
MiniBooNE is a placeholder tank generated locally by `sbn_loader`.

MicroBooNE uses the uboonecode production geometry
`microboonev12_nowires.gdml`: the cryostat, TPC, PMTs, CRT and the LArTF pit,
hall and local ground. The loader downloads it from a
[pinned SIREN-data revision](https://github.com/SIREN-Generator/SIREN-data/tree/df2d5a77fedfacafca0a913203609d17b8521f4e/detectors/SBN/v1/MicroBooNE),
which hosts the byte-identical file from `uboone/ubcore` `03c0bb06` with a
provenance README. It checks the SHA-256 (`a33e1d1d…d0215c`) on every load and
writes a copy, `microboonev12_nowires_siren.gdml`, without the LArSoft
`volVacuumSpace` placement: a 1.5 km vacuum box above grade that would
otherwise replace the composite's atmosphere. The copy is rebuilt on every
load, so an edited or stale one is never used.

The LArSoft world origin sits at BNB `(-1.24325, 0.0093, 463.363525)` m, the
inverse of the beam origin in MicroBooNE's own beam-to-detector transform
(`ubsim` `FluxReaderBNB.cxx`). This puts the TPC centre 468.55 m from the beam
origin, the published 468.5 m baseline (MICROBOONE-NOTE-1031). G4BNB's nominal
`(0, 0, 470)` m lies 1.45 m further downstream and is not used. As for ICARUS
and SBND, the detector origin is the centre of the active volume,
`volTPCActive`, which is offset `(-1.55, +0.97, 0)` cm from the TPC centre.

The transform is a pure translation, as in MicroBooNE's production flux
conversion (`BooNEtoGSimple.cxx`). `FluxReaderBNB.cxx` also applies a
mrad-scale rotation; it is not a beam-axis correction, shifts points by at most
3.6 cm across the TPC, and is omitted. JINST 12 (2017) P02017 and JINST 16
(2021) P04004 describe the detector and the building.
