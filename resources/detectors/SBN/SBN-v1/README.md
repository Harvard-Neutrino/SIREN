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
hall and local ground. The byte-identical upstream file from `uboone/ubcore`
`03c0bb06` is hosted with its provenance README in
[SIREN-data](https://github.com/SIREN-Generator/SIREN-data/tree/df2d5a77fedfacafca0a913203609d17b8521f4e/detectors/SBN/v1/MicroBooNE),
fetched from that pinned revision. Its SHA-256 `a33e1d1d…d0215c` is checked on
every load, not only on download, because the file is rewritten before it is
composed. The SIREN copy `microboonev12_nowires_siren.gdml` differs only by
dropping the LArSoft `volVacuumSpace` placement, a 1.5 km vacuum box above
grade that would otherwise replace the composite's atmosphere; it records the
source digest, and is rebuilt if that digest changes.

The LArSoft world origin sits at BNB `(-1.24325, 0.0093, 463.363525)` m,
inverted from the beam origin in MicroBooNE's own beam-to-detector transform
(`ubsim` `FluxReaderBNB.cxx`). The TPC centre is then 468.55 m from the beam
origin, the published 468.5 m baseline (MICROBOONE-NOTE-1031). The active
volume is the sector `volTPCActive`, offset `(-1.55, +0.97, 0)` cm from the
TPC-box centre; as for ICARUS and SBND the detector origin is its centre, so
the two share a z and the baseline is the same either way. G4BNB's nominal `(0, 0, 470)` m is 1.45 m
downstream of that centre and is not used. The edge is a pure translation, as
in MicroBooNE's production flux conversion (`BooNEtoGSimple.cxx`); the
mrad-scale rotation `FluxReaderBNB.cxx` also carries is not a beam-axis
correction and is omitted, at most 3.6 cm across the TPC. See JINST 12 (2017)
P02017 and JINST 16 (2021) P04004 for the detector and building.
