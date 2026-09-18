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
