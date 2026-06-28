# Charm-DIS production in SIREN

SIREN can produce charmed hadrons from neutrino deep-inelastic scattering and
decay them semileptonically -- the production chain behind the IceCube diffuse
multi-cascade tau/charm search. There are two interchangeable production paths;
both emit a `D` meson directly as a DIS secondary (no separate hadronization
vertex), which a `CharmMesonDecay` secondary process then decays.

## Two production paths

| | `PythiaDISCrossSection` | `QuarkDISFromSpline` |
|---|---|---|
| Final state | sampled live from **Pythia8** | analytic slow-rescaling `(xi, y)` |
| Rate | total `sigma(E)` spline (**required**) | total `sigma` spline |
| Differential | **optional** (`d2sigma/dx dy` spline) | required (`dsdxidy` spline) |
| Reweightable | unbiased-only by default; reweightable if a differential spline is supplied | fully reweightable |
| Current | **charged-current only** | CC and NC |
| Use it for | matching Pythia's full hadronization exactly | a fast, fully analytic, reweightable generator |

`PythiaDISCrossSection` is "Pythia as the integrated generator": `SampleFinalState`
is pure Pythia, and the splines are only the rate/weighting side. When no
differential spline is supplied, `FinalStateProbability` returns a constant that
cancels in the unbiased weight (only the total cross section matters); supply a
differential spline to get a real `dsigma/sigma` density for reweighting.

## Runtime requirements (`PythiaDISCrossSection`)

`SampleFinalState` needs Pythia8 + LHAPDF at runtime. The total-cross-section /
weighting path does **not** (it only reads the splines).

1. Build SIREN with `-DSIREN_WITH_PYTHIA8=ON`.
2. `export LHAPDF_DATA_PATH=<dir>` -- the parent of the installed PDF-set folders.
   The requested `pdf_set` (default `LHAPDF6:HERAPDF20_NLO_EIG`) **must be
   installed there**, or construction raises a clear error naming the set. Any
   installed LHAPDF set works (e.g. `LHAPDF6:CT18NLO`).
3. `export PYTHIA8DATA=<pythia>/share/Pythia8/xmldoc`.
4. macOS: `export DYLD_LIBRARY_PATH=$PREFIX/lib:<pythia>/lib:$DYLD_LIBRARY_PATH`.

NC (`interaction_type=2`) is rejected at construction: Z exchange does not force
charm in Pythia. Use `QuarkDISFromSpline` for NC charm.

## Generating the Pythia splines

The Pythia splines are not committed (they depend on the PDF set and Pythia
version). Regenerate them from the same Pythia config the sampler uses:

```bash
python generate_charm_pythia_splines.py --out ./splines \
    --pdf LHAPDF6:CT18NLO --emin 100 --emax 1e6 --npoints 17
```

This writes `pythia_charm_sigma.fits` (always) and `pythia_charm_dsdxdy.fits`
(unless `--no-differential`). The default grid (100 GeV - 1 PeV) covers the
diffuse analysis band. See `siren.pythia_charm_splines` for the library API.

The slow-rescaling `QuarkDISFromSpline` `dsdxidy_/sigma_nu-N-{cc,nc}-charm-*.fits`
splines are external analysis inputs (LHAPDF-derived, not produced by SIREN);
point `SIREN_CHARM_SPLINE_DIR` at them to run `DIS_IceCube_charm.py` and the
`test_quarkdis_slow_rescaling.py` tests.

## Charmed-hadron decay

`CharmMesonDecay` (the 2-body class used by the chain) handles `D0`, `D+`, `Ds`
and their anti-flavors: `D+/-` and `D0/D0bar` via `K`/`K*(892)` semileptonic
modes with a V-A matrix element, and `Ds` via `eta/eta'/phi` pure phase space.
`CharmMesonDecay3Body` is an alternate D0/D+-only implementation (it rejects
other species at construction). `Decay.TotalDecayLength(record)` returns the lab
decay length `beta*gamma*c*tau` **in meters** -- the cascade separation the
Taupede reconstruction targets.

Every D species a production cross section emits **must** have a registered
secondary decay process, or the injector silently drops it (see the comment in
`DIS_IceCube_charm.py`). The fragmentation fractions (`D0:D+/-:Ds = 0.60:0.23:0.15`,
renormalized /0.98) fold the unmodeled `Lambda_c` into the D mesons.

## Known limitations / deferred work

- **NC charm** is not supported by `PythiaDISCrossSection` (use `QuarkDISFromSpline`).
- **`Lambda_c`** (~0.09 raw fraction) is folded into the D mesons, not modeled
  as a distinct species (it has a different lifetime and semileptonic mix).
- **Per-event Pythia re-init** (~0.1-1 s/event): variable-energy mode is
  unsupported for `WeakBosonExchange`, so each event rebuilds Pythia at its beam
  energy. Fine for analysis-scale samples; a throughput item for very large runs.
