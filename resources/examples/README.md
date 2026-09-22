# SIREN examples

Every maintained script declares its chain with `siren.Vertex` objects: the
`distributions` of a vertex sample the injection, `physical` names the flux
and direction factors the weighter divides by, and `expand` rules say which
secondaries recurse into their own vertex (a terminal vertex declares
`expand=(siren.expand.depth_below(0),)`). `siren.Simulation` builds the
Layer-2 `Injector` and `Weighter` from the vertices and returns a
`siren.Results` whose `save(prefix)` writes HDF5, Parquet, and native
`.siren_events` files. Each script accepts `--events`, `--seed`, and
`--output`; the defaults reproduce the historical configurations.

| Directory | Scripts | Interface |
| --- | --- | --- |
| `example1/` | `DIS_IceCube.py`, `DIS_DUNE.py`, `DIS_ATLAS.py`: muon-neutrino CC DIS with the CSMS splines | `siren.Simulation` over one primary `Vertex` |
| `example1/` | `DIS_IceCube_charm.py`: charm DIS with per-D-meson energy-loss and decay vertices (needs external splines, see `README_charm.md`) | Layer 2: `siren.injection.Injector`, `siren.injection.Weighter`, `siren.generate` |
| `example2/` | `DipolePortal_{CCM,MiniBooNE,MINERvA,ND280UPGRD}.py`: DarkNews dipole-portal upscattering followed by the N4 decay | `siren.Simulation` with a primary and a secondary `Vertex` |
| `example3/` | `CC_MARLEY_CCM.py`: MARLEY CC interactions from pi-DAR neutrinos; `plot_energy_and_vertex.py` reads the Parquet output | `siren.Simulation` over one primary `Vertex` |

`docs/quickstart.md` and `docs/simulation.md` describe the interface; the
notebooks (`example1/PaperPlots.ipynb`, `example2/PaperPlots.ipynb`,
`additional_paper_plots/PaperPlots.ipynb`) are the paper's figure sources
and still use the deprecated `SIREN_Controller`. `legacy/` keeps the
`SIREN_Controller` versions of the scripts and `legacy_examples/` the
LeptonInjector-era scripts; neither is maintained.
