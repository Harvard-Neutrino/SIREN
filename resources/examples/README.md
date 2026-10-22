# Examples

- `authoring_vertex.py`: data-free explicit Python two-body sampler with an
  analytic accepted decay probability; see `../../docs/quickstart.md`.
- `example1/DIS_{ATLAS,DUNE,IceCube}.py`: native DIS with explicit physical and
  proposal distributions; requires the named detector, spline and flux data.
- `example2/DipolePortal_*.py`: generic DarkNews workflows, requiring DarkNews
  and their detector/flux inputs; CCM supports `--save-tables DIRECTORY`. These are not a validation of a complete BSM
  cascade or of any Dutta/Kim model.
- `example3/CC_MARLEY_CCM.py`: MARLEY charged-current table, CCM inputs, and native MARLEY build support.

The scripts expose `--help`; resource-dependent runs require the datasets.
Production generation requests strict failure and shortfall handling. Older
Controller examples are retained under `legacy/` for migration reference,
with the CCM ProcessBundle unpacking adapted to the new return contract.
Charm/beam and directed/tuning examples require later stack units and are not
part of the supported U05 workflow. Model-specific Dutta/Kim code belongs in
BSM-beam.
