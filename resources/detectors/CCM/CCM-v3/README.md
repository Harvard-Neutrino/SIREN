# CCM-v3: corrected facility with the simple CCM detector

This is an assembled transport geometry: the corrected target_sim Mark IV
target, moderators, shielding and floor, plus the six nested detector cylinders
from CCM-v2. `ccm-facility.gdml` contains both the facility and detector and can
be opened directly in a GDML viewer. It contains no CCM lead floor slab.

## Load and visualize

```python
import siren

model = siren.utilities.load_detector("CCM-v3")
fiducial = siren.utilities.get_fiducial_volume("CCM-v3")

from pathlib import Path
gdml = Path(siren.utilities.get_detector_model_path("CCM-v3")) / "ccm-facility.gdml"
siren.visualization.view(str(gdml))
```

Pass `model` to `SIREN_Controller` through its explicit
`detector_model=model` argument;
its legacy file-loading path cannot load this GDML resource.
`densities.dat` contains fiducial/frame metadata only, not a transport model.
No `materials.dat` is supplied, so legacy two-file resource loading fails
instead of silently loading an incomplete detector.

SIREN's normal latest-version selection makes `load_detector("CCM")` select
v3 when this resource is installed. Pin `"CCM-v2"` to reproduce the old
approximate facility; those source files are unchanged. This geometry update
does not migrate external beam samplers, event files or BSM-beam configurations.

## Coordinates and provisional placement

The GDML uses **target_sim world coordinates, in metres**. The loader sets
SIREN's detector origin and rotation so `DetectorPosition` queries, events and
the fiducial volume use **detector-local coordinates**:

```
p_facility = Rz(0.43017060669479124) p_detector + (10.4949, 19.86074, 0.2938)
```

The user selected the CAD cryostat axis for the horizontal position. The
simple outer cylinder's bottom is placed at the floor datum `z = -1.0162 m`.
This is an approximation for the retained simple model, not a surveyed active
detector elevation or a CAD vessel/support assembly. Its horizontal distance
from the target origin is about 22.46 m, replacing Darcy's old 23 m distance.
The detector's horizontal orientation *relative to the target direction*
retains Darcy's `0.654498 rad` convention. The installed PMT-column orientation
relative to the target has not been independently verified; a CAD-pattern fit
gives a different candidate and does not supersede this provisional reference.
The old absolute `-0.65 m` detector
height and the old facility axes must not be applied to this model.

Convert target_sim source vertices and directions explicitly:

```python
from siren.detector import GeometryPosition, GeometryDirection
from siren.math import Vector3D

vertex_detector = model.GeoPositionToDetPosition(GeometryPosition(Vector3D(x, y, z)))
direction_detector = model.GeoDirectionToDetDirection(GeometryDirection(Vector3D(ux, uy, uz)))
```

## Detector content and limits

All six CCM-v2 cylinder radii, full heights, bulk steel/argon/aluminium
densities and sector names are retained. The fiducial cylinder remains
`radius = 1.03385 m`, `full height = 1.2396 m`, centered at detector-local zero.
Two explicit material corrections make the exchange valid in Geant4:

- Vacuum is represented by `1e-25 g/cm3`, rather than exact zero.
- The old steel fractions summed to `1.0003`; each is divided by that sum.
  Steel's bulk density remains `7.83 g/cm3`; isotope target densities are
  consequently smaller by a factor of `1/1.0003` than in CCM-v2.

This version retains the approximate continuous aluminium PMT-frame cylinder.
Detailed PMTs, optical response, the CAD cryostat/support solids, and finite
steel trench covers are not included. The actual trench is covered; the
facility transport model still has an opening because cover footprints and
seams remain unresolved. The previous 38.1 mm thickness estimate is not a
measured dimension and has not been inserted as an assumed solid.

## Rebuild and provenance

`geometry.json` holds detector dimensions, materials and placement;
`facility.gdml` is the validated portable facility input. Rebuild the assembled
file and fiducial metadata with:

```bash
python build_geometry.py
```

There is no runtime generation, download or mutable shared cache. The builder
checks the facility input's SHA256 before composing the detector. Input hashes
also identify the unchanged CCM-v2 files used to construct this version.

The facility was exported with Geant4 11.3.2 from target_sim commit
`3e8a43685dce94e5ad9a4aae7afb356ff716693c` plus the local removal of the
nonexistent `CCMFloorLead` slab. The portable conversion preserves the exact
planar boundaries of concave extrusions, bakes one transformed first Boolean
operand and replaces a multi-union cutter with sequential subtractions.
The input hash is
`b543bc7431f92965d985d4e8ebd85c77ddf2ee73f6922c75767286cde65665e4`.
The original six-cylinder description comes from SIREN revision
`43edbbcfb33e0bf178ab0494489502cf9ddf68ba`.

Validation on 2026-09-17: 557 strict SIREN sectors; 498 closed, outward mesh
volumes with no null meshes; zero disagreements with native Geant4 at 37,440
material probes and 1,014 material-column paths; no sampled insert overlaps
at 100,000 surface points per volume and 1 micrometre tolerance. These checks
validate the exchange and assembly, not a measured detector survey.
