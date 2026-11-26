# Interactive geometry viewing

`siren.visualization.view` opens a GDML file or a SIREN `DetectorModel` in a
VTK window. It needs the `visualization-3d` extras (pyg4ometry, VTK) and a
display; pyvista is only required for models with `TriangularMesh` sectors.

```python
from siren import visualization

visualization.view("detector.gdml")
visualization.view(model, screenshot="detector.png")   # off-screen PNG, full detail
visualization.view("detector.gdml", timings=True)      # print phase times on return
visualization.view("detector.gdml", presets="views.json", section="Front")
```

The window opens immediately with a loading status. A worker process meshes a
coarse preview, then the requested detail, while the window stays interactive;
camera position, hidden materials and opacity survive refinement. Closing the
window cancels outstanding work, removes its temporary files and returns to
Python, where another `view()` can be opened at once. Each worker process pays
the pyg4ometry import, so on-demand loads take a few seconds to start.

## Options

| Option | Default | Effect |
| --- | --- | --- |
| `mesh_slices` | 48 | Full-detail slices/stacks for curved solids; `0`/`None` keeps the mesher defaults. |
| `preview_slices` | 12 | Preview resolution; `None` disables the preview. Skipped when full detail is cached. |
| `progressive` | `True` | `False` loads synchronously. Screenshots and `interactive=False` always wait for full detail. |
| `cache`, `cache_dir` | `True`, `~/.cache/siren/meshes` | Persistent mesh cache (or under `XDG_CACHE_HOME`); `cache=False` neither reads nor writes it. |
| `show_gas` | `False` | Volumes with density <= 0.05 g/cm3 start hidden; **g** loads and shows them. |
| `hidden_volumes` | `()` | Logical-volume name globs to defer; **v** loads and restores them. Their children keep their placements. |
| `instancing` | `True` | Repeated placements share GPU prototypes; reflections and shear use exact matrix actors. |
| `display` | `"full"` | `"exterior"` shows only outward faces of annotated region shells (below). |
| `presets` | `None` | Mapping or JSON file with colours, layers, sections and the opening view (below). |
| `section` | `None` | Section to start with: the label of a preset section, or a mapping with `origin` and `normal` in metres. |
| `panel` | `True` | Side panel in interactive windows; off for screenshots and `display="exterior"`. |
| `edges` | `True` | Dark lines along sharp edges. |
| `save_to` | presets file | Where **Save settings** writes; `./siren-view.json` without a presets file. |
| `clipper` | `False` | Start with a section through `clip_origin` (mm) with normal `clip_normal`, and a plane to drag. |
| `cutter` | `False` | Coloured outlines in the x, y and z planes. Sections and outlines force `display="full"`. |
| `coloured`, `axes`, `legend`, `bounding_box`, `picker`, `near_frac` | | Usual display controls. |

Cache files are disposable numeric NPZ files (no pickles), keyed by input and
included-file contents, resolution, backend and library versions. Corrupt
files are rebuilt, and new files are published atomically only when the inputs
were stable throughout the load. Volumes sharing a solid share one cached mesh.
GDML files with `<!DOCTYPE>`/`<!ENTITY>` declarations still load through
pyg4ometry's reader but are never cached (a warning says so). Files using
replicas, divisions or parameterised placements use pyg4ometry's eager
meshing, which limits deferral and cache savings. `<file>` references resolve
against the working directory, as in pyg4ometry.

## Sections

A section removes everything on the positive side of a plane. Each volume it
cuts is filled in its material's colour and the surfaces behind are paler. A
volume fills only the part that none of its daughters occupies, so hiding a
daughter leaves an opening rather than its mother's material. In a SIREN model,
where sectors overlap, the sector with the highest level is the one filled, as
in `DetectorModel`. The fill follows the display meshes, so round surfaces are
as fine as `mesh_slices`: the Earth drawn with 48 facets is up to 14 km off at
its surface. Gases are never filled. Sectors drawn directly as `TriangularMesh`
are clipped but not filled. Sectors that `to_gdml` writes as bounding boxes are
cut in their exact shape when pyvista is installed and left out of the cut, with
a warning, otherwise; their uncut surfaces are still drawn as boxes.

The panel lists the preset sections, a slider that cuts across x, y or z within
the part of the model in view, and **Other side**, which keeps the other half.
**k** shows a plane to drag and turn. Right-clicking a filled face identifies
the volume there.

## Colours and surfaces

The panel switches between custom colours and density colours, and makes every
surface solid or see-through. Rows for materials, layers and volumes show each
one's colour, and a click on a row shows or hides it, makes it see-through or
solid, or opens a palette, as the three buttons above the rows choose.
Right-clicking a volume lists its logical volume under **Volumes**, where it can
be changed apart from other volumes, in all its copies. A volume's own choice
beats its layer's, and a layer's beats its material's; the all-solid and
all-see-through buttons apply where none is made. Density colours ignore custom
colours but keep the other choices. Right-clicks pass through see-through
surfaces when something solid is behind them. **Save settings** writes the
presets with these choices, the current view first in the section list and the
camera as the opening view; a file it replaces is kept as `NAME.bak`.

## Presets

Presets describe a model for display only; nothing detector-specific is built
in: without presets the colours follow density, solids are slightly transparent
and the panel offers the axis cuts and the materials. Lengths are metres in the
geometry frame and every entry is optional.

```python
presets = {
    "title": "Detector hall", "note": "z is up.", "up": [0, 0, 1],
    "opacity": 1.0, "see_through": 0.25,    # of solid and of see-through surfaces
    "colour_mode": "custom", "surfaces": "solid",
    "home": {"centre": [0, 0, 0], "direction": [1, 1, 1], "width": 20, "height": 10},
    "hidden_volumes": ["Ground*"],          # not meshed until v is pressed
    "materials": [{"label": "Steel", "match": ["*Steel*", "G4_Fe"], "colour": "#b04030"}],
    "layers": [{"label": "Roof", "match": ["Roof*"], "shown": True, "edges": True,
                "opacity": "see-through"}],
    "volumes": [{"match": ["Cryostat*"], "colour": "#3f9bd1", "shown": True}],
    "sections": [{"label": "Front", "origin": [0, 0, 0], "normal": [0, -1, 0],
                  "centre": [0, 0, 2], "width": 12, "height": 6}],
}
```

`match` holds name globs; Geant4's `0x...` suffixes are ignored and the first
match wins. Materials, layers and volumes each take an optional `colour`,
`opacity` (`"solid"` or `"see-through"`) and `shown`; materials without a colour
keep their density colour. A layer matches placed or logical volume names, and
daughters belong to their mother's layer unless they match one themselves. A
volume rule matches the same names but not daughters; its `label` defaults to
its first glob. A section is viewed along
its normal; `up`, `centre`, `width` and `height` set the view, and the whole cut
is framed without them.

## Exterior display of region shells

`display="exterior"` keeps only the faces whose normals point radially outward
from the local z axis for volumes annotated as shells of revolution, removing
inner walls and the artificial end caps between regions. The mapping is
supplied by the caller; SIREN carries no detector-specific region names.

```python
regions = {
    "auxtype": "pmt_region",           # <auxiliary auxtype=...> key on logical volumes
    "styles": {                        # region -> legend label, colour, opacity
        "external_tpb": {"label": "PMT TPB", "colour": [0.95, 0.94, 0.82], "alpha": 1.0},
        "bare_transparent_glass": {"label": "PMT bare glass", "colour": [0.72, 0.88, 0.94], "alpha": 0.25},
    },
    "hidden": ["vacuum"],              # regions not displayed at all
}
visualization.view("detector-cad.gdml", display="exterior", regions=regions,
                   hidden_volumes=["CADCryostat*"])
```

Unannotated volumes are rendered normally. Colours are display choices, not
optical properties. These open surfaces never alter the GDML or simulation
solids; use `display="full"` for volume inspection.

## Controls

- Left drag orbits, scroll zooms; arrow keys pan and fly, shift-up/down moves vertically.
- Right-click identifies the placed volume (name, material, density; a supplied
  SIREN model adds its sector query); shift-right-click hides that material.
- In the panel, click to cut and to change rows (above); the wheel scrolls it.
- **g** gas volumes, **v** show everything hidden, **n/m** opacity of every surface, **c** section outlines,
  **k** plane to drag, **p** panel, **l** legend, **b** bounding box, **o** orientation cube,
  **h** help, **q**/**e** quit.

## Timings

Pass a dictionary as `timings` to receive phase times (import, export, parse,
mesh/cache, cache write, scene build, render), first-frame times
(`window_frame_seconds`, `first_geometry_frame_seconds`, `preview_frame_seconds`,
`detail_frame_seconds`), cache hits/misses per stage, shared-mesh reuses,
face and instance counts, and whether the load was cancelled. Times are
measured from the `view()` call. `direct_mesh_error` records why a model's
mesh sectors were not drawn (typically a missing pyvista), and
`request_errors` collects failures of on-demand loads, which leave the viewer
open with a warning.

The `backend="pyvista"` path (`view_pv`) keeps its own renderer and controls;
the options above apply to the pyg4ometry backend.
