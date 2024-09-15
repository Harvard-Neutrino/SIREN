# DetectorVisuals.wl

`DetectorVisuals.wl` is a Wolfram Language package that takes the detector dimension files in `/densities/` and material file in `/materials/` to create a 3D model in `MATHEMATICA`.

To call functions in the package, save the working notebook in `\visuals\` and add `\path\to\visuals` to `$Path` so `MATHEMATICA` can find the package file.

If the notebook directory is correctly saved in `\visuals\`, run

```
AppendTo[$Path, NotebookDirectory[]];
```

to append `\path\to\visuals\` to `MATHEMATICA` so it has access to the package. Next, run

```
<< DetectorVisuals`
```

to import the package. To create visualisation for an experiment, run

```
Visuals["YourExperiment", {Excluded Modules}]
```

to create the 3D model, excluding specified parts of the detector. For more examples, check out `Det_Visuals_Examples.nb`.
