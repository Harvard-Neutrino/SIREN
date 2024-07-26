# DetectorVisuals.wl

`DetectorVisuals.wl` is a Wolfram Language package that takes the detector dimension files in `/densities/` and the corresponding material file in `/materials/` to create a 3D visualisation in `MATHEMATICA`.

To use the package, first save your working notebook in the `\visuals\` directory, and add the `\path\to\visuals` to `$Path` so `MATHEMATICA` can find it.

If you correctly saved the working notebook in `\visuals\`, run 

```
AppendTo[$Path, NotebookDirectory[]];
```

so `MATHEMATICA` has access to the directory, then run 

```
<< DetectorVisuals`
```

To create visualisation for an experiment, run 

```
Visuals["YourExperiment", {Excluded Modules}]
```

to create the 3D model, excluding the specified parts of the detector. For more examples, check out `Det_Visuals_Examples.nb`.