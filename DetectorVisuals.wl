(* ::Package:: *)

BeginPackage["DetectorVisuals`"];

Visuals::usage = "Visuals[Experiment, Modules] displays the experiment as a 3D model without the specified Modules"

Begin["`Private`"];

Visuals[Experiment_, Modules_] :=
Module[{exp = Experiment, mods = Modules, dimensionFile, materialFile, keysDim, materials, color, boxes, polygons, cylinders, spheres,numd,visual},
  dimensionFile = Import[FileNames[All, FileNameJoin[{NotebookDirectory[], "/resources/Detectors/densities", exp}]][[-1]]];
  materialFile = Import[FileNames[All, FileNameJoin[{NotebookDirectory[], "/resources/Detectors/materials", exp}]][[-1]]];

  dimensionFile = Delete[dimensionFile, Position[dimensionFile, _?(#[[1]] == "#" &)]];
  dimensionFile = Delete[dimensionFile, Position[dimensionFile, {}]];
  dimensionFile = Table[If[MemberQ[d,"#"],Drop[d,{Position[d,"#"]//Flatten//First,-1}],d], {d, dimensionFile}];
  keysDim = Table[If[d[[1]]=="object",d[[Position[d,_?(StringQ[#]&)]//Flatten]][[3]],d[[1]]], {d, dimensionFile}];
  dimensionFile = dimensionFile[[Position[!MemberQ[mods, #] & /@ keysDim, True] // Flatten]];

  materialFile = Delete[materialFile, Position[materialFile, _?(#[[1]] == "#" &)]];
  materialFile = Delete[materialFile, Position[materialFile, {}]];
  materialFile = Table[If[MemberQ[d,"#"],Drop[d,{Position[d,"#"]//Flatten//First,-1}],d], {d, materialFile}];

  materialFile = materialFile[[Position[materialFile, _?(StringQ[#[[1]]] &)] // Flatten]];
  materials = Table[m[[1]], {m, materialFile}];

  color = RandomColor[Length[materials]];

  boxes = dimensionFile[[Position[dimensionFile, _?(#[[2]] == "box"&&#[[1]]=="object" &)] // Flatten]];
  polygons = dimensionFile[[Position[dimensionFile, _?(#[[2]] == "extr"&&#[[1]]=="object" &)] // Flatten]];
  cylinders = dimensionFile[[Position[dimensionFile, _?(#[[2]] == "cylinder"&&#[[1]]=="object" &)] // Flatten]];
  spheres = dimensionFile[[Position[dimensionFile, _?(#[[2]] == "sphere"&&#[[1]]=="object" &)] // Flatten]];

  numd = Table[{d[[Position[d,_?(NumberQ[#]&)]//Flatten]],color[[Position[materials,d[[-3]]]//Flatten]]}, {d, polygons}];
  visual = Evaluate[Flatten[{Table[Flatten[{{color[[Position[materials, b[[13]]] // Flatten]], Opacity[0.1]}, 
           Cuboid[{b[[3]], b[[4]], b[[5]]} - {b[[9]], b[[10]], b[[11]]}/2, {b[[3]], b[[4]], b[[5]]} + {b[[9]], b[[10]], b[[11]]}/2]}], {b, boxes}],
      Table[Flatten[{{color[[Position[materials, s[[11]]] // Flatten]], Opacity[0.1]},
           Sphere[{s[[3]], s[[4]], s[[5]]}, s[[9]]]}], {s, spheres}],
      Table[Flatten[{{color[[Position[materials, c[[13]]] // Flatten]], Opacity[0.1]}, 
           Cylinder[{{c[[3]], c[[4]], c[[5]]} - {0, 0, c[[11]]}/2, 
                {c[[3]], c[[4]], c[[5]]} + {0, 0, c[[11]]}/2}, c[[9]]]}], {c, cylinders}],
      Table[Flatten[{{n[[2]], Opacity[0.1]}, 
           Polyhedron[Table[Table[n[[1]][[1;;3]]+AppendTo[p,z[[1]]],
                  {p,ArrayReshape[n[[1]][[8;;-11]],{(Length[n[[1]]]-17)/2,2}]*z[[-1]]}], 
                  {z,ArrayReshape[n[[1]][[-9;;-2]],{2,4}]}]]}],{n,numd}]
  }]];
  Graphics3D[visual]
];

End[];

EndPackage[];
