# Units and coordinate conventions

Four-vectors are `(E, px, py, pz)` in GeV with metric `(+,-,-,-)`.
Positions in InteractionRecord are detector coordinates in metres. Times are
ns, target number densities cm^-3, cross sections cm^2, and decay widths GeV.
Use explicit conversions at external interfaces.

`geometry.Box(widths=[x, y, z], center=[x0, y0, z0])` uses full side lengths.
The scalar positional constructor `Box(x, y, z)` raises a configuration error;
Placement/copy constructors remain available. Other solids keep their documented
conventions, such as Cylinder radius and full length.

`Measure.SolidAngleRest()` and `Measure.SolidAngleLab()` name different angular
measures. Density and sampler must use the same frame. `check_closure` uses an
independent validation stream and reports unsupported references as incomplete;
a passing check at one energy and sample count is not general physics validation.

Injector archive version 2 stores expansion rules; SecondaryExpansion version 1
also stores a distinct primary policy. PhysicalProcess version 1 stores phase
spaces. Injector/random archives with current RNG state resume the stream;
older version-0 archives restart from their stored seed. Setting a seed or
calling reset explicitly starts a new stream. Results captures attempt counts
separately so later generation cannot change its stored normalization.
