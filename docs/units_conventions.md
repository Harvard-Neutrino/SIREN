# Units and conventions

## Units

Energies and masses are in GeV, lengths are in meters, and times are in
nanoseconds (ns). Internally, SIREN fixes `meter = 1` and `second = 1e9`, so a
time value of `1.0` is 1 ns; this is the convention HepMC3 export also uses
for vertex/particle times.

## Geometry arguments are FULL widths

`Box` and similar geometry constructors take FULL side lengths, not
half-widths:

```python
box = siren.geometry.Box(widths=(x, y, z), center=(cx, cy, cz))
```

Here `x`, `y`, `z` are the full extents of the box along each axis (a cube
with `widths=(2.0, 2.0, 2.0)` spans from -1.0 to +1.0 on each side of its
center). GDML `<box>` `x`/`y`/`z` attributes are also full lengths, not half
lengths.

This is a known footgun: reading a full-width value as a half-width silently
doubles a volume's size. Always use `widths=`/`center=` and confirm which
convention a given geometry constructor documents; do not assume it matches
GDML's `<box>` (some GDML solids, like `<trd>`/`<trap>`, are natively defined
in terms of half lengths, and SIREN's GDML parser halves those on read -- see
the parser's per-solid handling, not this page, for the authoritative list).

The positional form `Box(x, y, z)` would place the box at the ORIGIN,
silently dropping any intended center, so it is rejected: it unconditionally
raises `ConfigurationError`, naming the `widths=`/`center=` fix-it in its
message. Always use the explicit form:

```python
# Raises ConfigurationError: use widths=/center= instead.
box = siren.geometry.Box(1.0, 2.0, 3.0)

# Correct: explicit widths and center.
box = siren.geometry.Box(widths=(1.0, 2.0, 3.0), center=(0.0, 0.0, 5.0))
```

## Frames and coordinates

Positions stored in `InteractionRecord` are in DETECTOR coordinates (the
detector model's local frame), not a beamline or geocentric frame. Loaders
such as `siren.utilities.load_detector` handle the transform from a beamline
frame (e.g. BNB) into detector-local coordinates when composing geometry.

Phase-space measures name the frame a solid angle is defined in:
`Measure.SolidAngleRest()` is uniform solid angle in the parent's rest frame;
`Measure.SolidAngleLab()` is uniform solid angle in the lab frame. These are
not interchangeable -- a sampler built for one and checked (or weighted)
against the other will show a systematic angular bias. `check_closure`'s
frame heuristic flags exactly this: a declared `SolidAngleRest` measure whose
samples actually look isotropic in the lab frame (or vice versa).

## RNG and reproducibility (two-stream contract)

SIREN maintains two independent RNG streams:

- a GENERATION stream, which drives event sampling (the `Injector`'s seed);
- a VALIDATION stream, which drives gauges, probes, and diagnostics (for
  example `check_closure`'s `seed` argument, which defaults to `seed=0` and
  is entirely separate from any generation seed).

Adding or running a diagnostic only draws from the validation stream, so it
can never shift a generation-stream result -- running `check_closure` on a
model, or probing a phase-space mixture, does not perturb any event stream
produced by an `Injector` with the same seed. Fixing a seed makes a
generation run reproducible: same seed, same configuration, same SIREN
version implies an identical event stream.

Serialization does NOT resume an RNG stream mid-sequence. Neither the C++
archive (`save`/`load`) nor pickle carries the random engine, so a reloaded or
unpickled injector begins a fresh generation stream. Pickle re-seeds that
stream from the stored seed, so a fresh injector and an unpickled one produce
the same stream from the start; `load` starts unseeded. Either way, an injector
that had already generated some events continues differently after a round-trip
than it would have run uninterrupted -- the round-trip restarts the stream, it
does not resume it.

## Environment variables

`SIREN_STRICT=1` is the only supported environment toggle. Setting it
escalates SIREN's `DeprecationWarning`s (emitted for deprecated aliases like
`Simulation(n_events=...)`) to hard errors, which is useful for finding
remaining legacy call sites in a codebase before a deprecated form is
removed. This is implemented in `siren`'s package `__init__`. Positional
`Box(x, y, z)` is not a warning case: it raises `ConfigurationError`
unconditionally.

`SIREN_QUIET` is not a supported toggle: nothing in the library reads it.
Setting it has no effect; do not rely on it to suppress warnings or errors.
