# Beam-parent tables

`siren.dk2nu.dk2nu_to_primary_distribution` turns dk2nu-shaped arrays into
a native `PrimaryExternalDistribution`. Positions enter in centimetres,
momenta in GeV, and times in ns. Its `frame` argument rotates momenta and
polarization and transforms positions into detector coordinates, in metres.

The raw `dk2nu_data` arrays must be plain, unmasked NumPy arrays, as returned
by `read_dk2nu`. Masked kinematics, weights, and other raw fields are unsupported;
only `extra_columns` receives explicit selected-row mask validation.

Both this function and `dk2nu_to_csv` accept two optional keyword arguments:

- `masses`: a mapping from integer PDG codes to masses in GeV. These overrides
  take precedence. Otherwise the converters retain their conventional beam
  masses, including the neutral pion, then consult `siren.particles.mass` for
  other species. `particles.define` registrations are supported. An unresolved
  mass raises `ConfigurationError` instead of assuming a pion.
- `extra_columns`: a mapping from column names to numeric arrays aligned with
  the **full input**, before filtering. Values follow the selected rows into
  the event's `interaction_parameters` under the supplied names. Native
  kinematic, weight, time, polarization, and cached-density names are reserved.

For example, when the input carries source identities and inclusion probabilities:

```python
from siren import dk2nu

parents = dk2nu.dk2nu_to_primary_distribution(
    data,
    detector_model,
    parent_pdg=111,
    masses={111: 0.1349768},
    sampling_bias=row_bias,
    extra_columns={
        "beam_event_id": data["event_id"],
        "beam_file_id": data["file_id"],
        "source_sampling_probability": data["sampling_probability"],
    },
)
```

The native table stores numbers as doubles. Integer metadata must lie within
`[-2**53, 2**53]`; larger integer-typed values raise rather than lose source
identity, including integers in mixed numeric sequences. IDs already converted
to floating point cannot be checked for earlier precision loss; pass them as
integers. Other finite floating-point metadata may exceed this ID range. Split
an identifier into exact numeric components upstream when needed. Metadata values
must be finite and unmasked on selected rows; each column must be one-dimensional
and match the input row count. Missing entries on filtered-out rows do not
invalidate the retained rows. Strings and arbitrary Python objects are not
table columns.

Physical row weights are `nimpwt/POT`, evaluated in double precision even for
lower-precision input arrays. Inclusion probabilities are **metadata**:
the caller must already include any source-subsampling correction in `nimpwt`.
The converter does not apply another inverse-probability factor. Species and
invalid-importance-row filtering use the same mask for kinematics and metadata.

`sampling_bias` changes row selection separately from physical weights. Arrays
align with the full input; callables receive the selected original kinematics,
including the stored **production** energy. Every positive-weight row must
have positive sampling support. Biases must be finite and non-negative with
a finite positive sum. The native generation and physical densities account
for the resulting importance ratio.

CSV export writes **decay-point** on-shell energy, the raw `nimpwt`, normalized
`weight`, available time/polarization, and explicit extra columns with enough
digits to round-trip doubles. It requires finite positive POT. Use
`units_cm=False` before loading a CSV directly into SIREN, which expects metres.
`output_path` must be a filesystem path (string or `os.PathLike`); file-like
objects are unsupported. Output is always plain UTF-8 text with LF newlines,
including paths ending in `.gz` or `.bz2`, because the native CSV loader does
not decompress files.
The legacy `position_transform` callback transforms positions only; rotate the
input momentum and polarization arrays too if changing their spatial basis.
For a complete frame conversion, prefer the direct builder's `frame` argument.
CSV loading samples uniformly; a nonuniform proposal belongs in the direct
builder's `sampling_bias`, not an inert extra CSV column.

Tables and sampling weights persist through native injector archives. Direct
Python pickling of a native distribution is explicitly unsupported. No live
Python metadata or new event-record fields are introduced by these converters.

## Track-segment tables

Segment mode is an explicit opt-in: `PrimaryExternalDistribution(...,
segment_column="length")` or `SetSegmentColumn("length")` names the column
(metres) that turns each row into a straight track segment starting at
`x0/y0/z0` along the `px/py/pz` direction, such as a Geant4 step of a photon
inside a production target. Without the opt-in a column of any name, including
`length`, remains ordinary metadata, and archives written before segment mode
existed always load as point tables. Segment tables require `x0/y0/z0` and
`px/py/pz`, must not carry `x/y/z`, and every length must be finite and
positive; segments shorter than the detector model's 10 micrometre direction
threshold are integrated along the trajectory direction with the ordinary
sector integration, so micrometre steps, including ones crossing a material
boundary, weight to their closed form like longer ones. The
generation density is evaluated on the table's own segment to within a few
ulps of the coordinates, never a fixed absolute distance, so a nanometre
proposal does not claim records beyond its own end point. The other columns
keep their meaning: `E`, `m`, `t0`, `weight`
(physical primaries per exposure on the segment, for example the Geant4 step
weight per POT), and explicit sampling weights for biased row selection.
Unrecognised columns remain interaction parameters, and the length column is
written to each record's parameters so saved events reweight without the table.

The primary process is a scattering (a `CrossSection` on the nuclei the
detector model assigns along the segment) and must use
`weighting=siren.ExternalBounds()`; the Injector rejects `Fixed()` and
`Propagated()` for segment tables and rejects `ExternalBounds()` for point
tables. The interaction vertex is sampled uniformly along the segment
(generation density `1/length` per metre, declared as
`PrimaryPositionLongitudinal`), `InjectionBounds` are the segment end points,
and the weighter supplies the interaction probability and the normalised
position density integrated along the segment. The generation density is
evaluated on the table's own support: it is zero for a record whose start,
direction or vertex is not on the row's segment, so several injectors over
overlapping segments of different lengths pool correctly (an injector whose
proposal does not cover an event is skipped entirely; it is an error only when
no injector covers it). Lengths and bounds come from the table's own row
rather than from the record's column name, so pooled tables may use different
segment-column names. `SetSegmentColumn` is transactional: a rejected update
leaves the distribution exactly as it was.

Pooled tables must share a **row layout**. Records cache the sampled row
index, every table asked to evaluate a record reads that index into its own
rows, and each table's row densities are relative to the uniform measure over
its own rows. The `Weighter` therefore requires every
`PrimaryExternalDistribution` it sees (each injector's injection side and the
primary physical side) to list the same primaries at the same indices: the
same row count, the same columns other than `weight` and either table's
segment-length column, and equal values in those columns row by row
(`RowLayoutMismatch` reports the first difference). Segment lengths, physical
weights, sampling weights and the length column's name may differ, so a
biased sampler paired with a physically weighted evaluator, or overlapping
proposals of different lengths, pool as before. A reordered table, a subset
or superset of rows, or different kinematics at one index raises
`ConfigurationError` when the weighter is configured; pooling such tables
would count only the generating proposal (a factor-of-two bias for a reversed
copy) or mix row measures (an asymptotic 25% excess for a one-row subset of a
two-row table). Disjoint chunks of one table are weighted separately, each
against its own chunk as the physical table. The event weights sum to

    sum_i w_i (1 - exp(-int_i n sigma dl)),

the physical primaries per exposure times their interaction probability on
each segment, evaluated with the **detector model's** materials, not those of
the simulation that produced the table. Targets are chosen at the sampled
vertex from the detector model; a vertex drawn in material without any
configured target is an ordinary `NoTargetsOnPath` miss that counts as an
attempt and is retried, as for a ray missing the fiducial volume. The uniform
proposal is unbiased for any material profile along a segment; only the
variance depends on it. Energy or direction biasing of the injected primaries
is a sampling-weight choice, e.g. `s_i = length_i * b(E_i, angle_i)`; the
`weight` column returns to the physical ensemble automatically.
