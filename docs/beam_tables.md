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

## Complete files, exposure, and identity

Install the optional reader with `pip install 'siren[beam]'`. The reader accepts
filesystem paths and reads complete files. It sums all finite positive POT
entries in every input file, even when a species filter removes all of that
file's rows. Missing trees/exposure, duplicate content, malformed selected
ancestor chains, and negative or nonfinite importance weights raise. Entry
slicing is rejected: a prefix of a file is not an unbiased subsample with the
full file exposure. Perform validated upstream subsampling with recorded
inclusion probabilities and corrected importance weights when needed.

`source_file_id` is the first 52 bits of the file's SHA-256 digest; it is exact
in the native double-valued metadata. `source_entry` is the original entry
index, before filtering. Reversing input-file order changes row order but not
these identities. `source_files` retains each path, full digest, original POT,
and original/selected entry counts. Copies of the same file are rejected, and
52-bit collisions within a read are rejected. Preserve the full manifest
alongside exported results for cross-dataset identity checks.

Both converters automatically preserve `source_file_id`, `source_entry`, and
`source_sampling_probability`. `extra_columns` cannot overwrite them. The
reader sets inclusion probabilities to one; parent/channel selection is a
physical restriction, not random source subsampling.

Timing is independent per file. With `read_time=True`, unavailable timing is
represented by NaN and `time_available=False`; a later untimed file cannot erase
known times from earlier files. Converters reject selected rows with missing
times. Explicit `read_time=False` requests a timing-independent calculation.
Empty or inconsistent selected ancestor chains raise instead of losing rows.

`E` retains the production energy because it is a distinct recorded observable.
The reader does not alter the decay momentum. The converters, spectrum helper,
and summary explicitly construct decay energy from that momentum and resolved
mass. This is an interpretation of dk2nu's separate production/decay fields.
A native table with an explicit `E,px,py,pz,m` preserves all five supplied
components, including an off-shell vector; ingestion does not project it.

For validated complementary DIF/DAR productions, `combine_dif_dar` (or
`read_dif_dar`) partitions rows by kinetic energy **at decay** and rescales DAR
importance weights to DIF reference POT. Each row's physical weight per POT
is unchanged. Both samples must carry the same columns; one-sided metadata is
rejected. The full source manifests and original sample exposures are retained.
This helper requires a physically justified production split; it does not make
overlapping arbitrary beam samples independent.

## Weighting and replay

For N rows, proposal probability q_i and physical row weight w_i, the table's
relative generation density is N q_i and its physical density is N w_i/sum(w).
The absolute normalization is sum(w). An accepted event therefore contributes
w_i/q_i divided by the total generation attempts, multiplied by the remaining
physical factors. Failed attempts remain zero contributions to that estimator.

Pass the table explicitly to both `distributions=[parents]` and
`physical=[parents]` on the Vertex. An omitted physical declaration does not
apply the table's physical weights or absolute exposure normalization. The generic physical-density hook
prevents cancellation when the table has a `weight` column. Other distributions
retain their existing density semantics. Unweighted tables retain their legacy
semantics. Zero total physical weight is rejected.

Native pooling and reweighting require the same row layout: matching row count,
kinematics and metadata at every index. Column order may differ; physical and
proposal weights may differ. Reordered or subset tables raise during weighter
construction. Stable source IDs identify rows across readers; they do not make
arbitrarily reordered native tables share the same sampling measure.

Run the self-contained public interface example:

```sh
python resources/examples/external_sources.py --events 256 --output /tmp/siren-source-example
```

It uses synthetic parents and a generic two-body decay to demonstrate
`dk2nu_to_primary_distribution` → `Vertex` → `Simulation` → `Results.save`.
Its analytic total is 0.1 primaries per exposure. It is an API/normalization
example, not a physical beam or detector-rate prediction. Standard beam decays
and the beam/charm physics examples remain a separate U09 feature.
