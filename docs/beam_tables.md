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
