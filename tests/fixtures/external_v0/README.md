# External distribution version-0 archives

Written with published U05 `53c64638ee0f0ecf018d23e708668c6796448200`,
Python 3.13, macOS arm64, native cereal BinaryOutputArchive. These little-endian
64-bit fixtures contain two initial-position rows with inert metadata weights
`[2, 2]` or `[-1, 2]`, plus a sphere vertex sampler and HNLDipoleDecay.
They contain no model tables, event samples, or private inputs. The regression
checks load and resave without interpreting the old metadata as physical weights.
Generation script and matching runtime provenance are retained in
`artifacts/siren-stack-v2-u06-20261025/legacy_weight_archive.py` in SBN_PRISM.
