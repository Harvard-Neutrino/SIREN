# External distribution version-0 archives

Written with published U05 `53c64638ee0f0ecf018d23e708668c6796448200`,
Python 3.13, macOS arm64, native cereal BinaryOutputArchive. These little-endian
64-bit fixtures contain two initial-position rows with inert metadata weights
`[2, 2]` or `[-1, 2]`, plus a sphere vertex sampler and HNLDipoleDecay.
They contain no model tables, event samples, or private inputs. The regression
checks load and resave without interpreting the old metadata as physical weights.
The generator clears the injector's RNG before saving. Its native engine-state
text differs between libc++ and libstdc++, so embedding it would prevent the
table compatibility test from reaching the external distribution on Linux.
The loader uses the caller's RNG when the archived pointer is null. The test
checks that fallback, then resaves with the local engine. These fixtures cover
table-schema compatibility, not cross-standard-library RNG continuation.
Regenerate using the unmodified published U05 runtime, not the U06 runtime:

```sh
python tests/fixtures/external_v0/generate.py /tmp/siren-external-v0-new
```

`generate.py` writes the CSV input and native archive for both cases. Confirm
that the selected Python imports SIREN built from the U05 revision above before
running it. Current archives intentionally use a newer version. The fixture
regression loads these original bytes directly, then resaves with the current
version and verifies the old metadata semantics remain intact.
