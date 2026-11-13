# Random-engine archives

`SIREN_random` archives preserve the engine position so supported reloads resume
the stream. They do not promise identical event generation across platforms:
standard-library distributions and other numerical components may differ.

New archives use version 2 and identify the `std::mt19937_64` text encoding as
libc++ or libstdc++. Loading an archive from the other library fails explicitly;
use the originating library to resume that run. This applies to RNGs embedded
in injector archives and pickle snapshots as well as standalone RNG archives.
No fallback silently restarts the engine or substitutes the caller's seed.

Version 1 lacks an encoding identifier. Its strict state-word count distinguishes
the supported libc++ and libstdc++ layouts, so same-library archives remain
readable while the other layout is refused. States must contain classic-locale
unsigned decimal words, with no trailing data and a valid position where stored.
Malformed and unsupported archives leave the destination RNG unchanged.
Old version-0 archives contain only a seed and retain their historical
restart-from-seed behavior. New readers support versions 0, 1 and 2; older
readers explicitly reject version 2.

Cross-library continuation requires a separately specified portable state
format. The current format deliberately refuses that operation. Existing table
fixtures with no stored RNG use the injector's documented caller-engine fallback;
they do not establish portable RNG replay.
