"""Deprecated shim: the implementation lives in ``siren.tune``."""

import warnings

_EXPORTED = (
    "_kp_update",
    "group_directed_channels",
    "optimize_multichannel_weights",
    "optimize_chain_weights",
)


def __getattr__(name):
    if name in _EXPORTED:
        from . import tune as _tune
        warnings.warn(
            "siren.optimize is deprecated; use siren.tune",
            DeprecationWarning, stacklevel=2)
        return getattr(_tune, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
