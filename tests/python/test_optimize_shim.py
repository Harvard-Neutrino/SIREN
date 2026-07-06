"""siren.optimize is a deprecated shim re-exporting from siren.tune."""
import warnings

import pytest


def test_optimize_import_reexports_tune_symbols():
    import siren.optimize as _optimize
    import siren.tune as _tune

    for name in ("_kp_update", "group_directed_channels",
                 "optimize_multichannel_weights", "optimize_chain_weights"):
        assert getattr(_optimize, name) is getattr(_tune, name)


def test_optimize_attribute_access_warns():
    import siren.optimize as _optimize

    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        with pytest.raises(DeprecationWarning):
            _optimize.optimize_chain_weights


def test_from_optimize_import_still_works():
    from siren.optimize import optimize_multichannel_weights
    from siren.tune import optimize_multichannel_weights as tune_fn

    assert optimize_multichannel_weights is tune_fn
