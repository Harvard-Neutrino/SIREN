"""Supported native directed proposals retain Results identity and replay."""
import pickle
import numpy as np
import pytest
import siren
from test_authoring_creation_survival import facade


def run(tmp_path, seed, fraction=.7, target_z=80.):
    gen = facade(tmp_path)
    gen.seed = seed
    target = siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0, 0, target_z)), 5, 5, 5)
    gen.primary.kinematics = siren.Directed('Gamma', target, fraction=fraction)
    return siren.generate(gen, siren.injection.Weighter(gen), events=64,
                          on_failure='raise', on_shortfall='raise')


def test_directed_results_snapshot_pool_and_reweight(tmp_path):
    a, b = run(tmp_path, 11), run(tmp_path, 12)
    pooled = siren.Results.merge([a, b])
    assert pooled.attempts == a.attempts + b.attempts
    copied = pickle.loads(pickle.dumps(pooled))
    np.testing.assert_array_equal(copied.weights, pooled.weights)
    np.testing.assert_allclose(copied.reweight().weights, pooled.weights, rtol=2e-12, atol=0)
    pooled.save(tmp_path/'directed', snapshot=True)
    loaded = siren.Results.load(tmp_path/'directed')
    np.testing.assert_array_equal(loaded.weights, pooled.weights)
    np.testing.assert_allclose(loaded.reweight().weights, pooled.weights, rtol=2e-12, atol=0)


@pytest.mark.parametrize('change', [{'fraction': .4}, {'target_z': 100.}])
def test_directed_results_distinguish_changed_proposal(tmp_path, change):
    a, b = run(tmp_path, 11), run(tmp_path, 12, **change)
    with pytest.raises(siren.errors.ConfigurationError, match='differing config'):
        siren.Results.merge([a, b])
