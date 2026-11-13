import math
import pickle
import numpy as np
import pytest
import siren
from test_authoring_creation_survival import facade
from test_decay_channels import channel, P


def test_results_selected_survival_roundtrip_pooling_and_reweight(tmp_path):
    a, b = facade(tmp_path), facade(tmp_path)
    b.seed = 43
    ra = siren.generate(a,siren.injection.Weighter(a),events=8,on_failure='raise')
    rb = siren.generate(b,siren.injection.Weighter(b),events=8,on_failure='raise')
    expected = 0.2 * -math.expm1(-1.0) * math.exp(-5.0)
    assert sum(ra.weights) == pytest.approx(expected,rel=2e-12,abs=0)
    pooled=siren.Results.merge([ra,rb])
    restored=pickle.loads(pickle.dumps(pooled))
    assert sum(restored.weights) == pytest.approx(expected,rel=2e-12,abs=0)
    np.testing.assert_allclose(restored.reweight().weights,restored.weights,rtol=2e-12,atol=0)
    width=5/10*3*1.973e-16
    changed=restored.reweight(primary_interactions=[channel([P.NuMu,P.Gamma],.4*width,width),channel([P.NuE,P.Gamma],.3*width,width)])
    np.testing.assert_allclose(changed.weights,2*restored.weights,rtol=2e-12,atol=0)
    np.testing.assert_allclose(changed.reweight().weights,changed.weights,rtol=2e-12,atol=0)


@pytest.mark.parametrize('change',['survival','selection'])
def test_results_reject_different_survival_or_selected_channels(tmp_path,change):
    a,b=facade(tmp_path),facade(tmp_path)
    b.seed = 43
    if change=='survival':
        b.primary.weighting=siren.Propagated()
    else:
        b.primary.decay_channels=None
    ra=siren.generate(a,siren.injection.Weighter(a),events=8,on_failure='raise')
    rb=siren.generate(b,siren.injection.Weighter(b),events=8,on_failure='raise')
    with pytest.raises(siren.errors.ConfigurationError,match='differing config'):
        siren.Results.merge([ra,rb])
