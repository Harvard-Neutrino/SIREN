"""Tests for the Kleiss-Pittau channel-weight update rule (_kp_update).

Two things are checked:

1. The fixed-point difference between the two update rules, using the
   analytic surrogate W_i(alpha) = c_i / alpha_i (variance contribution
   falls as a channel gets more weight).  The canonical multiplicative
   rule converges to W_i = const (the variance minimum); the memoryless
   rule converges to alpha_i ~ sqrt(W_i), i.e. alpha_i ~ c_i^(1/3), which
   leaves the channels with unequal variance contributions.

2. Both rules leave an importance-sampling integral unbiased (Kleiss-Pittau
   holds for any valid alpha), verified with a small pure-Python two-channel
   Monte Carlo integral with a known answer.
"""

import numpy as np
import pytest

from siren.optimize import _kp_update


# ----------------------------------------------------------------------
# 1. Analytic fixed-point behaviour (deterministic, no RNG)
# ----------------------------------------------------------------------

def _iterate_surrogate(c, rule, n_iter=400):
    """Iterate _kp_update with the surrogate W_i(alpha) = c_i / alpha_i."""
    n = len(c)
    alpha = [1.0 / n] * n
    for _ in range(n_iter):
        W = [c[i] / alpha[i] for i in range(n)]
        new = _kp_update(alpha, W, rule, min_weight=0.0)
        assert new is not None
        alpha = new
    return alpha


def test_canonical_rule_converges_to_equal_variance():
    """alpha_i * sqrt(W_i): fixed point is W_i = const (alpha_i ~ c_i)."""
    c = [1.0, 8.0, 27.0]
    alpha = _iterate_surrogate(c, "alpha_sqrt_W")

    # Every channel ends up contributing equal variance.
    W = [c[i] / alpha[i] for i in range(len(c))]
    assert max(W) / min(W) == pytest.approx(1.0, abs=1e-6)

    # Equivalently, alpha_i is proportional to c_i.
    ratios = [alpha[i] / c[i] for i in range(len(c))]
    assert max(ratios) / min(ratios) == pytest.approx(1.0, rel=1e-4)


def test_memoryless_rule_converges_to_cube_root():
    """sqrt(W_i): fixed point is alpha_i ~ c_i^(1/3), W_i NOT equal."""
    c = [1.0, 8.0, 27.0]
    alpha = _iterate_surrogate(c, "sqrt_W")

    cube = [ci ** (1.0 / 3.0) for ci in c]
    ratios = [alpha[i] / cube[i] for i in range(len(c))]
    assert max(ratios) / min(ratios) == pytest.approx(1.0, rel=1e-3)

    # The memoryless rule leaves a real inter-channel variance imbalance.
    W = [c[i] / alpha[i] for i in range(len(c))]
    assert max(W) / min(W) > 2.0


def test_kp_update_validity_and_errors():
    """Returns normalized non-negative weights; rejects bad rules / inputs."""
    new = _kp_update([0.5, 0.5], [4.0, 1.0], "sqrt_W", min_weight=0.1)
    assert abs(sum(new) - 1.0) < 1e-12
    assert all(w >= 0.1 - 1e-12 for w in new)

    # All-zero contributions -> degenerate -> None (caller keeps old weights).
    assert _kp_update([0.5, 0.5], [0.0, 0.0], "sqrt_W") is None

    with pytest.raises(ValueError):
        _kp_update([0.5, 0.5], [1.0, 1.0], "not_a_rule")


# ----------------------------------------------------------------------
# 2. Empirical unbiasedness of the multichannel estimator
# ----------------------------------------------------------------------
#
# Target f(x) = 2x on [0,1] (integral = 1).  Two channels:
#   g1(x) = 1        (uniform; sample U)
#   g2(x) = 3 x^2    (sample U^(1/3))
# Mixture g(x) = a1 g1 + a2 g2.  The importance estimate mean(f/g) is
# unbiased for ANY valid (a1, a2); the KP update only changes its variance.

def _g(x, a):
    return a[0] * 1.0 + a[1] * 3.0 * x * x


def _draw(a, rng, n):
    pick2 = rng.random(n) < a[1]
    u = rng.random(n)
    return np.where(pick2, u ** (1.0 / 3.0), u)


def _estimate_and_W(a, rng, n):
    x = _draw(a, rng, n)
    f = 2.0 * x
    g = _g(x, a)
    w = f / g
    g1 = np.ones_like(x)
    g2 = 3.0 * x * x
    W = [float((w * w * g1 / g).mean()), float((w * w * g2 / g).mean())]
    return float(w.mean()), W


@pytest.mark.parametrize("rule", ["sqrt_W", "alpha_sqrt_W"])
def test_multichannel_estimator_stays_unbiased(rule):
    rng = np.random.default_rng(20240531)
    truth = 1.0
    a = [0.5, 0.5]
    for _ in range(8):
        est, W = _estimate_and_W(a, rng, n=80000)
        assert est == pytest.approx(truth, abs=0.02)
        new = _kp_update(a, W, rule, min_weight=1e-3)
        assert new is not None
        assert abs(sum(new) - 1.0) < 1e-12
        assert all(v >= 0.0 for v in new)
        a = new
    # Still unbiased after the weights have converged.
    est, _ = _estimate_and_W(a, rng, n=200000)
    assert est == pytest.approx(truth, abs=0.01)


# ----------------------------------------------------------------------
# 3. Turning off a redundant covering channel (directed-channel fallback)
# ----------------------------------------------------------------------
#
# A directed channel in its isotropic fallback samples 1/4pi -- a uniform
# proposal that covers the support but does NOT match a non-isotropic physical
# density f.  Analog: f(x) = 3 x^2 on [0,1]; a perfect physical channel
# (g_phys = f) and a uniform fallback channel (g_fall = 1).  The optimum is
# alpha = [1, 0] (g = f, zero variance).  The memoryless sqrt_W rule cannot get
# there -- near g=f every covering channel has W ~ 1, so alpha ~ sqrt(W) parks
# the fallback at a finite weight; the canonical alpha_sqrt_W rule decays it to
# the floor.  Regression for the fallback diluting the physical channel.

def _draw_phys_fallback(a, rng, n):
    pick_fallback = rng.random(n) < a[1]
    u = rng.random(n)
    # physical channel samples from f (inverse-CDF u^(1/3)); fallback uniform
    return np.where(pick_fallback, u, u ** (1.0 / 3.0))


def _w_and_W_phys_fallback(a, rng, n):
    x = _draw_phys_fallback(a, rng, n)
    f = 3.0 * x * x
    g_phys = 3.0 * x * x          # the perfect physical channel
    g_fall = np.ones_like(x)      # the isotropic-fallback analog
    g = a[0] * g_phys + a[1] * g_fall
    w = f / g
    W = [float((w * w * g_phys / g).mean()), float((w * w * g_fall / g).mean())]
    return w, W


def _optimize_phys_fallback(rule, n_iter=60, batch=40000, damping=0.5,
                            min_w=1e-3, seed=1):
    rng = np.random.default_rng(seed)
    a = [0.5, 0.5]   # [physical, fallback]
    for _ in range(n_iter):
        _, W = _w_and_W_phys_fallback(a, rng, batch)
        new = _kp_update(a, W, rule, min_weight=min_w)
        if new is not None:
            a = [damping * new[i] + (1.0 - damping) * a[i] for i in range(2)]
    w, _ = _w_and_W_phys_fallback(a, rng, batch)
    return a, float(w.std() / w.mean())


def test_canonical_rule_turns_off_redundant_fallback():
    a_mem, cv_mem = _optimize_phys_fallback("sqrt_W")
    a_can, cv_can = _optimize_phys_fallback("alpha_sqrt_W")

    # memoryless sqrt_W cannot turn the fallback off -- it parks it high
    assert a_mem[1] > 0.3
    # canonical alpha_sqrt_W decays the fallback toward the floor
    assert a_can[1] < 0.05
    assert a_can[0] > 0.95
    # ... recovering a much lower weight variance (the fallback was diluting f)
    assert cv_can < 0.3 * cv_mem
