"""Pilot-only channel optimization with an isolated production configuration."""
from dataclasses import dataclass, field
import copy
import math
import operator
import pickle

from .errors import ConfigurationError, WeightCalculationError


@dataclass(frozen=True)
class Plan:
    """Settings for Simulation.run(optimize=Plan(...)).

    events counts attempts per pilot round. Pilot and production seeds must
    differ. A positive min_weight preserves every component's support.
    """
    rounds: int = 5
    events: int = 2000
    seed: int = 1729
    rule: str = 'alpha_sqrt_W'
    damping: float = 0.5
    min_weight: float = 0.001
    metric: object = None
    recurse: bool = True

    def __post_init__(self):
        for name in ('rounds', 'events', 'seed'):
            value = getattr(self, name)
            if isinstance(value, bool):
                raise ValueError(name + ' must be an integer')
            try:
                value = operator.index(value)
            except TypeError as exc:
                raise ValueError(name + ' must be an integer') from exc
            if value < (0 if name == 'seed' else 1) or value > 2**32 - 1:
                raise ValueError(name + ' is outside its supported range')
        if self.rule not in ('alpha_sqrt_W', 'sqrt_W'):
            raise ValueError('Unknown tuning rule')
        if not math.isfinite(self.damping) or not 0 <= self.damping <= 1:
            raise ValueError('damping must be between zero and one')
        if not math.isfinite(self.min_weight) or not 0 < self.min_weight <= 1:
            raise ValueError('min_weight must be positive and at most one')
        if self.metric is not None and not callable(self.metric):
            raise TypeError('metric must be callable or None')
        if not isinstance(self.recurse, bool):
            raise TypeError('recurse must be boolean')


@dataclass
class OptimizationReport:
    """Pilot outcomes and the fresh production pair. ESS includes zero attempts.

    converged describes only the last two proposal updates, not a validated
    integral or statistical convergence. Input objects are never modified.
    """
    injector: object
    weighter: object
    seed: int
    ess_trajectory: list = field(default_factory=list)
    alpha_history: list = field(default_factory=list)
    attempts: list = field(default_factory=list)
    successes: list = field(default_factory=list)
    failures: list = field(default_factory=list)
    converged: bool = False


def _mixtures(engine):
    processes = [engine.GetPrimaryProcess()] + list(engine.GetSecondaryProcesses())
    roots, seen = [], set()
    for process in processes:
        for mixture in process.GetPhaseSpaceMap().values():
            if mixture not in seen:
                seen.add(mixture)
                roots.append(mixture)
    return roots


def _levels(roots, recurse=True):
    from .injection import NestedMixtureChannel
    out, seen = [], set()
    def visit(mc):
        if mc in seen:
            raise ConfigurationError('Tuning requires unshared acyclic mixtures')
        seen.add(mc)
        out.append(mc)
        if recurse:
            for ch in mc.channels:
                if isinstance(ch, NestedMixtureChannel):
                    visit(ch.mixture)
    for mc in roots:
        visit(mc)
    return out


def _ess(weights):
    if not weights:
        return 0.0
    scale = max(abs(w) for w in weights)
    if scale == 0:
        return 0.0
    scaled = [w / scale for w in weights]
    return math.fsum(scaled)**2 / math.fsum(w*w for w in scaled)


def _report(injector):
    from .dataclasses import InteractionTree
    report = injector.report()
    borrowed = report.last_failed_tree
    owned = InteractionTree()
    if borrowed is not None:
        owned.tree = list(borrowed.tree)
        owned.header = borrowed.header
    report.last_failed_tree = owned
    return report


def _clone_pair(injector, weighter):
    from .Injector import Injector
    from .Weighter import Weighter
    if type(injector) is not Injector or type(weighter) is not Weighter:
        raise TypeError('tune requires the public Injector and Weighter classes')
    sources = weighter.injectors or []
    if len(sources) != 1 or (sources[0] is not injector and sources[0] is not injector.engine):
        raise ConfigurationError('tune requires a weighter for the same single injector')
    # Native snapshots refuse Python model/callback state that cannot be retained.
    new_injector = pickle.loads(pickle.dumps(injector))
    source = copy.copy(weighter)
    source.event_factor = None
    frozen = pickle.loads(pickle.dumps(source))
    factor = pickle.loads(pickle.dumps(weighter.event_factor))
    new_weighter = Weighter(new_injector,
        primary_physical=frozen.primary_physical_distributions,
        secondary_physical=frozen.secondary_physical_distributions,
        overrides={'detector_model': frozen.detector_model,
                   'primary_interactions': frozen.primary_interactions,
                   'secondary_interactions': frozen.secondary_interactions},
        event_factor=factor)
    return new_injector, new_weighter


def tune(injector, weighter, *, events=2000, rounds=5, seed=1729,
         rule='alpha_sqrt_W', damping=0.5, min_weight=0.001,
         metric=None, recurse=True, on_failure='raise'):
    """Return a tuned, independent production pair and pilot diagnostics.

    The input pair must have fresh counters and native-serializable physics.
    Each round uses exactly events attempts. Failed attempts remain zero
    contributions; unexpected generation failures raise by default. Invalid
    weights and non-finite metrics always raise. Metrics may be signed.

    Generate production with report.injector/report.weighter, or use
    Simulation.run(optimize=Plan(...)). No pilot events enter Results.
    """
    Plan(rounds, events, seed, rule, damping, min_weight, metric, recurse)
    if on_failure not in ('raise', 'retry'):
        raise ValueError("on_failure must be 'raise' or 'retry'")
    if injector.injection_attempts or injector.injected_events:
        raise ConfigurationError('tune requires fresh counters; create a fresh injector')
    if injector.seed is None or injector.seed == seed:
        raise ConfigurationError('Set distinct explicit pilot and production seeds')
    pilot, target = _clone_pair(injector, weighter)
    original_rng = pilot.engine.GetRandom()
    pilot.seed = seed
    roots = _mixtures(pilot.engine)
    levels = _levels(roots, recurse)
    for mc in levels:
        if len(mc.weights) * min_weight > 1:
            raise ValueError('min_weight exceeds the mixture simplex')
        if any(w <= 0 for w in mc.weights):
            raise ConfigurationError('Tuning requires initially positive component weights')
    report = OptimizationReport(pilot, target, seed)
    for _ in range(rounds if roots else 0):
        pilot.reset(events)
        for mc in roots:
            mc.ResetAccumulators(recurse)
        trees = []
        for _ in range(events):
            tree = pilot.engine.GenerateEvent()
            if tree.tree:
                trees.append(tree)
            elif on_failure == 'raise':
                _report(pilot)._raise_for_failure()
        # Fix the denominator after the complete pilot batch, including misses.
        weights = []
        for tree in trees:
            w = target(tree)
            if metric is not None:
                w = metric(tree, w)
            try:
                w = float(w)
            except (ValueError, TypeError, OverflowError) as exc:
                raise WeightCalculationError('Metric must return a finite real scalar') from exc
            if not math.isfinite(w):
                raise WeightCalculationError('Metric must return a finite real scalar')
            weights.append(w)
            for datum in tree.tree:
                process = (pilot.engine.GetPrimaryProcess() if datum.is_root()
                           else pilot.engine.GetSecondaryProcessMap()[datum.record.signature.primary_type])
                if process.HasPhaseSpace(datum.record.signature):
                    process.GetPhaseSpace(datum.record.signature).Accumulate(
                        pilot.detector_model, datum.record, w, recurse)
        for mc in roots:
            mc.UpdateWeights(rule, damping, min_weight, recurse)
        report.ess_trajectory.append(_ess(weights))
        report.alpha_history.append([list(mc.weights) for mc in levels])
        report.attempts.append(pilot.injection_attempts)
        report.successes.append(len(trees))
        report.failures.append(_report(pilot))
    if len(report.alpha_history) >= 2:
        report.converged = all(abs(a-b) < 1e-3
            for before, after in zip(report.alpha_history[-2], report.alpha_history[-1])
            for a,b in zip(before, after))
    pilot.reset(injector.number_of_events)
    pilot.seed = injector.seed
    pilot.engine.SetRandom(original_rng)
    return report
