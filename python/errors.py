"""Typed error vocabulary for SIREN's high-level interface.

Re-exports the pybind-translated C++ exceptions (registered on
``siren.utilities`` with a ``RuntimeError`` base) under one namespace and adds
the pure-Python errors the Layer-2 surface raises. Every name here derives from
``RuntimeError`` (or ``Exception``) so existing ``pytest.raises(RuntimeError)``
call sites keep matching.
"""

from . import utilities as _utilities

# ---- C++ exceptions surfaced through pybind (siren.utilities) ----
ConfigurationError = _utilities.ConfigurationError
MeasureCompatibilityError = _utilities.MeasureCompatibilityError
WeightCalculationError = _utilities.WeightCalculationError
AddProcessFailure = _utilities.AddProcessFailure
SecondaryProcessFailure = _utilities.SecondaryProcessFailure
InjectionFailure = _utilities.InjectionFailure
PythonImplementationError = _utilities.PythonImplementationError


# ---- Pure-Python errors for the Layer-2 surface ----
class ClosureError(RuntimeError):
    """A sampler and its density disagree beyond tolerance.

    Raised by the closure gauge when a model's FinalStateProbability does not
    match the distribution its SampleFinalState draws from.
    """


class NotSerializableError(RuntimeError):
    """An injector or weighter cannot be round-tripped to disk.

    Raised by ``Injector.save`` when the configuration holds objects that do
    not survive serialization: Python trampoline-derived interactions or
    distributions, or a non-Propagated weighting mode / per-signature phase
    space (neither is archived by the C++ process, so a reloaded injector would
    silently carry the wrong physics).
    """

    def __init__(self, message, offenders=None):
        super().__init__(message)
        self.offenders = list(offenders) if offenders is not None else []


class InjectionShortfall(RuntimeError):
    """Generation produced fewer successes than requested.

    Carries the ``InjectionReport`` describing where attempts were lost so the
    caller can act on it. Raised (or, under ``on_shortfall='warn'``, attached
    to a warning) when the attempt budget is exhausted or efficiency falls
    below ``min_efficiency`` before the requested count is reached.
    """

    def __init__(self, message, report=None):
        super().__init__(message)
        self.report = report
