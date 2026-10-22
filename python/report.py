"""Generation outcome carried by strict failures and shortfall exceptions.

Reports read the native failure ledger without evaluating models.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional

from . import dataclasses as _dataclasses


def _reason_name(reason) -> str:
    return getattr(reason, "name", str(reason))


def _pdg_name(pdg: int) -> str:
    """Human name for a PDG code, falling back to the integer."""
    try:
        return _dataclasses.ParticleType(int(pdg)).name
    except (ValueError, TypeError):
        return str(pdg)


# One-line remedy per FailureReason. Keyed by the enum member name so a new
# reason without a hint still renders (its name plus a generic line).
_REASON_HINTS: Dict[str, str] = {
    "NoTargetsOnPath":
        "the sampled path crosses no interaction target; widen the vertex "
        "distribution's region or check the detector materials",
    "NoPathThroughVolume":
        "the primary direction never enters the injection volume; check the "
        "direction distribution and the volume placement",
    "NoColumnDepthSolution":
        "no column-depth solution along the path; the energy or geometry "
        "leaves the vertex sampler with no root",
    "KinematicallyForbidden":
        "the sampled kinematics are below threshold; check masses and the "
        "primary energy range",
    "UnregisteredSecondaryType":
        "a produced secondary has no registered Vertex; add a secondary "
        "Vertex for it or prune it in the expand list",
    "PrimaryVertexFailure":
        "the primary vertex could not be sampled; inspect the primary "
        "position/energy distributions",
    "TopLevelCatch":
        "an unclassified failure escaped to the top level; inspect the "
        "exemplar message",
    "Unspecified":
        "failure reason not tagged; inspect the exemplar message",
}


class VertexAttrition:
    """One (depth, particle, reason) bucket of injection failures."""

    __slots__ = ("depth", "pdg", "reason", "count", "exemplar")

    def __init__(self, depth: int, pdg: int, reason, count: int, exemplar: str):
        self.depth = depth
        self.pdg = pdg
        self.reason = reason
        self.count = count
        self.exemplar = exemplar

    @property
    def reason_name(self) -> str:
        return _reason_name(self.reason)

    @property
    def particle(self) -> str:
        return _pdg_name(self.pdg)

    @property
    def hint(self) -> str:
        return _REASON_HINTS.get(self.reason_name, "inspect the exemplar message")

    def __repr__(self):
        return ("VertexAttrition(depth={} particle={} reason={} count={})"
                .format(self.depth, self.particle, self.reason_name, self.count))


class InjectionReport:
    """Readable view of an injection run's outcome.

    attempts/successes are the realized counts; by_vertex is the attrition
    breakdown; last_failed_tree is the most recent partial tree on a primary
    failure (may be empty). dominant() names the reason that lost the most
    attempts.
    """

    __slots__ = ("attempts", "successes", "by_vertex", "last_failed_tree")

    def __init__(self, attempts: int, successes: int,
                 by_vertex: List[VertexAttrition],
                 last_failed_tree=None):
        self.attempts = attempts
        self.successes = successes
        self.by_vertex = list(by_vertex)
        self.last_failed_tree = last_failed_tree

    @classmethod
    def from_ledger(cls, ledger, *, attempts: int, successes: int,
                    last_failed_tree=None) -> "InjectionReport":
        """Build a report from an engine FailureLedger.

        ledger.entries() maps (depth, parent_pdg, FailureReason) to
        (count, exemplar).
        """
        buckets = []
        for (depth, pdg, reason), (count, exemplar) in ledger.entries().items():
            buckets.append(VertexAttrition(depth, pdg, reason, count, exemplar))
        buckets.sort(key=lambda b: b.count, reverse=True)
        return cls(attempts, successes, buckets, last_failed_tree)

    def _raise_for_failure(self):
        """Apply the strict policy after a failed generation attempt."""
        from .injection import FailureReason
        from .errors import GenerationFailure
        allowed = (FailureReason.NoPathThroughVolume, FailureReason.NoTargetsOnPath)
        failures = [bucket for bucket in self.by_vertex if bucket.reason not in allowed]
        if failures or not self.by_vertex:
            details = "\n".join(
                "{} at depth {}, parent {}: {}".format(
                    bucket.reason_name, bucket.depth, bucket.pdg, bucket.exemplar)
                for bucket in failures)
            raise GenerationFailure(
                "injection failed: " + (details or "no failure diagnosis"), report=self)

    @property
    def failures(self) -> int:
        return sum(b.count for b in self.by_vertex)

    @property
    def efficiency(self) -> float:
        if self.attempts <= 0:
            return float("nan")
        return self.successes / self.attempts

    def dominant(self) -> Optional[VertexAttrition]:
        """The attrition bucket that lost the most attempts, or None."""
        if not self.by_vertex:
            return None
        return max(self.by_vertex, key=lambda b: b.count)

    def __str__(self):
        lines = []
        eff = self.efficiency
        eff_str = "{:.1%}".format(eff) if math.isfinite(eff) else "n/a"
        lines.append("InjectionReport: {}/{} succeeded (efficiency {})".format(
            self.successes, self.attempts, eff_str))
        if not self.by_vertex:
            lines.append("  no failures recorded")
            return "\n".join(lines)
        lines.append("  {:<6} {:<14} {:<26} {:>8}".format(
            "depth", "particle", "reason", "count"))
        for b in self.by_vertex:
            lines.append("  {:<6} {:<14} {:<26} {:>8}".format(
                b.depth, b.particle, b.reason_name, b.count))
        dom = self.dominant()
        if dom is not None:
            lines.append("  dominant: {} ({}x) -- {}".format(
                dom.reason_name, dom.count, dom.hint))
        return "\n".join(lines)
