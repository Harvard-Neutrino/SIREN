"""Rendering layer over the engine's failure ledger and weight breakdown.

InjectionReport turns an Injector's FailureLedger into a readable attrition
table with per-reason remedies; WeightBreakdown turns a Weighter's
EventWeightBreakdown into a per-vertex weight decomposition that names the
culprit vertex behind an inf/zero/NaN weight. Neither computes physics; both
read structures the engine already produced.
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


class VertexWeightLine:
    """One vertex's factors within a weight breakdown."""

    __slots__ = ("depth", "pdg", "generation", "physical",
                 "channel_densities", "cancelled", "flags")

    def __init__(self, depth, pdg, generation, physical,
                 channel_densities, cancelled, flags):
        self.depth = depth
        self.pdg = pdg
        self.generation = generation
        self.physical = physical
        self.channel_densities = dict(channel_densities)
        self.cancelled = list(cancelled)
        self.flags = list(flags)

    @property
    def particle(self) -> str:
        return _pdg_name(self.pdg)

    @property
    def is_ok(self) -> bool:
        return (math.isfinite(self.physical) and self.physical > 0.0
                and math.isfinite(self.generation) and self.generation > 0.0)

    def __repr__(self):
        flag = "" if self.is_ok else " !!!"
        return ("VertexWeightLine(d={} {} gen={:.3e} phys={:.3e}{})".format(
            self.depth, self.particle, self.generation, self.physical, flag))


class WeightBreakdown:
    """Per-vertex decomposition of one event's weight.

    total is the event weight; vertices are the per-vertex factor lines.
    culprit() names the first vertex whose factors explain a non-finite or
    zero total.
    """

    __slots__ = ("total", "vertices")

    def __init__(self, total: float, vertices: List[VertexWeightLine]):
        self.total = total
        self.vertices = list(vertices)

    @classmethod
    def from_engine(cls, breakdown) -> "WeightBreakdown":
        """Build from an engine EventWeightBreakdown struct."""
        lines = []
        for v in breakdown.vertices:
            lines.append(VertexWeightLine(
                v.depth, v.vertex_pdg, v.generation, v.physical,
                v.channel_densities, v.cancelled, v.flags))
        return cls(breakdown.total, lines)

    def culprit(self) -> Optional[VertexWeightLine]:
        """The first vertex whose factors explain a bad total, or None."""
        if math.isfinite(self.total) and self.total > 0.0:
            return None
        for v in self.vertices:
            if not v.is_ok:
                return v
        return None

    def __str__(self):
        ok = math.isfinite(self.total) and self.total > 0.0
        lines = ["WeightBreakdown: total={:.6e}{}".format(
            self.total, "" if ok else " (unusable)")]
        for v in self.vertices:
            line = "  d={} {:<12} gen={:.3e} phys={:.3e}".format(
                v.depth, v.particle, v.generation, v.physical)
            extras = []
            if v.cancelled:
                extras.append("cancelled=" + ",".join(v.cancelled))
            if v.flags:
                extras.append("flags=" + ";".join(v.flags))
            if not v.is_ok:
                extras.append("!!!")
            if extras:
                line += "  [" + " ".join(extras) + "]"
            lines.append(line)
        culprit = self.culprit()
        if culprit is not None:
            lines.append("  culprit: d={} {} ({})".format(
                culprit.depth, culprit.particle,
                ";".join(culprit.flags) or "no positive support"))
        return "\n".join(lines)
