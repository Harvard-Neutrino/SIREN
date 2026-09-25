"""Frozen, full-support importance allocation over an identified source table.

Pilot means describe physical scored contribution conditional on a source row.
The native distribution owns the physical/proposal correction in production;
``sample`` returns that ratio explicitly for pilot/benchmark use. A zero pilot
estimate never removes a physically populated row. Pilot and production streams
must be independent; changing a sampler alone need not invalidate the means.
"""

import hashlib
import json
from pathlib import Path

import numpy as np


def canonical_metadata(metadata):
    """One JSON representation for construction, validation and persistence.

    NumPy scalars become Python scalars. Keys use JSON's string convention;
    reject collisions rather than silently lose a definition field.
    """
    if isinstance(metadata, np.generic):
        return canonical_metadata(metadata.item())
    if isinstance(metadata, dict):
        result = {}
        for key, value in metadata.items():
            if isinstance(key, np.generic):
                key = key.item()
            key = next(iter(json.loads(json.dumps({key: None}, allow_nan=False))))
            if key in result:
                raise ValueError("Source metadata keys collide after JSON conversion")
            result[key] = canonical_metadata(value)
        return result
    if isinstance(metadata, (tuple, list)):
        return [canonical_metadata(v) for v in metadata]
    return json.loads(json.dumps(metadata, allow_nan=False))


def fingerprint(metadata):
    return hashlib.sha256(
        json.dumps(
            canonical_metadata(metadata), sort_keys=True, allow_nan=False, separators=(",", ":")
        ).encode()
    ).hexdigest()


class SourceImportanceTable:
    def __init__(self, physical_weights, mean_scores, metadata, baseline_fraction=0.05):
        weights = np.asarray(physical_weights, dtype=float)
        means = np.asarray(mean_scores, dtype=float)
        if weights.ndim != 1 or not len(weights) or means.shape != weights.shape:
            raise ValueError(
                "Source weights and pilot means must be equal nonempty vectors"
            )
        if (
            not np.isfinite(weights).all()
            or not np.isfinite(means).all()
            or np.any(weights < 0)
            or np.any(means < 0)
        ):
            raise ValueError("Nonfinite or negative source weights/pilot scores")
        if (
            not 0 < baseline_fraction <= 1
            or not np.isfinite(weights.sum())
            or weights.sum() <= 0
        ):
            raise ValueError(
                "Positive source exposure and baseline fraction in (0,1] required"
            )
        self.metadata = canonical_metadata(metadata)
        for key in ["source", "physics", "geometry", "scoring"]:
            if key not in self.metadata:
                raise ValueError(f"Missing source-table definition: {key}")
        self.baseline_fraction = float(baseline_fraction)
        self.physical_weights = weights.copy()
        self.mean_scores = means.copy()
        self.physical = weights / weights.sum()
        # Scaling avoids overflow without changing the allocation.
        allocation = self.physical * (means / means.max() if means.max() > 0 else 0)
        importance = (
            allocation / allocation.sum() if allocation.sum() > 0 else self.physical
        )
        self.proposal = (
            self.baseline_fraction * self.physical
            + (1 - self.baseline_fraction) * importance
        )
        for a in [
            self.physical_weights,
            self.mean_scores,
            self.physical,
            self.proposal,
        ]:
            a.flags.writeable = False
        self.definition_hash = fingerprint(self.metadata)

    def sample(self, rng, size):
        indices = rng.choice(len(self.proposal), size=size, p=self.proposal)
        return indices, self.physical[indices] / self.proposal[indices]

    def to_distribution(self, keys, data, *, metadata):
        """Build the native source sampler and physical-density evaluator.

        ``data`` must contain the complete, ordered source rows and a ``weight``
        column in physical primaries per unit exposure. Its normalized weights
        must match the frozen pilot's physical ensemble; the absolute scale may
        differ (for example, counts in the pilot versus counts/POT here).

        Declare the returned distribution on both the injection and physical
        sides of the vertex. The native Weighter owns the normalization and p/q
        correction; do not multiply either into the resulting event weights.
        No row filtering or reordering is performed.
        """
        self.validate_definition(metadata)
        keys = list(keys)
        if len(keys) != len(set(keys)) or "weight" not in keys:
            raise ValueError("Unique source columns including physical 'weight' required")
        rows = np.asarray(data, dtype=float)
        if rows.shape != (len(self.physical), len(keys)) or not np.isfinite(rows).all():
            raise ValueError("Finite complete source rows must match the frozen table")
        weights = rows[:, keys.index("weight")]
        total = weights.sum()
        if np.any(weights < 0) or not np.isfinite(total) or total <= 0:
            raise ValueError("Positive finite physical source normalization required")
        if not np.allclose(weights / total, self.physical, rtol=1e-12, atol=0):
            raise ValueError("Physical row weights differ from the frozen source table")
        from .distributions import PrimaryExternalDistribution

        return PrimaryExternalDistribution(keys, rows.tolist(), self.proposal.tolist())

    def validate_definition(self, metadata):
        if fingerprint(metadata) != self.definition_hash:
            raise ValueError(
                "Source/physics/geometry/scoring definition changed; rebuild or revalidate table"
            )

    def save(self, path):
        """Write to a path (used exactly as given) or a writable binary file."""
        arrays = dict(
            physical_weights=self.physical_weights,
            mean_scores=self.mean_scores, baseline_fraction=self.baseline_fraction,
            metadata=json.dumps(self.metadata, sort_keys=True),
        )
        if hasattr(path, "write"):
            np.savez_compressed(path, **arrays)
            return
        # Opening the exact path avoids NumPy's implicit '.npz' suffix.
        with Path(path).open("wb") as stream:
            np.savez_compressed(stream, **arrays)

    @classmethod
    def load(cls, path, metadata):
        source = path if hasattr(path, "read") else Path(path)
        with np.load(source, allow_pickle=False) as data:
            result = cls(
                data["physical_weights"],
                data["mean_scores"],
                json.loads(str(data["metadata"])),
                float(data["baseline_fraction"]),
            )
        result.validate_definition(metadata)
        return result
