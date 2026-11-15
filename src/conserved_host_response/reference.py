"""Frozen-reference, discovery-signed mean-z scoring for the CRRM388 study.

The v0.1.0 score module is intentionally unchanged. This module implements the
current study with explicit reference fitting, projection and eligibility.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
from typing import Mapping, Sequence

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class CoveragePolicy:
    total_fraction: float = 0.70
    direction_fraction: float = 0.50
    subset_fraction: float = 0.50
    subset_minimum: int = 5


def _expression(x: pd.DataFrame) -> pd.DataFrame:
    if x.index.has_duplicates or x.columns.has_duplicates:
        raise ValueError("Gene and observation identifiers must be unique; collapse probes upstream")
    if not len(x.columns):
        raise ValueError("Expression has no observations")
    return x.astype(float)


class FrozenReference:
    """Store reference parameters once; projecting observations never refits them."""

    schema = "crrm388.signed-mean-z.reference.v1"

    def __init__(self, parameters: pd.DataFrame, metadata: dict):
        self.parameters = parameters.copy(deep=True)
        self.metadata = dict(metadata)

    @classmethod
    def fit(cls, expression: pd.DataFrame, reference_columns: Sequence[str], *,
            units: pd.Series | None = None, sd_threshold: float = 1e-8,
            nonzero_support: pd.Series | None = None) -> "FrozenReference":
        x = _expression(expression)
        columns = list(reference_columns)
        if len(columns) < 2 or len(columns) != len(set(columns)):
            raise ValueError("At least two distinct reference profiles are required")
        a = x.loc[:, columns]
        complete = np.isfinite(a.to_numpy()).all(axis=1)
        if units is None:
            means = a.mean(axis=1, skipna=False)
            sd = a.std(axis=1, ddof=1, skipna=False)
            mode = "equal-profile means and sample SD (ddof=1)"
            weights = np.repeat(1 / len(columns), len(columns))
        else:
            unit = units.reindex(columns)
            if unit.isna().any() or unit.nunique() < 2:
                raise ValueError("Cellular references require identified profiles from at least two biological units")
            counts = unit.value_counts()
            weights = np.array([1 / (unit.nunique() * counts[u]) for u in unit])
            values = a.to_numpy()
            means = pd.Series(values @ weights, index=a.index)
            sd = pd.Series(np.sqrt(((values - means.to_numpy()[:, None]) ** 2) @ weights), index=a.index)
            mode = "equal biological-unit total weight; equal profiles within unit; weighted population SD"
        supported = pd.Series(True, index=a.index)
        if nonzero_support is not None:
            supported = nonzero_support.reindex(a.index).fillna(False).astype(bool)
        eligible = complete & np.isfinite(sd) & sd.gt(sd_threshold) & supported
        reason = np.select([~complete, ~np.isfinite(sd), ~sd.gt(sd_threshold), ~supported],
                           ["nonfinite reference measurement", "nonfinite SD", "reference SD below threshold", "no nonzero count support"], default="scorable")
        params = pd.DataFrame({"mean": means, "sd": sd, "eligible": eligible, "reason": reason})
        params.index.name = "human_ensembl"
        return cls(params, {"reference_profiles": columns, "reference_n": len(columns),
                            "mode": mode, "sd_threshold": sd_threshold,
                            "profile_weights": dict(zip(columns, weights.tolist())),
                            "nonzero_count_support_required": nonzero_support is not None})

    def to_dict(self) -> dict:
        rows = self.parameters.reset_index().astype(object)
        rows = rows.where(pd.notna(rows), None)
        return {"schema": self.schema, "metadata": self.metadata, "genes": rows.to_dict("records")}

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), ensure_ascii=False, indent=2, allow_nan=False) + "\n")

    @classmethod
    def load(cls, path: str | Path) -> "FrozenReference":
        value = json.loads(Path(path).read_text())
        if value.get("schema") != cls.schema:
            raise ValueError("Not a CRRM388 frozen-reference artifact")
        return cls(pd.DataFrame(value["genes"]).set_index("human_ensembl"), value["metadata"])

    def fingerprint(self) -> str:
        content = json.dumps(self.to_dict(), sort_keys=True, allow_nan=False).encode()
        return hashlib.sha256(content).hexdigest()

    def project(self, expression: pd.DataFrame) -> pd.DataFrame:
        x = _expression(expression)
        ids = [g for g in self.parameters.index[self.parameters.eligible] if g in x.index]
        values = x.loc[ids]
        if not np.isfinite(values.to_numpy()).all():
            raise ValueError("Projection measurements must be finite; absent features are omitted, never filled with zero")
        p = self.parameters.loc[ids]
        return values.sub(p["mean"], axis=0).div(p["sd"], axis=0)


def score(expression: pd.DataFrame, reference: FrozenReference, genes: pd.DataFrame,
          subsets: Mapping[str, Sequence[str]] | None = None, *,
          omitted: Sequence[str] = (), policy: CoveragePolicy = CoveragePolicy()):
    """Return scores and explicit coverage/status, using unique core genes once.

    Subsets may overlap. Their means are never added to form the total. Missing
    or reference-invariant genes do not enter the actual scored denominator.
    """
    if genes.index.has_duplicates or not len(genes):
        raise ValueError("Core definition requires unique Ensembl identifiers")
    if not genes.frozen_direction.isin(["up", "down"]).all():
        raise ValueError("Discovery directions must be frozen up/down labels")
    subsets = dict(subsets or {})
    if any(len(v) != len(set(v)) for v in subsets.values()):
        raise ValueError("Each subset must contain unique identifiers")
    if any(not set(v).issubset(genes.index) for v in subsets.values()):
        raise ValueError("Annotation subsets must be members of the frozen core")
    z = reference.project(expression)
    eligible = [g for g in genes.index if g in z.index and g not in set(omitted)]
    signs = genes.frozen_direction.map({"up": 1.0, "down": -1.0})
    signed = z.loc[eligible].mul(signs.loc[eligible], axis=0)
    result = pd.DataFrame(index=expression.columns)
    result.index.name = "sample_id"
    frozen_up = genes.index[genes.frozen_direction.eq("up")]
    frozen_down = genes.index[genes.frozen_direction.eq("down")]
    up = [g for g in eligible if signs[g] == 1]
    down = [g for g in eligible if signs[g] == -1]
    passing = (len(eligible) / len(genes) >= policy.total_fraction and
               (not len(frozen_up) or len(up) / len(frozen_up) >= policy.direction_fraction) and
               (not len(frozen_down) or len(down) / len(frozen_down) >= policy.direction_fraction))
    result["CRRM388"] = signed.mean(axis=0) if passing else np.nan
    result["Up contribution"] = signed.loc[up].sum(axis=0) / len(eligible) if passing else np.nan
    result["Down contribution"] = signed.loc[down].sum(axis=0) / len(eligible) if passing else np.nan
    coverage = [{"set": "CRRM388", "nominal": len(genes), "scored": len(eligible),
                 "up_scored": len(up), "down_scored": len(down), "eligible": passing,
                 "status": "OK" if passing else "UNESTIMABLE_COVERAGE"}]
    for name, ids in subsets.items():
        hit = [g for g in ids if g in eligible]
        passing = len(hit) >= policy.subset_minimum and len(hit) / max(len(ids), 1) >= policy.subset_fraction
        result[name] = signed.loc[hit].mean(axis=0) if passing else np.nan
        coverage.append({"set": name, "nominal": len(ids), "scored": len(hit), "eligible": passing,
                         "status": "OK" if passing else "UNESTIMABLE_COVERAGE"})
    return result, pd.DataFrame(coverage), signed
