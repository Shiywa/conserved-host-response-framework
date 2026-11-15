"""Current three-model CRRM definition; independent of the legacy anchor CLI."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd

MODELS = ("A549", "mouse", "cotton_rat")

def reconstruct(evidence: pd.DataFrame, *, q_threshold=0.05, effect_threshold=0.58):
    if evidence.human_ensembl.duplicated().any():
        raise ValueError("Evidence requires one row per uniquely projected human gene")
    valid = np.ones(len(evidence), dtype=bool)
    directions = []
    for model in MODELS:
        tested = evidence[model + "_mapped_and_tested"].astype(str).str.lower().isin(["true", "1", "t"])
        fc = pd.to_numeric(evidence[model + "_log2FC"], errors="coerce")
        q = pd.to_numeric(evidence[model + "_FDR"], errors="coerce")
        valid &= tested & np.isfinite(fc) & np.isfinite(q) & q.lt(q_threshold) & fc.abs().gt(effect_threshold)
        directions.append(np.sign(fc))
    valid &= (directions[0] == directions[1]) & (directions[1] == directions[2])
    selected = evidence.loc[valid].copy()
    selected["frozen_direction"] = np.where(directions[0][valid] > 0, "up", "down")
    return selected.sort_values("human_ensembl").reset_index(drop=True)

def verify_definition(genes: pd.DataFrame, subsets: pd.DataFrame):
    if len(genes) != 388 or genes.human_ensembl.nunique() != 388:
        raise ValueError("Current definition must contain 388 unique genes")
    if genes.frozen_direction.value_counts().to_dict() != {"up": 280, "down": 108}:
        raise ValueError("Frozen directions differ from the current study")
    if not set(subsets.human_ensembl).issubset(genes.human_ensembl):
        raise ValueError("Subset contains a gene outside CRRM388")
    if subsets.duplicated(["set_name", "human_ensembl"]).any():
        raise ValueError("Duplicate within-subset membership")
    counts = subsets.groupby("set_name").size()
    for name, n in [("IFN-list overlap", 79), ("Outside IFN lists", 309), ("Neutrophil degranulation", 25)]:
        if counts.get(name) != n:
            raise ValueError("Incorrect nominal membership: " + name)
    return dict(members=388, up=280, down=108, nominal_subset_sizes=counts.to_dict())
