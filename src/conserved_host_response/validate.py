from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy import stats as _scipy_stats
except ImportError:
    _scipy_stats = None

from .io import read_numeric_table, write_table


def _bh(values):
    values = np.asarray(values, dtype=float)
    result = np.full_like(values, np.nan)
    valid = np.isfinite(values)
    p = values[valid]
    if not len(p):
        return result
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result[valid] = restored
    return result


def _normal_two_sided(z):
    return math.erfc(abs(float(z)) / math.sqrt(2))


def _exact_sign_test(delta):
    delta = np.asarray(delta)
    delta = delta[delta != 0]
    n = len(delta)
    if not n:
        return 1.0
    k = min(int(np.sum(delta > 0)), int(np.sum(delta < 0)))
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2**n)
    return min(1.0, 2 * tail)


def _independent_effect(case, control):
    case = np.asarray(case, dtype=float)
    control = np.asarray(control, dtype=float)
    case = case[np.isfinite(case)]
    control = control[np.isfinite(control)]
    n1, n0 = len(case), len(control)
    if min(n1, n0) < 2:
        return None
    s1, s0 = np.std(case, ddof=1), np.std(control, ddof=1)
    pooled = math.sqrt(((n1 - 1) * s1**2 + (n0 - 1) * s0**2) / (n1 + n0 - 2))
    d = (np.mean(case) - np.mean(control)) / pooled if pooled > 0 else np.nan
    correction = 1 - 3 / (4 * (n1 + n0) - 9)
    g = correction * d
    variance = (n1 + n0) / (n1 * n0) + g**2 / (2 * (n1 + n0 - 2))
    se = math.sqrt(variance)
    if _scipy_stats is not None:
        p_value = float(_scipy_stats.ttest_ind(case, control, equal_var=False, nan_policy="omit").pvalue)
        test_method = "Welch_t"
    else:
        standard_error = math.sqrt(s1**2 / n1 + s0**2 / n0)
        z_value = (np.mean(case) - np.mean(control)) / standard_error if standard_error > 0 else 0
        p_value = _normal_two_sided(z_value)
        test_method = "Welch_normal_approximation"
    return {
        "n_case": n1,
        "n_control": n0,
        "mean_case": float(np.mean(case)),
        "mean_control": float(np.mean(control)),
        "effect_type": "Hedges_g",
        "effect": float(g),
        "ci95_low": float(g - 1.96 * se),
        "ci95_high": float(g + 1.96 * se),
        "p_value": p_value,
        "test_method": test_method,
    }


def _paired_effect(frame, score_column, group_column, case_label, control_label, pair_column):
    subset = frame[[pair_column, group_column, score_column]].dropna()
    counts = subset.groupby([pair_column, group_column]).size()
    if (counts > 1).any():
        raise ValueError(f"More than one sample per pair/group for {score_column}")
    wide = subset.pivot(index=pair_column, columns=group_column, values=score_column)
    if case_label not in wide.columns or control_label not in wide.columns:
        return None
    delta = (wide[case_label] - wide[control_label]).dropna().to_numpy(dtype=float)
    if len(delta) < 2:
        return None
    mean_delta = float(np.mean(delta))
    sd_delta = float(np.std(delta, ddof=1))
    se = sd_delta / math.sqrt(len(delta))
    if _scipy_stats is not None:
        try:
            paired_p = float(_scipy_stats.wilcoxon(delta, zero_method="wilcox").pvalue)
        except ValueError:
            paired_p = 1.0
        test_method = "Wilcoxon_signed_rank"
    else:
        paired_p = _exact_sign_test(delta)
        test_method = "exact_sign_test"
    return {
        "n_case": len(delta),
        "n_control": len(delta),
        "mean_case": float(wide.loc[wide.index.intersection(wide.dropna().index), case_label].mean()),
        "mean_control": float(wide.loc[wide.index.intersection(wide.dropna().index), control_label].mean()),
        "effect_type": "paired_dz",
        "effect": mean_delta / sd_delta if sd_delta > 0 else np.nan,
        "ci95_low": mean_delta - 1.96 * se,
        "ci95_high": mean_delta + 1.96 * se,
        "p_value": paired_p,
        "test_method": test_method,
    }


def validate_scores(
    score_path: str | Path,
    metadata_path: str | Path,
    output_dir: str | Path,
    group_column: str,
    case_label: str,
    control_label: str,
    pair_column: str | None = None,
    score_columns: list[str] | None = None,
    filter_column: str | None = None,
    filter_value: str | None = None,
):
    scores = read_numeric_table(score_path)
    metadata = read_numeric_table(metadata_path)
    if "sample_id" not in scores.columns or "sample_id" not in metadata.columns:
        raise ValueError("Scores and metadata must both include sample_id")
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    scores["sample_id"] = scores["sample_id"].astype(str)
    frame = metadata.merge(scores, on="sample_id", how="inner", validate="one_to_one")
    if filter_column:
        if filter_column not in frame.columns:
            raise ValueError(f"Metadata lacks filter column: {filter_column}")
        expected = str(filter_value).strip().lower()
        frame = frame[frame[filter_column].astype(str).str.strip().str.lower() == expected].copy()
    if group_column not in frame.columns:
        raise ValueError(f"Metadata lacks group column: {group_column}")
    if score_columns is None:
        score_columns = [
            column for column in scores.columns if column != "sample_id" and pd.api.types.is_numeric_dtype(scores[column])
        ]
    pairing_audit = None
    if pair_column:
        if pair_column not in frame.columns:
            raise ValueError(f"Metadata lacks pair column: {pair_column}")
        pair_frame = frame[[pair_column, group_column]].dropna().copy()
        pair_frame[pair_column] = pair_frame[pair_column].astype(str)
        pair_frame[group_column] = pair_frame[group_column].astype(str)
        pair_rows = []
        for pair_id, group in pair_frame.groupby(pair_column, sort=True):
            n_case = int((group[group_column] == case_label).sum())
            n_control = int((group[group_column] == control_label).sum())
            if n_case > 1 or n_control > 1:
                status = "duplicated"
            elif n_case == 1 and n_control == 1:
                status = "complete"
            else:
                status = "incomplete"
            pair_rows.append(
                {
                    "pair_id": pair_id,
                    "n_case": n_case,
                    "n_control": n_control,
                    "pair_status": status,
                    "included_in_paired_effect": status == "complete",
                }
            )
        pairing_audit = pd.DataFrame(pair_rows)
        if (pairing_audit["pair_status"] == "duplicated").any():
            duplicated = pairing_audit.loc[pairing_audit["pair_status"] == "duplicated", "pair_id"].tolist()
            raise ValueError(f"Duplicated case/control observations within pairs: {', '.join(duplicated)}")

    rows = []
    for column in score_columns:
        if pair_column:
            result = _paired_effect(frame, column, group_column, case_label, control_label, pair_column)
        else:
            result = _independent_effect(
                frame.loc[frame[group_column].astype(str) == case_label, column],
                frame.loc[frame[group_column].astype(str) == control_label, column],
            )
        if result:
            rows.append({"score": column, "case": case_label, "control": control_label, **result})
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary["fdr"] = _bh(summary["p_value"])
        summary["direction"] = np.where(summary["effect"] > 0, "positive", np.where(summary["effect"] < 0, "negative", "zero"))
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_table(summary, output_dir / "validation_effects.tsv")
    if pairing_audit is not None:
        write_table(pairing_audit, output_dir / "pairing_audit.tsv")
        complete_pairs = int((pairing_audit["pair_status"] == "complete").sum())
        incomplete_pairs = int((pairing_audit["pair_status"] == "incomplete").sum())
    else:
        complete_pairs = 0
        incomplete_pairs = 0
    report = [
        "# Post-freeze validation report",
        "",
        f"- Comparison: {case_label} versus {control_label}",
        f"- Design: {'paired' if pair_column else 'independent'}",
        f"- Matched samples: {len(frame)}",
        f"- Population filter: {filter_column}={filter_value}" if filter_column else "- Population filter: none",
        f"- Exact complete pairs: {complete_pairs}; incomplete pairs excluded: {incomplete_pairs}" if pair_column else "- Exact complete pairs: not applicable",
        f"- Score components evaluated: {len(summary)}",
        "",
        "This report describes association after framework freezing. It does not establish diagnosis, prediction, calibration, decision utility, or causality.",
    ]
    (output_dir / "validation_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return summary
