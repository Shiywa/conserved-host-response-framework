"""Legacy v0.1.0 compatibility only; not evidence for current CRRM388 scoring."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd

from conserved_host_response.audit import audit_manifest
from conserved_host_response.construct import construct_framework
from conserved_host_response.freeze import freeze_framework
from conserved_host_response.score import score_framework
from conserved_host_response.validate import validate_scores


class LegacyV010ToolkitTests(unittest.TestCase):
    def test_audit_detects_shared_controls_and_study_family(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = pd.DataFrame(
                [
                    ["A", "family", "construction", "infection", "baseline", "human", "blood", "DOI:A"],
                    ["B", "family", "construction", "infection", "baseline", "human", "blood", "DOI:B"],
                    ["C", "other", "severity_association", "severity", "", "human", "blood", "DOI:C"],
                ],
                columns=[
                    "dataset_id", "study_family", "role", "biological_role",
                    "shared_control_group", "species", "tissue", "source_reference",
                ],
            )
            path = root / "manifest.tsv"
            manifest.to_csv(path, sep="\t", index=False)
            findings = audit_manifest(path, root / "audit")
            statuses = {row["check"]: row["status"] for row in findings}
            self.assertEqual(statuses["Construction study-family independence"], "WARN")
            self.assertEqual(statuses["Shared construction controls"], "WARN")
            self.assertEqual(statuses["Outcome and cellular data excluded from construction"], "PASS")

    def test_construct_freeze_score_validate(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            evidence = pd.DataFrame(
                [
                    ["G1", "A", "human_anchor", "up", True, True, True, "M1"],
                    ["G1", "B", "animal", "up", True, False, True, "M1"],
                    ["G2", "A", "human_anchor", "down", True, True, True, "M2"],
                    ["G2", "B", "animal", "up", True, False, True, "M2"],
                ],
                columns=["gene", "dataset_id", "system", "direction", "detected", "primary_deg", "supported", "module"],
            )
            evidence_path = root / "evidence.tsv"
            evidence.to_csv(evidence_path, sep="\t", index=False)
            config = {
                "framework": {"name": "toy", "disease": "toy", "version": "0.1.0"},
                "columns": {"gene": "gene", "module": "module", "direction": "direction"},
                "construction": {
                    "anchor_datasets": ["A"],
                    "conservation_systems": ["animal"],
                    "minimum_primary_anchor_datasets": 1,
                    "minimum_conservation_systems": 1,
                },
                "scoring": {
                    "method": "mean_gene_zscore",
                    "zscore_ddof": 1,
                    "invert_down_genes": False,
                    "missing_gene_policy": "available_gene_mean",
                },
            }
            config_path = root / "config.yaml"
            import json
            config_path.write_text(json.dumps(config), encoding="utf-8")
            _, selected = construct_framework(evidence_path, config_path, root / "construction")
            self.assertEqual(selected["gene"].tolist(), ["G1"])

            freeze_framework(root / "construction/selected_gene_table.tsv", config_path, root / "frozen")
            expression = pd.DataFrame(
                {"gene": ["G1", "G1", "OTHER"], "S1": [1, 3, 9], "S2": [3, 5, 9], "S3": [5, 7, 9], "S4": [7, 9, 9]}
            )
            expression_path = root / "expression.csv"
            expression.to_csv(expression_path, index=False)
            scores, coverage, _ = score_framework(
                expression_path,
                root / "frozen/frozen_gene_manifest.tsv",
                root / "frozen/frozen_scoring_spec.yaml",
                root / "scores",
            )
            self.assertEqual(len(scores), 4)
            self.assertEqual(coverage.loc[0, "genes_scored"], 1)
            self.assertTrue(np.isclose(scores["framework_total"].mean(), 0))

            metadata = pd.DataFrame(
                {"sample_id": ["S1", "S2", "S3", "S4"], "condition": ["control", "control", "case", "case"]}
            )
            metadata_path = root / "metadata.csv"
            metadata.to_csv(metadata_path, index=False)
            summary = validate_scores(
                root / "scores/sample_framework_scores.tsv",
                metadata_path,
                root / "validation",
                "condition",
                "case",
                "control",
            )
            effect = summary.loc[summary["score"] == "framework_total", "effect"].iloc[0]
            self.assertGreater(effect, 0)

    def test_freeze_rejects_duplicate_genes_and_is_content_deterministic(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            config = {
                "framework": {"name": "toy_non_rsv", "disease": "toy_non_rsv", "version": "1"},
                "columns": {"gene": "gene", "module": "module", "direction": "direction"},
                "scoring": {"method": "mean_gene_zscore", "invert_down_genes": False},
            }
            import json
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            duplicate = pd.DataFrame({"gene": ["g1", "G1"], "module": ["M", "M"], "direction": ["up", "up"]})
            duplicate_path = root / "duplicate.tsv"
            duplicate.to_csv(duplicate_path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "unique"):
                freeze_framework(duplicate_path, config_path, root / "duplicate_frozen")

            unique = pd.DataFrame({"gene": ["G1", "G2"], "module": ["M1", "M2"], "direction": ["up", "down"]})
            unique_path = root / "unique.tsv"
            unique.to_csv(unique_path, sep="\t", index=False)
            _, _, first = freeze_framework(unique_path, config_path, root / "frozen_a")
            _, _, second = freeze_framework(unique_path, config_path, root / "frozen_b")
            self.assertEqual(first["framework_sha256"], second["framework_sha256"])
            self.assertEqual((root / "frozen_a/SHA256SUMS").read_text(), (root / "frozen_b/SHA256SUMS").read_text())

    def test_score_reports_missing_and_invariant_frozen_genes_without_inversion(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = pd.DataFrame(
                {"gene": ["UP", "DOWN", "MISSING"], "module": ["M", "M", "M"], "direction": ["up", "down", "up"]}
            )
            manifest_path = root / "manifest.tsv"
            manifest.to_csv(manifest_path, sep="\t", index=False)
            spec = root / "spec.json"
            spec.write_text('{"invert_down_genes": false, "zscore_ddof": 1}', encoding="utf-8")
            expression = pd.DataFrame({"gene": ["UP", "DOWN", "INVARIANT"], "S1": [1, 1, 5], "S2": [3, 3, 5], "S3": [5, 5, 5]})
            expression_path = root / "expression.csv"
            expression.to_csv(expression_path, index=False)
            scores, coverage, z = score_framework(expression_path, manifest_path, spec, root / "scores")
            total = coverage.loc[coverage["component"] == "framework_total"].iloc[0]
            self.assertEqual(int(total["genes_detected"]), 2)
            self.assertEqual(int(total["genes_scored"]), 2)
            self.assertAlmostEqual(float(z.loc["DOWN", "S1"]), -1.0)
            self.assertTrue(np.isfinite(scores["framework_total"]).all())

    def test_anchor_presence_adapter_column_preserves_primary_deg_semantics(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            evidence = pd.DataFrame(
                [
                    ["G1", "A", "anchor", "up", True, False, True, "M", True],
                    ["G1", "B", "anchor", "up", True, False, True, "M", True],
                    ["G1", "C", "animal", "up", True, True, True, "M", True],
                ],
                columns=["gene", "dataset_id", "system", "direction", "detected", "primary_deg", "supported", "module", "evidence_positive"],
            )
            evidence_path = root / "evidence.tsv"
            evidence.to_csv(evidence_path, sep="\t", index=False)
            import json
            config_path = root / "config.json"
            config_path.write_text(json.dumps({
                "construction": {"anchor_datasets": ["A", "B"], "anchor_presence_column": "evidence_positive", "conservation_systems": ["animal"], "minimum_primary_anchor_datasets": 0, "minimum_conservation_systems": 1}
            }), encoding="utf-8")
            _, selected = construct_framework(evidence_path, config_path, root / "construction")
            self.assertEqual(selected["gene"].tolist(), ["G1"])

    def test_construct_excludes_declared_unassigned_modules(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            evidence = pd.DataFrame(
                [
                    ["G1", "A", "anchor", "up", True, True, True, "Unassigned"],
                    ["G1", "B", "animal", "up", True, True, True, "Unassigned"],
                ],
                columns=["gene", "dataset_id", "system", "direction", "detected", "primary_deg", "supported", "module"],
            )
            evidence_path = root / "evidence.tsv"
            evidence.to_csv(evidence_path, sep="\t", index=False)
            import json
            config_path = root / "config.json"
            config_path.write_text(json.dumps({
                "construction": {"anchor_datasets": ["A"], "conservation_systems": ["animal"], "minimum_primary_anchor_datasets": 1, "minimum_conservation_systems": 1, "excluded_modules": ["Unassigned"]}
            }), encoding="utf-8")
            audit, selected = construct_framework(evidence_path, config_path, root / "construction")
            self.assertTrue(selected.empty)
            self.assertIn("module_assignment", audit.loc[0, "failure_reasons"])

    def test_audit_blocks_clinical_construction_and_preserves_related_series_warning(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = pd.DataFrame(
                [
                    ["A", "family", "construction", "infection", "baseline", "human", "blood"],
                    ["B", "family", "construction", "infection", "baseline", "human", "blood"],
                    ["C", "other", "construction", "clinical outcome", "", "human", "blood"],
                ],
                columns=["dataset_id", "study_family", "role", "biological_role", "shared_control_group", "species", "tissue"],
            )
            path = root / "manifest.tsv"
            manifest.to_csv(path, sep="\t", index=False)
            statuses = {x["check"]: x["status"] for x in audit_manifest(path, root / "audit")}
            self.assertEqual(statuses["Construction study-family independence"], "WARN")
            self.assertEqual(statuses["Shared construction controls"], "WARN")
            self.assertEqual(statuses["Outcome and cellular data excluded from construction"], "FAIL")

    def test_paired_validation_audits_incomplete_pairs_and_rejects_duplicates(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            scores = pd.DataFrame({"sample_id": ["c1", "r1", "c2", "c3", "r3"], "framework_total": [0.0, 2.0, 0.3, 1.0, 2.5]})
            metadata = pd.DataFrame({"sample_id": ["c1", "r1", "c2", "c3", "r3"], "condition": ["control", "case", "control", "control", "case"], "pair": ["P1", "P1", "P2", "P3", "P3"]})
            scores_path, meta_path = root / "scores.tsv", root / "meta.tsv"
            scores.to_csv(scores_path, sep="\t", index=False)
            metadata.to_csv(meta_path, sep="\t", index=False)
            summary = validate_scores(scores_path, meta_path, root / "validation", "condition", "case", "control", "pair")
            self.assertEqual(int(summary.loc[0, "n_case"]), 2)
            audit = pd.read_csv(root / "validation/pairing_audit.tsv", sep="\t")
            self.assertEqual(set(audit["pair_status"]), {"complete", "incomplete"})

            duplicated = pd.concat([metadata, pd.DataFrame({"sample_id": ["r2"], "condition": ["case"], "pair": ["P1"]})], ignore_index=True)
            duplicated_scores = pd.concat([scores, pd.DataFrame({"sample_id": ["r2"], "framework_total": [1.0]})], ignore_index=True)
            duplicated.to_csv(meta_path, sep="\t", index=False)
            duplicated_scores.to_csv(scores_path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "Duplicated"):
                validate_scores(scores_path, meta_path, root / "duplicate_validation", "condition", "case", "control", "pair")


if __name__ == "__main__":
    unittest.main()
